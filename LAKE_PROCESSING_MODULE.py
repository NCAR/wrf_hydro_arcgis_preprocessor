#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      ksampson
#
# Created:     18/11/2019
# Copyright:   (c) ksampson 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------

'''
11/20/2019:
    This module is intended for use in pre-processing watebodies for inclusion
    into a WRF-Hydro domain.
'''

# Import Python core modules
import sys
sys.dont_write_bytecode = True

import os
import time
import csv

# Import 3rd-party modules
import numpy

# Pure python method of getting unique sums - from http://stackoverflow.com/questions/4373631/sum-array-by-number-in-numpy
from itertools import izip                                                      # Used in the group_min function
from operator import itemgetter                                                 # Used in the group_min function
import collections                                                              # Used in the group_min function

###################################################
# Globals
datestr = time.strftime("%Y_%m_%d")                                             # Date string to append to output files
FLID = 'COMID'                                                                  # Field name for the flowline IDs
LakeAssoc = 'WBAREACOMI'                                                        # Field name containing link-to-lake association
NoDownstream = [None, -1, 0]                                                    # A list of possible values indicating that the downstream segment is invalid (ocean, network endpoint, etc.)
LkNodata = 0                                                                    # Set the nodata value for the link-to-lake association
hydroSeq = 'HydroSeq'
save_Lake_Link_Type_arr = True                                                  # Switch for saving the Lake_Link_Type_arr array to CSV
##################################################

class TeeNoFile(object):
    '''
    Send print statements to a log file:
    http://web.archive.org/web/20141016185743/https://mail.python.org/pipermail/python-list/2007-May/460639.html
    https://stackoverflow.com/questions/11124093/redirect-python-print-output-to-logger/11124247
    '''
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self
    def close(self):
        if self.stdout is not None:
            sys.stdout = self.stdout
            self.stdout = None
        if self.file is not None:
            self.file.close()
            self.file = None
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)
    def flush(self):
        self.file.flush()
        self.stdout.flush()
    def __del__(self):
        self.close()

# Functions
def printMessages(arcpy, messages):
    '''provide a list of messages and they will be printed and returned as tool messages.'''
    for i in messages:
        print(i)
        #arcpy.AddMessage(i)

def group_min(l, g):
    '''Function for gathering minimum value from a set of groups'''
    groups = collections.defaultdict(int)
    for li, gi in izip(l, g):
        if li <= groups[gi] or groups[gi]==0:
            groups[gi] = li
    return groups

def Waterbody_SpatialJoin(arcpy, Flowline, Waterbody, fields, NJ=True, outDir=None):
    '''This function performs a simple intersection between the flowline network
    and the waterbodies. The output is a dictionary of waterbodies touched by
    each flowline.

    If optional parameter NJ is True, the spatial join will be performed. If False,
    the Waterbody is assumed to be the already-joined flowlines to waterbodies
    feature class.'''

    tic1 = time.time()
    printMessages(arcpy, ['      Fields requested: {0}, {1}'.format(fields[0], fields[1])])
    if NJ:
        printMessages(arcpy, ['      Starting to Intersect flowline network with waterbodies.'])
        if not outDir:
            OutFC = 'in_memory/IntersectFC'
        else:
            OutFC = os.path.join(outDir, 'IntersectFC.shp')

         # The Spatial Join method used below will create duplicate records of the Flowlines if more than one Waterbody is intersected by that flowline
        # Other match options: ["WITHIN","INTERSECT"]
        arcpy.SpatialJoin_analysis(target_features=Flowline, join_features=Waterbody, out_feature_class=OutFC, join_operation="JOIN_ONE_TO_MANY", join_type="KEEP_COMMON", match_option="HAVE_THEIR_CENTER_IN")
        intersectarr = arcpy.da.FeatureClassToNumPyArray(OutFC, fields)
        arcpy.Delete_management(OutFC)
    else:
        printMessages(arcpy, ['      Using previously intersected feature class of flowlines to waterbodies: {0}.'.format(Waterbody)])
        intersectarr = arcpy.da.FeatureClassToNumPyArray(Waterbody, fields)

        # Remove elements where there is no flowline:lake association (mirroring the "KEEP_COMMON" option in the Spatial Join step above).
        intersectarr = intersectarr[intersectarr[fields[-1]]!=LkNodata]

    # Populate dictionary of all flowline/lake intersections
    WaterbodyDict = {}
    print(intersectarr)
    for x in intersectarr:
        try:
            WaterbodyDict[x[fields[0]]] += [x[fields[1]]]
        except KeyError:
            WaterbodyDict[x[fields[0]]] = [x[fields[1]]]
    del intersectarr

    printMessages(arcpy, ['      Size of Dictionary: {0} flowlines comprising {1} lake intersections'.format(len(WaterbodyDict), sum([len(WaterbodyDict[i]) for i in WaterbodyDict]))])
    printMessages(arcpy, ['      Found {0} flowlines with connectivity to lakes.'.format(len(WaterbodyDict))])
    if NJ:
        printMessages(arcpy, ['      Finished intersecting flowline network with waterbodies in {0:3.2f}s'.format(time.time()-tic1)])
    else:
        printMessages(arcpy, ['      Finished using previously intersected feature class to determine flowline-to-waterbody associations in {0: 3.2f}s'.format(time.time()-tic1)])
    return WaterbodyDict

def set_problem(problem_lakes, LakeID, problemstr):
    '''This function is used to add elements to a dictionary where the values are
    a list and the keys may or may not exist. This speeds up adding elements to
    an existing dictionary.'''
    try:
        problem_lakes[LakeID] += [problemstr]
    except KeyError:
        problem_lakes[LakeID] = [problemstr]
    return problem_lakes

def get_inflow_segs(FLWBarr, localLk, FromComIDs, FromSegs):
    '''This function will list all inflow segments to a particular lake.
    Inputs:
        1) FLWBarr: Array containing flowlines and the associated lake
        2) localLk: The ID of the lake to examine
        3) FromComIDs: Dictionary of the From:To relationship for flowlines
        4) FromSegs: Dictionary of the upstream links for each flowline
    Globals:
        1) LakeAssoc: The field name used to associate flowlines with lakes.
    Output:
        inflows: Array of links that are inflows to the lake
    '''

    Lake_Links = FLWBarr[FLWBarr[LakeAssoc]==localLk]                           # Find all of the links inside the lake

    # Find all headwater segments as start points
    ToSeg_keys = Lake_Links[FLID]                                               # Array of all lake links for this local lake
    ToSeg_vals = numpy.array([FromComIDs.get(key) for key in ToSeg_keys])       # Array of all downstream links for all of the lake links
    ups = ToSeg_keys[~numpy.in1d(ToSeg_keys, ToSeg_vals)].tolist()              # List of all 'headwater' segments for this lake
    del ToSeg_keys, ToSeg_vals, Lake_Links

    # Gather all inflow segments from lake 'headwater' segments
    inflows = [FromSegs[up] for up in ups if up in FromSegs]                    # list of all contributing segments to the headwater segments for this lake
    inflows = numpy.array([item for sublist in inflows for item in sublist])    # Convert from list of lists to flat list
    del ups
    return inflows

def check_downstream(uplink, localLk, FromComIDs, Lake_LinkDict):
    '''This function will keep looking downstream to see if the link will flow
    back into the same lake and will stop if it encounteres a different lake.
    Inputs:
        1) uplink: ID of the flowline to examine
        2) localLk: ID of the lake under examination
        3) FromComIDs: Dictionary of the From:To relationship for flowlines.
        4) Lake_LinkDict: Dictionary of the Flowline:Lake relationships.
    Outputs:
        result: Boolean (True/False) of if there is a flow loop
        downlinks: List of links to re-associate with this lake
        counter: Number of flowlines downstream iterated over to find the loop
        '''

    # Setup initial values
    result = False                                                              # Default condition
    counter = 0                                                                 # Initiate counter
    downlinks = [uplink]                                                        # Initiate list of links including supplied link
    while not result:
        counter += 1                                                            # Advance the counter
        down = FromComIDs.get(uplink)                                           # Move down one segment
        if down in NoDownstream:
            break                                                               # Reached end of flow network. Break without returning True
        else:
            downlinks.append(down)                                              # Add this link to the list
            if down in Lake_LinkDict:
                if Lake_LinkDict.get(down) == localLk:                          # This means that the links flows eventually back into the same lake
                    result = True
                    break
                else:
                    break                                                       # Link is associated with another lake. Break without returning true.
        uplink=down                                                             # Make this segment the upstream segment
    #print '      Looked downstream %s links from %s in lake %s.' %(counter, uplink, localLk)
    return result, downlinks, counter

def get_lake_routing_info(FLWBarr, localLk, Lake_LinkDict, accum_val, FromComIDs, FromSegs):
    '''This function will list all inflow segments to a "lake" as well as additional
    datasets to assist in routing through the lake.

    Inputs:
        1) FLWBarr: Array containing flowlines and the associated lake
        2) localLk: The ID of the lake to be examined
        3) Lake_LinkDict: Dictionary of all link:lake associations in the domain (pass through to check_downstream function)
        4) accum_val: The value given to each segment before performing lake routed accumulations
        5) FromComIDs: Dictionary of the From:To relationship for flowlines
        6) FromSegs: Dictionary of the upstream links for each flowline
    Globals:
        1) LakeAssoc: A field describing which lake ComID a flowline is associated with
        2) FLID: The flowline ComID
    Outputs:
        1) Lake_LinksList: List of lake link ComIDs
        2) inflows: list of all contributing segments to the headwater segments for this lake
        3) SegVals: A dictionary initializing the accumulation of flow for each link in this lake
        4) SegVals2: A dictionary initializing the accumulation of flow for each link in this lake. This dict is used to find the outlet.
        5) ups: List of all 'headwater' segments for this lake
        6) newLakeLinks: Dictionary to store newly found flowline:lake associations

    If UseAll == True, the script will look outside of the link:lake associations
    to find segments that may exit the lake and then re-enter. If the flow re-enters
    the lake, it will include those segments.

    Input array FLWBarr must have fields [LakeAssoc, FLID], as defined in global variables.
    '''

    # Local variables
    UseAll = True                                                               # Switch to look outside of just the links that spatially intersect the lake

    # Use all links that are associated with lakes (either through spatial join or attribute join)
    Lake_Links = FLWBarr[FLWBarr[LakeAssoc]==localLk]                           # Find all of the links associated with the current lake
    Lake_LinksList = Lake_Links[FLID].tolist()                                  # List of lake link ComIDs

    # Find all headwater segments as start points
    ToSeg_keys = numpy.unique(Lake_Links[FLID])                                 # Array of all unique lake links for this local lake
    ToSeg_vals = numpy.unique(numpy.array([FromComIDs.get(key) for key in ToSeg_keys])) # Array of all unique downstream links for all of the lake links
    ups = ToSeg_keys[~numpy.in1d(ToSeg_keys, ToSeg_vals)].tolist()              # List of all 'headwater' segments for this lake

    # Find any non-assigned lake segments in the lake flow network by iterating down the network
    newLakeLinks = {}                                                           # Dictionary to store new flowline:lake associations
    counter = 0                                                                 # Initiate counter
    seen = []                                                                   # Initiate list of seen links
    if UseAll:
        # Search down the network of links in this lake to find the outlet
        changes = 1                                                             # Initiate the change detector in order to move down the network in the first iteration
        iteration = 0                                                           # Initiate the iteration counter
        while changes > 0:
            downs = list(set([FromComIDs.get(key) for key in ups]))             # Get unique list of downstream flowline ComIDs

            # Add any link:lake associations that were not already present. This is an attempt to eliminate flow looping out of lakes and back in.
            checklinks = list(set([item for item in downs if item not in Lake_LinksList]))  # Unique list of links that are downstream of lake links but not associated with that lake

            # Check these suspect links for downstream connectivity back to the lake
            if len(checklinks) > 0:
                for uplink in checklinks:
                    if uplink in seen:
                        continue
                    result, downlinks, counter1 = check_downstream(uplink, localLk, FromComIDs, Lake_LinkDict)
                    if result:
                        # If result == True, then this lake flows out and then back into itself
                        #print '      Looked downstream %s links from %s in lake %s.' %(counter1, uplink, localLk)
                        Lake_LinksList += downlinks                             # Add these links to the lake association
                        newLakeLinks.update({item:localLk for item in downlinks})   # Store these new flowline:lake associations
                    else:
                        # This downstream segment never returns to the lake. Disassociate.
                        Lake_LinksList = [item for item in Lake_LinksList if item not in downlinks[1:]] # Remove these segments from association with this lake
                    downs.remove(uplink)                                        # Remove this segment so that the loop can complete.
                    seen.append(uplink)

            downs = [item for item in downs if item in Lake_LinksList]          # Remove downstream segments that are not associated with this lake
            for key in downs:
                if key in NoDownstream:
                    continue                                                    # Downstream segment is not a valid ComID
            ups = downs                                                         # Change the upstream segments to examine in the next iteration to the downstream ones from the previous iteration
            changes = sum([1 for item in downs if item is not None])            # Number of changes in this iteration
            iteration += 1                                                      # Add one for each level
        counter += 1                                                            # Advance the counter

    Lake_LinksList = list(set(Lake_LinksList))                                  # Remove duplicate link IDs

    # Set local accumulation values to default accum_val (1)
    SegVals = {key:accum_val for key in Lake_LinksList}                         # Set all initial values to accum_val
    SegVals2 = {key:accum_val for key in Lake_LinksList}                        # Set all initial values to accum_val

    # Test to regenerate upstream inflows to account for newly added links (2/26/2018)
    ToSeg_keys = numpy.array(Lake_LinksList)                                    # Array of all unique lake links for this local lake
    ToSeg_vals = numpy.unique(numpy.array([FromComIDs.get(key) for key in ToSeg_keys])) # Array of all unique downstream links for all of the lake links

    # Re-gather the upstream segments for this lake
    ups = ToSeg_keys[~numpy.in1d(ToSeg_keys, ToSeg_vals)].tolist()              # List of all 'headwater' segments for this lake
    del ToSeg_keys, ToSeg_vals                                                  # Free up memory

    # Gather all inflow segments from lake 'headwater' segments. This may not account for oCONUS contributing links
    #inflows = [FromSegs[up] for up in ups if up in FromSegs]                    # list of all contributing segments to the headwater segments for this lake
    #inflows = numpy.unique(numpy.array([item for sublist in inflows for item in sublist]))  # Convert from list of lists to unique flat list

    # 2/27/2018: We can more closely replicate WRF-Hydro's TYPE=3 segments if we say inflows are all segments that flow into any lake segment that are not already accounted for
    inflows = []
    for link in Lake_LinksList:
        upsegs = FromSegs.get(link)                                             # All upstream segments for each lake link
        if upsegs is None:
            continue
        upsegs2 = [item for item in upsegs if item not in Lake_LinksList]       # All contributing links that are not already associated with the lake
        inflows += upsegs2
        del upsegs, upsegs2
    inflows = numpy.array(inflows)
    return Lake_LinksList, inflows, SegVals, SegVals2, ups, newLakeLinks

def Lake_Link_Type(arcpy, FLWBarr, FromComIDs, FLarr, subset=None):
    '''
    This function will assign a link type to each lake.
            3 = Lake Inflow Link
            2 = Internal Lake Link
            1 = Lake Outflow Link (based on accumulated flow in the lake)

    This function sorts the lakes by the minimum HydrSeq of all links associated
    with that lake. Thus, theoretically, lakes can be visited by increasing HydroSeq
    and no downstream searching for lakes should be necessary. ...

    This function will examine each lake and determine the outlets based on an
    accumulation of values through the topology of links within the lake. It will
    identify lakes with (potentially) multiple outlets, as well as internally
    draining lakes. Lakes with multiple outlets according to flow accumulation
    within the lake will have the minor outlets removed.

    2/16/2018: This function finally diagnoses all potential conflicts with WRF-Hydro.
        In the "for LakeID in LakeSeq[FLID]:" loop, any continue statements will
        eliminate that lake from being placed on the flow network. Conditions that
        prevent a lake from functioning in WRF-Hydro are:
            1) The lake has no Lake Link Type of 1 (outlet link). This seems to
               occur when a lake has only one interior link, or when all interior
               links flow directly to a nonexistent outlet link (endorheic).
            2) A lake flows directly into another lake.
            3) Multiple outlet links are defined
    '''

    # Options and defaults
    tic1 = time.time()                                                          # Initiate timer for this function
    #HydroSeq = False                                                            # Switch to use HydroSeq to find the lake outlet
    IterateList = True                                                          # Switch to iterate over lake links to find outlet(s)
    accum_val = 1                                                               # Values given to each segment before accumulation

    # Build default mapping files to store new lake mappings
    problem_lakes = {}                                                          # Problem dictionary
    Old_New_LakeComID = {}                                                      # Dictionary to store the 'old Lake ComID':'new Lake ComID' mapping
    ChainedLakes = {}                                                           # Dictionary to store number of lakes chained together
    Remove_Association = []                                                     # List used to remove a lake association from a flowline (for divergences in lakes)

    # Set output array dtype and field names
    dtype1 = dict(names=(FLID, 'LINK_TYPE', LakeAssoc, 'Accum1', 'Accum2'), formats=('<i4', '<i4', '<i4', '<i4', '<i4'))
    dtype2 = dict(names=(FLID, 'LINK_TYPE', LakeAssoc, 'Accum1', 'Accum2', 'Reason'), formats=('<i4', '<i4', '<i4', '<i4', '<i4', '<i4'))

    # Attempt to use lists instead of arrays
    Lake_Link_Type_COMID = []                                                   # Flowline ComID for flowlines associated with lakes
    Lake_Link_Type = []                                                         # The lake link type for this flowline
    Lake_Link_Type_WBAREACOMI = []                                              # The lake ComID associated with this flowline
    Lake_Link_Type_Accum1 = []                                                  # Added 2/17/2018 to keep track of lake local accumulation
    Lake_Link_Type_Accum2 = []                                                  # Added 2/17/2018 to keep track of lake local accumulation

    # Attempt to keep diagnostic information for lakes that are eliminated by this pre-processor
    Tossed_Lake_Link_Type_arr = numpy.empty(0, dtype=dtype2)                    # Generate new empty array to store tossed lake values

    # Build a dictionary here to pass into the get_lake_routing_info function
    Lake_LinkDict = {item[FLID]:item[LakeAssoc] for item in FLWBarr}            # Generate dictionary of link:lake associations

    # 1) Gather a list of all contributing segments for each segment in the system
    tic2 = time.time()
    FromSegs = {}
    for key,val in FromComIDs.items():
        try:
            FromSegs[val] += [key]
        except KeyError:
            FromSegs[val] = [key]
    printMessages(arcpy, ['      Completed FromSegs dictionary in {0:3.2f} seconds.'.format(time.time()-tic2)])

    # 2) Create a sorting of lakes that will start with the lake which has the lowest HydroSeq value in it's flowlines
    # Use indexing to grab elements common to both arrays while preserving order of one of the arrays (FLWBarr)
    tic2 = time.time()
    commons = FLarr[numpy.in1d(FLarr[FLID], FLWBarr[FLID])]                     # An array of all of the flowlines in the flowline array that are common to the flowline-waterbody association array
    xsorted = numpy.argsort(commons[FLID])                                      # Find the sorted order for the array that is to be sorted
    ypos = numpy.searchsorted(commons[xsorted][FLID], FLWBarr[numpy.in1d(FLWBarr[FLID], FLarr[FLID])][FLID]) # Search only the common values between arrays
    indices = xsorted[ypos]
    commons2 = commons[indices]                                                 # Re-order the input array to match the comparison array order
    del commons, xsorted, ypos, indices

    # Use the group_min function to find the minimum HydroSeq value for each group of WBAREACOMID values
    group_minDict = group_min(commons2[hydroSeq], FLWBarr[LakeAssoc])         # Find the minimum HydroSeq value for this lake
    del commons2                                                                # Free up memory

    # Construct an array to store the Lake COMID and Minimum Hydrosequence
    dtype = dict(names=(FLID, 'minHydroSeq'), formats=('<i4', '<f8'))
    LakeSeq = numpy.array(group_minDict.items(), dtype=dtype)                   # Create array of minimum HydroSeq values for each lake
    del group_minDict, dtype                                                    # Free up memory
    LakeSeq = LakeSeq[LakeSeq[FLID]>-9998]                                      # Clip off -9999, -9998
    LakeSeq.sort(order='minHydroSeq')                                           # Sort by minimum HydroSeq
    LakeSeqsize = LakeSeq.shape[0]                                              # Get the number of elements in the array
    printMessages(arcpy, ['      Completed sorting lakes by minimum HydroSeq in {0:3.2f} seconds.'.format(time.time()-tic2)])

    # 3)  Find outflow link and assign type 1 (but not if it flows into another lake)
    counter = 0                                                                 #
    counter3 = 0                                                                # Counter to keep track of multiple outlet lakes
    counter4 = 0                                                                # Counter to keep track of headwater lakes
    counter5 = 0                                                                # Counter to keep track of terminal lakes
    counter6 = 0                                                                # Counter to keep track of lakes immediately downstream
    counter7 = 0                                                                # Counter to keep track of lakes with outlet that drians to nowhere
    counter8 = 0                                                                # Counter to keep track of lakes with no Type=2 links
    tic2 = time.time()                                                          # Reset the counter to provide progressive print statements
    seen = []                                                                   # Lakes in this list have been 'seen', or reviewed
    if subset is not None:
        LakeSeq = LakeSeq[numpy.in1d(LakeSeq[FLID], subset[FLID])]              # Subset the list of lakes to a subset array if requested
        printMessages(arcpy, ['      Subsetted the lakes from {0} to {1} based on a provided subset array.'.format(LakeSeqsize, LakeSeq.shape[0])])
    for LakeID in LakeSeq[FLID]:
        if LakeID in seen:
            continue                                                            # This lake has already been examined

        # Get initial lake information, including adding looping outflows
        Lake_LinksList, inflows, SegVals, SegVals2, ups, newLakeLinks = get_lake_routing_info(FLWBarr, LakeID, Lake_LinkDict, accum_val, FromComIDs, FromSegs)

        # First, determine if this is a headwater lake. Must set 'continue statement to remove the lake
        if len(inflows) == 0:
            problem_lakes = set_problem(problem_lakes, LakeID, 'Headwater lake, no inflows')
            counter4 += 1                                                       # Advance the counter

        # See if any of the inflow links are on other lakes
        counter2 = 0                                                            # Counter to keep track of lakes immediately upstream for each lake
        num_uplakes = FLWBarr[numpy.in1d(FLWBarr[FLID], inflows)][LakeAssoc].tolist()   # This will be a list of lakes that are associated with this current lake's inflows
        if LakeID in num_uplakes:                                               # Check to make sure lake doesn't flow into itself
            # If the upstream lake is the same ID as the current lake, log the issue
            num_uplakes.remove(LakeID)                                          # If it does, remove that lake to elminate a flow loop
            problem_lakes = set_problem(problem_lakes, LakeID, 'Potential Flow Loop')
        to_alter = list(num_uplakes)                                            # Copy the list to a new list

        # If there are lakes immediately above this one, examine each one, then examine all lakes upstream of those, etc.
        # This will essentially merge any lakes that flow directly into another lake, giving the ID of the most downstream lake to all chained lakes above it.
        while len(num_uplakes) > 0:
            uplake = num_uplakes[0]                                             # Look at the first lake in the list
            if uplake in seen:                                                  # If this lake has already been examined, then skip it
                num_uplakes.remove(uplake)                                      # If this lake has been examined, remove it from this list
                continue
            Old_New_LakeComID[uplake] = LakeID                                  # This will give the upstream lake the ID of the downstream lake it flows directly into
            inflows2 = get_inflow_segs(FLWBarr, uplake, FromComIDs, FromSegs)   # Get all of the upstream lake's inflows and associate these flowlines with the new LakeID
            uplakes2 = FLWBarr[numpy.in1d(FLWBarr[FLID], inflows2)][LakeAssoc].tolist() # Get all the upstream lakes for this lake
            if uplake in uplakes2:
                uplakes2.remove(uplake)                                         # Remove any flow loops
            num_uplakes += uplakes2                                             # Add these new upstream lakes to the list of upstream lakes
            to_alter += uplakes2
            num_uplakes.remove(uplake)                                          # Remove this lake from the list so that iteration will stop when the list is empty
            counter2 += 1                                                       # Advance the counter
            seen.append(uplake)                                                 # Add this upstream lake to the list of examined lakes
            del inflows2, uplakes2                                              # Free up memory

        # Alter all Waterbody ComIDs at once to match the most downstream lake and regenerate lake information for the new, merged lake
        if len(to_alter) > 0:
            ChainedLakes[LakeID] = to_alter                                     # Add list of upstream lakes to this dictionary for this lake
            printMessages(arcpy, ['      Altering {0} lake COMID value(s) to {1}'.format(len(to_alter), LakeID)])
            for item in to_alter:
                problem_lakes = set_problem(problem_lakes, LakeID, 'Upstream segment for {0} belongs to another lake [{1}]'.format(LakeID, item))
            FLWBarr[LakeAssoc][numpy.in1d(FLWBarr[LakeAssoc], numpy.array(to_alter))] = LakeID  # Change all the IDs of the upstream lakes at once to the current lake ID
            Lake_LinksList, inflows, SegVals, SegVals2, ups, newLakeLinks = get_lake_routing_info(FLWBarr, LakeID, Lake_LinkDict, accum_val, FromComIDs, FromSegs)
        Lake_LinksList2 = Lake_LinksList + inflows.tolist()                     # Add all lake links to all inflows

        # Assign the initial LINK_TYPE values based on local lake links and inflow links
        #Lake_Link_Type_local = [3 if item in inflows.tolist() and item not in newLakeLinks else 2 for item in Lake_LinksList2]   # Append Default LINK_TYPE (2) for all items in the Lake_LinksList
        Lake_Link_Type_local = [3 if item in inflows.tolist() else 2 for item in Lake_LinksList2]   # Append Default LINK_TYPE (2) for all items in the Lake_LinksList

        # For each lake, search down the network of links in this lake to find the outlet flowline(s)
        changes = 1                                                             # Initiate the change detector in order to move down the network in the first iteration
        iteration = 0                                                           # Initiate the iteration counter
        while changes > 0:
            downs = list(set([FromComIDs.get(key) for key in ups]))             # Get unique list of downstream flowline ComIDs
            downs = [item for item in downs if item in Lake_LinksList2]         # Remove downstream segments that are not associated with this lake
            for key in downs:
                if key in NoDownstream:
                    # This might be where to catch an endorheic basin terminal lake
                    # This appears to never get triggered, probably because there are no 0 values in Lake_LinksList2
                    #problem_lakes = set_problem(problem_lakes, LakeID, 'Drains to ocean or internally')
                    counter5 += 1                                               # Advance the counter
                else:
                    # This downstream segment is a real flowline in the network
                    ups_local = FromSegs[key]                                   # Get all the upstream links for this link
                    vals = [SegVals.get(val) for val in ups_local if val in SegVals]    # Get the accumulated value for each upstream link (usually 1 for each upstream link)
                    newval = 1 + sum(vals)                                      # Add one to the sum of the accumulated values of all upstream links
                    SegVals[key] = newval                                       # Put the accumulated sum in the SegVals dictionary
                    SegVals2[key] = newval                                      # Put the accumulated sum in the SegVals2 dictionary
                    for up in ups_local:
                        SegVals2[up] = 0                                        # Set all upstream values to 0 in the SegVals2 dictionary
            ups = downs                                                         # Move down the network by one link
            changes = sum([1 for item in downs if item not in NoDownstream])    # Quantify the number of valid downstream links
            iteration += 1                                                      # Add one for each level
        seen.append(LakeID)                                                     # Add this lake to the list of lakes that have been seen already

        # Record the routed and unrouted accumulation values
        Accum1_local = [0 if item in inflows.tolist() else SegVals[item] for item in Lake_LinksList2]   # Added 2/17/2018 to keep track of lake local accumulation
        Accum2_local = [0 if item in inflows.tolist() else SegVals2[item] for item in Lake_LinksList2]  # Added 2/17/2018 to keep track of lake local accumulation

        # Are there multiple outflows? If so, remove all networks contributing to minor outlets
        # Note that this method only chooses the first link ID if multiple links have the same maximum accumulation value
        maxflow = SegVals2.keys()[SegVals2.values().index(max(SegVals2.values()))]  # Find the flowline COMID with the highest accumulation value
        outflows = [key for key,val in SegVals2.items() if val > 0]             # All links with accumulation still in them
        if len(outflows) > 1:
            problem_lakes = set_problem(problem_lakes, LakeID, 'Potentially multiple outlets. Secondary Outlets: {0}'.format(outflows))

            # Iterate upstream from the secondary outlet and remove associations.
            non_outlets = [item for item in outflows if item != maxflow]        # Create list of secondary outlets (outlets with fewer contributing flowlines than the maxflow)
            Remove_local = []                                                   # Initiate list of associations to remove
            ups = non_outlets
            while len(ups) > 0:
                Remove_local += ups                                             # Add the segments from the secondary outlet to the association removal list
                #ups = [FromSegs.get(key) for key in ups if key in Lake_LinksList2] # Move upstream in the flow network, only considering the local lake links
                ups = [FromSegs.get(key) for key in ups if key in Lake_LinksList]   # Move upstream in the flow network, only considering the local lake links (not inflows links)
                ups = [item for sublist in ups if sublist is not None for item in sublist]  # Flatten list including None values
            counter3 += len(Remove_local)                                       # Iterate counter to keep track of lake association removals
            Remove_Association += Remove_local                                  # Add these networks that drain to a secondary lake outlet to the association removal list

            # Before eliminating this lake, save the information in a separate table
            arrsize = Tossed_Lake_Link_Type_arr.shape[0]                        # Get current length of the array
            maskList = numpy.array([True if item in Remove_local else False for item in Lake_LinksList2])   # Create a mask list to mask other lists with
            Tossed_Lake_Link_Type_arr.resize(arrsize + maskList.sum(), refcheck=False)              # Add rows to the recarray to store new data
            Tossed_Lake_Link_Type_arr[FLID][arrsize:] = numpy.array(Lake_LinksList2)[maskList]      # Add Lake_LinksList to the list
            Tossed_Lake_Link_Type_arr['LINK_TYPE'][arrsize:] = numpy.array(Lake_Link_Type_local)[maskList] # Add Lake Link Types
            Tossed_Lake_Link_Type_arr[LakeAssoc][arrsize:] = [LakeID for item in Remove_local]      # Add WBAREACOMI lake associations
            Tossed_Lake_Link_Type_arr['Accum1'][arrsize:] = numpy.array(Accum1_local)[maskList]     # Added 2/17/2018 to keep track of lake local accumulation
            Tossed_Lake_Link_Type_arr['Accum2'][arrsize:] = numpy.array(Accum2_local)[maskList]     # Added 2/17/2018 to keep track of lake local accumulation
            Tossed_Lake_Link_Type_arr['Reason'][arrsize:] = 1                   # Reason: 'Potentially multiple outlets. Secondary Outlets: %s' %(outflows)

            # Use a reverse the masklist to remove these minor network elements draining to false outlets from the rest of the data
            Lake_LinksList2 = numpy.array(Lake_LinksList2)[~maskList].tolist()              # Subset the list to exclude the networks draining to a minor outlet
            Lake_Link_Type_local = numpy.array(Lake_Link_Type_local)[~maskList].tolist()    # Subset the list to exclude the networks draining to a minor outlet
            Accum1_local = numpy.array(Accum1_local)[~maskList].tolist()                    # Subset the list to exclude the networks draining to a minor outlet
            Accum2_local = numpy.array(Accum2_local)[~maskList].tolist()                    # Subset the list to exclude the networks draining to a minor outlet

        # Find if this has an outlet or not
        if FromComIDs.get(maxflow) in NoDownstream:
            # This is either a coastline waterbody or an internally draining (endorheic) lake
            # If these lakes are not given any type=1 outlet links, the lake will not be placed into RouteLink
            problem_lakes = set_problem(problem_lakes, LakeID, 'Link with maximum flow drains to ocean or internally')
            #Lake_Link_Type_local = [1 if com==maxflow else LT for com,LT in zip(Lake_LinksList2, Lake_Link_Type_local)]    # Test to see if lakes that flow to a nodata point can be kept in WRF-Hydro
            counter7 += 1                                                       # Advance the counter
        elif FLWBarr[FLWBarr[FLID]==FromComIDs.get(maxflow)].shape[0] > 0:      # This is a downstream lake
            # This is a real downstream segment
            problem_lakes = set_problem(problem_lakes, LakeID, 'Downstream of outlet segment is a lake segment')
            counter6 += 1                                                       # Advance the counter
        else:
            # This is a legitimate downstream flowline. Alter local lake link types to reflect the outlet link type
            Lake_Link_Type_local = [1 if com==maxflow else LT for com,LT in zip(Lake_LinksList2, Lake_Link_Type_local)]

        # Additional check: Does the lake have only Type=3 and Type=1 links in it? If so, elminate. Added 2/15/2018
        if 1 not in list(set(Lake_Link_Type_local)):
            problem_lakes = set_problem(problem_lakes, LakeID, 'This lake has no type=1 (outlet) link in it. Eliminating...')
            counter8 += 1                                                       # Advance the counter

            # Before eliminating this lake, save the information in a separate table
            arrsize = Tossed_Lake_Link_Type_arr.shape[0]                        # Get current length of the array
            Tossed_Lake_Link_Type_arr.resize(arrsize + len(Lake_LinksList2), refcheck=False)    # Add rows to the recarray to store new data
            Tossed_Lake_Link_Type_arr[FLID][arrsize:] = Lake_LinksList2         # Add Lake_LinksList to the list
            Tossed_Lake_Link_Type_arr['LINK_TYPE'][arrsize:] = Lake_Link_Type_local         # Add Lake Link Types
            Tossed_Lake_Link_Type_arr[LakeAssoc][arrsize:] = [LakeID for item in Lake_LinksList2]         # Add WBAREACOMI lake associations
            Tossed_Lake_Link_Type_arr['Accum1'][arrsize:] = Accum1_local        # Added 2/17/2018 to keep track of lake local accumulation
            Tossed_Lake_Link_Type_arr['Accum2'][arrsize:] = Accum2_local        # Added 2/17/2018 to keep track of lake local accumulation
            Tossed_Lake_Link_Type_arr['Reason'][arrsize:] = 2                   # Reason: 'This lake has no type=1 (outlet) link in it. Eliminating...'
            continue                                                            # Added 2/15/2018 as a test to eliminate these lakes from RouteLink

        # Iterate the lake counter and print a statement every so often
        counter += 1                                                            # Advance the main counter
        if counter % 1000 == 0:
            printMessages(arcpy, ['        {0} lakes processed in {1:3.2f} seconds.'.format(counter, time.time()-tic2)])
            tic2 = time.time()

        # Store link type information in lists
        Lake_Link_Type_COMID += Lake_LinksList2                                 # Add Lake_LinksList to the list
        Lake_Link_Type_WBAREACOMI += [LakeID for item in Lake_LinksList2]       # Append LakeID for all items in the Lake_LinksList
        Lake_Link_Type += Lake_Link_Type_local
        Lake_Link_Type_Accum1 += Accum1_local                                   # Added 2/17/2018 to keep track of lake local accumulation
        Lake_Link_Type_Accum2 += Accum2_local                                   # Added 2/17/2018 to keep track of lake local accumulation
        del Lake_Link_Type_local

    # Subset Lake_Link_Type_arr to just the links in XXX
    Lake_Link_Type_arr = numpy.zeros(len(Lake_Link_Type_COMID), dtype=dtype1)
    Lake_Link_Type_arr[FLID] = Lake_Link_Type_COMID                             # Populate with lake COMIDs
    Lake_Link_Type_arr['LINK_TYPE'] = Lake_Link_Type                            # Add Lake Link Types
    Lake_Link_Type_arr[LakeAssoc] = Lake_Link_Type_WBAREACOMI                   # Add WBAREACOMI lake associations
    Lake_Link_Type_arr['Accum1'] = Lake_Link_Type_Accum1                        # Add unrouted lake accumulation
    Lake_Link_Type_arr['Accum2'] = Lake_Link_Type_Accum2                        # Add routed lake accumulation
    del Lake_Link_Type_COMID, Lake_Link_Type, Lake_Link_Type_WBAREACOMI, dtype1, dtype2, Lake_LinkDict

    # Remove WBAREACOMI association for the multiple outlets (other than that with max flow)
    FLWBarr = FLWBarr[~numpy.in1d(FLWBarr[FLID], numpy.array(Remove_Association))]   # Remove from the array any items that need the association removed
    Lake_Link_Type_arr[LakeAssoc][numpy.in1d(Lake_Link_Type_arr[FLID], numpy.array(Remove_Association))] = 0  # Old Way (left alot of zeros) Remove WBAREACOMI lake associations
    Lake_Link_Type_arr = Lake_Link_Type_arr[Lake_Link_Type_arr[LakeAssoc]!= 0]  # Now remove all lake associations with lake ID = 0 (added 2/14/2018)

    # Clean up and return
    printMessages(arcpy, ['        {0} lake associations eliminated due to multiple outlets.'.format(counter3)])
    printMessages(arcpy, ['        {0} lakes are a headwater lake (no inflows).'.format(counter4)])
    printMessages(arcpy, ['        {0} lakes have an endorheic lake lake condition.'.format(counter5)])
    printMessages(arcpy, ['        {0} lakes eliminated due to having a lake immediately downstream.'.format(counter6)])
    printMessages(arcpy, ['        {0} lakes have an outlet that drains to nowhere.'.format(counter7)])
    printMessages(arcpy, ['        {0} lakes eliminated due to having no type 1 (outlet) segments associated with it.'.format(counter8)])
    printMessages(arcpy, ['      Examined {0} lakes in {1:3.2f} seconds.'.format(len(seen), time.time()-tic1)])
    return Lake_Link_Type_arr, problem_lakes, seen, ChainedLakes, Old_New_LakeComID, FLWBarr, Remove_Association, Tossed_Lake_Link_Type_arr

def main(arcpy, outDir, Flowline, Waterbody, FromComIDs, order, fields, Subset_arr=None, NJ=False, datestr=datestr):
    '''

    '''

    # Setup Logging
    LakeDiagnosticFile = os.path.join(outDir, 'Lake_Preprocssing_Info.txt')
    tee = TeeNoFile(LakeDiagnosticFile, 'w')

    # Set up logging
    tic1 = time.time()
    printMessages(arcpy, ['    Lake module initiated on {0}'.format(time.ctime())])

    if Flowline is not None:
        # This is the normal case. The user wishes to evaluate flowline:waterbody associations for either a flowline feature class
        # which has the associations specified, or perform a spatial join between those flowlines and a Waterbody feature class.
        outDir = outDir     # Will save spatial joined feature class to disk
        #outDir = None       # Will not save spatial joined feature class to disk
        WaterbodyDict = Waterbody_SpatialJoin(arcpy, Flowline, Waterbody, fields, NJ, outDir=outDir)     # Build table of Lake Connectivity for each reach based on simple intersection
    else:
        # This is the case where the user is only wishing to submit a dictionary of flowline:waterbody associations for evaluation
        # Thus, choose Flowline=None and Waterbody=WaterbodyDict in the main() function arguments.
        # Populate dictionary of all flowline/lake intersections
        WaterbodyDict = {item[0]:[item[1]] for item in Waterbody.items()}   # Convert to lists

    # Create an array of all flowlines associated with all lakes from WaterbodyDict
    dtype = dict(names=(FLID,LakeAssoc), formats=('<i4', '<i4'))
    FLWBarr = numpy.array([(item[0], item[1][0]) for item in WaterbodyDict.items()], dtype=dtype)   # Grab the first lake association for any flowline
    printMessages(arcpy, ['    Found {0} unique lake ComIDs from flowline association'.format(numpy.unique(FLWBarr[LakeAssoc]).shape[0])])

    # Gather all link lake types
    Lake_Link_Type_arr, problem_lakes, seen, ChainedLakes, Old_New_LakeComID, FLWBarr, Remove_Association, Tossed_Lake_Link_Type_arr = Lake_Link_Type(arcpy, FLWBarr, FromComIDs, order, subset=Subset_arr)
    unique_lakes = numpy.unique(Lake_Link_Type_arr[LakeAssoc]).shape[0]
    printMessages(arcpy, ['    Found {0} unique lake comID values.'.format(unique_lakes)])
    printMessages(arcpy, ['    Found {0} outlet flowlines.'.format(Lake_Link_Type_arr[Lake_Link_Type_arr['LINK_TYPE']==1].shape[0])])
    printMessages(arcpy, ['    Found {0} internal flowlines.'.format(Lake_Link_Type_arr[Lake_Link_Type_arr['LINK_TYPE']==2].shape[0])])
    printMessages(arcpy, ['    Found {0} contributing flowlines.'.format(Lake_Link_Type_arr[Lake_Link_Type_arr['LINK_TYPE']==3].shape[0])])
    printMessages(arcpy, ['    Problems: {0}'.format(len(problem_lakes))])
    printMessages(arcpy, ['    Chained Lakes: {0}'.format(ChainedLakes.keys())])

    # Write the lake problem file to a CSV format file
    LakeProblemFile = os.path.join(outDir, 'Lake_Problems.csv')
    printMessages(arcpy, ['    Writing Lake_Problems dictionary to file: {0}'.format(LakeProblemFile)])
    with open(LakeProblemFile,'wb') as f:
        w = csv.writer(f)
        w.writerows(problem_lakes.items())

    # Write the dictionary to disk as CSV file that shows which lakes have been merged (added 1/16/2017)
    Old_New_LakeComIDFile = os.path.join(outDir, 'Old_New_LakeComIDs.csv')
    printMessages(arcpy, ['    Writing Lake merging dictionary to file: {0}'.format(Old_New_LakeComIDFile)])
    with open(Old_New_LakeComIDFile,'wb') as f:
        w = csv.writer(f)
        w.writerows(Old_New_LakeComID.items())

    # Save the Lake_Link_Type array
    if save_Lake_Link_Type_arr:
        Lake_Link_Type_File = os.path.join(outDir, 'Lake_Link_Types.csv')
        printMessages(arcpy, ['    Writing Lake_Link_Type array to file: {0}'.format(Lake_Link_Type_File)])
        numpy.savetxt(Lake_Link_Type_File, Lake_Link_Type_arr, fmt='%i', delimiter=",")

    # Save the tossed Lake_Link_Type array
    Lake_Link_Type_File2 = os.path.join(outDir, 'Tossed_Lake_Link_Types.csv')
    printMessages(arcpy, ['    Writing Tossed Lake_Link_Type array to file: {0}'.format(Lake_Link_Type_File2)])
    numpy.savetxt(Lake_Link_Type_File2, Tossed_Lake_Link_Type_arr, fmt='%i', delimiter=",")

    # Re-assign the array to a dictionary type
    #WaterbodyDict = {item[FLID]:item[LakeAssoc] for item in FLWBarr}      # Output to dictionary
    WaterbodyDict = {item[FLID]:item[LakeAssoc] for item in Lake_Link_Type_arr[Lake_Link_Type_arr['LINK_TYPE']<3]}      # Output to dictionary from new Lake Link Type Array (excluding inflow links)

    # Clean up and return
    #del Subset_arr, dtype, problem_lakes, seen, ChainedLakes, Remove_Association, FLWBarr
    printMessages(arcpy, ['    Finished building Lake Association and flowline connectivity tables.  Time elapsed: {0:3.2f} seconds.'.format(time.time()-tic1)])
    tee.close()
    del tee
    return WaterbodyDict, Lake_Link_Type_arr, Old_New_LakeComID

if __name__ == '__main__':
    pass