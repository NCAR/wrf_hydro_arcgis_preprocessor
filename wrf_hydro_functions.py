# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# ** Copyright UCAR (c) 2016
# ** University Corporation for Atmospheric Research(UCAR)
# ** National Center for Atmospheric Research(NCAR)
# ** Research Applications Laboratory(RAL)
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# ** 2016/4/27 10:00:00
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

# --- Import Modules --- #
import sys
sys.dont_write_bytecode = True                                                  # Do not write .pyc files
import os
import csv
import time
import math
import zipfile
from zipfile import ZipFile, ZipInfo
from arcpy.sa import *
import numpy
import netCDF4                                                                  # Part of ArcGIS 10.3, added to this script 6/11/2015
import shutil                                                                   # Added 07/10/2015 for ExmineOutputs
from collections import defaultdict                                             # Added 09/03/2015Needed for topological sorting algorthm
from itertools import takewhile, count                                          # Added 09/03/2015Needed for topological sorting algorthm
# --- End Import Modules --- #

# --- Module Configurations --- #
# --- End Module Configurations --- #

# --- Globals --- #

# Initiate dictionaries of GEOGRID projections and parameters
#   See http://www.mmm.ucar.edu/wrf/users/docs/user_guide_V3/users_guide_chap3.htm#_Description_of_the_1
projdict = {1: 'Lambert Conformal Conic',
            2: 'Polar Stereographic',
            3: 'Mercator',
            6: 'Cylindrical Equidistant'}
CF_projdict = {1: "lambert_conformal_conic",
                2: "polar_stereographic",
                3: "mercator",
                6: "latitude_longitude"}

# Output file types
outNCType = 'NETCDF4_CLASSIC'                                                   # Define the output netCDF version for RouteLink.nc and LAKEPARM.nc
#outNCType = 'NETCDF3_64BIT'                                                     # Set output netCDF format for spatial metdata files
RTfmt = ["NC"]                                                                  # Set output RouteLink table format here ["CSV", "NC"], ["NC"], or ["CSV"]
LKfmt = ["NC", "TBL"]                                                                  # Set output LAKEPARM table format here ["NC", "TBL"], ["NC"], or ["TBL"]

# Output file names
LK_nc = 'LAKEPARM.nc'                                                           # Default Lake parameter table name
RT_nc = 'Route_Link.nc'                                                         # Default Route Link parameter table name
FullDom = 'Fulldom_hires.nc'                                                    # Default Full Domain routing grid nc file
LDASFile = 'GEOGRID_LDASOUT_Spatial_Metadata.nc'                                # Defualt LDASOUT domain grid nc file
GW_ASCII = 'gw_basns_geogrid.txt'                                               # Default Groundwater Basins ASCII grid output

# Other Global Variables
NoDataVal = -9999                                                               # Default NoData value for gridded variables
Threshold = 0.75                                                                # Minimum Lake Size threshold  (sq. km)
walker = 3                                                                      # Number of cells to walk downstream before gaged catchment delineation
LK_walker = 3                                                                   # Number of cells to walk downstream to get minimum lake elevation
z_limit = 1000.0                                                                # Maximum fill depth (z-limit) between a sink and it's pour point
lksatfac_val = 1000.0                                                           # Default LKSATFAC value

# Channel Routing default parameters
Qi = 0                              # Initial Flow in link (cms)
MusK = 3600                         # Muskingum routing time (s)
MusX = 0.2                          # Muskingum weighting coefficient
n = 0.035                           # Manning's roughness
ChSlp = 0.05                        # Channel Side Slope (%; drop/length)
BtmWdth = 5                         # Bottom Width of Channel (m)
Kc = 0                              # channel conductivity (mm/hour)
#type_ = 0                           # Link Type (0 = stream)

#Default Lake Routing parameters
OrificeC = 0.1
OrificA = 1.0
WeirC = 0.4
#WeirL = 0.0                                                                    # Old default weir length (0m)
WeirL = 10.0                                                                    # New default prescribed by D. Yates 5/11/2017 (10m default weir length)

#Custom Geotransformations for spheroid-to-sphere translation
geoTransfmName = "GeoTransform_Null_WRFHydro"                                   # Custom Geotransformation Name
customGeoTransfm = "GEOGTRAN[METHOD['Null']]"                                   # Null Transformation
#customGeoTransfm = "GEOGTRAN[METHOD['Geocentric_Translation'],PARAMETER['X_Axis_Translation',''],PARAMETER['Y_Axis_Translation',''],PARAMETER['Z_Axis_Translation','']]"   # Zero-parameter Geocentric Translation

# --- End Globals --- #

# --- Classes --- #
class ZipCompat(ZipFile):
    def __init__(self, *args, **kwargs):
        ZipFile.__init__(self, *args, **kwargs)

    def extract(self, member, path=None):
        if not isinstance(member, ZipInfo):
            member = self.getinfo(member)
        if path is None:
            path = os.getcwd()
        return self._extract_member(member, path)

    def extractall(self, path=None, members=None, pwd=None):
        if members is None:
            members = self.namelist()
        for zipinfo in members:
            self.extract(zipinfo, path)

    def _extract_member(self, member, targetpath):
        if (targetpath[-1:] in (os.path.sep, os.path.altsep)
            and len(os.path.splitdrive(targetpath)[1]) > 1):
            targetpath = targetpath[:-1]
        if member.filename[0] == '/':
            targetpath = os.path.join(targetpath, member.filename[1:])
        else:
            targetpath = os.path.join(targetpath, member.filename)
        targetpath = os.path.normpath(targetpath)
        upperdirs = os.path.dirname(targetpath)
        if upperdirs and not os.path.exists(upperdirs):
            os.makedirs(upperdirs)
        if member.filename[-1] == '/':
            if not os.path.isdir(targetpath):
                os.mkdir(targetpath)
            return targetpath
        target = file(targetpath, "wb")
        try:
            target.write(self.read(member.filename))
        finally:
            target.close()
        return targetpath
# --- End Classes --- #

# --- Functions --- #
def makeoutncfile(arcpy, raster, outname, outvar, projdir, loglines):
    outfile = os.path.join(projdir, outname)
    arcpy.RasterToNetCDF_md(raster, outfile, outvar, "", "x", "y")
    loglines.append('    Process: %s completed without error' %outname)
    arcpy.AddMessage(loglines[-1])
    loglines.append('         Output File: %s' %outfile)
    arcpy.AddMessage(loglines[-1])
    return loglines

def getxy(inraster, projdir, loglines):
    import arcgisscripting
    gp = arcgisscripting.create()
    gp.OverWriteOutput = 1

    # Check out any necessary licenses
    gp.CheckOutExtension("spatial")
    loglines.append('    Starting Process: Converting raster to XMap/YMap')
    gp.ScratchWorkspace = projdir

    OutLyr = 'gplayer'
    gp.MakeRasterLayer_management(inraster, OutLyr)

    # Set environments
    gp.outputCoordinateSystem = OutLyr
    gp.snapRaster= OutLyr
    gp.extent = OutLyr
    gp.cellSize = OutLyr

    # Perform map algebra
    xmap = os.path.join(projdir, 'xmap')
    ymap = os.path.join(projdir, 'ymap')
    result1 = gp.SingleOutputMapAlgebra("$$XMap", xmap)
    result2 = gp.SingleOutputMapAlgebra("$$YMap", ymap)

    # Clean up
    gp.delete_management(OutLyr)
    del OutLyr

    loglines.append('    Conversion of input raster to XMap/YMap completed without error.')
    return xmap, ymap, loglines

def zipws(arcpy, zipfile, path, zip, keep, nclist):
    path = os.path.normpath(path)
    for dirpath, dirnames, filenames in os.walk(path):
        for file in filenames:
            if file in nclist:
                try:
                    if keep:
                        zip.write(os.path.join(dirpath, file), os.path.join(os.sep + os.path.join(dirpath, file)[len(path) + len(os.sep):]))
                except Exception as e:
                    arcpy.AddMessage(e)
                    arcpy.AddWarning((file, e[0]))

def zipUpFolder(arcpy, folder, outZipFile, nclist):
    try:
        zip = zipfile.ZipFile(outZipFile, 'w', zipfile.ZIP_DEFLATED, allowZip64=True)
        zipws(arcpy, zipfile, str(folder), zip, 'CONTENTS_ONLY', nclist)
        zip.close()
    except RuntimeError:
        pass

def Examine_Outputs(arcpy, in_zip, out_folder, skipfiles=[]):
    """This tool takes the output zip file from the ProcessGeogrid script and creates a raster
     from each output NetCDF file. The Input should be a .zip file that was created using the
     WRF Hydro pre-processing tools. The tool will create the folder which will contain the
     results (out_folder), if that folder does not already exist."""

    # Set environments and create scratch workspace
    arcpy.env.overwriteOutput = True
    out_sfolder = arcpy.CreateScratchName("temp", data_type="Folder", workspace=arcpy.env.scratchFolder)
    os.mkdir(out_sfolder)

    # Unzip to a known location (make sure no other nc files live here)
    ZipCompat(in_zip).extractall(out_sfolder)

    dellist = []                                                                # Keep a directory of files to delete

    # Iterate through unzipped files and copy to output directory as necessary
    for dirpath, dirnames, filenames in os.walk(out_sfolder):
        for file in filenames:
            infile = os.path.join(dirpath, file)
            dellist.append(infile)

            # Copy skipped files over to new directory
            if file in skipfiles:
                shutil.copy2(infile, out_folder)
                arcpy.AddMessage('  File Copied: %s' %file)
                del file

            if file.endswith('.nc'):

                # Trap to eliminate Parameter tables in NC format from this extraction
                if file.endswith(LK_nc) or file.endswith(RT_nc):
                    shutil.copy2(infile, out_folder)
                    arcpy.AddMessage('  File Created: %s' %file)
                    del file
                    continue

                # Establish an object for reading the input NetCDF file
                rootgrp = netCDF4.Dataset(infile, 'r')

                ### Find 2D variables (y,x)
                ##ncVariableNames_2D = [variable for variable, ncvar in rootgrp.variables.iteritems() if ncvar.dimensions==('y', 'x')]
                ##rootgrp.close()
                ##
                ### Loop through global variables in NetCDF file to gather projection information
                ##for variablename in ncVariableNames_2D:
                ##    outRasterLayer = variablename
                ##    arcpy.MakeNetCDFRasterLayer_md(infile, variablename, 'x', 'y', outRasterLayer, "", "", "BY_VALUE")
                ##    RasterObj = arcpy.Raster(outRasterLayer)
                ##    RasterObj.save(os.path.join(out_folder, outRasterLayer))
                ##    #arcpy.Raster(outRasterLayer).save(outRasterLayer + '.tif')
                ##    arcpy.AddMessage('  File Created: %s' %outRasterLayer)
                ##    arcpy.Delete_management(outRasterLayer)
                ##    del RasterObj, variablename, outRasterLayer
                ##del file, infile, rootgrp, ncVariableNames_2D, variable, ncvar

                # Using netCDF4 library - still need point2, DX, DY
                GT = rootgrp.variables['ProjectionCoordinateSystem'].GeoTransform.split(" ")
                PE_string = rootgrp.variables['ProjectionCoordinateSystem'].esri_pe_string
                arcpy.AddMessage('  GeoTransform: %s' %GT)
                DX = float(GT[1])
                DY = abs(float(GT[5]))                                          # In GeoTransform, Y is usually negative
                arcpy.AddMessage('  DX: %s' %DX)
                arcpy.AddMessage('  DY: %s' %DY)
                sr = arcpy.SpatialReference()
                sr.loadFromString(PE_string.replace('"', "'"))
                point = arcpy.Point(float(GT[0]), float(GT[3]) - float(DY*len(rootgrp.dimensions['y'])))    # Calculate LLCorner value from GeoTransform (ULCorner)
                arcpy.env.outputCoordinateSystem = sr
                for variablename, ncvar in rootgrp.variables.iteritems():
                    if ncvar.dimensions==('y', 'x'):
                        outRasterLayer = variablename
                        nc_raster = arcpy.NumPyArrayToRaster(ncvar[:], point, DX, DY)
                        arcpy.CalculateStatistics_management(nc_raster)
                        arcpy.DefineProjection_management(nc_raster, sr)
                        nc_raster.save(os.path.join(out_folder, outRasterLayer))
                        arcpy.AddMessage('  File Created: %s' %outRasterLayer)
                        del nc_raster, variablename, outRasterLayer
                rootgrp.close()
                del file, infile, rootgrp, ncvar
                continue

            if file.endswith('.txt'):
                outRasterLayer = os.path.basename(infile)[:10] + 'txt'
                outRaster = os.path.join(out_folder, outRasterLayer)
                #outRaster = os.path.join(out_folder, outRasterLayer + '.tif')
                arcpy.ASCIIToRaster_conversion(infile, outRaster, 'INTEGER')
                arcpy.AddMessage('  File Created: %s' %file)                      # %outRaster
                del file, infile, outRaster, outRasterLayer
                continue

            if file.endswith('.shp'):
                newshp = os.path.join(out_folder, file)
                arcpy.CopyFeatures_management(infile, newshp)
                arcpy.AddMessage('  File Created: %s' %str(file))
                del file, infile, newshp
                continue

            if file.endswith('.csv') or file.endswith('.TBL'):
                shutil.copy2(infile, out_folder)
                arcpy.AddMessage('  File Created: %s' %file)
                del file, infile
                continue
            else:
                continue
    del dirpath, dirnames, filenames

    # Remove each file from the temporary extraction directory
    for infile in dellist:
        os.remove(infile)
    arcpy.AddMessage('Extraction of WRF routing grids completed.')
    return out_sfolder

def georeference_geogrid_file(arcpy, in_nc, Variable):
    """The first step of the process chain, in which the input NetCDF file gets
    georeferenced and projection files are created"""

    # First step: Import and georeference NetCDF file
    loglines = ['Step 1: NetCDF Conversion initiated... (%s)'%Variable]         # Initiate log list for this process
    arcpy.AddMessage(loglines[-1])

    # Establish an object for reading the input NetCDF file
    ncFP = arcpy.NetCDFFileProperties(in_nc)

    # Loop through global variables in NetCDF file to gather projection information
    ncAttributeNames = ncFP.getAttributeNames("")

    # Find out which projection this GEOGRID file is in
    map_pro = ncFP.getAttributeValue("", 'MAP_PROJ')
    loglines.append('    Map Projection: %s' %projdict[map_pro])
    arcpy.AddMessage(loglines[-1])

    # Collect grid corner XY and DX DY for creating ascii raster later
    if 'corner_lats' in ncAttributeNames:
        corner_lat = ncFP.getAttributeValue("", 'corner_lats')      # Note: The values returned are corner points of the mass grid
    if 'corner_lons' in ncAttributeNames:
        corner_lon = ncFP.getAttributeValue("", 'corner_lons')      # Note: The values returned are corner points of the mass grid
    if 'DX' in ncAttributeNames:
        DX = ncFP.getAttributeValue("", 'DX')
    if 'DY' in ncAttributeNames:
        DY = ncFP.getAttributeValue("", 'DY')

    # Collect necessary information to put together the projection file
    if 'TRUELAT1' in ncAttributeNames:
        standard_parallel_1 = ncFP.getAttributeValue("", 'TRUELAT1')
    if 'TRUELAT2' in ncAttributeNames:
        standard_parallel_2 = ncFP.getAttributeValue("", 'TRUELAT2')
    if 'STAND_LON' in ncAttributeNames:
        central_meridian = ncFP.getAttributeValue("", 'STAND_LON')
    if 'POLE_LAT' in ncAttributeNames:
        pole_latitude = ncFP.getAttributeValue("", 'POLE_LAT')
    if 'POLE_LON' in ncAttributeNames:
        pole_longitude = ncFP.getAttributeValue("", 'POLE_LON')
    if 'MOAD_CEN_LAT' in ncAttributeNames:
        loglines.append('    Using MOAD_CEN_LAT for latitude of origin.')
        arcpy.AddMessage(loglines[-1])
        latitude_of_origin = ncFP.getAttributeValue("", 'MOAD_CEN_LAT')         # Added 2/26/2017 by KMS
    elif 'CEN_LAT' in ncAttributeNames:
        loglines.append('    Using CEN_LAT for latitude of origin.')
        arcpy.AddMessage(loglines[-1])
        latitude_of_origin = ncFP.getAttributeValue("", 'CEN_LAT')

    # Process: Make NetCDF Raster Layer
    Dimensions = ncFP.getDimensionsByVariable(Variable)
    Y_Dimension = Dimensions[1]
    X_Dimension = Dimensions[2]
    NC_Raster_Layer = "NC_RASTER"
    arcpy.MakeNetCDFRasterLayer_md(in_nc, Variable, X_Dimension, Y_Dimension, NC_Raster_Layer, "", "", "BY_VALUE")

    # Process: Raster to Numpy Array
    data = arcpy.RasterToNumPyArray(NC_Raster_Layer)
    arcpy.Delete_management(NC_Raster_Layer)
    del NC_Raster_Layer

    # Create Projection file with information from NetCDF global attributes
    sr2 = arcpy.SpatialReference()
    if map_pro == 1:
        # Lambert Conformal Conic
        if 'standard_parallel_2' in locals():
            projname = 'Lambert_Conformal_Conic_2SP'
            loglines.append('    Using Standard Parallel 2 in Lambert Conformal Conic map projection.')
            arcpy.AddMessage(loglines[-1])
        else:
            # According to http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Lambert_Conformal_Conic
            projname = 'Lambert_Conformal_Conic_1SP'
            standard_parallel_2 = standard_parallel_1
            latitude_of_origin = standard_parallel_1
        Projection_String = ('PROJCS["Lambert_Conformal_Conic",'
                             'GEOGCS["GCS_Sphere",'
                             'DATUM["D_Sphere",'
                             'SPHEROID["Sphere",6370000.0,0.0]],'
                             'PRIMEM["Greenwich",0.0],'
                             'UNIT["Degree",0.0174532925199433]],'
                             'PROJECTION["' + projname + '"],'
                             'PARAMETER["false_easting",0.0],'
                             'PARAMETER["false_northing",0.0],'
                             'PARAMETER["central_meridian",' + str(central_meridian) + '],'
                             'PARAMETER["standard_parallel_1",' + str(standard_parallel_1) + '],'
                             'PARAMETER["standard_parallel_2",' + str(standard_parallel_2) + '],'
                             'PARAMETER["latitude_of_origin",' + str(latitude_of_origin) + '],'
                             'UNIT["Meter",1.0]]')

    elif map_pro == 2:
        # Polar Stereographic

        # Set up pole latitude
        phi1 = float(standard_parallel_1)

        ### Back out the central_scale_factor (minimum scale factor?) using formula below using Snyder 1987 p.157 (USGS Paper 1395)
        ##phi = math.copysign(float(pole_latitude), float(latitude_of_origin))    # Get the sign right for the pole using sign of CEN_LAT (latitude_of_origin)
        ##central_scale_factor = (1 + (math.sin(math.radians(phi1))*math.sin(math.radians(phi))) + (math.cos(math.radians(float(phi1)))*math.cos(math.radians(phi))))/2

        # Method where central scale factor is k0, Derivation from C. Rollins 2011, equation 1: http://earth-info.nga.mil/GandG/coordsys/polar_stereographic/Polar_Stereo_phi1_from_k0_memo.pdf
        # Using Rollins 2011 to perform central scale factor calculations. For a sphere, the equation collapses to be much  more compact (e=0, k90=1)
        central_scale_factor = (1 + math.sin(math.radians(abs(phi1))))/2        # Equation for k0, assumes k90 = 1, e=0. This is a sphere, so no flattening

        # Hardcode in the pole as the latitude of origin (correct assumption?) Added 4/4/2017 by KMS
        standard_parallel_1 = pole_latitude

        loglines.append('      Central Scale Factor: %s' %central_scale_factor)
        arcpy.AddMessage(loglines[-1])
        Projection_String = ('PROJCS["Sphere_Stereographic",'
                             'GEOGCS["GCS_Sphere",'
                             'DATUM["D_Sphere",'
                             'SPHEROID["Sphere",6370000.0,0.0]],'
                             'PRIMEM["Greenwich",0.0],'
                             'UNIT["Degree",0.0174532925199433]],'
                             'PROJECTION["Stereographic"],'
                             'PARAMETER["False_Easting",0.0],'
                             'PARAMETER["False_Northing",0.0],'
                             'PARAMETER["Central_Meridian",' + str(central_meridian) + '],'
                             'PARAMETER["Scale_Factor",' + str(central_scale_factor) + '],'
                             'PARAMETER["Latitude_Of_Origin",' + str(standard_parallel_1) + '],'
                             'UNIT["Meter",1.0]]')

    elif map_pro == 3:
        # Mercator Projection
        Projection_String = ('PROJCS["Sphere_Mercator",'
                             'GEOGCS["GCS_Sphere",'
                             'DATUM["D_Sphere",'
                             'SPHEROID["Sphere",6370000.0,0.0]],'
                             'PRIMEM["Greenwich",0.0],'
                             'UNIT["Degree",0.0174532925199433]],'
                             'PROJECTION["Mercator"],'
                             'PARAMETER["False_Easting",0.0],'
                             'PARAMETER["False_Northing",0.0],'
                             'PARAMETER["Central_Meridian",' + str(central_meridian) + '],'
                             'PARAMETER["Standard_Parallel_1",' + str(standard_parallel_1) + '],'
                             'UNIT["Meter",1.0]]')

    elif map_pro == 6:
        # Cylindrical Equidistant (or Rotated Pole)
        if pole_latitude != float(90) or pole_longitude != float(0):
            # if pole_latitude, pole_longitude, or stand_lon are changed from thier default values, the pole is 'rotated'.
            loglines.append('    Cylindrical Equidistant projection with a rotated pole is not currently supported.')
            arcpy.AddMessage(loglines[-1])
            sys.exit(1)
        else:
            # Check units (linear unit not used in this projection).  GCS?
            Projection_String = ('PROJCS["Sphere_Equidistant_Cylindrical",'
                                 'GEOGCS["GCS_Sphere",'
                                 'DATUM["D_Sphere",'
                                 'SPHEROID["Sphere",6370000.0,0.0]],'
                                 'PRIMEM["Greenwich",0.0],'
                                 'UNIT["Degree",0.0174532925199433]],'
                                 'PROJECTION["Equidistant_Cylindrical"],'
                                 'PARAMETER["False_Easting",0.0],'
                                 'PARAMETER["False_Northing",0.0],'
                                 'PARAMETER["Central_Meridian",' + str(central_meridian) + '],'
                                 'PARAMETER["Standard_Parallel_1",' + str(standard_parallel_1) + '],'
                                 'UNIT["Meter",1.0]]')                          # 'UNIT["Degree", 1.0]]') # ?? For lat-lon grid?

    sr2.loadFromString(Projection_String)

    ##    # Create a custom geotransformation for this conversion (WGS84 to WRF Sphere)
    ##    geoTransfmName = 'GCS_WGS_1984_To_WRF_Sphere'
    ##    #sr_WGS84 = arcpy.SpatialReference(4326)                                    # WGS 1984
    ##    customGeoTransfm = "GEOGTRAN[METHOD['Geocentric_Translation'],PARAMETER['X_Axis_Translation',0.0],PARAMETER['Y_Axis_Translation',0.0],PARAMETER['Z_Axis_Translation',0.0]]"
    ##    arcpy.CreateCustomGeoTransformation_management(geoTransfmName, sr_WGS84, sr2, customGeoTransfm)

    # Create a point geometry object from gathered corner point data
    sr1 = arcpy.SpatialReference(104128)                                        # Using EMEP Sphere (6370000m)
    point = arcpy.Point()
    point.X = corner_lon
    point.Y = corner_lat
    pointGeometry = arcpy.PointGeometry(point, sr1)
    projpoint = pointGeometry.projectAs(sr2)                                    # Optionally add transformation method:

    # Get projected X and Y of the point geometry and adjust so lower left center becomes lower left corner
    point2 = arcpy.Point((projpoint.firstPoint.X - (DX/2)),(projpoint.firstPoint.Y - (DY/2)))   # Adjust by half a grid cell
    x00 = point2.X                                                              # X value doesn't change from LLcorner to UL corner
    y00 = point2.Y + float(DY*data.shape[0])                                    # Adjust Y from LL to UL

    # Process: Numpy Array to Raster
    nc_raster = arcpy.NumPyArrayToRaster(data, point2, DX, DY, NoDataVal)

    # GeoTransform
    GeoTransformStr = '%s %s %s %s %s %s ' %(x00, DX, 0, y00, 0, -DY)
    loglines.append('    GeoTransform: ' + GeoTransformStr)
    arcpy.AddMessage(loglines[-1])
    del data, ncFP

    # Process: Define Projection
    arcpy.DefineProjection_management(nc_raster, sr2)
    loglines.append('    Step 1 completed without error.')
    arcpy.AddMessage(loglines[-1])
    return nc_raster, sr2, Projection_String, map_pro, GeoTransformStr, loglines

def create_CF_NetCDF(arcpy, in_raster, rootgrp, map_pro, projdir, DXDY_dict, sr, GeoTransformStr, addLatLon=False, notes='', loglines=[], addVars=[]):
    """This function will create the netCDF file with CF conventions for the grid
    description. Valid output formats are 'GEOGRID', 'ROUTING_GRID', and 'POINT'.
    The output NetCDF will have the XMAP/YMAP created for the x and y variables
    and the LATITUDE and LONGITUDE variables populated from the XLAT_M and XLONG_M
    variables in the GEOGRID file or in the case of the routing grid, populated
    using the getxy function."""

    tic1 = time.time()
    loglines.append('Creating CF-netCDF File.')
    arcpy.AddMessage(loglines[-1])

    # Local variables
    copydims = []
    copyvars = []
    copyatts = []
    ##    copyatts = ['CEN_LAT','CEN_LON', 'TRUELAT1', 'TRUELAT2', 'MOAD_CEN_LAT', 'STAND_LON',
    ##                'POLE_LAT', 'POLE_LON', 'corner_lats', 'corner_lons', 'MAP_PROJ']       # only keep these global attributes from GEOGRID file

    # For copying dimensions, variables, attributes from input GEOGRID file to output spatial metadata file
    ##    # Use only necessary GEOGRID dimensions, variables, and attributes
    ##    if format_out == 'LDASOUT':
    ##
    ##        # Copy certain dimensions from the GEOGRID
    ##        copydims = ['Time', 'west_east', 'south_north']
    ##        copyvars = ['XLAT_M', 'XLONG_M']

    ##    # Copy dimensions if necessary
    ##    for dname, dimension in rootgrp2.dimensions.iteritems():
    ##        if dname in copydims:
    ##            rootgrp.createDimension(dname, len(dimension) if not dimension.isunlimited() else None)
    ##
    ##    # Copy certain spatial variables and any coordinate variables
    ##    for v_name, variable in rootgrp2.variables.iteritems():
    ##        if v_name in copyvars or v_name in copydims:
    ##            var_vals = variable[:]                                              # Read variable values
    ##            varatts = {k: variable.getncattr(k) for k in variable.ncattrs()}    # Read variable attributes
    ##            var_dtype = variable.datatype                                       # Read variable data type
    ##            var_dims = variable.dimensions                                      # Read variable dimensions
    ##            outVar = rootgrp.createVariable(v_name, var_dtype, var_dims)        # Create variable
    ##            outVar.setncatts(varatts)                                           # Set variable attributes
    ##            outVar[:] = var_vals                                                # Set variable values
    ##
    ##    # Copy global attributes (remove strings as this may be causing a problem)
    ##    inglobalatts = rootgrp2.__dict__                                            # Global attributes from GEOGRID file
    ##    outglobalatts = {key:value for key,value in inglobalatts.iteritems() if key in copyatts}
    ##    outglobalatts.update(DXDY_dict)                                             # Add in the new DX/DY values
    ##    rootgrp.setncatts(outglobalatts)
    ##    rootgrp2.close()                                                            # Close the GEOGRID file

    # Gather projection information from input raster projection
    descData = arcpy.Describe(in_raster)
    dim1size = descData.width
    dim2size = descData.height
    sr = descData.SpatialReference
    srType = sr.type
    PE_string = sr.exportToString()
    PE_string = PE_string.replace("'", '"')                                     # Replace ' with " so Esri can read the PE String properly when running NetCDFtoRaster
    loglines.append('    Esri PE String: %s' %PE_string)
    arcpy.AddMessage(loglines[-1])

    # Find name for the grid mapping
    if CF_projdict.get(map_pro) is not None:
        grid_mapping = CF_projdict[map_pro]
    else:
        grid_mapping = sr.name
        loglines.append('    Map Projection of input raster (not a WRF projection): %s' %grid_mapping)
        arcpy.AddMessage(loglines[-1])

    # Must handle difference between ProjectionCoordinateSystem and LatLonCoordinateSystem
    if srType == 'Geographic':
        CoordSysVarName = "LatLonCoordinateSystem"
    elif srType == 'Projected':
        CoordSysVarName = "ProjectionCoordinateSystem"
        #proj_units = sr.linearUnitName.lower()                                  # sr.projectionName wouldn't work for a GEOGCS
        proj_units = 'm'                                                        # Change made 11/3/2016 by request of NWC

    # Create Dimensions
    dim_y = rootgrp.createDimension('y', dim2size)
    dim_x = rootgrp.createDimension('x', dim1size)
    print '    Dimensions created after {0: 8.2f} seconds.'.format(time.time()-tic1)

    # Create coordinate variables
    var_y = rootgrp.createVariable('y', 'f8', 'y')                              # (64-bit floating point)
    var_x = rootgrp.createVariable('x', 'f8', 'x')                              # (64-bit floating point)
    if srType == 'Geographic':
        #var_y.standard_name = ''
        #var_x.standard_name = ''
        var_y.long_name = "latitude coordinate"
        var_x.long_name = "longitude coordinate"
        var_y.units = "degrees_north"
        var_x.units = "degrees_east"
        var_y._CoordinateAxisType = "Lat"
        var_x._CoordinateAxisType = "Lon"
    elif srType == 'Projected':
        # Set variable attributes
        var_y.standard_name = 'projection_y_coordinate'
        var_x.standard_name = 'projection_x_coordinate'
        var_y.long_name = 'y coordinate of projection'
        var_x.long_name = 'x coordinate of projection'
        var_y.units = proj_units                                                # was 'meter', now 'm'
        var_x.units = proj_units                                                # was 'meter', now 'm'
        var_y._CoordinateAxisType = "GeoY"                                      # Use GeoX and GeoY for projected coordinate systems only
        var_x._CoordinateAxisType = "GeoX"                                      # Use GeoX and GeoY for projected coordinate systems only
        var_y.resolution = float(DXDY_dict['DY'])                               # Added 11/3/2016 by request of NWC
        var_x.resolution = float(DXDY_dict['DX'])                               # Added 11/3/2016 by request of NWC

        # Scalar projection variable - http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/StandardCoordinateTransforms.html
        proj_var = rootgrp.createVariable(CoordSysVarName, 'S1')                # (Scalar Char variable)
        proj_var._CoordinateAxes = 'y x'                                        # Coordinate systems variables always have a _CoordinateAxes attribute
        proj_var._CoordinateTransformType = "Projection"
        proj_var.transform_name = grid_mapping                                  # grid_mapping
        proj_var.grid_mapping_name = grid_mapping                               # for CF compatibility
        proj_var._CoordinateAxes = 'y x'                                        # Optional for dealing with implicit coordinate systems
        proj_var.esri_pe_string = PE_string                                     # For ArcGIS
        proj_var.spatial_ref = PE_string                                        # For GDAl
        proj_var.GeoTransform = GeoTransformStr                                 # For GDAl - GeoTransform array

        # Projection specific parameters - http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/StandardCoordinateTransforms.html
        if map_pro == 1:
            # Lambert Conformal Conic

            # Required transform variables
            proj_var.standard_parallel = sr.standardParallel1, sr.standardParallel2     # Double
            proj_var.longitude_of_central_meridian = float(sr.centralMeridianInDegrees) # Double
            proj_var.latitude_of_projection_origin = float(sr.latitudeOfOrigin)         # Double

            # Optional tansform variable attributes
            proj_var.false_easting = float(sr.falseEasting)                             # Double  Always in the units of the x and y projection coordinates
            proj_var.false_northing = float(sr.falseNorthing)                           # Double  Always in the units of the x and y projection coordinates
            proj_var.earth_radius = 6370000.0                                           # OPTIONAL. Parameter not read by Esri. Default CF sphere: 6371.229 km
            #proj_var.semi_major_axis =  6370000.0                                      # Double - optional Lambert Conformal Conic parameter
            #proj_var.semi_minor_axis =  6370000.0                                      # Double - optional Lambert Conformal Conic parameter
            #proj_var.inverse_flattening =   0.0                                        # Double - optional Lambert Conformal Conic parameter

        elif map_pro == 2:
            # Polar Stereographic

            # Required transform variables
            proj_var.longitude_of_projection_origin = float(sr.longitudeOfOrigin)       # Double - proj_var.straight_vertical_longitude_from_pole = ''
            proj_var.latitude_of_projection_origin = float(sr.latitudeOfOrigin)         # Double
            proj_var.scale_factor_at_projection_origin = float(sr.scaleFactor)          # Double

            # Optional tansform variable attributes
            proj_var.false_easting = float(sr.falseEasting)                             # Double  Always in the units of the x and y projection coordinates
            proj_var.false_northing = float(sr.falseNorthing)                           # Double  Always in the units of the x and y projection coordinates
            #proj_var.semi_major_axis =  6370000.0                                      # Double - optional Lambert Conformal Conic parameter
            #proj_var.semi_minor_axis =  6370000.0                                      # Double - optional Lambert Conformal Conic parameter
            #proj_var.inverse_flattening =   0.0                                        # Double - optional Lambert Conformal Conic parameter

        elif map_pro == 3:
            # Mercator

            # Required transform variables
            proj_var.longitude_of_projection_origin = float(sr.longitudeOfOrigin)       # Double
            proj_var.latitude_of_projection_origin = float(sr.latitudeOfOrigin)         # Double
            proj_var.standard_parallel = float(sr.standardParallel1)                    # Double

        elif map_pro == 6:
            # Cylindrical Equidistant or rotated pole

            #http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#appendix-grid-mappings
            # Required transform variables
            #proj_var.grid_mapping_name = "latitude_longitude"                  # or "rotated_latitude_longitude"

            #loglines.append('        Cylindrical Equidistant projection not supported.')
            #arcpy.AddMessage(loglines[-1])
            #raise SystemExit
            pass                                                                # No extra parameters needed for latitude_longitude

    # For prefilling additional variables and attributes on the same 2D grid, given as a list [[<varname>, <vardtype>, <long_name>],]
    for varinfo in addVars:
        ncvar = rootgrp.createVariable(varinfo[0], varinfo[1], ('y', 'x'))
        ncvar.esri_pe_string = PE_string
        ncvar.grid_mapping = CoordSysVarName
        #ncvar.spatial_ref = PE_string                                          # For GDAl
        #ncvar.GeoTransform = GeoTransformStr                                   # For GDAl - GeoTransform array
        #ncvar.long_name = varinfo[2]
        #ncvar.units = varinfo[3]

    # Get x and y variables for the netCDF file
    xmap, ymap, loglines2 = getxy(in_raster, projdir, [])
    loglines += loglines2
    ymaparr = arcpy.RasterToNumPyArray(ymap)
    xmaparr = arcpy.RasterToNumPyArray(xmap)

    # Assumes even spacing in y across domain
    var_y[:] = ymaparr[:,0]
    var_x[:] = xmaparr[0,:]
    arcpy.Delete_management(xmap)
    arcpy.Delete_management(ymap)
    del xmap, ymap, loglines2

    loglines.append('    Coordinate variables and variable attributes set after {0: 8.2f} seconds.'.format(time.time()-tic1))
    arcpy.AddMessage(loglines[-1])

    if addLatLon == True:

        loglines.append('    Proceeding to add LATITUDE and LONGITUDE variables after {0: 8.2f} seconds.'.format(time.time()-tic1))
        arcpy.AddMessage(loglines[-1])

        # Populate this file with 2D latitude and longitude variables
        # Latitude and Longitude variables (WRF)
        lat_WRF = rootgrp.createVariable('LATITUDE', 'f4', ('y', 'x'))          # (32-bit floating point)
        lon_WRF = rootgrp.createVariable('LONGITUDE', 'f4', ('y', 'x'))         # (32-bit floating point)
        lat_WRF.long_name = 'latitude coordinate'                               # 'LATITUDE on the WRF Sphere'
        lon_WRF.long_name = 'longitude coordinate'                              # 'LONGITUDE on the WRF Sphere'
        lat_WRF.units = "degrees_north"
        lon_WRF.units = "degrees_east"
        lat_WRF._CoordinateAxisType = "Lat"
        lon_WRF._CoordinateAxisType = "Lon"
        lat_WRF.grid_mapping = CoordSysVarName                                  # This attribute appears to be important to Esri
        lon_WRF.grid_mapping = CoordSysVarName                                  # This attribute appears to be important to Esri
        lat_WRF.esri_pe_string = PE_string
        lon_WRF.esri_pe_string = PE_string
        #lat_WRF.spatial_ref = PE_string                                         # For GDAl
        #lon_WRF.spatial_ref = PE_string                                         # For GDAl
        #lat_WRF.GeoTransform = GeoTransformStr                                  # For GDAl - GeoTransform array
        #lon_WRF.GeoTransform = GeoTransformStr                                  # For GDAl - GeoTransform array

        # Missing value attribute not needed yet
        #missing_val = numpy.finfo(numpy.float32).min                            # Define missing data variable based on numpy
        #lat_WRF.missing_value = missing_val                                     # Float sys.float_info.min?
        #lon_WRF.missing_value = missing_val                                     # Float sys.float_info.min?

        '''Adding the Esri PE String in addition to the CF grid mapping attributes
        is very useful. Esri will prefer the PE string over other CF attributes,
        allowing a spherical datum to be defined. Esri can interpret the coordinate
        system variable alone, but will assume the datum is WGS84. This cannot be
        changed except when using an Esri PE String.'''

        ##    # Create a new coordinate system variable
        ##    LatLonCoordSysVarName = "LatLonCoordinateSystem"
        ##    latlon_var = rootgrp.createVariable(LatLonCoordSysVarName, 'S1')            # (Scalar Char variable)
        ##    latlon_var._CoordinateAxes = 'LATITUDE LONGITUDE'                           # Coordinate systems variables always have a _CoordinateAxes attribute

        # Data variables need _CoodinateSystems attribute
        lat_WRF._CoordinateAxisType = "Lat"
        lon_WRF._CoordinateAxisType = "Lon"
        lat_WRF._CoordinateSystems = CoordSysVarName
        lon_WRF._CoordinateSystems = CoordSysVarName
        ##    lat_WRF._CoordinateSystems = "%s %s" %(CoordSysVarName, LatLonCoordSysVarName)        # For specifying more than one coordinate system
        ##    lon_WRF._CoordinateSystems = "%s %s" %(CoordSysVarName, LatLonCoordSysVarName)        # For specifying more than one coordinate system

        # Create latitude and longitude rasters
        try:
            # Try to trap any errors in this try statement

            # Get lat and lon grids on WRF Sphere
            wkid = 104128                                                       # Using EMEP Sphere (6370000m)
            loglines2, xout, yout, xmap, ymap = create_lat_lon_rasters(arcpy, projdir, in_raster, wkid)
            loglines += loglines2
            arcpy.Delete_management(xmap)
            arcpy.Delete_management(ymap)
            del xmap, ymap, loglines2

            # Populate netCDF variables using converted numpy arrays
            youtarr = arcpy.RasterToNumPyArray(yout)
            xoutarr = arcpy.RasterToNumPyArray(xout)
            loglines.append('        youtarr.shape : %s, %s' %(youtarr.shape[0], youtarr.shape[1]))
            arcpy.AddMessage(loglines[-1])
            loglines.append('        xoutarr.shape : %s, %s' %(xoutarr.shape[0], xoutarr.shape[1]))
            arcpy.AddMessage(loglines[-1])
            loglines.append('        lat_WRF.shape : %s, %s' %(lat_WRF.shape[0], lat_WRF.shape[1]))
            arcpy.AddMessage(loglines[-1])
            loglines.append('        lon_WRF.shape : %s, %s' %(lon_WRF.shape[0], lon_WRF.shape[1]))
            arcpy.AddMessage(loglines[-1])
            lat_WRF[:] = youtarr
            lon_WRF[:] = xoutarr
            del xout, yout, youtarr, xoutarr, in_raster

            loglines.append('    Variables populated after {0: 8.2f} seconds.'.format(time.time()-tic1))
            arcpy.AddMessage(loglines[-1])
            loglines.append('    Process completed without error.')
            arcpy.AddMessage(loglines[-1])

            loglines.append('    LATITUDE and LONGITUDE variables and variable attributes set after {0: 8.2f} seconds.'.format(time.time()-tic1))
            arcpy.AddMessage(loglines[-1])

        except Exception as e:
            loglines.append('    Process did not complete. Error: %s' %e)
            arcpy.AddMessage(loglines[-1])

    # Global attributes
    rootgrp.GDAL_DataType = 'Generic'
    rootgrp.Conventions = 'CF-1.5'                                              # Maybe 1.0 is enough?
    rootgrp.Source_Software = 'WRF-Hydro GIS Pre-processor'
    rootgrp.history = 'Created %s' %time.ctime()
    rootgrp.processing_notes = notes
    loglines.append('    netCDF global attributes set after {0: 8.2f} seconds.'.format(time.time()-tic1))
    arcpy.AddMessage(loglines[-1])

    # Return the netCDF file to the calling script
    return rootgrp, grid_mapping, loglines

def domain_shapefile(arcpy, in_raster, out_shp, sr2):
    """This process creates a shapefile that bounds the GEOGRID file domain. This
    requires the Spatial Analyst extension."""

    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
        from arcpy.sa import *

    arcpy.AddMessage('Step 2: Build constant raster and convert to polygon...')

    # Set environments
    arcpy.env.overwriteOutput = True
    arcpy.env.outputCoordinateSystem = sr2
    descData = arcpy.Describe(in_raster)
    extent = descData.Extent                                                    # Overkill?
    arcpy.env.snapRaster = in_raster
    cellsize = descData.children[0].meanCellHeight
    arcpy.env.cellSize = in_raster

    # Build constant raster
    outConstRaster = CreateConstantRaster(1, "INTEGER", cellsize, extent)

    # Raster to Polygon conversion
    arcpy.RasterToPolygon_conversion(outConstRaster, out_shp, "NO_SIMPLIFY")    #, "VALUE")

    # Clean up
    arcpy.AddMessage('    Finished building GEOGRID domain boundary shapefile: %s.' %out_shp)

def create_high_res_topogaphy(arcpy, in_raster, hgt_m_raster, cellsize, sr2, projdir):
    """The second step creates a high resolution topography raster using a hydrologically-
    corrected elevation dataset (currently either HydroSHEDS or NHDPlusv2)."""

    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
        from arcpy.sa import *

    # Second part of the process
    loglines = ['Step 2 initiated...']                                          # Initiate log list for this process
    arcpy.AddMessage(loglines[-1])

    #Get the extent information from raster object
    arcpy.MakeRasterLayer_management(hgt_m_raster, 'hgt_m_Layer')
    descData = arcpy.Describe('hgt_m_Layer')
    extent = descData.Extent
    arcpy.env.snapRaster = 'hgt_m_Layer'                                            # Does this work or does it need to be hgt_m_raster?

    # Test to make sure hgt_m is an integer multiple of supplied output resolution
    cellsize1 = descData.children[0].meanCellHeight
    loglines.append('    The GEOGRID File resolution is %sm' %str(cellsize1))
    arcpy.AddMessage(loglines[-1])
    cellsize2 = (cellsize1/cellsize)
    loglines.append('    The High-resolution dataset will be %sm' %str(cellsize2))
    arcpy.AddMessage(loglines[-1])

    # List of coordinates from extent and create a polygon geometry object using an array object
    boundaryPolygon = extent.polygon                                            # Added 2016-02-11 to simplify code
    extent1 = extent                                                            # Added 2016-02-11 to simplify code

    # Now project the polygon object to the raster catalog spatial reference
    sr3 = arcpy.Describe(in_raster).spatialReference
    arcpy.CreateCustomGeoTransformation_management(geoTransfmName, sr2, sr3, customGeoTransfm)
    loglines.append('    Tranformation: %s' %geoTransfmName)
    arcpy.AddMessage(loglines[-1])
    projpoly = boundaryPolygon.projectAs(sr3, geoTransfmName)                     # Should be: u'NAD_1983_To_WGS_1984_1'

    # Create raster layer from input raster or mosaic dataset
    MosaicLayer = "MosaicLayer"
    arcpy.MakeRasterLayer_management(in_raster, MosaicLayer, "#", projpoly.extent)
    loglines.append('    MakeRasterLayer process completed without error.')
    arcpy.AddMessage(loglines[-1])
    loglines.append('    The coarse grid has %s rows and %s columns.' %(arcpy.Describe(hgt_m_raster).height, arcpy.Describe(hgt_m_raster).width))
    arcpy.AddMessage(loglines[-1])
    loglines.append('    The input elevation grid (before reprojection) has %s rows and %s columns.' %(arcpy.Describe(MosaicLayer).height, arcpy.Describe(MosaicLayer).width))
    arcpy.AddMessage(loglines[-1])

    # Set environments to force creation of high-res raster to have exact extent and cellsize needed
    arcpy.env.extent = extent1                                                      # using extent directly doesn't work.
    arcpy.env.outputCoordinateSystem = sr2
    arcpy.env.cellSize = cellsize2
    arcpy.env.snapRaster = hgt_m_raster                                             # Redundant?

    # Now project the polygon object to the raster catalog spatial reference
    mosprj = os.path.join(projdir, 'mosaicprj')
    descData = arcpy.Describe('hgt_m_Layer')
    extent = descData.Extent
    loglines.append('    Projecting input elevation data to WRF coordinate system.')
    arcpy.AddMessage(loglines[-1])
    arcpy.ProjectRaster_management(MosaicLayer, mosprj, sr2, "NEAREST", cellsize2, geoTransfmName)
    loglines.append('      Finished projecting input elevation data to WRF coordinate system.')
    arcpy.AddMessage(loglines[-1])
    loglines.append('    The fine grid (before ExtractByMask) has %s rows and %s columns.' %(arcpy.Describe(mosprj).height, arcpy.Describe(mosprj).width))
    arcpy.AddMessage(loglines[-1])

    # Extract By Mask
    arcpy.env.cellSize = cellsize2
    mosprj2 = ExtractByMask(mosprj, hgt_m_raster)                               # Thin the raster down from the projected raster.
    arcpy.Delete_management(mosprj)
    mosprj2.save(os.path.join(projdir, 'mosaicprj'))

    # Check that the number of rows and columns are correct
    loglines.append('    Fine Grid has %s rows and %s columns.' %(arcpy.Describe(mosprj2).height, arcpy.Describe(mosprj2).width))
    arcpy.AddMessage(loglines[-1])

    # Clean up
    arcpy.Delete_management("MosaicLayer")
    del MosaicLayer

    # Finish
    loglines.append('    Step 2 completed without error.')
    arcpy.AddMessage(loglines[-1])
    return mosprj, cellsize1, cellsize2, loglines

def create_lat_lon_rasters(arcpy, projdir, mosprj, wkid=104128):
    """The third function in the process is to create the latitude and longitude
    rasters that are necessary for running wrf-hydro.  The latitude and longitude
    that are given by the cell values in the resulting output grids are the latitude
    and longitude of the cell center."""

    # Set buffer distance so that data is still present when converting between Sphere and WGS84
    # Added 2016-02-11 for conversion to WGS84 instead of sphere
    #buffdist = 0.2        # Distance in degrees that ~22000m (maximum sphere-spheroid error) represents (really 0.199)

    # Initiate logging
    loglines = ['    Latitude and Longitude grid generation initiated...']
    arcpy.AddMessage(loglines[-1])

    # Create Constant Raster in-memory so we don't have to send the hgt_m layer
    arcpy.env.outputCoordinateSystem = mosprj
    arcpy.env.snapRaster = mosprj
    arcpy.env.extent = mosprj
    arcpy.env.cellSize = mosprj
    OutRas = CreateConstantRaster(1)                                            # Create constant raster with value of 1

    # Set coordinate reference systems
    sr1 = arcpy.SpatialReference(wkid)                                          # Use provided wkid for CRS or default = EMEP Sphere (6370000m)
    sr2 = arcpy.Describe(mosprj).spatialReference

    # Project constant raster to geocentric coordinates system
    projraster = os.path.join(projdir, 'mosprj2')
    arcpy.ResetEnvironments()
    arcpy.env.outputCoordinateSystem = sr1                                      #arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(104128)           # EMEP sphere (same as WRF sphere)
    arcpy.CopyRaster_management(OutRas, projraster)

    ##    # Expand raster to accomodate datum shift of up to 23km (north-south direction only) Added 2016-02-11
    ##    extent1 = arcpy.Describe(projraster).extent
    ##    newextent = arcpy.Extent(extent1.XMin, extent1.YMin-buffdist, extent1.XMax, extent1.YMax+buffdist)  #Only buffer in the Y direction
    ##    arcpy.ResetEnvironments()
    ##    arcpy.env.outputCoordinateSystem = sr1
    ##    arcpy.env.snapRaster = projraster
    ##    arcpy.env.extent = newextent
    ##    arcpy.env.cellSize = projraster
    ##    OutRas2 = CreateConstantRaster(1)                                           # Create constant raster with value of 1
    ##    projraster2 = os.path.join(projdir, 'mosprj3')
    ##    arcpy.CopyRaster_management(OutRas2, projraster2)
    ##    #del extent1, OutRas2

    # Create xmap/ymap grids
    #xmap, ymap, loglines = getxy(projraster2, projdir, loglines)
    xmap, ymap, loglines = getxy(projraster, projdir, loglines)
    loglines.append('        XMAP and YMAP rasters created.')
    arcpy.AddMessage(loglines[-1])

    # Set environments
    arcpy.ResetEnvironments()                                                   # Added 2016-02-11
    arcpy.env.outputCoordinateSystem = mosprj
    arcpy.env.snapRaster = mosprj
    arcpy.env.extent = mosprj
    arcpy.env.cellSize = mosprj

    # Project Raster back to original projection
    CellSize = arcpy.Describe(mosprj).MeanCellHeight
    xout = os.path.join(projdir, 'xoutput')
    yout = os.path.join(projdir, 'youtput')

    # Project the input CRS to the output CRS
    # Change to "BILINEAR" or "CUBIC" interpolation method to reduce the issue James found of redundancy in some coordinate values.
    arcpy.ProjectRaster_management(xmap, xout, sr2, "BILINEAR", CellSize, "#", "#", sr1)        # "NEAREST"
    arcpy.ProjectRaster_management(ymap, yout, sr2, "BILINEAR", CellSize, "#", "#", sr1)        # "NEAREST"

    # Test - extract by  mask
    xout2 = ExtractByMask(xout, OutRas)     # Added 8/2/2014
    yout2 = ExtractByMask(yout, OutRas)     # Added 8/2/2014
    arcpy.Delete_management(xout)
    arcpy.Delete_management(yout)
    arcpy.Delete_management(projraster)

    loglines.append('        Latitude and Longitude NetCDF Files created.')
    arcpy.AddMessage(loglines[-1])
    loglines.append('        Latitude and Longitude grid generation completed without error.')
    arcpy.AddMessage(loglines[-1])
    return loglines, xout2, yout2, xmap, ymap

def Routing_Table(arcpy, projdir, sr2, channelgrid, fdir, Elev, Strahler, loglines):
    """If "Create reach-based routing files?" is selected, this function will create
    the Route_Link.nc table and Streams.shp shapefiles in the output directory."""

    # Stackless topological sort algorithm, adapted from: http://stackoverflow.com/questions/15038876/topological-sort-python
    def sort_topologically_stackless(graph):

        '''This function will navigate through the list of segments until all are accounted
        for. The result is a sorted list of which stream segments should be listed
        first. Simply provide a topology dictionary {Fromnode:[ToNode,...]} and a sorted list
        is produced that will provide the order for navigating downstream. This version
        is "stackless", meaning it will not hit the recursion limit of 1000.'''

        levels_by_name = {}
        names_by_level = defaultdict(set)

        def add_level_to_name(name, level):
            levels_by_name[name] = level
            names_by_level[level].add(name)

        def walk_depth_first(name):
            stack = [name]
            while(stack):
                name = stack.pop()
                if name in levels_by_name:
                    continue

                if name not in graph or not graph[name]:
                    level = 0
                    add_level_to_name(name, level)
                    continue

                children = graph[name]

                children_not_calculated = [child for child in children if child not in levels_by_name]
                if children_not_calculated:
                    stack.append(name)
                    stack.extend(children_not_calculated)
                    continue

                level = 1 + max(levels_by_name[lname] for lname in children)
                add_level_to_name(name, level)

        for name in graph:
            walk_depth_first(name)

        list1 = list(takewhile(lambda x: x is not None, (names_by_level.get(i, None) for i in count())))
        list2 = [item for sublist in list1 for item in sublist][::-1]               # Added by KMS 9/2/2015 to reverse sort the list
        list3 = [x for x in list2 if x is not None]                                 # Remove None values from list
        return list3

    loglines.append('    Routing table will be created...')
    arcpy.AddMessage(loglines[-1])

    # Get grid information from channelgrid
    descdata = arcpy.Describe(channelgrid)
    extent = descdata.extent
    cellsizeY = descdata.meanCellHeight
    cellsizeX = descdata.meanCellWidth
    del descdata

    # Set output coordinate system environment
    sr1 = arcpy.SpatialReference(104128)                                        # Using EMEP/WRF Sphere (6370000m)
    arcpy.env.outputCoordinateSystem = sr2

    # Output files
    outStreams = os.path.join(projdir, 'Streams.shp')
    RoutingCSV = os.path.join(projdir, 'Route_Link.csv')
    OutFC = os.path.join("in_memory", "Nodes")
    OutFC2 = os.path.join("in_memory", "NodeElev")
    OutFC3 = os.path.join("in_memory", "NodeOrder")
    outRaster = "in_memory/LINKID"

    # Build Stream Features shapefile
    StreamToFeature(channelgrid, fdir, outStreams, "NO_SIMPLIFY")
    loglines.append('        Stream to features step complete.')
    arcpy.AddMessage(loglines[-1])

    # Set environments based on channelgrid
    arcpy.env.snapRaster = channelgrid
    arcpy.env.extent =  channelgrid

    # Create a raster based on the feature IDs that were created in StreamToFeature
    arcpy.FeatureToRaster_conversion(outStreams, 'ARCID', outRaster, channelgrid)                   # Must do this to get "ARCID" field into the raster
    maxValue = arcpy.SearchCursor(outStreams, "", "", "", 'ARCID' + " D").next().getValue('ARCID')  # Gather highest "ARCID" value from field of segment IDs
    maxRasterValue = arcpy.GetRasterProperties_management(outRaster, "MAXIMUM")                     # Gather maximum "ARCID" value from raster
    if int(maxRasterValue[0]) > maxValue:
        loglines.append('        Setting linkid values of %s to Null.' %maxRasterValue[0])
        arcpy.AddMessage(loglines[-1])
        whereClause = "VALUE = %s" %int(maxRasterValue[0])
        outRaster = SetNull(outRaster, outRaster, whereClause)                  # This should eliminate the discrepency between the numbers of features in outStreams and outRaster
    outRaster = Con(IsNull(outRaster) == 1, NoDataVal, outRaster)               # Set Null values to -9999 in LINKID raster

    # Create new Feature Class and begin populating it
    arcpy.CreateFeatureclass_management("in_memory", "Nodes", "POINT")
    arcpy.AddField_management(OutFC, "NODE", "LONG")

    # Initiate dictionaries for storing topology information
    From_To = {}                                                                # From_Node/To_Node information
    Nodes = {}                                                                  # Node firstpoint/lastpoint XY information
    NodesLL = {}                                                                # Stores the projected node Lat/Lon information in EMEP Sphere GCS
    Lengths = {}                                                                # Gather the stream feature length
    StrOrder = {}                                                               # Store stream order for each node

    # Enter for loop for each feature/row to gather length and endpoints of stream segments                                                       # Create an array and point object needed to create features
    point = arcpy.Point()
    with arcpy.da.SearchCursor(outStreams, ['SHAPE@', 'ARCID', 'FROM_NODE', 'TO_NODE']) as rows:                       # Start SearchCursor to look through all linesegments
        for row in rows:
            ID = row[1]                                                         # Get Basin ARCID
            From_To[ID] = [row[2], row[3]]                                      # Store From/To node for each segment
            feat = row[0]                                                       # Create the geometry object 'feat'

            if feat.isMultipart == False:                                       # Make sure that each part is a single part feature
                firstpoint = feat.firstPoint                                    # First point feature geometry
                lastpoint = feat.lastPoint                                      # Last point feature geometry
                Lengths[ID] = feat.length                                       # Store length of the stream segment

                # Gather the X and Y of the top and bottom ends
                for i in firstpoint,lastpoint:                                  # Now put this geometry into new feature class
                    point.X = i.X
                    point.Y = i.Y
                    pointGeometry = arcpy.PointGeometry(point, sr2)
                    projpoint = pointGeometry.projectAs(sr1)                    # Convert to latitude/longitude on the sphere
                    projpoint1 = projpoint.firstPoint
                    if i == firstpoint:                                         # Top Point
                        if row[2] in Nodes:
                            continue                                            # Skip entry if record already exists in dictionary
                        Nodes[row[2]] = (i.X, i.Y)
                        NodesLL[row[2]] = (projpoint1.X, projpoint1.Y)
                    elif i == lastpoint:                                        # Bottom Point
                        if row[3] in Nodes:
                            continue                                            # Skip entry if record already exists in dictionary
                        Nodes[row[3]] = (i.X, i.Y)
                        NodesLL[row[3]] = (projpoint1.X, projpoint1.Y)
            else:
                loglines.append('This is a multipart line feature and cannot be handled.')
                arcpy.AddMessage(loglines[-1])
    #rows.reset()                                                                # Not sure why this is necessary
    del row, rows, point, ID, feat, firstpoint, lastpoint, pointGeometry, projpoint, projpoint1, i, sr1, sr2
    loglines.append('        Done reading streams layer.')
    arcpy.AddMessage(loglines[-1])

    # Make a point feature class out of the nodes
    IC = arcpy.da.InsertCursor(OutFC, ['SHAPE@', 'NODE'])
    NodesXY = Nodes                                                             # Copy the Nodes dictionary before it gets modified
    for node in Nodes.keys():                                                   # Now we have to adjust the points that fall outside of the raster edge

        # Adjust X
        if Nodes[node][0] <= extent.XMin:
            Nodes[node] = ((Nodes[node][0] + (cellsizeX/2)), Nodes[node][1])
        elif Nodes[node][0] >= extent.XMax:
            Nodes[node] = ((Nodes[node][0] - (cellsizeX/2)), Nodes[node][1])

        # Adjust Y
        if Nodes[node][1] <= extent.YMin:
            Nodes[node] = (Nodes[node][0], (Nodes[node][1] + (cellsizeY/2)))
        elif Nodes[node][1] >= extent.YMax:
            Nodes[node] = (Nodes[node][0], (Nodes[node][1] - (cellsizeY/2)))

        IC.insertRow([Nodes[node], node])                                       # Insert point and ID information into the point feature class

    del IC, node, Nodes, extent, cellsizeY, cellsizeX
    loglines.append('        Done building Nodes layer with adjustments.')
    arcpy.AddMessage(loglines[-1])

    # Get the elevation values for the nodes feature class
    ExtractValuesToPoints(OutFC, Elev, OutFC2, "NONE", "VALUE_ONLY")
    loglines.append('        Done extracting elevations to points.')
    arcpy.AddMessage(loglines[-1])
    NodeElev = {row[0]: row[1] for row in arcpy.da.SearchCursor(OutFC2, ['NODE', 'RASTERVALU'])}
    loglines.append('        Done reading node elevations.')
    arcpy.AddMessage(loglines[-1])
    arcpy.Delete_management(OutFC2)                                             # Clean up

    # Incorporate Strahler Order
    ExtractValuesToPoints(OutFC, Strahler, OutFC3, "NONE", "VALUE_ONLY")
    with arcpy.da.SearchCursor(OutFC3, ['NODE', 'RASTERVALU']) as rows:
        for row in rows:
            if row[1] <= 0:                                                     # Reclass -9999 values to 1
                order = 1
            else:
                order = row[1]
            StrOrder[row[0]] = order
    loglines.append('        Done reading Strahler stream orders.')
    arcpy.AddMessage(loglines[-1])

    arcpy.Delete_management(OutFC)                                              # Clean up
    arcpy.Delete_management(OutFC3)                                             # Clean up

    # Add stream order into the streams shapefile
    arcpy.AddField_management(outStreams, "Order_", "SHORT")                    # Add field for "Order_"
    arcpy.AddField_management(outStreams, "index", "LONG")                      # Add field for "index" that gives the topologically sorted order of streams in the out nc file
    with arcpy.da.UpdateCursor(outStreams, ['ARCID', 'Order_']) as rows:        # Start UpdateCursor to add the stream order information
        for row in rows:
            row[1] = StrOrder[From_To[row[0]][0]]                               # This gets the stream order of the upstream node
            rows.updateRow(row)

    # Deconstruct from Node space to segment space
    Arc_From = {x:From_To[x][0] for x in From_To}                                       # Build ARCID-keyed topology of The ARCID:FromNode
    From_Arc = {From_To[x][0]:x for x in From_To}                                       # Build Node-keyed topology of The FromNode:ARCID
    From_To2 = {From_To[x][0]:From_To[x][1] for x in From_To}                           # Build Node-keyed topology of The FromNode:ToNode
    Arc_From_To = {item:From_Arc.get(From_To2[Arc_From[item]]) for item in Arc_From}    # Build ARCID-keyed topology of the ARCID:ToARCID   ".get()" allows None to be placed in dictionary
    To_From = {From_To[x][1]:From_To[x][0] for x in From_To}                            # Build Node-keyed topology of The ToNode:FromNode
    Arc_To_From = {item:From_Arc.get(To_From.get(Arc_From[item])) for item in Arc_From} # Build ARCID-keyed topology of the ARCID:FromARCID. ".get()" allows None to be placed in dictionary

    # Get the order of segments according to a simple topological sort
    whereclause = "%s > 1" %arcpy.AddFieldDelimiters(outStreams, "Order_")
    Straglers = [row[0] for row in arcpy.da.SearchCursor(outStreams, 'ARCID', whereclause) if Arc_To_From.get(row[0]) is None]    # These are not picked up by the other method
    tic2 = time.time()
    order = sort_topologically_stackless({item:[From_Arc.get(From_To2[Arc_From[item]])] for item in Arc_From})    # 'order' variable is a list of LINK IDs that have been reordered according to a simple topological sort
    loglines.append('        Time elapsed for sorting: %ss' %(time.time()-tic2))
    arcpy.AddMessage(loglines[-1])

    # Fix Streams shapefile from "FROM_NODE" and "TO_NODE" to "FROM_ARCID" and "TO_ARCID"
    arcpy.AddField_management(outStreams, "From_ArcID", "LONG", "#", "#", "#", "#", "NULLABLE")     # Add field for "From_ArcID"
    arcpy.AddField_management(outStreams, "To_ArcID", "LONG", "#", "#", "#", "#", "NULLABLE")       # Add field for "To_ArcID"
    with arcpy.da.UpdateCursor(outStreams, ("ARCID", "From_ArcID", "To_ArcID", "Order_", "index")) as cursor:
        for row in cursor:
            arcid = row[0]
            row[4] = order.index(arcid)                                         # Assign field 'index' with the topologically sorted order (adds a bit of time to the process)
            if arcid in Straglers:
                row[3] = 1                                                      # Deal with issue of some segments being assigned higher orders than they should.
            if Arc_To_From[arcid] is not None:
                row[1] = Arc_To_From.get(arcid)
            if Arc_From_To[arcid] is not None:
                row[2] = Arc_From_To.get(arcid)
            cursor.updateRow(row)
    arcpy.DeleteField_management(outStreams, ['FROM_NODE', 'TO_NODE'])          # Delete node-based fields

    if "NC" in RTfmt:
        # To create a netCDF parameter file
        rootgrp = netCDF4.Dataset(RoutingCSV.replace('.csv', '.nc'), 'w', format=outNCType)     # 'NETCDF4'

        # Create dimensions and set other attribute information
        dim1 = 'linkDim'
        dim = rootgrp.createDimension(dim1, len(order))

        # Create fixed-length variables
        ids = rootgrp.createVariable('link', 'i4', (dim1))                      # Variable (32-bit signed integer)
        froms = rootgrp.createVariable('from','i4',(dim1))                      # Variable (32-bit signed integer)
        tos = rootgrp.createVariable('to','i4',(dim1))                          # Variable (32-bit signed integer)
        slons = rootgrp.createVariable('lon', 'f4', (dim1))                     # Variable (32-bit floating point)
        slats = rootgrp.createVariable('lat', 'f4', (dim1))                     # Variable (32-bit floating point)
        selevs = rootgrp.createVariable('alt', 'f4', (dim1))                    # Variable (32-bit floating point)
        #types = rootgrp.createVariable('type','i4',(dim1))                      # Variable (32-bit signed integer)
        orders = rootgrp.createVariable('order','i4',(dim1))                    # Variable (32-bit signed integer)
        Qis = rootgrp.createVariable('Qi', 'f4', (dim1))                        # Variable (32-bit floating point)
        MusKs = rootgrp.createVariable('MusK','f4',(dim1))                      # Variable (32-bit floating point)
        MusXs = rootgrp.createVariable('MusX', 'f4', (dim1))                    # Variable (32-bit floating point)
        Lengthsnc = rootgrp.createVariable('Length', 'f4', (dim1))              # Variable (32-bit floating point)
        ns = rootgrp.createVariable('n', 'f4', (dim1))                          # Variable (32-bit floating point)
        Sos = rootgrp.createVariable('So', 'f4', (dim1))                        # Variable (32-bit floating point)
        ChSlps = rootgrp.createVariable('ChSlp', 'f4', (dim1))                  # Variable (32-bit floating point)
        BtmWdths = rootgrp.createVariable('BtmWdth','f4',(dim1))                # Variable (32-bit floating point)
        Times = rootgrp.createVariable('time', 'f4')                            # Scalar Variable (32-bit floating point)
        geo_x = rootgrp.createVariable('x', 'f4', (dim1))                       # Variable (32-bit floating point)
        geo_y = rootgrp.createVariable('y', 'f4', (dim1))                       # Variable (32-bit floating point)
        Kcs = rootgrp.createVariable('Kchan', 'i2', (dim1))                        # Variable (16-bit signed integer)

        # Set variable descriptions
        ids.long_name = 'Link ID'
        froms.long_name = 'From Link ID'
        tos.long_name = 'To Link ID'
        slons.long_name = 'longitude of the start node'
        slats.long_name = 'latitude of the start node'
        selevs.long_name = 'Elevation in meters at start node'
        #types.long_name = 'Link type'
        orders.long_name = 'Stream order (Strahler)'
        Qis.long_name = 'Initial flow in link (CMS)'
        MusKs.long_name = 'Muskingum routing time (s)'
        MusXs.long_name = 'Muskingum weighting coefficient'
        Lengthsnc.long_name = 'Stream length (m)'
        ns.long_name = "Manning's roughness"
        Sos.long_name = 'Slope (%; drop/length)'
        ChSlps.long_name = 'Channel side slope (%; drop/length)'
        BtmWdths.long_name = 'Bottom width of channel'
        geo_x.long_name = "x coordinate of projection"
        geo_y.long_name = "y coordinate of projection"
        Kcs.long_name = "channel conductivity"

        # Variable attributes for CF compliance
        slons.units = 'degrees_east'                                            # For compliance
        slats.units = 'degrees_north'                                           # For compliance
        slons.standard_name = 'longitude'                                       # For compliance
        slats.standard_name = 'latitude'                                        # For compliance
        Times.standard_name = 'time'                                            # For compliance
        Times.long_name = 'time of measurement'                                 # For compliance
        Times.units = 'days since 2000-01-01 00:00:00'                          # For compliance
        selevs.standard_name = "height"                                         # For compliance
        selevs.units = "m"                                                      # For compliance
        selevs.positive = "up"                                                  # For compliance
        selevs.axis = "Z"                                                       # For compliance
        ids.cf_role = "timeseries_id"                                           # For compliance
        geo_x.standard_name = "projection_x_coordinate"
        geo_y.standard_name = "projection_y_coordinate"
        geo_x.units = "m"
        geo_y.units = "m"
        Kcs.units = "mm h-2"
        slons.standard_name = 'longitude'                                       # For compliance with NCO
        slats.standard_name = 'latitude'                                        # For compliance with NCO

        # Coordinates for CF-compliance
        froms.coordinates = 'lat lon'                                           # For compliance
        tos.coordinates = 'lat lon'                                             # For compliance
        #types.coordinates = 'lat lon'                                           # For compliance
        orders.coordinates = 'lat lon'                                          # For compliance
        Qis.coordinates = 'lat lon'                                             # For compliance
        MusKs.coordinates = 'lat lon'                                           # For compliance
        MusXs.coordinates = 'lat lon'                                           # For compliance
        Lengthsnc.coordinates = 'lat lon'                                       # For compliance
        ns.coordinates = 'lat lon'                                              # For compliance
        Sos.coordinates = 'lat lon'                                             # For compliance
        ChSlps.coordinates = 'lat lon'                                          # For compliance
        BtmWdths.coordinates = 'lat lon'                                        # For compliance
        geo_x.coordinates = 'lat lon'                                           # For compliance
        geo_y.coordinates = 'lat lon'                                           # For compliance
        Kcs.coordinates = 'lat lon'                                             # For compliance

        # Fill in global attributes
        rootgrp.featureType = 'timeSeries'                                      # For compliance
        rootgrp.history = 'Created %s' %time.ctime()

        loglines.append('        Starting to fill in routing table NC file.')
        arcpy.AddMessage(loglines[-1])

        ids[:] = numpy.array(order)                                                             # Fill in id field information
        fromnodes = [From_To[arcid][0] for arcid in order]                                      # The FROM node from the streams shapefile (used later as a key))
        tonodes = [From_To[arcid][1] for arcid in order]                                        # The TO node from the streams shapefile (used later as a key)
        drops = [NodeElev[fromnode]-NodeElev[tonode] for fromnode, tonode in zip(fromnodes, tonodes)]
        drops = [x if x>0 else 0 for x in drops]                                                # Replace negative values with 0

        # Set variable value arrays
        fromlist = [Arc_To_From[arcid] for arcid in order]                                      # List containes 'None' values, must convert to numpy.nan
        tolist = [Arc_From_To[arcid] for arcid in order]                                        # List containes 'None' values, must convert to numpy.nan

        # Change None values to 0.  Could alternatively use numpy.nan
        froms[:] = numpy.array(map(lambda x: 0 if x==None else x, fromlist))                    # Note that the From in this case is the ARCID of any of the immediately upstream contributing segments
        tos[:] = numpy.array(map(lambda x: 0 if x==None else x, tolist))                        # Note that the To in this case is the ARCID of the immediately downstream segment

        # Fill in other variables
        slons[:] = numpy.array([NodesLL[fromnode][0] for fromnode in fromnodes])
        slats[:] = numpy.array([NodesLL[fromnode][1] for fromnode in fromnodes])
        geo_x[:] = numpy.array([NodesXY[fromnode][0] for fromnode in fromnodes])
        geo_y[:] = numpy.array([NodesXY[fromnode][1] for fromnode in fromnodes])
        selevs[:] = numpy.array([round(NodeElev[fromnode], 3) for fromnode in fromnodes])       # Round to 3 digits
        Lengthsnc[:] = numpy.array([round(Lengths[arcid], 1) for arcid in order])               # Round to 1 digit

        # Modify order and slope arrays
        order_ = [1 if arcid in Straglers else StrOrder[From_To[arcid][0]] for arcid in order]  # Deal with issue of some segments being assigned higher orders than they should.
        orders[:] = numpy.array(order_)
        Sos_ = numpy.round(numpy.array(drops).astype(float)/Lengthsnc[:], 3)    # Must convert list to float to result in floats
        numpy.place(Sos_, Sos_==0, [0.005])                                     # Set minimum slope to be 0.005
        Sos[:] = Sos_[:]

        # Set default arrays
        #types[:] = type_
        Qis[:] = Qi
        MusKs[:] = MusK
        MusXs[:] = MusX
        ns[:] = n
        ChSlps[:] = ChSlp
        BtmWdths[:] = BtmWdth
        Times[:] = 0
        Kcs[:] = Kc

        loglines.append('        Done writing NC file to disk.')
        arcpy.AddMessage(loglines[-1])

        # Close file
        rootgrp.close()

    # To create an ASCII CSV file:
    if "CSV" in RTfmt:
        with open(RoutingCSV, 'wb') as fp:
            a = csv.writer(fp, delimiter=',')
            data = [['link', 'from', 'to', 'start_lon', 'start_lat', 'start_elev', 'order', 'Qi', 'MusK', 'MusX', 'Length', 'n', 'So', 'ChSlp', 'BtmWdth', 'Kc']]     # 'type', 'LkHZArea', 'LkMxH', 'WeirC', 'WeirL', 'OrificeC', 'OrificeA', 'OrificeE']]
            for arcid in order:
                fromnode = From_To[arcid][0]                                    # The FROM node from the streams shapefile (used later as a key)
                tonode = From_To[arcid][1]                                      # The TO node from the streams shapefile (used later as a key)
                #from_ = Arc_To_From[arcid]                                     # Note that the From in this case is the ARCID of ANY of the immediately upstream contributing segments
                from_ = 0                                                       # I think we set this to 0 because From doesn't really tell us anything (9/3/2015)
                to = Arc_From_To.get(arcid, 0)                                  # Note that the To in this case is the ARCID of the immediately downstream segment
                start_lon = NodesLL[fromnode][0]
                start_lat = NodesLL[fromnode][1]
                start_elev = round(NodeElev[fromnode], 3)                       # Round to 3 digits
                Length = round(Lengths[arcid], 1)                               # Round to 1 digit
                drop = NodeElev[fromnode]-NodeElev[tonode]

                # Deal with issue of some segments being assigned higher orders than they should.
                if arcid in Straglers:
                    order_ = 1
                else:
                    order_ = StrOrder[fromnode]
                if drop <= 0:                                                   # Set negative slopes to 0
                    drop = 0
                So = round(drop/Length, 3)                                      # Round to 3 digits
                if So < 0.005:
                    So == 0.005                                                 # Set minimum slope to be 0.005
                data.append([arcid, from_, to, start_lon, start_lat, start_elev, order_, Qi, MusK, MusX, Length, n, So, ChSlp, BtmWdth, Kc])     # type_, LkHZArea, LkMxH, WeirC, WeirL, OrificeC, OrificeA, OrificeE])
            a.writerows(data)
        loglines.append('        Routing Table:\n            %s Lines\n            %s Nodes.' %(len(From_To.keys()), len(NodesLL.keys())))
        arcpy.AddMessage(loglines[-1])
        loglines.append('        Done writing CSV table to disk.')
        arcpy.AddMessage(loglines[-1])

    loglines.append('    Routing table created without error.')
    arcpy.AddMessage(loglines[-1])
    return outRaster, loglines

def add_reservoirs(arcpy, channelgrid, in_lakes, flac, projdir, fill2, cellsize, sr2, loglines, lakeIDfield=None, Threshold=0):
    """This function is intended to add reservoirs into the model grid stack, such
    that the channelgrid and lake grids are modified to accomodate reservoirs and
    lakes."""

    # Get information about the input domain and set environments
    arcpy.env.cellSize = cellsize
    arcpy.env.extent =  channelgrid
    arcpy.env.outputCoordinateSystem = sr2
    arcpy.env.snapRaster = channelgrid

    # Use extent of channelgrid raster to add a feature layer of lake polygons
    outshp = projdir + os.path.sep + 'in_lakes_clip.shp'
    arcpy.CopyFeatures_management(in_lakes, outshp)

    # Create new area field
    Field1 = 'AREASQKM'
    if Field1 not in [field.name for field in arcpy.ListFields(outshp)]:
        arcpy.AddField_management(outshp, Field1, "FLOAT")
        arcpy.CalculateField_management(outshp, Field1, '!shape.area@squarekilometers!', "PYTHON_9.3")
    field1delim = arcpy.AddFieldDelimiters(outshp, Field1)

    # Use extent of channelgrid raster to add a feature layer of lake polygons
    Field2 = 'FTYPE'
    if Field2 in [field.name for field in arcpy.ListFields(outshp)]:
        Values = ["'LakePond'", "'Reservoir'"]
        field2delim = arcpy.AddFieldDelimiters(outshp, Field2)
        ftypesql = ' AND (%s = %s OR %s = %s)' %(field2delim, Values[0], field2delim, Values[1])
    else:
        ftypesql = ''
    where_clause = """%s >= %s%s""" %(field1delim, Threshold, ftypesql)
    arcpy.MakeFeatureLayer_management(outshp, "Lakeslyr", where_clause)

    if lakeIDfield is None:
        loglines.append('    Adding auto-incremented lake ID field (1...n)')
        arcpy.AddMessage(loglines[-1])

        # Renumber lakes 1-n (addfield, calculate field)
        lakeID = "NEWID"
        arcpy.AddField_management("Lakeslyr", lakeID, "LONG")
        expression = 'autoIncrement()'
        code_block = """rec = 0\ndef autoIncrement():\n    global rec\n    pStart = 1\n    pInterval = 1\n    if (rec == 0):\n        rec = pStart\n    else:\n        rec = rec + pInterval\n    return rec"""
        arcpy.CalculateField_management("Lakeslyr", lakeID, expression, "PYTHON_9.3", code_block)
    else:
        loglines.append('    Using provided lake ID field: %s' %lakeIDfield)
        arcpy.AddMessage(loglines[-1])
        lakeID = lakeIDfield                                                    # Use existing field specified by 'lakeIDfield' parameter

    # Create a raster from the lake polygons that matches the channelgrid layer
    outRastername = os.path.join(projdir, "Lakesras")
    outfeatures = os.path.join(projdir, 'lakes.shp')
    arcpy.CopyFeatures_management("Lakeslyr", outfeatures)
    arcpy.PolygonToRaster_conversion(outfeatures, lakeID, outRastername, "MAXIMUM_AREA")      # This tool requires ArcGIS for Desktop Advanced OR Spatial Analyst Extension
    #arcpy.FeatureToRaster_conversion(outfeatures, lakeID, outRastername)        # This tool requires only ArcGIS for Desktop Basic, but does not allow a priority field

    # Hack to convert Lakesras to 16bit integer
    outRaster1 = (channelgrid * 0) + Raster(outRastername)                      # Create 16bit ratser for Lakesras out of channelgrid
    outRaster = Con(IsNull(outRaster1)==1, NoDataVal, outRaster1)               # Convert Null or NoData to -9999

    # Con statement for lakes where channelgrid = -9999
    zonstat = ZonalStatistics(outRastername, "Value", flac, "MAXIMUM")          # Get maximum flow accumulation value for each lake
    lakeacc = Con(outRastername, flac)                                          # Flow accumulation over lakes only
    TestCon = Con(lakeacc == zonstat, outRastername, NoDataVal)                 # Bottom of lake channelgrid pixel = lake number
    NewChannelgrid = Con(IsNull(TestCon)==1, channelgrid, TestCon)              # This is the new lake-routed channelgrid layer

    # Now march down a set number of pixels to get minimum lake elevation
    tolerance = int(float(arcpy.GetRasterProperties_management(channelgrid, 'CELLSIZEX').getOutput(0)) * LK_walker)
    SnapPour = SnapPourPoint(SetNull(TestCon, TestCon, "VALUE = %s" %NoDataVal), flac, tolerance)                          # Snap lake outlet pixel to FlowAccumulation with tolerance
    del flac
    outTable2 = r'in_memory/zonalstat2'
    Sample(fill2, SnapPour, outTable2, "NEAREST")
    min_elevs = {row[1]: row[-1] for row in arcpy.da.SearchCursor(outTable2, "*")}

    # Convert lakes raster to polygons and gather size & elevations
    loglines.append('    Gathering lake parameter information.')
    arcpy.AddMessage(loglines[-1])
    outTable = r'in_memory/zonalstat'
    zontable = ZonalStatisticsAsTable(outRastername, 'VALUE', fill2, outTable, "DATA", "MIN_MAX_MEAN")
    areas = {row[0]: row[1] for row in arcpy.da.SearchCursor(outTable, ['VALUE', 'AREA'])}                      # Searchcursor on zonal stats table
    max_elevs = {row[0]: row[1] for row in arcpy.da.SearchCursor(outTable, ['VALUE', 'MAX'])}                   # Searchcursor on zonal stats table

    OrificEs = {x:(min_elevs[x] + ((max_elevs[x] - min_elevs[x])/3)) for x in min_elevs.keys()}             # OrificElevation is 1/3 between the low elevation and max lake elevation
    Elevation = {x:(min_elevs[x] + (((max_elevs[x] - min_elevs[x])/3) * 2)) for x in min_elevs.keys()}      # Elevation is 2/3 between the low elevation and max lake elevation
    WeirH_vals = {x:(min_elevs[x] + ((max_elevs[x] - min_elevs[x]) * 0.9)) for x in min_elevs.keys()}       # WierH is 0.9 of the distance between the low elevation and max lake elevation

    #  Gather centroid lat/lons
    out_lake_raster = os.path.join(projdir, "out_lake_raster.shp")
    out_lake_raster_dis = os.path.join(projdir, "out_lake_raster_dissolve.shp")
    arcpy.RasterToPolygon_conversion(outRastername, out_lake_raster, "NO_SIMPLIFY", "VALUE")
    arcpy.Dissolve_management(out_lake_raster, out_lake_raster_dis, "GRIDCODE", "", "MULTI_PART")               # Dissolve to eliminate multipart features
    arcpy.Delete_management(out_lake_raster)                                                                    # Added 9/4/2015

    # Create a point geometry object from gathered lake centroid points
    loglines.append('    Starting to gather lake centroid information.')
    arcpy.AddMessage(loglines[-1])
    sr1 = arcpy.SpatialReference(4326)                                          # GCS_WGS_1984
    point = arcpy.Point()
    cen_lats = {}
    cen_lons = {}
    shapes = {row[0]: row[1] for row in arcpy.da.SearchCursor(out_lake_raster_dis, ['GRIDCODE', 'SHAPE@XY'])}
    arcpy.Delete_management(out_lake_raster_dis)                                                                # Added 9/4/2015
    for shape in shapes:
        point.X = shapes[shape][0]
        point.Y = shapes[shape][1]
        pointGeometry = arcpy.PointGeometry(point, sr2)
        projpoint = pointGeometry.projectAs(sr1)                                # Optionally add transformation method:
        cen_lats[shape] = projpoint.firstPoint.Y
        cen_lons[shape] = projpoint.firstPoint.X
    loglines.append('    Done gathering lake centroid information.')
    arcpy.AddMessage(loglines[-1])

    # Create Lake parameter file
    LakeTBL = os.path.join(projdir, 'LAKEPARM.TBL')
    loglines.append('    Starting to create lake parameter table.')
    arcpy.AddMessage(loglines[-1])

    loglines.append('        Lakes Table: %s Lakes' %len(areas.keys()))
    arcpy.AddMessage(loglines[-1])

    # Create NetCDF output table
    if "NC" in LKfmt:

        # To create a netCDF parameter file
        rootgrp = netCDF4.Dataset(LakeTBL.replace('.TBL', '.nc'), 'w', format=outNCType)     # 'NETCDF4'

        # Create dimensions and set other attribute information
        dim1 = 'nlakes'
        #dim2 = 'DOY'
        dim = rootgrp.createDimension(dim1, len(min_elevs))
        #DOY = rootgrp.createDimension(dim2, 365)

        # Create coordinate variables
        ids = rootgrp.createVariable('lake_id','i4',(dim1))                          # Variable (32-bit signed integer)
        ids[:] = numpy.array(min_elevs.keys())                                  # Variable (32-bit signed integer)

        # Create fixed-length variables
        LkAreas = rootgrp.createVariable('LkArea','f8',(dim1))                  # Variable (64-bit floating point)
        LkMxHs = rootgrp.createVariable('LkMxH', 'f8', (dim1))                  # Variable (64-bit floating point)
        WeirCs = rootgrp.createVariable('WeirC', 'f8', (dim1))                  # Variable (64-bit floating point)
        #WeirCs = rootgrp.createVariable('WeirC', 'f8', (dim1, dim2))            # Variable (64-bit floating point)
        WeirLs = rootgrp.createVariable('WeirL', 'f8', (dim1))                  # Variable (64-bit floating point)
        OrificeCs = rootgrp.createVariable('OrificeC', 'f8', (dim1))            # Variable (64-bit floating point)
        #OrificeCs = rootgrp.createVariable('OrificeC', 'f8', (dim1, dim2))      # Variable (64-bit floating point)
        OrificeAs = rootgrp.createVariable('OrificeA', 'f8', (dim1))            # Variable (64-bit floating point)
        OrificeEs = rootgrp.createVariable('OrificeE', 'f8', (dim1))            # Variable (64-bit floating point)
        lats = rootgrp.createVariable('lat', 'f4', (dim1))                      # Variable (32-bit floating point)
        longs = rootgrp.createVariable('lon', 'f4', (dim1))                     # Variable (32-bit floating point)
        elevations = rootgrp.createVariable('alt', 'f4', (dim1))                # Variable (32-bit floating point)
        Times = rootgrp.createVariable('time', 'f8', (dim1))                    # Variable (64-bit floating point)
        Discharges = rootgrp.createVariable('Discharge', 'f8', (dim1))          # Variable (64-bit floating point)
        WeirHs = rootgrp.createVariable('WeirH', 'f8', (dim1))                  # Variable (64-bit floating point)
        AscendOrder = rootgrp.createVariable('ascendingIndex', 'i4', (dim1))    # Variable (32-bit signed integer)

        # Set variable descriptions
        ids.long_name = 'Lake ID'
        LkAreas.long_name = 'Gridded lake area (sq. km)'
        LkMxHs.long_name = 'Maximum lake elevation (m ASL)'
        WeirCs.long_name = 'Weir coefficient'
        WeirLs.long_name = 'Weir length (m)'
        OrificeCs.long_name = 'Orifice coefficient'
        OrificeAs.long_name = 'Orifice cross-sectional area (sq. m)'
        OrificeEs.long_name = 'Orifice elevation (m ASL)'
        WeirHs.long_name = 'Weir Height (m ASL)'
        lats.long_name = 'latitude of the lake centroid'
        longs.long_name = 'longitude of the lake centroid'
        elevations.long_name = 'vertical distance above mean sea level in (m ASL)'
        Discharges.long_name = 'Default Discharge'                              # 'Climatological Discharge' if specifying a non-default or time-varying discharge
        longs.units = 'degrees_east'                                            # For compliance
        lats.units = 'degrees_north'                                            # For compliance
        longs.standard_name = 'longitude'                                       # For compliance
        lats.standard_name = 'latitude'                                         # For compliance
        Times.standard_name = 'time'                                            # For compliance
        Times.long_name = 'time of measurement'                                 # For compliance
        Times.units = 'days since 2000-01-01 00:00:00'                          # For compliance
        elevations.standard_name = 'height'                                     # For compliance
        elevations.units = 'm'                                                  # For compliance
        elevations.positive = 'up'                                              # For compliance
        elevations.axis = 'Z'                                                   # For compliance
        Discharges.units = 'CMS'
        WeirHs.units = 'm'
        ids.cf_role = "timeseries_id"                                           # For compliance
        AscendOrder.long_name = 'Index to use for sorting IDs (ascending)'

        # Coordinates for CF-compliance
        LkAreas.coordinates = 'lat lon'                                         # For compliance
        LkMxHs.coordinates = 'lat lon'                                          # For compliance
        WeirCs.coordinates = 'lat lon'                                          # For compliance
        WeirLs.coordinates = 'lat lon'                                          # For compliance
        OrificeCs.coordinates = 'lat lon'                                       # For compliance
        OrificeAs.coordinates = 'lat lon'                                       # For compliance
        OrificeEs.coordinates = 'lat lon'                                       # For compliance
        Discharges.coordinates = 'lat lon'                                      # For compliance
        WeirHs.coordinates = 'lat lon'                                          # For compliance

        # Fill in global attributes
        rootgrp.featureType = 'timeSeries'                                      # For compliance
        rootgrp.history = 'Created %s' %time.ctime()

        loglines.append('        Starting to fill in lake parameter table NC file.')
        arcpy.AddMessage(loglines[-1])

        AscendOrder[:] = numpy.argsort(ids[:])                                  # Use argsort to give the ascending sort order for IDs. Added by KMS 4/4/2017
        LkAreas[:] = numpy.array([float(areas[lkid])/float(1000000) for lkid in min_elevs.keys()])  # Divide by 1M for kilometers^2
        LkMxHs[:] = numpy.array([max_elevs[lkid] for lkid in min_elevs.keys()])
        WeirCs[:] = WeirC
        WeirLs[:] = WeirL
        OrificeCs[:] = OrificeC
        OrificeAs[:] = OrificA
        Times[:] = 0
        OrificeEs[:] = numpy.array([OrificEs[lkid] for lkid in min_elevs.keys()])           # Orifice Elevation is 1/3 between 'min' and max lake elevation.
        lats[:] = numpy.array([cen_lats[lkid] for lkid in min_elevs.keys()])
        longs[:] = numpy.array([cen_lons[lkid] for lkid in min_elevs.keys()])
        elevations[:] = numpy.array([Elevation[lkid] for lkid in min_elevs.keys()])         # Base Elevation is 2/3 betwen 'min' and max lake elevation.
        WeirHs[:] = numpy.array([WeirH_vals[lkid] for lkid in min_elevs.keys()])            # WierH is 0.9 of the distance between the low elevation and max lake elevation
        #Discharges[:] = 0

        loglines.append('        Done writing LAKEPARM.nc table to disk.')
        arcpy.AddMessage(loglines[-1])

        # Close file
        rootgrp.close()

    # Create .TBL output table (ASCII)
    if "TBL" in LKfmt:

        with open(LakeTBL, 'wb') as fp:
            a = csv.writer(fp, dialect='excel-tab', quoting=csv.QUOTE_NONE)
            #a.writerow(['lake', 'LkArea', 'LkMxH', 'WeirC', 'WeirL', 'OrificeC', 'OrificeA', 'OrificeE', 'lat', 'long', 'elevation', 'WeirH'])
            for lkid in min_elevs.keys():
                lkarea = float(areas[lkid])/float(1000000)                      # Divide by 1M for kilometers^2
                lkmaxelev = max_elevs[lkid]
                baseelev = Elevation[lkid]                                      # Base Elevation is 2/3 betwen 'min' and max lake elevation.
                OrificeE = OrificEs[lkid]                                       # Orifice Elevation is 1/3 between 'min' and max lake elevation.
                cen_lat = cen_lats[lkid]
                cen_lon = cen_lons[lkid]
                WeirH = WeirH_vals[lkid]
                a.writerow([lkid, lkarea, lkmaxelev, WeirC, WeirL, OrificeC, OrificA, OrificeE, cen_lat, cen_lon, baseelev, WeirH])   #COMID?

        loglines.append('        Done writing LAKEPARM.TBL table to disk.')
        arcpy.AddMessage(loglines[-1])

    # Process: Output Channelgrid
    channelgrid_arr = arcpy.RasterToNumPyArray(NewChannelgrid)
    outRaster_arr = arcpy.RasterToNumPyArray(outRaster)

    # Clean up by deletion
    loglines.append('    Lake parameter table created without error.')
    arcpy.AddMessage(loglines[-1])

    # Clean up
    del sr1, point, outRaster1
    arcpy.Delete_management("Lakeslyr")
    arcpy.Delete_management(outshp)
    #arcpy.Delete_management(outfeatures)
    arcpy.Delete_management('in_memory')

    #outRaster.save(outRastername)
    del outRaster                                                               # Added 9/4/2015
    del channelgrid                                                             # Added 9/6/2015
    del projdir, fill2, cellsize, sr2, in_lakes
    return arcpy, channelgrid_arr, outRaster_arr, loglines

def add_reservoirs_nogrid(arcpy, channelgrid, in_lakes, flac, projdir, fill2, cellsize, sr2, loglines, lakeIDfield=None):
    """
    This function is intended to add reservoirs into the model grid stack, such
    that the channelgrid and lake grids are modified to accomodate reservoirs and
    lakes.

    This version does not attempt to subset the lakes by a size threshold, nor
    does it filter based on FTYPE.
    3/23/2017
    """

    # Get information about the input domain and set environments
    arcpy.env.cellSize = cellsize
    arcpy.env.extent =  channelgrid
    arcpy.env.outputCoordinateSystem = sr2
    arcpy.env.snapRaster = channelgrid

    # Use extent of channelgrid raster to add a feature layer of lake polygons
    outshp = projdir + os.path.sep + 'in_lakes_clip.shp'
    arcpy.CopyFeatures_management(in_lakes, outshp)

    # Create new area field
    Field1 = 'AREASQKM'
    if Field1 not in [field.name for field in arcpy.ListFields(outshp)]:
        arcpy.AddField_management(outshp, Field1, "FLOAT")
        arcpy.CalculateField_management(outshp, Field1, '!shape.area@squarekilometers!', "PYTHON_9.3")
    field1delim = arcpy.AddFieldDelimiters(outshp, Field1)
    arcpy.MakeFeatureLayer_management(outshp, "Lakeslyr")

    if lakeIDfield is None:
        loglines.append('    Adding auto-incremented lake ID field (1...n)')
        arcpy.AddMessage(loglines[-1])

        # Renumber lakes 1-n (addfield, calculate field)
        lakeID = "NEWID"
        arcpy.AddField_management("Lakeslyr", lakeID, "LONG")
        expression = 'autoIncrement()'
        code_block = """rec = 0\ndef autoIncrement():\n    global rec\n    pStart = 1\n    pInterval = 1\n    if (rec == 0):\n        rec = pStart\n    else:\n        rec = rec + pInterval\n    return rec"""
        arcpy.CalculateField_management("Lakeslyr", lakeID, expression, "PYTHON_9.3", code_block)
    else:
        loglines.append('    Using provided lake ID field: %s' %lakeIDfield)
        arcpy.AddMessage(loglines[-1])
        lakeID = lakeIDfield                                                    # Use existing field specified by 'lakeIDfield' parameter

    # Create a raster from the lake polygons that matches the channelgrid layer
    outRastername = os.path.join(projdir, "Lakesras")
    outfeatures = os.path.join(projdir, 'lakes.shp')
    arcpy.CopyFeatures_management("Lakeslyr", outfeatures)
    arcpy.PolygonToRaster_conversion(outfeatures, lakeID, outRastername, "MAXIMUM_AREA")      # This tool requires ArcGIS for Desktop Advanced OR Spatial Analyst Extension
    #arcpy.FeatureToRaster_conversion(outfeatures, lakeID, outRastername)        # This tool requires only ArcGIS for Desktop Basic, but does not allow a priority field

    # Hack to convert Lakesras to 16bit integer
    outRaster1 = (channelgrid * 0) + Raster(outRastername)                      # Create 16bit ratser for Lakesras out of channelgrid
    outRaster = Con(IsNull(outRaster1)==1, NoDataVal, outRaster1)               # Convert Null or NoData to -9999

    # Con statement for lakes where channelgrid = -9999
    zonstat = ZonalStatistics(outRastername, "Value", flac, "MAXIMUM")          # Get maximum flow accumulation value for each lake
    lakeacc = Con(outRastername, flac)                                          # Flow accumulation over lakes only
    TestCon = Con(lakeacc == zonstat, outRastername, NoDataVal)                 # Bottom of lake channelgrid pixel = lake number
    NewChannelgrid = Con(IsNull(TestCon)==1, channelgrid, TestCon)              # This is the new lake-routed channelgrid layer

    # Now march down a set number of pixels to get minimum lake elevation
    tolerance = int(float(arcpy.GetRasterProperties_management(channelgrid, 'CELLSIZEX').getOutput(0)) * LK_walker)
    SnapPour = SnapPourPoint(SetNull(TestCon, TestCon, "VALUE = %s" %NoDataVal), flac, tolerance)                          # Snap lake outlet pixel to FlowAccumulation with tolerance
    del flac
    outTable2 = r'in_memory/zonalstat2'
    Sample(fill2, SnapPour, outTable2, "NEAREST")
    min_elevs = {row[1]: row[-1] for row in arcpy.da.SearchCursor(outTable2, "*")}

    # Convert lakes raster to polygons and gather size & elevations
    loglines.append('    Gathering lake parameter information.')
    arcpy.AddMessage(loglines[-1])
    outTable = r'in_memory/zonalstat'
    zontable = ZonalStatisticsAsTable(outRastername, 'VALUE', fill2, outTable, "DATA", "MIN_MAX_MEAN")
    areas = {row[0]: row[1] for row in arcpy.da.SearchCursor(outTable, ['VALUE', 'AREA'])}                      # Searchcursor on zonal stats table
    max_elevs = {row[0]: row[1] for row in arcpy.da.SearchCursor(outTable, ['VALUE', 'MAX'])}                   # Searchcursor on zonal stats table

    OrificEs = {x:(min_elevs[x] + ((max_elevs[x] - min_elevs[x])/3)) for x in min_elevs.keys()}             # OrificElevation is 1/3 between the low elevation and max lake elevation
    Elevation = {x:(min_elevs[x] + (((max_elevs[x] - min_elevs[x])/3) * 2)) for x in min_elevs.keys()}      # Elevation is 2/3 between the low elevation and max lake elevation
    WeirH_vals = {x:(min_elevs[x] + ((max_elevs[x] - min_elevs[x]) * 0.9)) for x in min_elevs.keys()}       # WierH is 0.9 of the distance between the low elevation and max lake elevation

    #  Gather centroid lat/lons
    out_lake_raster = os.path.join(projdir, "out_lake_raster.shp")
    out_lake_raster_dis = os.path.join(projdir, "out_lake_raster_dissolve.shp")
    arcpy.RasterToPolygon_conversion(outRastername, out_lake_raster, "NO_SIMPLIFY", "VALUE")
    arcpy.Dissolve_management(out_lake_raster, out_lake_raster_dis, "GRIDCODE", "", "MULTI_PART")               # Dissolve to eliminate multipart features
    arcpy.Delete_management(out_lake_raster)                                                                    # Added 9/4/2015

    # Create a point geometry object from gathered lake centroid points
    loglines.append('    Starting to gather lake centroid information.')
    arcpy.AddMessage(loglines[-1])
    sr1 = arcpy.SpatialReference(4326)                                          # GCS_WGS_1984
    point = arcpy.Point()
    cen_lats = {}
    cen_lons = {}
    shapes = {row[0]: row[1] for row in arcpy.da.SearchCursor(out_lake_raster_dis, ['GRIDCODE', 'SHAPE@XY'])}
    arcpy.Delete_management(out_lake_raster_dis)                                                                # Added 9/4/2015
    for shape in shapes:
        point.X = shapes[shape][0]
        point.Y = shapes[shape][1]
        pointGeometry = arcpy.PointGeometry(point, sr2)
        projpoint = pointGeometry.projectAs(sr1)                                # Optionally add transformation method:
        cen_lats[shape] = projpoint.firstPoint.Y
        cen_lons[shape] = projpoint.firstPoint.X
    loglines.append('    Done gathering lake centroid information.')
    arcpy.AddMessage(loglines[-1])

    # Create Lake parameter file
    LakeTBL = os.path.join(projdir, 'LAKEPARM.TBL')
    loglines.append('    Starting to create lake parameter table.')
    arcpy.AddMessage(loglines[-1])

    loglines.append('        Lakes Table: %s Lakes' %len(areas.keys()))
    arcpy.AddMessage(loglines[-1])

    # Create NetCDF output table
    if "NC" in LKfmt:

        # To create a netCDF parameter file
        rootgrp = netCDF4.Dataset(LakeTBL.replace('.TBL', '.nc'), 'w', format=outNCType)     # 'NETCDF4'

        # Create dimensions and set other attribute information
        dim1 = 'nlakes'
        #dim2 = 'DOY'
        dim = rootgrp.createDimension(dim1, len(min_elevs))
        #DOY = rootgrp.createDimension(dim2, 365)

        # Create coordinate variables
        ids = rootgrp.createVariable('lake_id','i4',(dim1))                          # Variable (32-bit signed integer)
        ids[:] = numpy.array(min_elevs.keys())                                  # Variable (32-bit signed integer)

        # Create fixed-length variables
        LkAreas = rootgrp.createVariable('LkArea','f8',(dim1))                  # Variable (64-bit floating point)
        LkMxHs = rootgrp.createVariable('LkMxH', 'f8', (dim1))                  # Variable (64-bit floating point)
        WeirCs = rootgrp.createVariable('WeirC', 'f8', (dim1))                  # Variable (64-bit floating point)
        #WeirCs = rootgrp.createVariable('WeirC', 'f8', (dim1, dim2))            # Variable (64-bit floating point)
        WeirLs = rootgrp.createVariable('WeirL', 'f8', (dim1))                  # Variable (64-bit floating point)
        OrificeCs = rootgrp.createVariable('OrificeC', 'f8', (dim1))            # Variable (64-bit floating point)
        #OrificeCs = rootgrp.createVariable('OrificeC', 'f8', (dim1, dim2))      # Variable (64-bit floating point)
        OrificeAs = rootgrp.createVariable('OrificeA', 'f8', (dim1))            # Variable (64-bit floating point)
        OrificeEs = rootgrp.createVariable('OrificeE', 'f8', (dim1))            # Variable (64-bit floating point)
        lats = rootgrp.createVariable('lat', 'f4', (dim1))                      # Variable (32-bit floating point)
        longs = rootgrp.createVariable('lon', 'f4', (dim1))                     # Variable (32-bit floating point)
        elevations = rootgrp.createVariable('alt', 'f4', (dim1))                # Variable (32-bit floating point)
        Times = rootgrp.createVariable('time', 'f8', (dim1))                    # Variable (64-bit floating point)
        Discharges = rootgrp.createVariable('Discharge', 'f8', (dim1))          # Variable (64-bit floating point)
        WeirHs = rootgrp.createVariable('WeirH', 'f8', (dim1))                  # Variable (64-bit floating point)
        AscendOrder = rootgrp.createVariable('ascendingIndex', 'i4', (dim1))    # Variable (32-bit signed integer)

        # Set variable descriptions
        ids.long_name = 'Lake ID'
        LkAreas.long_name = 'Gridded lake area (sq. km)'
        LkMxHs.long_name = 'Maximum lake elevation (m ASL)'
        WeirCs.long_name = 'Weir coefficient'
        WeirLs.long_name = 'Weir length (m)'
        OrificeCs.long_name = 'Orifice coefficient'
        OrificeAs.long_name = 'Orifice cross-sectional area (sq. m)'
        OrificeEs.long_name = 'Orifice elevation (m ASL)'
        WeirHs.long_name = 'Weir Height (m ASL)'
        lats.long_name = 'latitude of the lake centroid'
        longs.long_name = 'longitude of the lake centroid'
        elevations.long_name = 'vertical distance above mean sea level in (m ASL)'
        Discharges.long_name = 'Default Discharge'                              # 'Climatological Discharge' if specifying a non-default or time-varying discharge
        longs.units = 'degrees_east'                                            # For compliance
        lats.units = 'degrees_north'                                            # For compliance
        longs.standard_name = 'longitude'                                       # For compliance
        lats.standard_name = 'latitude'                                         # For compliance
        Times.standard_name = 'time'                                            # For compliance
        Times.long_name = 'time of measurement'                                 # For compliance
        Times.units = 'days since 2000-01-01 00:00:00'                          # For compliance
        elevations.standard_name = 'height'                                     # For compliance
        elevations.units = 'm'                                                  # For compliance
        elevations.positive = 'up'                                              # For compliance
        elevations.axis = 'Z'                                                   # For compliance
        Discharges.units = 'CMS'
        WeirHs.units = 'm'
        ids.cf_role = "timeseries_id"                                           # For compliance
        AscendOrder.long_name = 'Index to use for sorting IDs (ascending)'

        # Coordinates for CF-compliance
        LkAreas.coordinates = 'lat lon'                                         # For compliance
        LkMxHs.coordinates = 'lat lon'                                          # For compliance
        WeirCs.coordinates = 'lat lon'                                          # For compliance
        WeirLs.coordinates = 'lat lon'                                          # For compliance
        OrificeCs.coordinates = 'lat lon'                                       # For compliance
        OrificeAs.coordinates = 'lat lon'                                       # For compliance
        OrificeEs.coordinates = 'lat lon'                                       # For compliance
        Discharges.coordinates = 'lat lon'                                      # For compliance
        WeirHs.coordinates = 'lat lon'                                          # For compliance

        # Fill in global attributes
        rootgrp.featureType = 'timeSeries'                                      # For compliance
        rootgrp.history = 'Created %s' %time.ctime()

        loglines.append('        Starting to fill in lake parameter table NC file.')
        arcpy.AddMessage(loglines[-1])

        AscendOrder[:] = numpy.argsort(ids[:])                                  # Use argsort to give the ascending sort order for IDs. Added by KMS 4/4/2017
        LkAreas[:] = numpy.array([float(areas[lkid])/float(1000000) for lkid in min_elevs.keys()])  # Divide by 1M for kilometers^2
        LkMxHs[:] = numpy.array([max_elevs[lkid] for lkid in min_elevs.keys()])
        WeirCs[:] = WeirC
        WeirLs[:] = WeirL
        OrificeCs[:] = OrificeC
        OrificeAs[:] = OrificA
        Times[:] = 0
        OrificeEs[:] = numpy.array([OrificEs[lkid] for lkid in min_elevs.keys()])           # Orifice Elevation is 1/3 between 'min' and max lake elevation.
        lats[:] = numpy.array([cen_lats[lkid] for lkid in min_elevs.keys()])
        longs[:] = numpy.array([cen_lons[lkid] for lkid in min_elevs.keys()])
        elevations[:] = numpy.array([Elevation[lkid] for lkid in min_elevs.keys()])         # Base Elevation is 2/3 betwen 'min' and max lake elevation.
        WeirHs[:] = numpy.array([WeirH_vals[lkid] for lkid in min_elevs.keys()])            # WierH is 0.9 of the distance between the low elevation and max lake elevation
        #Discharges[:] = 0

        loglines.append('        Done writing LAKEPARM.nc table to disk.')
        arcpy.AddMessage(loglines[-1])

        # Close file
        rootgrp.close()

    # Create .TBL output table (ASCII)
    if "TBL" in LKfmt:

        with open(LakeTBL, 'wb') as fp:
            a = csv.writer(fp, dialect='excel-tab', quoting=csv.QUOTE_NONE)
            #a.writerow(['lake', 'LkArea', 'LkMxH', 'WeirC', 'WeirL', 'OrificeC', 'OrificeA', 'OrificeE', 'lat', 'long', 'elevation', 'WeirH'])
            for lkid in min_elevs.keys():
                lkarea = float(areas[lkid])/float(1000000)                      # Divide by 1M for kilometers^2
                lkmaxelev = max_elevs[lkid]
                baseelev = Elevation[lkid]                                      # Base Elevation is 2/3 betwen 'min' and max lake elevation.
                OrificeE = OrificEs[lkid]                                       # Orifice Elevation is 1/3 between 'min' and max lake elevation.
                cen_lat = cen_lats[lkid]
                cen_lon = cen_lons[lkid]
                WeirH = WeirH_vals[lkid]
                a.writerow([lkid, lkarea, lkmaxelev, WeirC, WeirL, OrificeC, OrificA, OrificeE, cen_lat, cen_lon, baseelev, WeirH])   #COMID?

        loglines.append('        Done writing LAKEPARM.TBL table to disk.')
        arcpy.AddMessage(loglines[-1])

    # Process: Output Channelgrid
    channelgrid_arr = arcpy.RasterToNumPyArray(NewChannelgrid)
    outRaster_arr = arcpy.RasterToNumPyArray(outRaster)

    # Clean up by deletion
    loglines.append('    Lake parameter table created without error.')
    arcpy.AddMessage(loglines[-1])

    # Clean up
    del sr1, point, outRaster1
    arcpy.Delete_management("Lakeslyr")
    arcpy.Delete_management(outshp)
    #arcpy.Delete_management(outfeatures)
    arcpy.Delete_management('in_memory')

    #outRaster.save(outRastername)
    del outRaster                                                               # Added 9/4/2015
    del channelgrid                                                             # Added 9/6/2015
    del projdir, fill2, cellsize, sr2, in_lakes
    return arcpy, channelgrid_arr, outRaster_arr, loglines

def adjust_to_landmask(arcpy, in_raster, LANDMASK, sr2, projdir, inunits):
    """Only to be used if the domain is surrounded by water (LANDMASK=0), and there
    is an issue with WRF Hydro due to 0 releif and 'water' and 'land' having the
    same elevation."""

    """This function is an attempt to alter the elevation data such that values
    above 0 match exactly the land areas from the coarse grid LANDMASK. Values
    that are defined as water in the LANDMASK layer are then set to an elevation
    of 0. Values where LANDMASK = 1 (land) will be at least 0.5m above sea level.
    This is done because the model cannot route water through what it thinks is
    ocean/water based on the LANDMASK."""

    # Start logging
    loglines = ['Step 4(a) initiated...']
    arcpy.AddMessage(loglines[-1])

    mosprj2 = os.path.join(projdir, 'mosaicprj2')

    # Set environments
    arcpy.env.outputCoordinateSystem = sr2
    descData = arcpy.Describe(in_raster)
    extent = descData.Extent                                                    # Overkill?
    arcpy.env.snapRaster = in_raster
    cellsize = descData.children[0].meanCellHeight
    loglines.append('    Cellsize: %s.' %cellsize)
    arcpy.AddMessage(loglines[-1])
    arcpy.env.cellSize = in_raster
    loglines.append('    Environments set.')
    arcpy.AddMessage(loglines[-1])

    # Adjust height to raise all land areas by depending on input z units
    if inunits == 'm':
        above = 0.5
    elif inunits == 'cm':                                                       # trap for fixing cm to m conversion
        above = 50

    # Create landmask on the fine grid
    loglines.append('    Creating resampled LANDMASK layer from GEOGRID file.')
    arcpy.AddMessage(loglines[-1])
    LANDMASK2 = os.path.join(projdir, "LANDMASK")
    arcpy.Resample_management(LANDMASK, LANDMASK2, cellsize, "NEAREST")
    loglines.append('    Finished creating resampled LANDMASK layer from GEOGRID file.')
    arcpy.AddMessage(loglines[-1])

    # Conditional statement for creating elevation increase
    LANDMASK3 = Raster(LANDMASK2)                                               # Redundant?
    RAISEDRASTER = Con(LANDMASK3, above, "", "VALUE = 1")

    # The below equation assumes that any NoData value in the input raster means ocean, and will be given a value of 0.
    RAISEDRASTER2 = RAISEDRASTER + Con(IsNull(in_raster), 0, in_raster, "VALUE = 1")  # Add the raisedraster to the topography layer, only where the landmask isn't NoData
    del RAISEDRASTER, LANDMASK3

    # Clean up
    RAISEDRASTER2.save(os.path.join(projdir, 'mosaicprj2'))
    arcpy.Delete_management(LANDMASK2)
    del RAISEDRASTER2, LANDMASK2

    # Report and return objects
    loglines.append('    Topography corrected to match coarse grid LANDMASK.')
    arcpy.AddMessage(loglines[-1])
    return mosprj2, loglines

def build_groundwater_buckets(build_option=1):
    '''
    5/17/2017: This function will build the groundwater bucket grids and parameter
               tables.

    1) A set of basins must be provided. This is a grid of watershed pixels, in
       which each value on the grid corresponds to a basin.

    Build Options:
        1) Build option 1 will biuld the groundwater buckets from ...

    NOTES:
       * Groundwater buckets are currently resolved on the LSM (coarse/GEOGRID)
         grid. In the future this may change.
       * The ID values for groundwater buckets must be numbered 1...n, and will
         not directly reflect the COMID values of individual basins or pour points.
         Much like the LAKEPARM and LAKEGRID values, a mapping must be made between
         input basins/lakes and the outputs using the shapefiles/parameter tables
         output by these tools.
    '''
    return

def sa_functions(arcpy, rootgrp, basin_mask, mosprj, ovroughrtfac_val, retdeprtfac_val, projdir, in_csv, threshold, inunits, LU_INDEX, cellsize1, cellsize2, routing, in_lakes, lakeIDfield=None):
    """The last major function in the processing chain is to perform the spatial
    analyst functions to hydrologically process the input raster datasets."""

    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
        from arcpy.sa import *
    arcpy.env.overwriteOutput = True

    # Fourth part of the process
    loglines = ['Step 4 initiated...']
    arcpy.AddMessage(loglines[-1])

    # Set Basin mask attribute to boolean from ArcGIS text
    if basin_mask == 'true':
        bsn_msk = True
        loglines.append('    Channelgrid will be masked to basins.')
        arcpy.AddMessage(loglines[-1])
    elif basin_mask == 'false':
        bsn_msk = False
        loglines.append('    Channelgrid will not be masked to basins.')
        arcpy.AddMessage(loglines[-1])

    if routing:
        loglines.append('    Reach-based routing files will be created.')
        arcpy.AddMessage(loglines[-1])
    else:
        loglines.append('    Reach-based routing files will not be created.')
        arcpy.AddMessage(loglines[-1])

    # Set environments
    arcpy.MakeRasterLayer_management(mosprj, 'mosaicprj')
    sr2 = arcpy.Describe('mosaicprj').spatialReference
    arcpy.env.outputCoordinateSystem = sr2

    # Process: Fill DEM
    #fill = Fill(mosprj)                                                        # Sink-filling without z-limit
    fill = Fill(mosprj, z_limit)                                                # Sink-filling using z-limit (global)
    if inunits == 'm':
        fill1 = Float(fill)
        del fill
    elif inunits == 'cm':                                                       # trap for fixing cm to m conversion
        fill1 = fill/float(100)
        del fill

    # Process: Flow Direction
    fdir = FlowDirection(fill1)
    fdir_var = rootgrp.variables['FLOWDIRECTION']
    fdir_arr = arcpy.RasterToNumPyArray(fdir)
    fdir_var[:] = fdir_arr
    loglines.append('    Process: FLOWDIRECTION written to output netCDF.')
    arcpy.AddMessage(loglines[-1])
    del fdir_arr

    # Process: Flow Accumulation (intermediate
    flac = FlowAccumulation(fdir, '#', 'FLOAT')
    flac_var = rootgrp.variables['FLOWACC']
    flac_arr = arcpy.RasterToNumPyArray(flac)
    flac_var[:] = flac_arr
    loglines.append('    Process: FLOWACC written to output netCDF.')
    arcpy.AddMessage(loglines[-1])
    del flac_arr

    # Set NoData Value for topography
    fill2 = Con(IsNull(fill1) == 1, NoDataVal, fill1)
    del fill1
    nodataval = flac.noDataValue
    if nodataval is None:
        nodataval = NoDataVal                                                   # Fix to help this portion run in ArcMap
    arcpy.SetRasterProperties_management(fill2, "ELEVATION", "#", "#", [[1, nodataval]])     # Must set NoData away from -Infinityf
    loglines.append('    Diagnosis: Finished setting rasterproperties for Fill grid.')
    arcpy.AddMessage(loglines[-1])
    fill2_var = rootgrp.variables['TOPOGRAPHY']
    fill2_arr = arcpy.RasterToNumPyArray(fill2)
    fill2_var[:] = fill2_arr
    loglines.append('    Process: TOPOGRAPHY written to output netCDF.')
    arcpy.AddMessage(loglines[-1])
    del fill2_arr

    # Create stream channel raster according to threshold
    strm = SetNull(flac, '1', 'VALUE < %s' % threshold)
    channelgrid = Con(IsNull(strm)==0, 0, NoDataVal)

    # Uncomment these if you want to hardcode the script to use a certain channelgrid raster
    #strm = arcpy.Raster(r'E:\Projects\domains\URG\URGRB_2014_10_30\Experiments\newstrm1')
    #channelgrid = arcpy.Raster(r'E:\Projects\domains\URG\URGRB_2014_10_30\Experiments\newchnlgrd3')

    # Create initial constant raster of -9999
    constraster = CreateConstantRaster(NoDataVal, "INTEGER")
    constraster_arr = arcpy.RasterToNumPyArray(constraster)

    # Create initial constant raster of value retdeprtfac_val
    retdeprtfac_value = retdeprtfac_val
    inraster2 = CreateConstantRaster(retdeprtfac_value, "FLOAT")
    inraster2_var = rootgrp.variables['RETDEPRTFAC']
    inraster2_arr = arcpy.RasterToNumPyArray(inraster2)
    inraster2_var[:] = inraster2_arr
    loglines.append('    Process: RETDEPRTFAC written to output netCDF.')
    arcpy.AddMessage(loglines[-1])
    del retdeprtfac_value, inraster2, inraster2_arr

    # Create initial constant raster of ovroughrtfac_val
    ovroughrtfac_value = ovroughrtfac_val
    inraster3 = CreateConstantRaster(ovroughrtfac_value, "FLOAT")
    inraster3_var = rootgrp.variables['OVROUGHRTFAC']
    inraster3_arr = arcpy.RasterToNumPyArray(inraster3)
    inraster3_var[:] = inraster3_arr
    loglines.append('    Process: OVROUGHRTFAC written to output netCDF.')
    arcpy.AddMessage(loglines[-1])
    del ovroughrtfac_value, inraster3, inraster3_arr

    # Create initial constant raster of LKSATFAC - added 3/30/2017
    inraster4 = CreateConstantRaster(lksatfac_val, "FLOAT")
    inraster4_var = rootgrp.variables['LKSATFAC']
    inraster4_arr = arcpy.RasterToNumPyArray(inraster4)
    inraster4_var[:] = inraster4_arr
    loglines.append('    Process: LKSATFAC written to output netCDF.')
    arcpy.AddMessage(loglines[-1])
    del inraster4, inraster4_arr

    # Process: Stream Order
    order = StreamOrder(strm, fdir)         # Default = "STRAHLER"
    order2 = Con(IsNull(order) == 1, NoDataVal, order)
    order2_var = rootgrp.variables['STREAMORDER']
    order2_arr = arcpy.RasterToNumPyArray(order2)
    order2_var[:] = order2_arr
    loglines.append('    Process: STREAMORDER written to output netCDF.')
    arcpy.AddMessage(loglines[-1])
    del order2_arr

    if routing:
        linkid, loglines = Routing_Table(arcpy, projdir, sr2, strm, fdir, fill2, order2, loglines)
        linkid_var = rootgrp.variables['LINKID']
        linkid_arr = arcpy.RasterToNumPyArray(linkid)
        linkid_var[:] = linkid_arr
        loglines.append('    Process: LINKID written to output netCDF.')
        arcpy.AddMessage(loglines[-1])
        del linkid_arr
    del order, order2

    # Find out if forecast points are chosen, then set mask for them
    if in_csv is not None:
        # Make feature layer from CSV
        loglines.append('    Forecast points provided and basins being delineated.')
        arcpy.AddMessage(loglines[-1])
        frxst_layer = 'frxst_layer'
        sr1 = arcpy.SpatialReference(4326)                                      # GCS_WGS_1984
        arcpy.MakeXYEventLayer_management(in_csv, 'LON', 'LAT', frxst_layer, sr1)
        tolerance = int(float(arcpy.GetRasterProperties_management('mosaicprj', 'CELLSIZEX').getOutput(0)) * walker)
        tolerance1 = int(float(arcpy.GetRasterProperties_management('mosaicprj', 'CELLSIZEX').getOutput(0)))
        frxst_raster = SnapPourPoint(frxst_layer, flac, tolerance1, 'FID')
        frxst_raster2 = Con(IsNull(frxst_raster) == 0, frxst_raster, NoDataVal)
        frxst_raster2_var = rootgrp.variables['frxst_pts']
        frxst_raster2_arr = arcpy.RasterToNumPyArray(frxst_raster2)
        frxst_raster2_var[:] = frxst_raster2_arr
        loglines.append('    Process: frxst_pts written to output netCDF.')
        arcpy.AddMessage(loglines[-1])
        del frxst_raster2_arr

        SnapPour = SnapPourPoint(frxst_layer, flac, tolerance)
        arcpy.Delete_management(frxst_layer)
        arcpy.Delete_management(frxst_raster2)

        # Delineate above points
        outWatershed = Watershed(fdir, SnapPour, 'VALUE')
        watershedraster = os.path.join(projdir, 'watersheds')
        outWatershed.save(watershedraster)
        outWatershed2 = Con(IsNull(outWatershed) == 0, outWatershed, NoDataVal)
        outWatershed2_arr = arcpy.RasterToNumPyArray(outWatershed2)

        # Groundwater basins are for now the same as the forecast basins
        gw_basns_var = rootgrp.variables['basn_msk']
        gw_basns_var[:] = outWatershed2_arr
        loglines.append('    Process: gw_basns written to output netCDF.')
        arcpy.AddMessage(loglines[-1])
        del outWatershed2_arr

        # Process: Resample gw_basns grid to a lower resolution
        gw_basns_geogrid = os.path.join(projdir, 'watersheds2')
        outASCII = os.path.join(projdir, GW_ASCII)
        arcpy.env.cellSize = cellsize1                                          # Set cellsize environment to coarse grid
        arcpy.Resample_management(outWatershed2, gw_basns_geogrid, cellsize1, "MAJORITY")
        arcpy.RasterToASCII_conversion(gw_basns_geogrid, outASCII)
        arcpy.env.cellSize = cellsize2                                          # Set cellsize environment back to fine grid
        loglines.append('    Process: %s completed without error' %GW_ASCII)
        arcpy.AddMessage(loglines[-1])

        # Set mask for future raster output
        if bsn_msk == True:
            channelgrid2 = Con(outWatershed2 >= 0, IsNull(strm), Con(IsNull(strm), 1, -1))  # Converts channelgrid values inside basins to 0, outside to -1
            channelgrid = Con(channelgrid2==1, NoDataVal, channelgrid2)
            del channelgrid2
        del outWatershed2

    else:
        # Handle the change if input forecast points are not provided
        constraster_var = rootgrp.variables['frxst_pts']
        constraster_var[:] = constraster_arr
        loglines.append('    Process: frxst_pts was empty. Constant value raster created.')
        arcpy.AddMessage(loglines[-1])

        constraster_var2 = rootgrp.variables['basn_msk']
        constraster_var2[:] = constraster_arr
        loglines.append('    Process: gw_basns was empty. Constant value raster created.')
        arcpy.AddMessage(loglines[-1])

    if in_lakes is None:
        # Process: Output LAKEGRID
        constraster_var3 = rootgrp.variables['LAKEGRID']
        constraster_var3[:] = constraster_arr
        loglines.append('    Process: LAKEGRID written to output netCDF.')
        arcpy.AddMessage(loglines[-1])

        # Process: Output Channelgrid
        channelgrid_var = rootgrp.variables['CHANNELGRID']
        channelgrid_arr = arcpy.RasterToNumPyArray(channelgrid)
        channelgrid_var[:] = channelgrid_arr
        loglines.append('    Process: CHANNELGRID written to output netCDF.')
        arcpy.AddMessage(loglines[-1])

    else:
        # Alter Channelgrid for reservoirs
        arcpy, channelgrid_arr, outRaster_arr, loglines = add_reservoirs(arcpy, channelgrid, in_lakes, flac, projdir, fill2, cellsize2, sr2, loglines, lakeIDfield)

        # Process: Output LAKEGRID
        outRaster_var = rootgrp.variables['LAKEGRID']
        outRaster_var[:] = outRaster_arr
        loglines.append('    Process: LAKEGRID written to output netCDF.')
        arcpy.AddMessage(loglines[-1])

        # Process: Output Channelgrid
        channelgrid_var = rootgrp.variables['CHANNELGRID']
        channelgrid_var[:] = channelgrid_arr
        loglines.append('    Process: CHANNELGRID written to output netCDF.')
        arcpy.AddMessage(loglines[-1])

    del constraster, fdir, flac, strm, channelgrid, fill2, channelgrid_arr, constraster_arr

    # Process: Resample LU_INDEX grid to a higher resolution
    LU_INDEX2 = os.path.join(projdir, "LU_INDEX")
    arcpy.Resample_management(LU_INDEX, LU_INDEX2, cellsize2, "NEAREST")
    LU_INDEX2_var = rootgrp.variables['landuse']
    LU_INDEX2_arr = arcpy.RasterToNumPyArray(LU_INDEX2)
    LU_INDEX2_var[:] = LU_INDEX2_arr
    loglines.append('    Process: landuse written to output netCDF.')
    arcpy.AddMessage(loglines[-1])
    del LU_INDEX2_arr

    # Clean up
    arcpy.Delete_management('mosaicprj')
    arcpy.Delete_management(mosprj)
    loglines.append('    Step 4 completed without error.')
    arcpy.AddMessage(loglines[-1])
    return rootgrp, loglines

if __name__ == '__main__':
    '''Protect this script against import from another script.'''
    pass
# --- End Functions --- #