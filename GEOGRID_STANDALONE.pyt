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
import shutil
import time
import numpy
import arcpy
import netCDF4                                                                  # Packaged with ArcGIS 10.3 and higher
import re                                                                       # Added 10/11/2016 for string matching in netCDF global attributes
import importlib
import copy                                                                     # Added 11/19/2019 to allow copying of class objects
from distutils.version import StrictVersion, LooseVersion

# Test current version of Python's ability to reload a module
# https://stackoverflow.com/questions/961162/reloading-module-giving-nameerror-name-reload-is-not-defined
try:
    reload  # Python 2.7
except NameError:
    try:
        from importlib import reload    # Python 3.4+
    except ImportError:
        from imp import reload          # Python 3.0 - 3.3

### Test if the packaging module is installed:
##if sys.version_info[0] >= 3:
##    from packaging import version

# Specify import path and append to PATH
configfile = '~/wrf_hydro_functions.py'
sys.path.insert(1,os.path.dirname(os.path.expanduser(configfile)))
import wrf_hydro_functions as wrfh                                              # Function script packaged with this toolbox
reload(wrfh)                                                                    # Re-load the function script in case of script changes
# --- End Import Modules --- #

# --- Module Configurations --- #
arcpy.env.overwriteOutput = True                                                # Allow overwriting of outputs
arcpy.env.parallelProcessingFactor = "100%"                                     # Use all of the cores on the machine for tools that respect this environment.
if arcpy.CheckExtension("Spatial") == "Available":                              # Check if the Spatial Analyst extention is available
    arcpy.CheckOutExtension("Spatial")                                          # Check out Spatial Analyst extention license
    from arcpy.sa import *                                                      # Add Spatial Analyst functions into namespace
# --- End Module Configurations --- #

# Find out if we are in 32-bit or 64-bit
if sys.maxsize > 2**32:
    bit64 = True
else:
    bit64 = False

# --- Globals --- #
#outNCType = 'NETCDF3_64BIT'                                                     # Set output netCDF format for spatial metdata files. This was the default before 7/31/2018
outNCType = 'NETCDF4_CLASSIC'                                                   # Define the output netCDF version for RouteLink.nc and LAKEPARM.nc

# Processing Notes to insert into netCDF global attributes
# Processing notes for Spatial Metdata files
processing_notes_SM = '''Created: %s''' %time.ctime()

# Processing notes for the FULLDOM (Routing Grid) file
processing_notesFD = '''Created: %s''' %time.ctime()

# List of all possible routing-stack files to keep between the working directory and output .zip files
nclist = [wrfh.LDASFile,
            wrfh.FullDom,
            'gw_basns.nc',
            wrfh.GW_ASCII,
            'gw_basns_geogrid.prj',
            wrfh.RT_nc,
            'Route_Link.csv',
            wrfh.LK_tbl,
            wrfh.LK_nc,
            'streams.shp', 'streams.shx', 'streams.shp.xml', 'streams.sbx', 'streams.sbn', 'streams.prj', 'streams.dbf',
            'lakes.shp', 'lakes.shx', 'lakes.shp.xml', 'lakes.sbx', 'lakes.sbn', 'lakes.prj', 'lakes.dbf',
            wrfh.GW_nc,
            wrfh.GWGRID_nc,
            wrfh.GW_TBL,
            'Lake_Problems.csv',
            'Old_New_LakeComIDs.csv',
            'Lake_Link_Types.csv',
            'Tossed_Lake_Link_Types.csv',
            'Lake_Preprocssing_Info.txt',
            'Lakes_with_minimum_depth.csv']

# Groundwater input options
GW_with_Stack = True                                                            # Switch for building default groundwater inputs with any routing stack
defaultGWmethod = 'FullDom LINKID local basins'                                 # Provide the default groundwater basin generation method. Options ['FullDom basn_msk variable', 'FullDom LINKID local basins', 'Polygon Shapefile or Feature Class']
in_GWPolys = None                                                               # The polygon shapefile to use if defaultGWmethod == 'Polygon Shapefile or Feature Class'

# Methods test switches
coordMethod1 = True                                                             # Interpolate GEOGRID latitude and longitude coordinate arrays
coordMethod2 = False                                                            # Transform coordinate pairs at each grid cell from projected to geocentric

# Spatial Metadata

# Adding global attributes from the GEOGRID file to the Spatial Metadata files and
# to the Fulldom file can be useful, because then those files can be used as input
# to other GIS tools, without requiring a WPS GEOGRID file.
Globals_to_Add = ['MAP_PROJ', 'corner_lats', 'corner_lons', 'TRUELAT1', 'TRUELAT2', 'STAND_LON', 'POLE_LAT', 'POLE_LON', 'MOAD_CEN_LAT', 'CEN_LAT']
# --- End Globals --- #

# --- Toolbox Classes --- #
class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "WRFHydro_GIS_Pre-Processor"
        self.alias = ""
        self.description = "This is a standalone ArcGIS geoprocessing toolbox for WRF-Hydro."

        # List of tool classes associated with this toolbox
        self.tools = [ProcessGeogridFile,
                        ExportGrid,
                        ExamineOutputs,
                        ExportPRJ,
                        GenerateLatLon,
                        SpatialMetadataFile,
                        DomainShapefile,
                        Reach_Based_Routing_Addition,
                        Lake_Parameter_Addition,
                        GWBUCKPARM]

class ProcessGeogridFile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Process GEOGRID File"
        self.description = "This tool takes an input WRF GEOGRID file in NetCDF format" + \
                           " and uses the HGT_M grid and an input high-resolution elevation grid" + \
                           "to produce a high-resolution hydrologically processed output."
        self.canRunInBackground = True                                          #self.canRunInBackground = False
        self.category = "Processing"

    def getParameterInfo(self):
        """Define parameter definitions"""

        in_nc = arcpy.Parameter(
            displayName="Input GEOGRID File",
            name="in_nc",
            datatype="File",
            parameterType="Required",
            direction="Input")

        in_csv = arcpy.Parameter(
            displayName="Forecast Points (CSV)",
            name="in_csv",
            datatype="File",
            parameterType="Optional",
            direction="Input")
        # To define a file filter that includes .csv and .txt extensions, set the filter list to a list of file extension names
        in_csv.filter.list = ['csv']

        basin_mask = arcpy.Parameter(
            displayName="Mask CHANNELGRID variable to forecast basins?",
            name="basin_mask",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        basin_mask.value = False

        RB_routing = arcpy.Parameter(
            displayName="Create reach-based routing (RouteLink) files?",
            name="RB_routing",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        RB_routing.value = False

        Lake_routing = arcpy.Parameter(
            displayName="Create lake parameter (LAKEPARM) file?",
            name="Lake_routing",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        Lake_routing.value = False

        in_reservoirs = arcpy.Parameter(
            displayName="Reservoirs Shapefile or Feature Class",
            name="in_reservoirs",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Input")

        in_raster = arcpy.Parameter(
            displayName="Input Elevation Raster",
            name="in_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        cellsize = arcpy.Parameter(
            displayName="Regridding (nest) Factor",
            name="cellsize",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        cellsize.value = 10

        threshold = arcpy.Parameter(
            displayName="Number of routing grid cells to define stream",
            name="threshold",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        threshold.value = 200

        ovroughrtfac_val = arcpy.Parameter(
            displayName="OVROUGHRTFAC Value",
            name="ovroughrtfac_val",
            datatype="Any Value",
            parameterType="Required",
            direction="Input")
        ovroughrtfac_val.value = 1.0
        ovroughrtfac_val.category = "Parameter Values"

        retdeprtfac_val = arcpy.Parameter(
            displayName="RETDEPRTFAC Value",
            name="retdeprtfac_val",
            datatype="Any Value",
            parameterType="Required",
            direction="Input")
        retdeprtfac_val.value = 1.0
        retdeprtfac_val.category = "Parameter Values"

        out_zip = arcpy.Parameter(
            displayName="Output ZIP File",
            name="out_zip",
            datatype="File",
            parameterType="Required",
            direction="Output")
        out_zip.value = 'WRF_Hydro_routing_grids.zip'

        parameters = [in_nc, in_csv, basin_mask, RB_routing, Lake_routing, in_reservoirs, in_raster, cellsize, threshold, ovroughrtfac_val, retdeprtfac_val, out_zip]   #, in_LakeIDField
        return parameters

    def isLicensed(self):
        """Allow the tool to execute, only if the ArcGIS Spatial Analyst extension
        is available."""
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False                                                        # tool cannot be executed
        return True                                                             # tool can be executed

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        # Only activate masking if a CSV file has been input
        if not parameters[1].altered:
          parameters[2].enabled = False
        else:
          parameters[2].enabled = True

        # Only activate Lake input parameter and lake ID field parameter if requesting LAKEPARM file
        if parameters[4].value == True:
            parameters[5].enabled = True
        else:
            parameters[5].enabled = False
            parameters[5].value = ''

        ##        # Populate lake ID field combo box with list of Integer type fields from reservoirs shapefile
        ##        if parameters[5].altered:
        ##            in_lakes_file = parameters[5].valueAsText
        ##            FieldNames = [field.name for field in arcpy.ListFields(in_lakes_file, "", "Integer")]
        ##            parameters[6].filter.list = FieldNames

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""

        # Ensure that if a LAKEPARM table is requested, a lakes feature class is provided
        if parameters[4].value == True and parameters[5].value is None:
            parameters[5].setErrorMessage("Must specify Reservoirs Shapefile or Feature Class.")
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        reload(wrfh)                                                            # Reload in case code changes have been made
        tic = time.time()                                                       # Initiate timer
        isError = False                                                         # Starting Error condition (no errors)

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        in_csv = parameters[1].valueAsText
        basin_mask = parameters[2].value
        routing = parameters[3].value
        Lake_routing = parameters[4].value
        in_reservoir = parameters[5].valueAsText
        in_raster = parameters[6].valueAsText
        cellsize = parameters[7].value
        threshold = parameters[8].value
        ovroughrtfac_val = parameters[9].value
        retdeprtfac_val = parameters[10].value
        out_zip = parameters[11].valueAsText

        # Prepare output log file
        outtable = os.path.join(os.path.dirname(out_zip), os.path.basename(out_zip) + '.log')
        tee = wrfh.TeeNoFile(outtable, 'w')
        wrfh.printMessages(arcpy, ['Begining processing on {0}'.format(time.ctime())])
        wrfh.printMessages(arcpy, ['64-bit: {0}'.format(bit64)])
        wrfh.printMessages(arcpy, ['Input parameters:'])
        for param in parameters:
            wrfh.printMessages(arcpy, ['    Parameter: {0}: {1}'.format(param.displayName, param.valueAsText)])

        # Interpret the input for reservoir routing
        if Lake_routing:
            in_lakes = in_reservoir
        else:
            in_lakes = None
        wrfh.printMessages(arcpy, ['{0}'.format(in_lakes)])

        # Create scratch directory for temporary outputs
        projdir = os.path.join(os.path.dirname(out_zip), 'scratchdir')
        if os.path.exists(projdir):
            shutil.rmtree(projdir)
        os.makedirs(projdir)
        arcpy.env.overwriteOutput = True
        arcpy.env.workspace = projdir
        arcpy.env.scratchWorkspace = projdir

        '''Pre-defining the variables and populating variable attributes is
        a much faster strategry than creating and populating each variable
        sequentially, especially for netCDF3 versions. Also, unsigned integer
        types are only allowed in NETCDF4.'''
        # List of variables to create [<varname>, <vardtype>, <long_name>]
        varList2D = [['CHANNELGRID', 'i4', ''],
                    ['FLOWDIRECTION', 'i2', ''],
                    ['FLOWACC', 'i4', ''],
                    ['TOPOGRAPHY', 'f4', ''],
                    ['RETDEPRTFAC', 'f4', ''],
                    ['OVROUGHRTFAC', 'f4', ''],
                    ['STREAMORDER', 'i1', ''],
                    ['frxst_pts', 'i4', ''],
                    ['basn_msk', 'i4', ''],
                    ['LAKEGRID', 'i4', ''],
                    ['landuse', 'f4', ''],
                    ['LKSATFAC', 'f4', '']]
        # Add variables depending on the input options
        if routing:
            varList2D.append(['LINKID', 'i4', ''])

        # Step 1 - Georeference geogrid file
        rootgrp = netCDF4.Dataset(in_nc, 'r')                                   # Establish an object for reading the input NetCDF files
        globalAtts = rootgrp.__dict__                                           # Read all global attributes into a dictionary
        if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
            rootgrp.set_auto_mask(False)                                        # Change masked arrays to old default (numpy arrays always returned)
        coarse_grid = wrfh.WRF_Hydro_Grid(arcpy, rootgrp)                       # Instantiate a grid object
        fine_grid = copy.copy(coarse_grid)                                      # Copy the grid object for modification
        fine_grid.regrid(cellsize)                                              # Regrid to the coarse grid
        GeoTransform1 = coarse_grid.GeoTransformStr()
        wrfh.printMessages(arcpy, ['    The GEOGRID File resolution is {0}sm'.format(str(coarse_grid.DX))])
        wrfh.printMessages(arcpy, ['    Proj4: {0}'.format(coarse_grid.proj4)])     # Print Proj.4 string to screen
        wrfh.printMessages(arcpy, ['    GeoTransform: {0}'.format(GeoTransform1)])  # Print affine transformation to screen.
        wrfh.printMessages(arcpy, ['    Created projection definition from input NetCDF GEOGRID file.'])

        # Build output raster from numpy array of the GEOGRID variable requested. This will be used as a template later on
        LU_INDEX = coarse_grid.numpy_to_Raster(arcpy, wrfh.flip_grid(rootgrp.variables['LU_INDEX'][0]))
        out_nc1 = os.path.join(projdir, wrfh.LDASFile)
        rootgrp1 = netCDF4.Dataset(out_nc1, 'w', format=outNCType)              # wrfh.outNCType)
        rootgrp1 = wrfh.create_CF_NetCDF(arcpy, coarse_grid, rootgrp1, addLatLon=False, notes=processing_notes_SM)
        rootgrp1.proj4 = coarse_grid.proj4                                      # Add proj.4 string as a global attribute
        for item in Globals_to_Add + ['DX', 'DY']:
            if item in globalAtts:
                rootgrp1.setncattr(item, globalAtts[item])
        rootgrp1.close()
        del rootgrp1

        # Step 2 - Create high resolution topography layers
        mosprj = wrfh.create_high_res_topogaphy(arcpy, in_raster, LU_INDEX, cellsize, fine_grid.proj, projdir)

        # Build latitude and longitude arrays for Fulldom_hires netCDF file
        if coordMethod1:
            latArr, lonArr = wrfh.coordMethod1(arcpy, coarse_grid, fine_grid, rootgrp, projdir=projdir)
        elif coordMethod2:
            latArr, lonArr = wrfh.coordMethod2(arcpy, fine_grid)

        # Create FULLDOM file
        out_nc2 = os.path.join(projdir, wrfh.FullDom)
        rootgrp2 = netCDF4.Dataset(out_nc2, 'w', format=outNCType)              # wrfh.outNCType)
        rootgrp2 = wrfh.create_CF_NetCDF(arcpy, fine_grid, rootgrp2, addLatLon=True,
            notes=processing_notesFD, addVars=varList2D, latArr=latArr, lonArr=lonArr)
        del latArr, lonArr

        # Add some global attribute metadata to the Fulldom file, including relevant WPS attributes for defining the model coordinate system
        rootgrp2.geogrid_used = in_nc                                           # Paste path of geogrid file to the Fulldom global attributes
        rootgrp2.proj4 = fine_grid.proj4                                        # Add proj.4 string as a global attribute
        rootgrp2.DX = fine_grid.DX                                              # Add X resolution as a global attribute
        rootgrp2.DY = fine_grid.DY                                              # Add Y resolution as a global attribute
        for item in Globals_to_Add:
            if item in globalAtts:
                rootgrp2.setncattr(item, globalAtts[item])
        rootgrp.close()
        del item, globalAtts

        # Process: Resample LU_INDEX grid to the routing grid resolution
        LU_INDEX2 = os.path.join("in_memory", "LU_INDEX")
        LU_INDEX_fine = fine_grid.project_to_model_grid(arcpy, LU_INDEX, LU_INDEX2, resampling="NEAREST") # Regrid from GEOGRID resolution to routing grid resolution
        LU_INDEX2_arr = arcpy.RasterToNumPyArray(LU_INDEX2)
        rootgrp2.variables['landuse'][:] = LU_INDEX2_arr
        wrfh.printMessages(arcpy, ['    Process: landuse written to output netCDF.'])
        arcpy.Delete_management(LU_INDEX)                                       # Added 4/19/2017 to allow zipws to complete
        arcpy.Delete_management(LU_INDEX2)
        del LU_INDEX2, LU_INDEX2_arr, LU_INDEX_fine

        try:
            # Step 4 - Hyrdo processing functions
            rootgrp2 = wrfh.sa_functions(arcpy, rootgrp2, basin_mask, mosprj,
                ovroughrtfac_val, retdeprtfac_val, projdir, in_csv, threshold,
                routing, in_lakes=in_lakes)
            rootgrp2.close()
            del rootgrp2
        except Exception as e:
            wrfh.printMessages(arcpy, ['Exception: {0}'.format(e)])
            rootgrp2.close()
            isError = True

        if not isError:
            try:
                if GW_with_Stack:
                    # Build groundwater files
                    GWBasns, GWBasns_arr = wrfh.build_GW_Basin_Raster(arcpy, out_nc2,
                        projdir, defaultGWmethod, fine_grid, in_Polys=in_GWPolys)
                    wrfh.build_GW_buckets(arcpy, projdir, GWBasns, GWBasns_arr, coarse_grid, Grid=True)
                    del GWBasns, GWBasns_arr
                del GeoTransform1, LU_INDEX
            except Exception as e:
                wrfh.printMessages(arcpy, ['Exception: {0}'.format(e)])
                isError = True

        # Clean up and give finishing message
        if isError:
            wrfh.printMessages(arcpy, ['Error encountered after {0} seconds.'.format(time.time()-tic)])
            arcpy.env.workspace = projdir
            for infile in arcpy.ListDatasets():
                arcpy.Delete_management(infile)
            arcpy.Delete_management(projdir)
            arcpy.AddError("ERROR")
            raise SystemExit
        else:
            # zip the folder
            zipper = wrfh.zipUpFolder(arcpy, projdir, out_zip, nclist)
            wrfh.printMessages(arcpy, ['Completed without error in {0} seconds.'.format(time.time()-tic)])
            arcpy.env.workspace = projdir
            for infile in arcpy.ListDatasets():
                arcpy.Delete_management(infile)
            arcpy.Delete_management(projdir)
        tee.close()
        del tee
        return

class ExportGrid(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Export grid from GEOGRID file"
        self.description = "This tool takes an input WRF Geogrid file in NetCDF format" + \
                           " and uses the specified variable's grid to produce a raster."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_nc = arcpy.Parameter(
            displayName="Input GEOGRID File",
            name="in_nc",
            datatype="File",
            parameterType="Required",
            direction="Input")

        # Fourth parameter
        var_name = arcpy.Parameter(
            displayName="Variable Name",
            name="var_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        var_name.filter.type = "ValueList"

        # Third parameter
        out_raster = arcpy.Parameter(
            displayName="Output Raster",
            name="out_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        parameters = [in_nc, var_name, out_raster]
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        if parameters[0].altered:
            in_nc_file = parameters[0].valueAsText

            # Establish an object for reading the input NetCDF file
            rootgrp = netCDF4.Dataset(in_nc_file, 'r')

            # Loop through global variables in NetCDF file to gather projection information
            ncMassgridNames = [varName for varName,ncvar in rootgrp.variables.items() if 'south_north' and 'west_east' in ncvar.dimensions]
            parameters[1].filter.list = ncMassgridNames
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        reload(wrfh)                                                            # Reload in case code changes have been made

        # Set environments
        arcpy.env.overwriteOutput = True

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        Variable = parameters[1].valueAsText
        out_raster = parameters[2].valueAsText

        # Step 1 - Georeference geogrid file
        rootgrp = netCDF4.Dataset(in_nc, 'r')                                   # Establish an object for reading the input NetCDF files
        if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
            rootgrp.set_auto_mask(False)                                        # Change masked arrays to old default (numpy arrays always returned)
        coarse_grid = wrfh.WRF_Hydro_Grid(arcpy, rootgrp)                       # Instantiate a grid object
        nc_raster = coarse_grid.numpy_to_Raster(arcpy, wrfh.flip_grid(rootgrp.variables[Variable][0]))

        # Set environments and save
        nc_raster.save(out_raster)
        del nc_raster

        arcpy.AddMessage('    Process completed without error.')
        arcpy.AddMessage('    Output Raster: %s' %out_raster)
        return

class ExamineOutputs(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Examine Outputs of GIS Preprocessor"
        self.description = "This tool takes the output zip file from the ProcessGeogrid script" + \
                           "and creates a raster from each 2D variable in each NetCDF file." + \
                           "" + \
                           "The Input should be a .zip file that was created using the WRF Hydro pre-" + \
                           "processing tools.  The Output Folder parameter should be set to a non-existent " +\
                           "folder location.  The tool will create the folder which will contain the results."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_zip = arcpy.Parameter(
            displayName="Input ZIP File",
            name="in_zip",
            datatype="File",
            parameterType="Required",
            direction="Input")

        # Output parameter
        out_folder = arcpy.Parameter(
            displayName="Output Folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Output")
        out_folder.defaultEnvironmentName = "workspace"

        parameters = [in_zip, out_folder]
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        reload(wrfh)                                                            # Reload in case code changes have been made

        # Initiate input and output parameters
        in_zip = parameters[0].valueAsText
        out_folder = parameters[1].valueAsText

        # Create output directory
        os.mkdir(out_folder)

        # Use wrfh to perform process
        out_sfolder = wrfh.Examine_Outputs(arcpy, in_zip, out_folder, skipfiles=[])
        return

class ExportPRJ(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Export ESRI projection file (PRJ) from GEOGRID file"
        self.description = "This tool takes an input WRF Geogrid file in NetCDF format" + \
                           " and uses the specified variable's projection parameters" + \
                           " to produce a projection file."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_nc = arcpy.Parameter(
            displayName="Input GEOGRID File",
            name="in_nc",
            datatype="File",
            parameterType="Required",
            direction="Input")

        # Second parameter
        out_prj = arcpy.Parameter(
            displayName="Output Projection File (.prj)",
            name="out_prj",
            datatype="DEPrjFile",
            parameterType="Required",
            direction="Output")

        parameters = [in_nc, out_prj]
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        reload(wrfh)                                                            # Reload in case code changes have been made

        # Set environments
        arcpy.env.overwriteOutput = True

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        out_prj = parameters[1].valueAsText
        Variable = 'HGT_M'

        # Step 1 - Georeference geogrid file
        rootgrp = netCDF4.Dataset(in_nc, 'r')                                   # Establish an object for reading the input NetCDF files
        coarse_grid = wrfh.WRF_Hydro_Grid(arcpy, rootgrp)                       # Instantiate a grid object

        # Optional save .prj file
        with open(out_prj, 'w') as prj_file:
             prj_file.write(coarse_grid.WKT)
        rootgrp.close()
        del rootgrp, coarse_grid
        arcpy.AddMessage('    Process completed without error.')
        arcpy.AddMessage('    Output Projection File: %s' %out_prj)
        return

class GenerateLatLon(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Generate Latitude and Longitude Rasters"
        self.description = "This tool takes an input raster (most likely produced" + \
                           " using the ExportGrid tool) and uses that grid to produce" + \
                           " latitude and longitude ESRI GRID rasters."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_raster = arcpy.Parameter(
            displayName="Input Raster",
            name="in_raster",
            datatype="Raster Dataset",
            parameterType="Required",
            direction="Input")
        # datatype="DERasterDataset" keyword can be used at 10.1 SP1

        # Input parameter
        in_method = arcpy.Parameter(
            displayName="Method for deriving geocentric coordinates",
            name="in_method",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        in_method.filter.type = "ValueList"
        in_method.filter.list = ['1: Interpolate from re-projected geocentric coordinates (faster, more error)',
                                 '2: Direct conversion of each cell centroid coordinate pairs (slower, more accurate)]']

        # Output parameter
        out_folder = arcpy.Parameter(
            displayName="Output Folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")
        out_folder.defaultEnvironmentName = "workspace"

        parameters = [in_raster, in_method, out_folder]
        return parameters

    def isLicensed(self):
        """Allow the tool to execute, only if the ArcGIS Spatial Analyst extension
        is available."""
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False                                                            # tool cannot be executed
        return True                                                                 # tool can be executed

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        # Set environments
        reload(wrfh)                                                            # Reload in case code changes have been made
        arcpy.env.overwriteOutput = True

        # Gather all necessary parameters
        inraster = parameters[0].valueAsText
        method = parameters[1].valueAsText
        projdir = parameters[2].valueAsText
        wrfh.printMessages(arcpy, ['Input Raster Dataset: {0}'.format(inraster)])
        wrfh.printMessages(arcpy, ['Method selected: {0}'.format(method)])
        wrfh.printMessages(arcpy, ['Directory to be used for outputs: {0}'.format(projdir)])

        in_raster = arcpy.Raster(inraster)
        descData = arcpy.Describe(in_raster)
        sr_in = descData.spatialReference
        extent = descData.Extent

        # Method 1 - Create latitude and longitude rasters
        if method.startswith('1'):
            xout2, yout2, xmap, ymap = wrfh.create_lat_lon_rasters(arcpy, projdir, in_raster, wrfh.wkt_text)
            arcpy.Delete_management(xmap)
            arcpy.Delete_management(ymap)
            del xmap, ymap

        # # Method 2: Transform each point from projected coordinates to geocentric coordinates
        elif method.startswith('2'):
            wrfh.printMessages(arcpy, ['  Deriving geocentric coordinates from direct transformation of input coordinates.'])
            DY = descData.meanCellHeight
            DX = descData.meanCellWidth
            wgs84_proj = arcpy.SpatialReference()                               # Projectipm coordinate system is specified by wkt_text in globals
            wgs84_proj.loadFromString(wrfh.wkt_text)                            # Load the Sphere datum CRS using WKT
            xmap, ymap = wrfh.getxy(arcpy, in_raster, projdir)
            xmap_arr = arcpy.RasterToNumPyArray(xmap)                            # Read channel grid array
            ymap_arr = arcpy.RasterToNumPyArray(ymap)                            # Read channel grid array
            lonArr2, latArr2 = wrfh.ReprojectCoords(arcpy, xmap_arr, ymap_arr, sr_in, wgs84_proj)  # Transform coordinate arrays
            arcpy.Delete_management(xmap)
            arcpy.Delete_management(ymap)
            xout2 = arcpy.NumPyArrayToRaster(lonArr2, extent.lowerLeft, DX, DY)
            yout2 = arcpy.NumPyArrayToRaster(latArr2, extent.lowerLeft, DX, DY)
            del xmap, ymap, wgs84_proj, latArr2, lonArr2, DY, DX, xmap_arr, ymap_arr

        # Write to disk
        xout2.save(os.path.join(projdir, 'longitude'))
        yout2.save(os.path.join(projdir, 'latitude'))
        wrfh.printMessages(arcpy, ['    Process completed without error.'])
        del sr_in, in_raster
        return

class SpatialMetadataFile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Build Spatial Metadata File"
        self.description = "This tool takes an input GEOGRID and uses that grid " + \
                           " information to produce spatial metadata files against " + \
                           " the multiple resolutions of WRF Hydro output files."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_nc = arcpy.Parameter(
            displayName="Input GEOGRID File",
            name="in_nc",
            datatype="File",
            parameterType="Optional",
            direction="Input")

        format_out = arcpy.Parameter(
            displayName="Output Grid Resolution",
            name="format_out",
            datatype="String",
            parameterType="Required",
            direction="Input")
        # Set a value list for input raster catalog
        format_out.filter.type = "ValueList"
        format_out.filter.list = ["LDASOUT", "RTOUT"]
        format_out.value = "LDASOUT"

        # Input parameter
        factor = arcpy.Parameter(
            displayName="Regridding Factor",
            name="factor",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")
        factor.value = 1

        # Output parameter
        out_nc = arcpy.Parameter(
            displayName="Output netCDF File",
            name="out_nc",
            datatype="File",
            parameterType="Required",
            direction="Output")

        # Output parameter
        latlon_vars = arcpy.Parameter(
            displayName="Include LATITUDE and LONGITUDE 2D variables?",
            name="latlon_vars",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        latlon_vars.value = False

        parameters = [in_nc, format_out, factor, out_nc, latlon_vars]
        return parameters

    def isLicensed(self):
        """Allow the tool to execute, only if the ArcGIS Spatial Analyst extension
        is available."""
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False                                                            # tool cannot be executed
        return True                                                                 # tool can be executed

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        # Only activate Regridding Factor parameter if requesting CHRTOUT resolution
        if parameters[1].altered:
            if parameters[1].value == "LDASOUT":
                parameters[2].enabled = False
                parameters[2].value = 1
            elif parameters[1].value == "RTOUT":
                parameters[2].enabled = True
            else:
                parameters[2].enabled = False
                parameters[2].value = 1
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        reload(wrfh)                                                            # Reload in case code changes have been made

        # Coordinate Attribute Conventions from: http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/tutorial/CoordinateAttributes.html
        # Esri's use of CF conventions discussed here: http://pro.arcgis.com/en/pro-app/help/data/multidimensional/spatial-reference-for-netcdf-data.htm

        # Set environments
        tic1 = time.time()
        arcpy.env.overwriteOutput = True

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        format_out = parameters[1].valueAsText
        factor = float(parameters[2].valueAsText)
        out_nc = parameters[3].valueAsText
        latlon_vars = parameters[4].value
        projdir = os.path.dirname(out_nc)

        # Prepare output log file
        outtable = os.path.join(projdir, os.path.basename(out_nc) + '.log')
        tee = wrfh.TeeNoFile(outtable, 'w')

        # Print informational messages
        wrfh.printMessages(arcpy, ['Begining processing on {0}'.format(time.ctime())])
        wrfh.printMessages(arcpy, ['64-bit: {0}'.format(bit64)])
        wrfh.printMessages(arcpy, ['Input parameters:'])
        for param in parameters:
            wrfh.printMessages(arcpy, ['    Parameter: {0}: {1}'.format(param.displayName, param.valueAsText)])
        wrfh.printMessages(arcpy, ['Directory to be used for outputs: {0}'.format(projdir)])

        # Georeference geogrid file
        rootgrp = netCDF4.Dataset(in_nc, 'r')                                   # Establish an object for reading the input NetCDF files
        globalAtts = rootgrp.__dict__                                           # Read all global attributes into a dictionary
        if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
            rootgrp.set_auto_mask(False)                                        # Change masked arrays to old default (numpy arrays always returned)
        coarse_grid = wrfh.WRF_Hydro_Grid(arcpy, rootgrp)                       # Instantiate a grid object

        # Record GEOGRID MAP_PROJ attribute
        wrfh.printMessages(arcpy, ['    Map Projection of GEOGRID: {0}'.format(wrfh.projdict[coarse_grid.map_pro])])
        wrfh.printMessages(arcpy, ['    Esri PE String: {0}'.format(coarse_grid.WKT)])

        # Create high resolution raster for RTOUT output using Spatial Analyst CreateConstantRaster function
        if format_out == "RTOUT":
            grid_obj = copy.copy(coarse_grid)                                      # Copy the grid object for modification
            grid_obj.regrid(factor)                                                # Regrid to the coarse grid
            wrfh.printMessages(arcpy, ['    New Resolution: {0} {1}'.format(grid_obj.DX, grid_obj.DY)])
        else:
            grid_obj = coarse_grid

        if latlon_vars:
            if coordMethod1:
                if format_out == "RTOUT":
                    latArr, lonArr = wrfh.coordMethod1(arcpy, coarse_grid, grid_obj, rootgrp, projdir=projdir)
                else:
                    latArr = wrfh.flip_grid(rootgrp.variables['XLAT_M'][0])     # Extract array of GEOGRID latitude values
                    lonArr = wrfh.flip_grid(rootgrp.variables['XLONG_M'][0])    # Extract array of GEOGRID longitude values
            elif coordMethod2:
                latArr, lonArr = wrfh.coordMethod2(arcpy, fine_grid)
        else:
            latArr = lonArr = None

        # Create the netCDF file with spatial metadata
        rootgrp = netCDF4.Dataset(out_nc, 'w', format=outNCType)
        rootgrp = wrfh.create_CF_NetCDF(arcpy, grid_obj, rootgrp, addLatLon=latlon_vars,
            notes=processing_notes_SM, latArr=latArr, lonArr=lonArr)
        for item in Globals_to_Add + ['DX', 'DY']:
            if item in globalAtts:
                rootgrp.setncattr(item, globalAtts[item])
        rootgrp.close()
        wrfh.printMessages(arcpy, ['Completed without error in {0} seconds.'.format(time.time()-tic1)])
        tee.close()
        del tee, latArr, lonArr
        return

class DomainShapefile(object):

    """This function taks the WRF GEOGRID file and outputs a polygon shapefile using the bounding coordinate information."""

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Domain Boundary Shapefile"
        self.description = "This tool takes an WRF Geogrid file and creates a single" + \
                           " polygon shapefile that makes up the boundary of the domain" + \
                           " of the M-grid (HGT_M, for example)."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_nc = arcpy.Parameter(
            displayName="Input GEOGRID File",
            name="in_nc",
            datatype="File",
            parameterType="Required",
            direction="Input")

        # Output parameter
        out_shp = arcpy.Parameter(
            displayName="Output Shapefile",
            name="out_shp",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        parameters = [in_nc, out_shp]
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False                                                            # tool cannot be executed
        return True                                                                 # tool can be executed

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        reload(wrfh)                                                            # Reload in case code changes have been made

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        out_shp = parameters[1].valueAsText

        # Create scratch directory for temporary outputs
        projdir = os.path.dirname(in_nc)
        arcpy.env.overwriteOutput = True
        arcpy.env.workspace = projdir
        arcpy.env.scratchWorkspace = projdir

        # Georeference geogrid file
        rootgrp = netCDF4.Dataset(in_nc, 'r')                                   # Establish an object for reading the input NetCDF files
        coarse_grid = wrfh.WRF_Hydro_Grid(arcpy, rootgrp)                       # Instantiate a grid object
        ConstRaster = coarse_grid.numpy_to_Raster(arcpy, numpy.ones((coarse_grid.nrows, coarse_grid.ncols)))
        coarse_grid.domain_shapefile(arcpy, ConstRaster, out_shp)
        return

class Reach_Based_Routing_Addition(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Add reach-based routing"
        self.description = "This tool takes an input Routing Stack .ZIP file," + \
                           " typically created by the Process GEOGRID FILE tool," + \
                           " and adds in the reach-based routing table and grids."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_zip = arcpy.Parameter(
            displayName="Routing Stack ZIP File",
            name="in_zip",
            datatype="File",
            parameterType="Required",
            direction="Input")
        in_zip.filter.list = ['zip']

        # Output parameter
        out_zip = arcpy.Parameter(
            displayName="Output ZIP File",
            name="out_zip",
            datatype="File",
            parameterType="Required",
            direction="Output")

        parameters = [in_zip, out_zip]
        return parameters

    def isLicensed(self):
        """Allow the tool to execute, only if the ArcGIS Spatial Analyst extension
        is available."""
        try:
            if not arcpy.CheckExtension("Spatial") == "Available":
                raise Exception
        except Exception:
            return False                                                            # tool cannot be executed
        return True                                                                 # tool can be executed

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        reload(wrfh)                                                            # Reload in case code changes have been made

        # Gather all necessary parameters
        in_zip = parameters[0].valueAsText
        out_zip = parameters[1].valueAsText

        # Create scratch directory for temporary outputs
        projdir = arcpy.CreateScratchName("temp", data_type="Folder", workspace=arcpy.env.scratchFolder)
        os.mkdir(projdir)

        # Set environments
        arcpy.env.overwriteOutput = True
        arcpy.env.workspace = projdir
        arcpy.env.scratchWorkspace = projdir

        # Unzip
        wrfh.printMessages(arcpy, ['Begining processing on {0}'.format(time.ctime())])
        wrfh.printMessages(arcpy, ['Beginning to extract WRF routing grids...'])

        # Unzip to a known location (make sure no other nc files live here)
        FullDom = wrfh.FullDom
        GWBasins = wrfh.GWGRID_nc
        out_sfolder = wrfh.Examine_Outputs(arcpy, in_zip, projdir, skipfiles=[FullDom, GWBasins])

        # Add a check for lakes in the routing stack. Terminate if lakes are found
        if os.path.isfile(os.path.join(projdir, wrfh.LK_nc)):
            # This means that LAKEPARM is in the input routing stack. Terminate
            msg ='The input routing stack already has lakes in it. it is not a good idea ' \
                'to add reaches because network connectivity may be compromised. Exiting...'
            messages.addErrorMessage(msg)
            raise SystemExit

        # Prepare other rasters for the Routing Table function
        fdir = arcpy.Raster(os.path.join(projdir, 'flowdirection'))
        fill2 = arcpy.Raster(os.path.join(projdir, 'topography'))
        order2 = arcpy.Raster(os.path.join(projdir, 'streamorder'))
        frxst_raster = arcpy.Raster(os.path.join(projdir, 'frxst_pts'))         # Added 08/23/2018 by KMS to include forecast points in reach-based routing file
        channelgrid = arcpy.Raster(os.path.join(projdir, 'channelgrid'))        # Added 06/03/2019 by KMS

        # Get georeference inforormation from unzipped raster layer
        sr = arcpy.Describe(fdir).spatialReference
        WKT = sr.exportToString().replace("'", '"')
        arcpy.env.outputCoordinateSystem = sr

        # Change CHANNELGRID so that it has values of 1 and Null
        strm = SetNull(channelgrid, 1, "VALUE = %s" %wrfh.NoDataVal)

        # Use topography raster to create reach-based routing files
        if frxst_raster.maximum == float(wrfh.NoDataVal):                       # Added 08/23/2018 by KMS to include forecast points in reach-based routing file
            frxst_raster = None                                                 # Default is no forecast points for reach-based routing file
        else:
            frxst_raster = SetNull(frxst_raster, frxst_raster, "VALUE = %s" %wrfh.NoDataVal)
        linkid = wrfh.Routing_Table(arcpy, projdir, sr, strm, fdir, fill2, order2, gages=frxst_raster)
        linkid_arr = arcpy.RasterToNumPyArray(linkid)

        # Add new LINKID grid to the FullDom file
        rootgrp = netCDF4.Dataset(os.path.join(projdir, FullDom), 'r+')         # Read+ object on old FullDom file
        linkvar = rootgrp.createVariable('LINKID', 'i4', ('y','x'))
        linkvar.setncatts(rootgrp.variables['TOPOGRAPHY'].__dict__)             # Steal variable attributes from a known variable
        linkvar[:] = linkid_arr                                                 # Populate variable using array
        rootgrp.close()

        # Perform cleanup before zipping up output
        del linkvar, linkid_arr, linkid, strm, fdir, fill2, order2, sr, rootgrp, WKT

        # Zip everything back up (all possible files)
        zipper = wrfh.zipUpFolder(arcpy, projdir, out_zip, nclist)

        arcpy.Delete_management('in_memory')
        arcpy.AddMessage('Process completed without error.')
        arcpy.AddMessage('Output ZIP File: %s' %out_zip)
        del out_zip

        try:
            arcpy.env.workspace = projdir
            for infile in arcpy.ListDatasets():
                arcpy.Delete_management(infile)
            arcpy.Delete_management(projdir)
        except:
            arcpy.AddMessage('Could not delete scratch folder: %s' %projdir)
            arcpy.AddMessage('You will have to delete this yourself after closing ArcGIS applications.')
        del projdir
        return

class Lake_Parameter_Addition(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Add Lake Parameters"
        self.description = "This tool takes an input Routing Stack .ZIP file," + \
                           " typically created by the Process GEOGRID File tool," + \
                           " and adds in the reservoir parameter file and grids." + \
                           "Only use this tool if the existing routing stack was " + \
                           " built without reservoirs in the first place."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_zip = arcpy.Parameter(
            displayName="Routing Stack ZIP File",
            name="in_zip",
            datatype="File",
            parameterType="Required",
            direction="Input")
        in_zip.filter.list = ['zip']

        # Input parameter
        in_reservoirs = arcpy.Parameter(
            displayName="Reservoirs Shapefile or Feature Class",
            name="in_reservoirs",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input")

        # Output parameter
        out_zip = arcpy.Parameter(
            displayName="Output ZIP File",
            name="out_zip",
            datatype="File",
            parameterType="Required",
            direction="Output")

        parameters = [in_zip, in_reservoirs, out_zip]   # , in_LakeIDField
        return parameters

    def isLicensed(self):
        """Allow the tool to execute, only if the ArcGIS Spatial Analyst extension
        is available."""
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False                                                            # tool cannot be executed
        return True                                                                 # tool can be executed

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        reload(wrfh)                                                            # Reload in case code changes have been made
        import arcpy                                                            # Not at all sure why this must be here, but it must

        # Gather all necessary parameters
        in_zip = parameters[0].valueAsText
        in_lakes = parameters[1].valueAsText
        out_zip = parameters[2].valueAsText

        wrfh.printMessages(arcpy, ['Begining processing on {0}'.format(time.ctime())])
        wrfh.printMessages(arcpy, ['Beginning to extract WRF-Hydro routing grids...'])

        # Create scratch directory for temporary outputs
        projdir = arcpy.CreateScratchName("temp", data_type="Folder", workspace=arcpy.env.scratchFolder)
        os.mkdir(projdir)
        arcpy.env.overwriteOutput = True
        arcpy.env.workspace = projdir
        arcpy.env.scratchWorkspace = projdir

        # Unzip to a known location (make sure no other nc files live here)
        FullDom = wrfh.FullDom
        GWBasins = wrfh.GWGRID_nc
        out_sfolder = wrfh.Examine_Outputs(arcpy, in_zip, projdir, skipfiles=[FullDom, GWBasins])

        # Prepare other rasters for the Lake Routing function
        channelgrid = arcpy.Raster(os.path.join(projdir, 'CHANNELGRID'))
        flac = arcpy.Raster(os.path.join(projdir, 'flowacc'))
        fill2 = arcpy.Raster(os.path.join(projdir, 'topography'))

        # Check to see if this was a reach-based routing domain
        if os.path.exists(os.path.join(projdir, 'LINKID')):
            Gridded = False
        else:
            Gridded = True

        # Get georeference inforormation from unzipped raster layer
        sr = arcpy.Describe(channelgrid).spatialReference
        arcpy.env.outputCoordinateSystem = sr
        cellsize = channelgrid.meanCellHeight

        # Run the lake addition function
        arcpy, channelgrid_arr, lakegrid_arr = wrfh.add_reservoirs(arcpy, channelgrid,
            in_lakes, flac, projdir, fill2, cellsize, sr, lakeIDfield=None, Gridded=Gridded)
        del flac, fill2, channelgrid, sr, cellsize, in_lakes, in_zip

        # Add new LINKID grid to the FullDom file
        rootgrp = netCDF4.Dataset(os.path.join(projdir, FullDom), 'r+')         # Read+ object on old FullDom file

        # Save the gridded lake array to the Fulldom file
        if Gridded:
            rootgrp.variables['LAKEGRID'][:] = lakegrid_arr                     # Write array to output netCDF file
            rootgrp.variables['CHANNELGRID'][:] = channelgrid_arr               # Populate variable using array
        else:
            rootgrp.variables['LAKEGRID'][:] = wrfh.NoDataVal                   # No need to populate LAKEGRID if not using gridded lakes

        rootgrp.close()
        del rootgrp, lakegrid_arr, channelgrid_arr

        # Zip everything back up (all possible files)
        zipper = wrfh.zipUpFolder(arcpy, projdir, out_zip, nclist)
        del zipper

        arcpy.Delete_management('in_memory')
        arcpy.AddMessage('Process completed without error.')
        arcpy.AddMessage('Output ZIP File: %s' %out_zip)
        del out_zip

        try:
            shutil.rmtree(projdir)
        except:
            arcpy.AddMessage('Could not delete scratch folder: %s' %projdir)
            arcpy.AddMessage('You will have to delete this yourself after closing ArcGIS applications.')
        del projdir
        return

class GWBUCKPARM(object):

    """This function will build a GWBUCKPARM table out of a variety of inputs."""

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Build Groundwater Inputs"
        self.description = "This tool takes a FullDom file and will use the basn_msk " + \
                           " grid to build a 1D GWBUCKPARM Table and a 2D GWBASINS grid."
        self.canRunInBackground = True
        self.category = "Utilities"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter
        in_nc = arcpy.Parameter(
            displayName="Input FullDom File",
            name="in_nc",
            datatype="File",
            parameterType="Required",
            direction="Input")

        # Input parameter
        in_geo = arcpy.Parameter(
            displayName="Input Geogrid File",
            name="in_geo",
            datatype="File",
            parameterType="Required",
            direction="Input")

        # Input parameter
        in_method = arcpy.Parameter(
            displayName="Method for deriving groundwater basins",
            name="in_method",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        in_method.filter.type = "ValueList"
        in_method.filter.list = ['FullDom basn_msk variable', 'FullDom LINKID local basins', 'Polygon Shapefile or Feature Class'] #

        # Input parameter
        in_Polys = arcpy.Parameter(
            displayName="Polygon feature class to define groundwater basins",
            name="in_Polygons",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Input")
        in_Polys.enabled = False

        # Output parameter
        out_dir = arcpy.Parameter(
            displayName="Output Directory",
            name="out_dir",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

        parameters = [in_nc, in_geo, in_method, in_Polys, out_dir]
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        try:
            if not arcpy.CheckExtension("Spatial") == "Available":
                raise Exception
        except Exception:
            return False                                                        # tool cannot be executed
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        # Only activate masking if a CSV file has been input
        if parameters[2].altered == True:
            if parameters[2].valueAsText == 'Polygon Shapefile or Feature Class':
                parameters[3].enabled = True
            else:
                parameters[3].enabled = False
        else:
            parameters[3].enabled = False

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        reload(wrfh)                                                            # Reload in case code changes have been made

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        in_geo = parameters[1].valueAsText
        in_method = parameters[2].valueAsText
        in_Polys = parameters[3].valueAsText
        out_dir = parameters[4].valueAsText

        outtable = os.path.join(out_dir, 'GroundwaterBasins.log')
        tee = wrfh.TeeNoFile(outtable, 'w')

        # Print inputs to screen
        wrfh.printMessages(arcpy, ['Begining processing on {0}'.format(time.ctime())])
        wrfh.printMessages(arcpy, ['Input parameters:'])
        for param in parameters:
            wrfh.printMessages(arcpy, ['    Parameter: {0}: {1}'.format(param.displayName, param.valueAsText)])

        # Create scratch directory for temporary outputs
        projdir = os.path.join(out_dir, 'gw_scratchdir')
        if os.path.exists(projdir):
            wrfh.printMessages(arcpy, ['  Using existing scratch directory'])
        else:
            wrfh.printMessages(arcpy, ['  Creating scratch directory {0}'.format(projdir)])
            os.makedirs(projdir)
        arcpy.env.overwriteOutput = True
        arcpy.env.workspace = projdir
        arcpy.env.scratchWorkspace = projdir

        # Georeference the Fulldom file
        rootgrp = netCDF4.Dataset(in_nc, 'r')
        fine_grid = wrfh.WRF_Hydro_Grid(arcpy, rootgrp)                         # Instantiate fine grid object
        rootgrp2 = netCDF4.Dataset(in_geo, 'r')
        coarse_grid = wrfh.WRF_Hydro_Grid(arcpy, rootgrp2)                      # Instantiate coarse grid object
        rootgrp.close()
        rootgrp2.close()

        # Use the fine grid file to create the basin inputs
        wrfh.printMessages(arcpy, ['    GeoTransform: {0}'.format(fine_grid.GeoTransformStr())])    # Print affine transformation to screen.
        GWBasns, GWBasns_arr = wrfh.build_GW_Basin_Raster(arcpy, in_nc, projdir, in_method, fine_grid, in_Polys=in_Polys)

        # Resample to coarse grid
        wrfh.build_GW_buckets(arcpy, out_dir, GWBasns, GWBasns_arr, coarse_grid, Grid=True)
        del fine_grid, coarse_grid, rootgrp, rootgrp2

        # Clean up
        shutil.rmtree(projdir)
        tee.close()
        del projdir, GWBasns, GWBasns_arr, in_geo, in_Polys, in_method, tee
        return

# --- End Toolbox Classes --- #