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
import wrf_hydro_functions                                                      # Function script packaged with this toolbox
reload(wrf_hydro_functions)                                                     # Re-load the function script in case of script changes
# --- End Import Modules --- #

# --- Module Configurations --- #
arcpy.env.overwriteOutput = True                                                # Allow overwriting of outputs
# --- End Module Configurations --- #

# Find out if we are in 32-bit or 64-bit
if sys.maxsize > 2**32:
    bit64 = True
else:
    bit64 = False

# --- Globals --- #
inunits = 'm'                                                                   # Set units for input Elevation raster dataset: 'm' or 'cm'
outNCType = 'NETCDF3_64BIT'                                                     # Set output netCDF format for spatial metdata files
LK_nc = wrf_hydro_functions.LK_nc                                               # Grab default lake parameter table name from function script
RT_nc = wrf_hydro_functions.RT_nc                                               # Grab default Route Link parameter table name from function script

# Processing Notes to insert into netCDF global attributes
# Processing notes for Spatial Metdata files
processing_notes_SM = '''Created: %s''' %time.ctime()

# Processing notes for the FULLDOM (Routing Grid) file
processing_notesFD = '''Created: %s''' %time.ctime()
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
                        Lake_Parameter_Addition]

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

        in_LakeIDField = arcpy.Parameter(
            displayName="ID field (Integer) for identifying lakes",
            name="in_IDfield",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        in_LakeIDField.filter.type = "ValueList"

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

        parameters = [in_nc, in_csv, basin_mask, RB_routing, Lake_routing, in_reservoirs, in_LakeIDField, in_raster, cellsize, threshold, ovroughrtfac_val, retdeprtfac_val, out_zip]
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
            parameters[6].enabled = True
        else:
            parameters[5].enabled = False
            parameters[6].enabled = False

        # Populate lake ID field combo box with list of Integer type fields from reservoirs shapefile
        if parameters[5].altered:
            in_lakes_file = parameters[5].valueAsText
            FieldNames = [field.name for field in arcpy.ListFields(in_lakes_file, "", "Integer")]
            parameters[6].filter.list = FieldNames

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""

        # Ensure that if a LAKEPARM table is requested, a lakes feature class is provided
        if parameters[4].value == True and parameters[5].value is None:
            parameters[5].setErrorMessage("Must specify Reservoirs Shapefile or Feature Class.")
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        tic = time.time()
        reload(wrf_hydro_functions)                                             # Reload in case code changes have been made
        isError = False                                                         # Starting Error condition (no errors)

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        in_csv = parameters[1].valueAsText
        basin_mask = parameters[2].valueAsText
        routing = parameters[3].value
        Lake_routing = parameters[4].value
        in_reservoir = parameters[5].valueAsText
        lakeIDfield = parameters[6].valueAsText
        in_raster = parameters[7].valueAsText
        cellsize = parameters[8].value
        threshold = parameters[9].value
        ovroughrtfac_val = parameters[10].value
        retdeprtfac_val = parameters[11].value
        out_zip = parameters[12].valueAsText

        # Prepare output log file
        outtable = open(os.path.join(os.path.dirname(out_zip), os.path.basename(out_zip) + '.log'), "w")
        loglines = ['Begining processing on %s' %time.ctime()]
        loglines.append('64-bit: %s' %bit64)
        loglines.append('Input parameters:')
        for param in parameters:
            loglines.append('    Parameter: %s: %s' %(param.displayName, param.valueAsText))
        outtable.writelines("\n".join(loglines) + "\n")

        if lakeIDfield is None:
            loglines = ['  No Lake ID Field provided']
            arcpy.AddMessage(loglines[-1])

        # Interpret the input for reservoir routing
        if Lake_routing is False:
            in_lakes = None
        else:
            in_lakes = in_reservoir

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
        sequentially, especially for netCDF3 versions.'''
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
                    ['landuse', 'f4', '']]
        # Add variables depending on the input options
        if routing:
            varList2D.append(['LINKID', 'i4', ''])

        # Step 1 - Georeference geogrid file
        LU_INDEX, sr2, Projection_String, map_pro, GeoTransform, loglines = wrf_hydro_functions.georeference_geogrid_file(arcpy, in_nc, 'LU_INDEX')            # Process: Generate LU Index grid
        outtable.writelines("\n".join(loglines) + "\n")
        hgt_m_raster, sr2, Projection_String, map_pro, GeoTransform, loglines = wrf_hydro_functions.georeference_geogrid_file(arcpy, in_nc, 'HGT_M')
        outtable.writelines("\n".join(loglines) + "\n")
        ##LANDMASK, sr2, Projection_String, map_pro, GeoTransform, loglines = wrf_hydro_functions.georeference_geogrid_file(arcpy, in_nc, 'LANDMASK')            # Process: Generate LANDMASK grid (1=Land, 0=Water)
        ##outtable.writelines("\n".join(loglines) + "\n")

        # Create spatial metadata file for GEOGRID/LDASOUT grids
        descData = arcpy.Describe(hgt_m_raster)
        DXDY_dict = {u'DX': float(descData.meanCellWidth), u'DY': float(descData.meanCellHeight)}
        out_nc1 = os.path.join(projdir, wrf_hydro_functions.LDASFile)
        rootgrp1 = netCDF4.Dataset(out_nc1, 'w', format=outNCType)              # wrf_hydro_functions.outNCType)
        rootgrp1, grid_mapping, loglines = wrf_hydro_functions.create_CF_NetCDF(arcpy, hgt_m_raster, rootgrp1, map_pro, projdir, DXDY_dict,
                sr2, GeoTransform, addLatLon=False, notes=processing_notes_SM, loglines=loglines)
        rootgrp1.close()
        del descData

        # Step 2 - Create high resolution topography layers
        mosprj, cellsize1, cellsize2, loglines = wrf_hydro_functions.create_high_res_topogaphy(arcpy, in_raster, hgt_m_raster, cellsize, sr2, projdir)
        outtable.writelines("\n".join(loglines) + "\n")

        # Create FULLDOM file
        descData2 = arcpy.Describe(mosprj)
        DXDY_dict = {u'DX': float(descData2.meanCellWidth), u'DY': float(descData2.meanCellHeight)}
        GT_bits = GeoTransform.split(" ")                                       # Split up GeoTransform string for replacing DX and DY
        GeoTransform = '%s %s %s %s %s %s ' %(GT_bits[0], DXDY_dict['DX'], GT_bits[2], GT_bits[3], GT_bits[4], -DXDY_dict['DY'])    # Alter DX/DY
        out_nc2 = os.path.join(projdir, wrf_hydro_functions.FullDom)
        rootgrp2 = netCDF4.Dataset(out_nc2, 'w', format=outNCType)              # wrf_hydro_functions.outNCType)
        rootgrp2, grid_mapping, loglines = wrf_hydro_functions.create_CF_NetCDF(arcpy, mosprj, rootgrp2, map_pro, projdir, DXDY_dict,
                sr2, GeoTransform, addLatLon=True, notes=processing_notesFD, loglines=loglines, addVars=varList2D)
        del descData2

        ### Step X(a) - Test to match LANDMASK - Only used for areas surrounded by water (LANDMASK=0)
        ##mosprj2, loglines = wrf_hydro_functions.adjust_to_landmask(arcpy, mosprj, LANDMASK, sr2, projdir, inunits)
        ##outtable.writelines("\n".join(loglines) + "\n")
        ##del LANDMASK

        try:

            # Step 4 - Hyrdo processing functions
            rootgrp2, loglines = wrf_hydro_functions.sa_functions(arcpy, rootgrp2, basin_mask, mosprj, ovroughrtfac_val, retdeprtfac_val, projdir, in_csv, threshold, inunits, LU_INDEX, cellsize1, cellsize2, routing, in_lakes, lakeIDfield) # , mosprj2,
            outtable.writelines("\n".join(loglines) + "\n")
            rootgrp2.close()

        except Exception as e:
            loglines.append('Exception: %s' %e)
            arcpy.AddMessage(loglines[-1])
            rootgrp2.close()
            isError = True

        # zip the folder
        nclist = ['GEOGRID_LDASOUT_Spatial_Metadata.nc',
                    wrf_hydro_functions.FullDom,
                    'gw_basns.nc',
                    'gw_basns_geogrid.txt',
                    'gw_basns_geogrid.prj',
                    RT_nc,
                    'Route_Link.csv',
                    'LAKEPARM.TBL',
                    LK_nc,
                    'streams.shp', 'streams.shx', 'streams.shp.xml', 'streams.sbx', 'streams.sbn', 'streams.prj', 'streams.dbf',
                    'lakes.shp', 'lakes.shx', 'lakes.shp.xml', 'lakes.sbx', 'lakes.sbn', 'lakes.prj', 'lakes.dbf']
        zipper = wrf_hydro_functions.zipUpFolder(arcpy, projdir, out_zip, nclist)

        # Clean up and give finishing message
        del LU_INDEX, hgt_m_raster
        if isError:
            loglines = ['Error encountered after %s seconds.' %(time.time()-tic)]
            arcpy.AddMessage(loglines[-1])
            raise SystemExit
        else:
            loglines = ['Completed without error in %s seconds.' %(time.time()-tic)]
            arcpy.AddMessage(loglines[-1])
            shutil.rmtree(projdir)
        outtable.write(loglines[-1])
        outtable.close()
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
            ncFP = arcpy.NetCDFFileProperties(in_nc_file)

            # Loop through global variables in NetCDF file to gather projection information
            ncVarNames = ncFP.getVariablesByDimension('west_east')
            ncMassgridNames = []
            for x in ncVarNames:
                mgridvar = ncFP.getAttributeValue(x, 'stagger')                 # Only use variables on Massgrid for now ('M')
                if mgridvar == 'M':
                    ncMassgridNames.append(x)
            parameters[1].filter.list = ncMassgridNames
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        reload(wrf_hydro_functions)                                             # Reload in case code changes have been made

        # Set environments
        arcpy.env.overwriteOutput = True

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        Variable = parameters[1].valueAsText
        out_raster = parameters[2].valueAsText

        # Use wrf_hydro_functions to perform process
        nc_raster, sr2, Projection_String, map_pro, GeoTransform, loglines = wrf_hydro_functions.georeference_geogrid_file(arcpy, in_nc, Variable)

        # Set environments and save
        arcpy.env.outputCoordinateSystem = sr2
        nc_raster.save(out_raster)
        arcpy.DefineProjection_management(out_raster, sr2)
        del nc_raster

        arcpy.AddMessage('    Process completed without error.')
        arcpy.AddMessage('    Output Raster: %s' %out_raster)
        return

class ExamineOutputs(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Examine Outputs of GIS Preprocessor"
        self.description = "This tool takes the output zip file from the ProcessGeogrid script" + \
                           "and creates a raster from each output NetCDF file." + \
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

        reload(wrf_hydro_functions)                                             # Reload in case code changes have been made

        # Initiate input and output parameters
        in_zip = parameters[0].valueAsText
        out_folder = parameters[1].valueAsText

        # Create output directory
        os.mkdir(out_folder)

        # Use wrf_hydro_functions to perform process
        out_sfolder = wrf_hydro_functions.Examine_Outputs(arcpy, in_zip, out_folder)
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

        reload(wrf_hydro_functions)                                             # Reload in case code changes have been made

        # Set environments
        arcpy.env.overwriteOutput = True

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        out_prj = parameters[1].valueAsText
        Variable = 'HGT_M'

        # Use wrf_hydro_functions to perform process
        nc_raster, sr2, Projection_String, map_pro, GeoTransform, loglines = wrf_hydro_functions.georeference_geogrid_file(arcpy, in_nc, Variable)
        del nc_raster, sr2

        # Optional save .prj file
        with open(out_prj, 'w') as prj_file:
             prj_file.write(Projection_String)

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

        # Output parameter
        out_folder = arcpy.Parameter(
            displayName="Output Folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")
        out_folder.defaultEnvironmentName = "workspace"

        parameters = [in_raster, out_folder]
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
        arcpy.env.overwriteOutput = True

        # Gather all necessary parameters
        inraster = parameters[0].valueAsText
        projdir = parameters[1].valueAsText
        arcpy.AddMessage('Input Raster Dataset: %s' %inraster)
        arcpy.AddMessage('Directory to be used for outputs: %s' %projdir)

        # Create latitude and longitude rasters
        in_raster = arcpy.Raster(inraster)
        loglines, xout2, yout2, xmap, ymap = wrf_hydro_functions.create_lat_lon_rasters(arcpy, projdir, in_raster)
        arcpy.Delete_management(xmap)
        arcpy.Delete_management(ymap)
        del xmap, ymap, loglines

        # Write to disk
        xout2.save(os.path.join(projdir, 'longitude'))
        yout2.save(os.path.join(projdir, 'latitude'))

        arcpy.AddMessage('    Process completed without error.')
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
        in_nc.filter.list = ['nc']

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

        reload(wrf_hydro_functions)                                             # Reload in case code changes have been made

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

        # Print informational messages
        arcpy.AddMessage('Input Raster Dataset: %s' %in_nc)
        arcpy.AddMessage('Output Grid Resolution: %s' %format_out)
        arcpy.AddMessage('Output Regridding Factor: %s' %factor)
        arcpy.AddMessage('Directory to be used for outputs: %s' %projdir)
        arcpy.AddMessage('Output netCDF File: %s' %out_nc)

        # Prepare output log file
        outtable = open(os.path.join(projdir, os.path.basename(out_nc) + '.log'), "w")
        loglines = ['Begining processing on %s' %time.ctime()]
        loglines.append('64-bit: %s' %bit64)
        loglines.append('Input parameters:')
        for param in parameters:
            loglines.append('    Parameter: %s: %s' %(param.displayName, param.valueAsText))

        # Georeference geogrid file
        in_raster, sr, Projection_String, map_pro, GeoTransformStr, loglines2 = wrf_hydro_functions.georeference_geogrid_file(arcpy, in_nc, 'HGT_M')
        loglines += loglines2
        descData = arcpy.Describe(in_raster)
        DXDY_dict = {u'DX': descData.meanCellWidth/factor, u'DY': descData.meanCellHeight/factor}
        loglines.append('    New Resolution: %s %s' %(DXDY_dict[u'DX'], DXDY_dict[u'DY']))
        arcpy.AddMessage(loglines[-1])

        # Record GEOGRID MAP_PROJ attribute
        loglines.append('    Map Projection of GEOGRID: %s' %wrf_hydro_functions.projdict[map_pro])
        arcpy.AddMessage(loglines[-1])
        loglines.append('    Esri PE String: %s' %Projection_String)
        arcpy.AddMessage(loglines[-1])

        # Create high resolution raster for RTOUT output using Spatial Analyst CreateConstantRaster function
        if format_out == "RTOUT":
            if arcpy.CheckExtension("Spatial") == "Available":
                arcpy.CheckOutExtension("Spatial")
                from arcpy.sa import *
            arcpy.env.snapRaster = in_raster
            arcpy.env.outputCoordinateSystem = sr
            in_raster = CreateConstantRaster(1, "INTEGER", DXDY_dict[u'DX'], descData.Extent)

        # Create the netCDF file with spatial metadata
        rootgrp = netCDF4.Dataset(out_nc, 'w', format=outNCType)
        rootgrp, grid_mapping, loglines = wrf_hydro_functions.create_CF_NetCDF(arcpy, in_raster, rootgrp, map_pro, projdir,
                DXDY_dict, sr, GeoTransformStr, latlon_vars, notes=processing_notes_SM, loglines=loglines)

        rootgrp.close()
        del in_raster
        loglines += ['Completed without error in %s seconds.' %(time.time()-tic1)]

        # Clean up and give finishing message
        arcpy.AddMessage(loglines[-1])
        outtable.writelines("\n".join(loglines) + "\n")
        outtable.close()
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

        reload(wrf_hydro_functions)                                             # Reload in case code changes have been made

        # Gather all necessary parameters
        in_nc = parameters[0].valueAsText
        out_shp = parameters[1].valueAsText

        # Create scratch directory for temporary outputs
        projdir = os.path.dirname(in_nc)
        arcpy.env.overwriteOutput = True
        arcpy.env.workspace = projdir
        arcpy.env.scratchWorkspace = projdir

        # Get coarse grid raster
        hgt_m_raster, sr2, Projection_String, map_pro, GeoTransform, loglines = wrf_hydro_functions.georeference_geogrid_file(arcpy, in_nc, 'HGT_M')
        del Projection_String, loglines

        # Use coarse grid raster to create shapefile
        wrf_hydro_functions.domain_shapefile(arcpy, hgt_m_raster, out_shp, sr2)

        # Delete intermediate files in scratch dir
        del hgt_m_raster, sr2

        return parameters

class Reach_Based_Routing_Addition(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Add reach-based routing"
        self.description = "This tool takes an input Routing Stack .ZIP file," + \
                           " typically created by the ProcessGeogridFile tool," + \
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

        reload(wrf_hydro_functions)                                             # Reload in case code changes have been made
        from arcpy.sa import *

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
        loglines = ['Begining processing on %s' %time.ctime()]
        arcpy.AddMessage('Beginning to extract WRF routing grids...')

        # Unzip to a known location (make sure no other nc files live here)
        FullDom = wrf_hydro_functions.FullDom
        out_sfolder = wrf_hydro_functions.Examine_Outputs(arcpy, in_zip, projdir, skipfiles=[FullDom])

        # Prepare other rasters for the Routing Table function
        fdir = arcpy.Raster(os.path.join(projdir, 'flowdirection'))
        fill2 = arcpy.Raster(os.path.join(projdir, 'topography'))
        order2 = arcpy.Raster(os.path.join(projdir, 'streamorder'))

        # Get georeference inforormation from unzipped raster layer
        sr = arcpy.Describe(fdir).spatialReference
        arcpy.env.outputCoordinateSystem = sr

        # Change CHANNELGRID so that it has values of 1 and Null
        strm = SetNull(arcpy.Raster('channelgrid'), 1, "VALUE = -9999")

        # Use topography raster to create reach-based routing files
        linkid, loglines = wrf_hydro_functions.Routing_Table(arcpy, projdir, sr, strm, fdir, fill2, order2, loglines)
        linkid_arr = arcpy.RasterToNumPyArray(linkid)

        # Add new LINKID grid to the FullDom file
        rootgrp = netCDF4.Dataset(os.path.join(projdir, FullDom), 'r+')     # Read+ object on old FullDom file
        linkvar = rootgrp.createVariable('LINKID', 'i4', ('y','x'))
        linkvar.setncatts(rootgrp.variables['TOPOGRAPHY'].__dict__)             # Steal variable attributes from a known variable
        linkvar[:] = linkid_arr                                                 # Populate variable using array
        rootgrp.close()

        # Perform cleanup before zipping up output
        del linkvar, linkid_arr, linkid, strm, fdir, fill2, order2, sr, rootgrp

        # Attempt to write a new FullDom file instead of adding a variable to the existing one, which can be slow
        ##        # Write new netCDF file to output directory (raster than adding a new variable to existing NC file, for netCDF3)
        ##        old_NC = os.rename(os.path.join(projdir, FullDom), os.path.join(projdir, 'OldNC.nc'))   # Rename old fulldom File
        ##        rootgrp_old = netCDF4.Dataset(old_NC, 'r')                              # Read object on old FullDom file
        ##        rootgrp_new = netCDF4.Dataset(os.path.join(projdir, FullDom), 'w')  # Write object on old FullDom file
        ##
        ##        # Copy dimensions from old FullDom file to new FullDom file
        ##        for dimname, dim in rootgrp_old.dimensions.iteritems():
        ##            rootgrp_new.createDimension(dimname, len(dim))                      # Copy other dimensions from the WRF-Hydro output file
        ##
        ##        # Copy variables from WRF-Hydro output file, adding variable attributes as necessary
        ##        for varname, ncvar in rootgrp_old.variables.iteritems():
        ##            var = rootgrp_new.createVariable(varname, ncvar.dtype, ncvar.dimensions)    # Create variable and define dimensions
        ##            var.setncatts(ncvar.__dict__)                                       # Set variable attributes from old FullDom file
        ##        rootgrp_new.createVariable('LINKID', 'i4', ('y','x'))
        ##        linkvar.setncatts(rootgrp_new.variables['TOPOGRAPHY'].__dict__)         # Steal variable attributes from a known variable
        ##        rootgrp_new.setncatts(rootgrp_old.__dict__)                             # Copy global attributes from old FullDom file
        ##
        ##        # Add variable values last (makes the script run faster)
        ##        for varname, ncvar in rootgrp_old.variables.iteritems():
        ##            var = rootgrp_new.variables[varname]
        ##            var[:] = ncvar[:]                                                   # Copy the variable data into the newly created variable
        ##        rootgrp_new.variables['LINKID'][:] = linkid_arr                         # Populate variable using array
        ##        rootgrp_old.close()
        ##        rootgrp_new.close()
        ##        os.remove(os.path.join(projdir, FullDom))                           # Remove the old FullDom file

        # Zip everything back up (all possible files)
        nclist = ['GEOGRID_LDASOUT_Spatial_Metadata.nc',
                    FullDom,
                    'gw_basns.nc',
                    'gw_basns_geogrid.txt',
                    'gw_basns_geogrid.prj',
                    RT_nc,
                    'Route_Link.csv',
                    'LAKEPARM.TBL',
                    LK_nc,
                    'streams.shp', 'streams.shx', 'streams.shp.xml', 'streams.sbx', 'streams.sbn', 'streams.prj', 'streams.dbf',
                    'lakes.shp', 'lakes.shx', 'lakes.shp.xml', 'lakes.sbx', 'lakes.sbn', 'lakes.prj', 'lakes.dbf']
        zipper = wrf_hydro_functions.zipUpFolder(arcpy, projdir, out_zip, nclist)

        arcpy.Delete_management('in_memory')
        arcpy.AddMessage('Process completed without error.')
        arcpy.AddMessage('Output ZIP File: %s' %out_zip)
        del out_zip, nclist

        try:
            shutil.rmtree(projdir)
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
                           " and adds in the lake parameter file and grids."
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

        # Input parameter
        in_LakeIDField = arcpy.Parameter(
            displayName="ID field (INTEGER) for identifying lakes",
            name="in_IDfield",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        in_LakeIDField.filter.type = "ValueList"

        # Output parameter
        out_zip = arcpy.Parameter(
            displayName="Output ZIP File",
            name="out_zip",
            datatype="File",
            parameterType="Required",
            direction="Output")

        parameters = [in_zip, in_reservoirs, in_LakeIDField, out_zip]
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
        if parameters[1].altered:
            in_lakes_file = parameters[1].valueAsText
            FieldNames = [field.name for field in arcpy.ListFields(in_lakes_file, "", "Integer")]
            parameters[2].filter.list = FieldNames
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        reload(wrf_hydro_functions)                                             # Reload in case code changes have been made
        import arcpy                                                            # Not at all sure why this must be here, but it must

        # Gather all necessary parameters
        in_zip = parameters[0].valueAsText
        in_lakes = parameters[1].valueAsText
        lakeIDfield = parameters[2].valueAsText
        out_zip = parameters[3].valueAsText

        loglines = ['Begining processing on %s' %time.ctime()]
        arcpy.AddMessage('Beginning to extract WRF-Hydro routing grids...')

        if lakeIDfield is None:
            loglines = ['  No Lake ID Field provided']
            arcpy.AddMessage(loglines[-1])

        # Create scratch directory for temporary outputs
        projdir = arcpy.CreateScratchName("temp", data_type="Folder", workspace=arcpy.env.scratchFolder)
        os.mkdir(projdir)
        arcpy.env.overwriteOutput = True
        arcpy.env.workspace = projdir
        arcpy.env.scratchWorkspace = projdir

        # Unzip to a known location (make sure no other nc files live here)
        FullDom = wrf_hydro_functions.FullDom
        out_sfolder = wrf_hydro_functions.Examine_Outputs(arcpy, in_zip, projdir, skipfiles=[FullDom])

        # Prepare other rasters for the Lake Routing function
        channelgrid = arcpy.Raster(os.path.join(projdir, 'CHANNELGRID'))
        flac = arcpy.Raster(os.path.join(projdir, 'flowacc'))
        fill2 = arcpy.Raster(os.path.join(projdir, 'topography'))

        # Get georeference inforormation from unzipped raster layer
        sr = arcpy.Describe(channelgrid).spatialReference
        arcpy.env.outputCoordinateSystem = sr
        cellsize = channelgrid.meanCellHeight

        # Run the lake addition function
        arcpy, channelgrid_arr, lakegrid_arr, loglines = wrf_hydro_functions.add_reservoirs(arcpy, channelgrid, in_lakes, flac, projdir, fill2, cellsize, sr, loglines, lakeIDfield)
        del flac, fill2, channelgrid, sr, cellsize, in_lakes, in_zip

        # Add new LINKID grid to the FullDom file
        rootgrp = netCDF4.Dataset(os.path.join(projdir, FullDom), 'r+')     # Read+ object on old FullDom file
        rootgrp.variables['LAKEGRID'][:] = lakegrid_arr                         # Populate variable using array
        rootgrp.variables['CHANNELGRID'][:] = channelgrid_arr                   # Populate variable using array
        rootgrp.close()
        del rootgrp, lakegrid_arr, channelgrid_arr

        # Zip everything back up (all possible files)
        nclist = ['GEOGRID_LDASOUT_Spatial_Metadata.nc',
                    FullDom,
                    'gw_basns.nc',
                    'gw_basns_geogrid.txt',
                    'gw_basns_geogrid.prj',
                    RT_nc,
                    'Route_Link.csv',
                    'LAKEPARM.TBL',
                    LK_nc,
                    'streams.shp', 'streams.shx', 'streams.shp.xml', 'streams.sbx', 'streams.sbn', 'streams.prj', 'streams.dbf',
                    'lakes.shp', 'lakes.shx', 'lakes.shp.xml', 'lakes.sbx', 'lakes.sbn', 'lakes.prj', 'lakes.dbf']
        zipper = wrf_hydro_functions.zipUpFolder(arcpy, projdir, out_zip, nclist)
        del zipper, nclist

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
# --- End Toolbox Classes --- #