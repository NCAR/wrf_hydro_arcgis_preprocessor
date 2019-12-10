# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# ** Copyright UCAR (c) 2018
# ** University Corporation for Atmospheric Research(UCAR)
# ** National Center for Atmospheric Research(NCAR)
# ** Research Applications Laboratory(RAL)
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# ** 2016/4/27 10:00:00
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

# Written by Kevin Sampson, NCAR
# Modified 2019/05/29

# --- Import Modules --- #
import sys
import os
import csv                                                                      # For reading input forecast points in CSV format
import time
import math
import zipfile
from zipfile import ZipFile, ZipInfo
import shutil
from collections import defaultdict                                             # Added 09/03/2015Needed for topological sorting algorthm
from itertools import takewhile, count                                          # Added 09/03/2015Needed for topological sorting algorthm                                                                  # Added 07/10/2015 for ExmineOutputs
import netCDF4
import numpy
from arcpy.sa import *                                                                 # Part of ArcGIS 10.3, added to this script 6/11/2015

try:
    import LAKE_PROCESSING_MODULE as LK_Routing
except:
    pass
# --- End Import Modules --- #

# --- Module Configurations --- #
sys.dont_write_bytecode = True                                                  # Do not write compiled (.pyc) files
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
               6: "latitude_longitude",
               0: "crs"}

# Output file types
outNCType = 'NETCDF4_CLASSIC'                                                   # Define the output netCDF version for RouteLink.nc and LAKEPARM.nc

# Default output file names
FullDom = 'Fulldom_hires.nc'                                                    # Default Full Domain routing grid nc file
LDASFile = 'GEOGRID_LDASOUT_Spatial_Metadata.nc'                                # Defualt LDASOUT domain grid nc file
LK_nc = 'LAKEPARM.nc'                                                           # Default Lake parameter table name [.nc]
LK_tbl = 'LAKEPARM.TBL'                                                         # Default Lake parameter table name [.TBL]
RT_nc = 'Route_Link.nc'                                                         # Default Route Link parameter table name
GW_nc = 'GWBUCKPARM.nc'                                                         # Default groundwater bucket parameter table name
GWGRID_nc = 'GWBASINS.nc'
GW_ASCII = 'gw_basns_geogrid.txt'                                               # Default Groundwater Basins ASCII grid output
GW_TBL = 'GWBUCKPARM.TBL'
StreamSHP = 'streams.shp'                                                       # Default streams shapefile name
TempLakeFIle = 'in_lakes_clip.shp'                                              # Temporary output lake file (clipped to domain)
LakesShp = 'lakes.shp'                                                          # Default lake shapefile name

# Options
maskRL = False                                                                  # Allow masking of channels in RouteLink file. May cause WRF-Hydro to crash if True
PpVersion = 'v5.1.1 (12/2019)'                                                  # WRF-Hydro ArcGIS Pre-processor version to add to FullDom metadata
CFConv = 'CF-1.5'                                                               # CF-Conventions version to place in the 'Conventions' attribute of RouteLink files

# Other Global Variables
NoDataVal = -9999                                                               # Default NoData value for gridded variables
walker = 3                                                                      # Number of cells to walk downstream before gaged catchment delineation
LK_walker = 3                                                                   # Number of cells to walk downstream to get minimum lake elevation
z_limit = 1000.0                                                                # Maximum fill depth (z-limit) between a sink and it's pour point
lksatfac_val = 1000.0                                                           # Default LKSATFAC value (unitless coefficient)
minDepth = 1.0                                                                  # Minimum active lake depth for lakes with no elevation variation

# Global attributes for altering the sphere radius used in computations. Do not alter sphere_radius for standard WRF-Hydro simulations
sphere_radius = 6370000.0                                                       # Radius of sphere to use (WRF Default = 6370000.0m)
wkt_text = "GEOGCS['GCS_Sphere_CUSTOM',DATUM['D_Sphere',SPHEROID['Sphere',%s,0.0]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];-400 -400 1000000000;-100000 10000;-100000 10000;8.99462786704589E-09;0.001;0.001;IsHighPrecision" %sphere_radius

# Unify all coordinate system variables to have the same name ("crs"). Ths makes it easier for WRF-Hydro output routines to identify the variable and transpose it to output files
crsVarname = True                                                               # Switch to make all coordinate system variables = "crs" instead of related to the coordinate system name
crsVar = CF_projdict[0]                                                         # Expose this as a global for other functions in other scripts to use
#crsVar = 'ProjectionCoordinateSystem'

# Point time-series CF-netCDF file coordinate system
'''Note that the point netCDF files are handled using a separate coordinate system than the grids.
This is because input data are usually in WGS84 or some other global spheroidal datum. We treat
these coordinates as though there is no difference between a sphere and a spheroid with respect
to latitude. Thus, we can choose an output coordinate system for the points, although no
transformation is performed. Properly transforming the points back and forth betwen sphere and
spheroid dramatically increases the runtime of the tools, with clear obvious benefit.'''
pointCF = True                                                                  # Switch to turn on CF-netCDF point time-series metadata attributes
pointSR = 4326                                                                  # The spatial reference system of the point time-series netCDF files (RouteLink, LAKEPARM). NAD83=4269, WGS84=4326

# Channel Routing default parameters
Qi = 0                                                                          # Initial Flow in link (cms)
MusK = 3600                                                                     # Muskingum routing time (s)
MusX = 0.2                                                                      # Muskingum weighting coefficient
n = 0.035                                                                       # Manning's roughness
ChSlp = 0.05                                                                    # Channel Side Slope (%; drop/length)
BtmWdth = 5                                                                     # Bottom Width of Channel (m)
Kc = 0                                                                          # channel conductivity (mm/hour)

#Default Lake Routing parameters
OrificeC = 0.1
OrificA = 1.0
WeirC = 0.4
WeirL = 10.0                                                                    # New default prescribed by D. Yates 5/11/2017 (10m default weir length). Old default weir length (0.0m).
ifd_Val = 0.90                                                                  # Default initial fraction water depth (90%)
out_LKtype = ['nc']                                                             # Default output lake parameter file format ['nc', 'ascii']
defaultLakeID = "NEWID"                                                         # Default new LakeID field for lakes, numbered 1..n

# Default groundwater bucket (GWBUCKPARM) parameters
coeff = 1.0000                                                                  # Bucket model coefficient
expon = 3.000                                                                   # Bucket model exponent
zmax = 50.00                                                                    # Conceptual maximum depth of the bucket
zinit = 10.0000                                                                 # Initial depth of water in the bucket model
out_2Dtype = ['nc']                                                             # Default output 2D groundwater bucket grid format ['nc', 'ascii']

#Custom Geotransformations for spheroid-to-sphere translation
geoTransfmName = "GeoTransform_Null_WRFHydro"                                   # Custom Geotransformation Name
customGeoTransfm = "GEOGTRAN[METHOD['Null']]"                                   # Null Transformation
#customGeoTransfm = "GEOGTRAN[METHOD['Geocentric_Translation'],PARAMETER['X_Axis_Translation',''],PARAMETER['Y_Axis_Translation',''],PARAMETER['Z_Axis_Translation','']]"   # Zero-parameter Geocentric Translation
#customGeoTransfm = "GEOGTRAN[METHOD['Geocentric_Translation'],PARAMETER['X_Axis_Translation',0.0],PARAMETER['Y_Axis_Translation',0.0],PARAMETER['Z_Axis_Translation',0.0]]"

# Elevation resampling method. Options: ["NEAREST", "BILINEAR", "CUBIC", "MAJORITY"]
ElevResampleMethod = "BILINEAR"                                                 # Method to be used for resampling high-resolution elevation to the model grid

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
        if (targetpath[-1:] in (os.path.sep, os.path.altsep) and len(os.path.splitdrive(targetpath)[1]) > 1):
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
        target = open(targetpath, "wb")                                         # 5/31/2019: Supporting Python3
        try:
            target.write(self.read(member.filename))
        finally:
            target.close()
        return targetpath

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

class WRF_Hydro_Grid:
    '''
    Class with which to create the WRF-Hydro grid representation. Provide grid
    information to initiate the class, and use getgrid() to generate a grid mesh
    and index information about the intersecting cells.

    Note:  The i,j index begins with (1,1) in the upper-left corner.
    '''
    def __init__(self, arcpy, rootgrp):
        """
        The first step of the process chain, in which the input NetCDF file gets
        georeferenced and projection files are created

        9/27/2017: This function was altered to use netCDF4-python library instead
        of arcpy.NetCDFFileProperties. Also changed the corner_lat and corner_lon
        index to reflect the lower left corner point rather than lower left center as
        the known point.

        10/6/2017: Proj4 string generation was added, with definitions adapted from
        Adapted from https://github.com/NCAR/wrf-python/blob/develop/src/wrf/projection.py

        5/29/2019: Removed MakeNetCDFRasterLayer_md function call to avoid bugs in ArcGIS
        10.7 (BUG-000122122, BUG-000122125).
        """
        corner_index = 12                                                       # 12 = Lower Left of the Unstaggered grid

        # First step: Import and georeference NetCDF file
        printMessages(arcpy, ['Step 1: NetCDF Conversion initiated...'])        # Initiate log list for this process

        # Read input WPS GEOGRID file
        # Loop through global variables in NetCDF file to gather projection information
        dimensions = rootgrp.dimensions
        globalAtts = rootgrp.__dict__                                           # Read all global attributes into a dictionary
        map_pro = globalAtts['MAP_PROJ']                                        # Find out which projection this GEOGRID file is in
        printMessages(arcpy, ['    Map Projection: %s' %projdict[map_pro]])

        # Collect grid corner XY and DX DY for creating ascii raster later
        if 'corner_lats' in globalAtts:
            corner_lat = float(globalAtts['corner_lats'][corner_index])         # Note: The values returned are corner points of the mass grid
        if 'corner_lons' in globalAtts:
            corner_lon = float(globalAtts['corner_lons'][corner_index])         # Note: The values returned are corner points of the mass grid
        if 'DX' in globalAtts:
            self.DX = float(globalAtts['DX'])
        if 'DY' in globalAtts:
            self.DY = float(globalAtts['DY'])

        # Collect necessary information to put together the projection file
        if 'TRUELAT1' in globalAtts:
            standard_parallel_1 = globalAtts['TRUELAT1']
        if 'TRUELAT2' in globalAtts:
            standard_parallel_2 = globalAtts['TRUELAT2']
        if 'STAND_LON' in globalAtts:
            central_meridian = globalAtts['STAND_LON']
        if 'POLE_LAT' in globalAtts:
            pole_latitude = globalAtts['POLE_LAT']
        if 'POLE_LON' in globalAtts:
            pole_longitude = globalAtts['POLE_LON']
        if 'MOAD_CEN_LAT' in globalAtts:
            printMessages(arcpy, ['    Using MOAD_CEN_LAT for latitude of origin.'])
            latitude_of_origin = globalAtts['MOAD_CEN_LAT']         # Added 2/26/2017 by KMS
        elif 'CEN_LAT' in globalAtts:
            printMessages(arcpy, ['    Using CEN_LAT for latitude of origin.'])
            latitude_of_origin = globalAtts['CEN_LAT']

        # Check to see if the input netCDF is a WPS-Generated Geogrid file.
        if 'TITLE' in globalAtts and 'GEOGRID' in globalAtts['TITLE']:
            self.isGeogrid = True
        else:
            self.isGeogrid = False
        del globalAtts

        # Handle expected dimension names from either Geogrid or Fulldom
        if 'south_north' in dimensions:
            self.nrows = len(dimensions['south_north'])
        elif 'y' in dimensions:
            self.nrows = len(dimensions['y'])
        if 'west_east' in dimensions:
            self.ncols = len(dimensions['west_east'])
        elif 'x' in dimensions:
            self.ncols = len(dimensions['x'])
        del dimensions

        # Create Projection file with information from NetCDF global attributes
        sr2 = arcpy.SpatialReference()
        if map_pro == 1:
            # Lambert Conformal Conic
            if 'standard_parallel_2' in locals():
                projname = 'Lambert_Conformal_Conic_2SP'
                printMessages(arcpy, ['    Using Standard Parallel 2 in Lambert Conformal Conic map projection.'])
            else:
                # According to http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Lambert_Conformal_Conic
                projname = 'Lambert_Conformal_Conic_1SP'
                standard_parallel_2 = standard_parallel_1
            Projection_String = ('PROJCS["Sphere_Lambert_Conformal_Conic",'
                                 'GEOGCS["GCS_Sphere",'
                                 'DATUM["D_Sphere",'
                                 'SPHEROID["Sphere",' + str(sphere_radius) + ',0.0]],'
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
            proj4 = ("+proj=lcc +units=m +a={} +b={} +lat_1={} +lat_2={} +lat_0={} +lon_0={} +x_0=0 +y_0=0 +k_0=1.0 +nadgrids=@null +wktext  +no_defs ".format(
                str(sphere_radius),
                str(sphere_radius),
                str(standard_parallel_1),
                str(standard_parallel_2),
                str(latitude_of_origin),
                str(central_meridian)))

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

            printMessages(arcpy, ['      Central Scale Factor: %s' %central_scale_factor])
            Projection_String = ('PROJCS["Sphere_Stereographic",'
                                 'GEOGCS["GCS_Sphere",'
                                 'DATUM["D_Sphere",'
                                 'SPHEROID["Sphere",' + str(sphere_radius) + ',0.0]],'
                                 'PRIMEM["Greenwich",0.0],'
                                 'UNIT["Degree",0.0174532925199433]],'
                                 'PROJECTION["Stereographic"],'
                                 'PARAMETER["False_Easting",0.0],'
                                 'PARAMETER["False_Northing",0.0],'
                                 'PARAMETER["Central_Meridian",' + str(central_meridian) + '],'
                                 'PARAMETER["Scale_Factor",' + str(central_scale_factor) + '],'
                                 'PARAMETER["Latitude_Of_Origin",' + str(standard_parallel_1) + '],'
                                 'UNIT["Meter",1.0]]')
            proj4 = ("+proj=stere +units=meters +a={} +b={} +lat0={} +lon_0={} +lat_ts={}".format(
                str(sphere_radius),
                str(sphere_radius),
				str(pole_latitude),
				str(central_meridian),
				str(standard_parallel_1)))

        elif map_pro == 3:
            # Mercator Projection
            Projection_String = ('PROJCS["Sphere_Mercator",'
                                 'GEOGCS["GCS_Sphere",'
                                 'DATUM["D_Sphere",'
                                 'SPHEROID["Sphere",' + str(sphere_radius) + ',0.0]],'
                                 'PRIMEM["Greenwich",0.0],'
                                 'UNIT["Degree",0.0174532925199433]],'
                                 'PROJECTION["Mercator"],'
                                 'PARAMETER["False_Easting",0.0],'
                                 'PARAMETER["False_Northing",0.0],'
                                 'PARAMETER["Central_Meridian",' + str(central_meridian) + '],'
                                 'PARAMETER["Standard_Parallel_1",' + str(standard_parallel_1) + '],'
                                 'UNIT["Meter",1.0]]')
            proj4 = ("+proj=merc +units=meters +a={} +b={} +lon_0={} +lat_ts={}".format(
                str(sphere_radius),
				str(sphere_radius),
				str(central_meridian),
				str(standard_parallel_1)))

        elif map_pro == 6:
            # Cylindrical Equidistant (or Rotated Pole)
            if pole_latitude != float(90) or pole_longitude != float(0):
                # if pole_latitude, pole_longitude, or stand_lon are changed from thier default values, the pole is 'rotated'.
                printMessages(arcpy, ['    Cylindrical Equidistant projection with a rotated pole is not currently supported.'])
                ##proj4 = ("+proj=ob_tran +o_proj=latlon +a={} +b={} +to_meter={} +o_lon_p={} +o_lat_p={} +lon_0={}".format(
                ##                        sphere_radius,
                ##                        sphere_radius,
                ##                        math.radians(1),
                ##                        180.0 - pole_longitude,
                ##                        pole_latitude,
                ##                        180.0 + pole_longitude))
                sys.exit(1)
            else:
                # Check units (linear unit not used in this projection).  GCS?
                Projection_String = ('PROJCS["Sphere_Equidistant_Cylindrical",'
                                     'GEOGCS["GCS_Sphere",'
                                     'DATUM["D_Sphere",'
                                     'SPHEROID["Sphere",' + str(sphere_radius) + ',0.0]],'
                                     'PRIMEM["Greenwich",0.0],'
                                     'UNIT["Degree",0.0174532925199433]],'
                                     'PROJECTION["Equidistant_Cylindrical"],'
                                     'PARAMETER["False_Easting",0.0],'
                                     'PARAMETER["False_Northing",0.0],'
                                     'PARAMETER["Central_Meridian",' + str(central_meridian) + '],'
                                     'PARAMETER["Standard_Parallel_1",' + str(standard_parallel_1) + '],'
                                     'UNIT["Meter",1.0]]')                          # 'UNIT["Degree", 1.0]]') # ?? For lat-lon grid?
                proj4 = ("+proj=eqc +units=meters +a={} +b={} +lon_0={}".format(str(sphere_radius), str(sphere_radius), str(central_meridian)))

        sr2.loadFromString(Projection_String)

        # Create a point geometry object from gathered corner point data
        sr1 = arcpy.SpatialReference()
        sr1.loadFromString(wkt_text)                                                # Load the Sphere datum CRS using WKT
        point = arcpy.Point()
        point.X = corner_lon
        point.Y = corner_lat
        pointGeometry = arcpy.PointGeometry(point, sr1)
        projpoint = pointGeometry.projectAs(sr2)                                    # Optionally add transformation method:

        # Get projected X and Y of the point geometry
        point2 = arcpy.Point(projpoint.firstPoint.X, projpoint.firstPoint.Y)    # Lower Left Corner Point

        # Populate class attributes
        self.x00 = point2.X                                                     # X value doesn't change from LLcorner to UL corner
        self.y00 = point2.Y + float(abs(self.DY)*self.nrows)                    # Adjust Y from LL to UL
        self.proj = sr2
        self.WKT = sr2.exportToString().replace("'", '"')                       # Replace ' with " so Esri can read the PE String properly when running NetCDFtoRaster
        self.proj4 = proj4
        self.point = point2
        self.map_pro = map_pro
        printMessages(arcpy, ['    Step 1 completed without error.'])
        del point, sr1, projpoint, pointGeometry

    def regrid(self, regrid_factor):
        '''
        Change the grid cell spacing while keeping all other grid parameters
        the same.
        '''
        self.DX = float(self.DX)/float(regrid_factor)
        self.DY = float(self.DY)/float(regrid_factor)
        self.nrows = int(self.nrows*regrid_factor)
        self.ncols = int(self.ncols*regrid_factor)
        print('  New grid spacing: dx={0}, dy={1}'.format(self.DX, self.DY))
        print('  New dimensions: rows={0}, cols={1}'.format(self.nrows, self.ncols))
        return self

    def GeoTransform(self):
        '''
        Return the affine transformation for this grid. Assumes a 0 rotation grid.
        (top left x, w-e resolution, 0=North up, top left y, 0 = North up, n-s pixel resolution (negative value))
        '''
        return (self.x00, self.DX, 0, self.y00, 0, -self.DY)

    def GeoTransformStr(self):
        return ' '.join([str(item) for item in self.GeoTransform()])

    def getxy(self):
        """
        This function will use the affine transformation (GeoTransform) to produce an
        array of X and Y 1D arrays. Note that the GDAL affine transformation provides
        the grid cell coordinates from the upper left corner. This is typical in GIS
        applications. However, WRF uses a south_north ordering, where the arrays are
        written from the bottom to the top.

        The input raster object will be used as a template for the output rasters.
        """
        print('    Starting Process: Building to XMap/YMap')

        # Build i,j arrays
        j = numpy.arange(self.nrows) + float(0.5)                               # Add 0.5 to estimate coordinate of grid cell centers
        i = numpy.arange(self.ncols) + float(0.5)                               # Add 0.5 to estimate coordinate of grid cell centers

        # col, row to x, y   From https://www.perrygeo.com/python-affine-transforms.html
        x = (i * self.DX) + self.x00
        y = (j * -self.DY) + self.y00
        del i, j

        # Create 2D arrays from 1D
        xmap = numpy.repeat(x[numpy.newaxis, :], y.shape, 0)
        ymap = numpy.repeat(y[:, numpy.newaxis], x.shape, 1)
        del x, y
        print('    Conversion of input raster to XMap/YMap completed without error.')
        return xmap, ymap

    def grid_extent(self, arcpy):
        '''
        Return the grid bounding extent [xMin, yMin, xMax, yMax]
        '''
        xMax = self.x00 + (float(self.ncols)*self.DX)
        yMin = self.y00 + (float(self.nrows)*-self.DY)
        extent_obj = arcpy.Extent(self.x00, yMin, xMax, self.y00)
        return extent_obj

    def numpy_to_Raster(self, arcpy, in_arr, value_to_nodata=NoDataVal):
        '''This funciton takes in an input netCDF file, a variable name, the ouput
        raster name, and the projection definition and writes the grid to the output
        raster. This is useful, for example, if you have a FullDom netCDF file and
        the GEOGRID that defines the domain. You can output any of the FullDom variables
        to raster.
        Only use 2D arrays as input.'''

        # Process: Numpy Array to Raster
        nc_raster = arcpy.NumPyArrayToRaster(in_arr, self.point, self.DX, self.DY, value_to_nodata)

        # Process: Define Projection
        arcpy.DefineProjection_management(nc_raster, self.proj)
        return nc_raster

    def domain_shapefile(self, arcpy, in_raster, outputFile):
        """This process creates a shapefile that bounds the GEOGRID file domain. This
        requires the Spatial Analyst extension."""

        printMessages(arcpy, ['Step 2: Build constant raster and convert to polygon...'])

        # Set environments
        arcpy.env.overwriteOutput = True
        arcpy.env.outputCoordinateSystem = self.proj

        # Build constant raster
        RasterLayer = "RasterLayer"
        arcpy.MakeRasterLayer_management(in_raster, RasterLayer)
        outConstRaster = Con(IsNull(RasterLayer) == 1, 0, 0)                    # New method 10/29/2018: Use Con to eliminate NoData cells if there are any.

        # Raster to Polygon conversion
        arcpy.RasterToPolygon_conversion(outConstRaster, outputFile, "NO_SIMPLIFY")    #, "VALUE")

        # Clean up
        arcpy.Delete_management(outConstRaster)
        printMessages(arcpy, ['    Finished building GEOGRID domain boundary shapefile: {0}.'.format(outputFile)])

    def xy_to_grid_ij(self, x, y):
        '''
        This function converts a coordinate in (x,y) to the correct row and column
        on a grid. Code from: https://www.perrygeo.com/python-affine-transforms.html

        Grid indices are 0-based.
        '''
        # x,y to col,row.
        col = int((x - self.x00) / self.DX)
        row = abs(int((y - self.y00) / self.DY))
        return row, col

    def grid_ij_to_xy(self, col, row):
        '''
        This function converts a 2D grid index (i,j) the grid cell center coordinate
        (x,y) in the grid coordinate system.
        Code from: https://www.perrygeo.com/python-affine-transforms.html

        Grid indices are 0-based.
        '''
        # col, row to x, y
        x = (col * self.DX) + self.x00 + self.DX/2.0
        y = (row * -self.DY) + self.y00 - self.DY/2.0
        return x, y

    def project_to_model_grid(self, arcpy, in_raster, out_raster, resampling=ElevResampleMethod):
        '''
        Project an input raster to the output cellsize and extent of the grid object.
        '''

        # Set environments
        arcpy.ResetEnvironments()                                               # Added 2016-02-11
        arcpy.env.outputCoordinateSystem = self.proj
        arcpy.env.extent = self.grid_extent(arcpy)
        arcpy.env.cellSize = self.DX

        # Project the input CRS to the output CRS
        arcpy.ProjectRaster_management(in_raster, out_raster, self.proj, resampling, self.DX, "#", "#", in_raster.spatialReference)

    def getgrid(self, arcpy, x, y):
        '''
        There are nine positions relative to a grid cell center (inclusive).
           [center, left center, right center, bottom center, top center, LL, UL, UR, LR]

        Method with which to create the coordinates relative to a grid cell center
        coordinate in unprojected cartesian space. Use getgrid() to generate a
        set of points on the staggered grid around the center point. Used to
        re-generate corner_lats and corner_lons in a geogrid file.
        '''

        # Store all points in a list
        points = []                                                             # Initiate list to store point geometries

        # Use multipliers of DX/DY to get the relative coordinates based on the center
        multipliersX = [0, -0.5, 0.5, 0, 0, -0.5, -0.5, 0.5, 0.5]               # DX multipliers
        multipliersY = [0, 0, 0, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5]               # DY multipliers
        for multiX, multiY in zip(multipliersX, multipliersY):
            point = arcpy.Point(float(x+(abs(self.DX)*float(multiX))), float(y+(abs(self.DY)*float(multiY))))
            points.append(arcpy.PointGeometry(point, self.proj))
        del point, multipliersX, multipliersY, multiX, multiY
        return points

# --- End Classes --- #

# --- Functions --- #
def printMessages(arcpy, messages):
    '''provide a list of messages and they will be printed and returned as tool messages.'''
    for i in messages:
        print(i)
        arcpy.AddMessage(i)

def makeoutncfile(arcpy, raster, outname, outvar, projdir):
    outfile = os.path.join(projdir, outname)
    arcpy.RasterToNetCDF_md(raster, outfile, outvar, "", "x", "y")
    printMessages(arcpy, ['    Process: {0} completed without error'.format(outname)])
    printMessages(arcpy, ['         Output File: {0}'.format(outfile)])

def ReprojectCoords(arcpy, xcoords, ycoords, src_srs, tgt_srs):
    '''
    Adapted from:
        https://gis.stackexchange.com/questions/57834/how-to-get-raster-corner-coordinates-using-python-gdal-bindings
     Reproject a list of x,y coordinates.
    '''
    tic1 = time.time()

    ravel_x = numpy.ravel(xcoords)
    ravel_y = numpy.ravel(ycoords)
    trans_x = numpy.zeros(ravel_x.shape, ravel_x.dtype)
    trans_y = numpy.zeros(ravel_y.shape, ravel_y.dtype)

    point = arcpy.Point()
    for num, (x, y) in enumerate(zip(ravel_x, ravel_y)):
        point.X = x
        point.Y = y
        pointGeometry = arcpy.PointGeometry(point, src_srs)
        projpoint = pointGeometry.projectAs(tgt_srs)                            # Optionally add transformation method:
        x1 = projpoint.firstPoint.X
        y1 = projpoint.firstPoint.Y
        trans_x[num] = x1
        trans_y[num] = y1

    # reshape transformed coordinate arrays of the same shape as input coordinate arrays
    trans_x = trans_x.reshape(*xcoords.shape)
    trans_y = trans_y.reshape(*ycoords.shape)
    print('Completed transforming coordinate pairs [{0}] in {1: 3.2f} seconds.'.format(num, time.time()-tic1))
    return trans_x, trans_y

def zipws(arcpy, path, zip, keep, nclist):
    path = os.path.normpath(path)
    for dirpath, dirnames, filenames in os.walk(path):
        for file in filenames:
            if file in nclist:
                try:
                    if keep:
                        zip.write(os.path.join(dirpath, file), os.path.join(os.sep + os.path.join(dirpath, file)[len(path) + len(os.sep):]))
                except Exception as e:
                    printMessages(arcpy, ['Exception: {0}'.format(e)])
                    arcpy.AddWarning((file, e))

def zipUpFolder(arcpy, folder, outZipFile, nclist):
    try:
        zip = zipfile.ZipFile(outZipFile, 'w', zipfile.ZIP_DEFLATED, allowZip64=True)
        zipws(arcpy, str(folder), zip, 'CONTENTS_ONLY', nclist)
        zip.close()
    except RuntimeError:
        pass

def flip_grid(array):
    '''This function takes a three dimensional array and flips it up-down to
    correct for the netCDF storage of these grids.'''
    array = array[:, ::-1]                                                     # Flip 2+D grid up-down
    return array

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
                printMessages(arcpy, ['  File Copied: {0}'.format(file)])

            # Ignore the GEOGRID LDASOUT spatial metadata file
            if file == LDASFile:
                shutil.copy2(infile, out_folder)

                printMessages(arcpy, ['  File Copied: {0}'.format(file)])
                continue

            if file.endswith('.nc'):

                # Trap to eliminate Parameter tables in NC format from this extraction
                if file.endswith(LK_nc) or file.endswith(RT_nc) or file.endswith(GW_nc):
                    shutil.copy2(infile, out_folder)
                    printMessages(arcpy, ['  File Created: {0}'.format(file)])
                    del file
                    continue

                # Establish an object for reading the input NetCDF file
                rootgrp = netCDF4.Dataset(infile, 'r')

                # Using netCDF4 library - still need point2, DX, DY
                GT = rootgrp.variables[crsVar].GeoTransform.split(" ")
                if 'esri_pe_string' in rootgrp.variables[crsVar].__dict__:
                    PE_string = rootgrp.variables[crsVar].esri_pe_string
                elif 'spatial_ref' in rootgrp.variables[crsVar].__dict__:
                    PE_string = rootgrp.variables[crsVar].spatial_ref
                printMessages(arcpy, ['  GeoTransform: {0}'.format(GT)])
                DX = float(GT[1])
                DY = abs(float(GT[5]))                                          # In GeoTransform, Y is usually negative
                printMessages(arcpy, ['  DX: {0}'.format(DX)])
                printMessages(arcpy, ['  DY: {0}'.format(DY)])
                sr = arcpy.SpatialReference()
                sr.loadFromString(PE_string.replace('"', "'"))
                point = arcpy.Point(float(GT[0]), float(GT[3]) - float(DY*len(rootgrp.dimensions['y'])))    # Calculate LLCorner value from GeoTransform (ULCorner)
                arcpy.env.outputCoordinateSystem = sr
                for variablename, ncvar in rootgrp.variables.items():
                    if ncvar.dimensions == ('y', 'x'):
                        outRasterLayer = variablename
                        outArr = numpy.array(ncvar[:])
                        nc_raster = arcpy.NumPyArrayToRaster(outArr, point, DX, DY)
                        arcpy.CalculateStatistics_management(nc_raster)
                        arcpy.DefineProjection_management(nc_raster, sr)
                        nc_raster.save(os.path.join(out_folder, outRasterLayer))
                        printMessages(arcpy, ['  File Created: {0}'.format(outRasterLayer)])
                        del nc_raster, variablename, outRasterLayer
                rootgrp.close()
                del file, infile, rootgrp, ncvar
                continue

            if file.endswith('.shp'):
                newshp = os.path.join(out_folder, file)
                arcpy.CopyFeatures_management(infile, newshp)
                printMessages(arcpy, ['  File Created: {0}'.format(str(file))])
                del file, infile, newshp
                continue

            # These file formats are legacy output files, but allow the tool to work with older routing stacks
            if file.endswith('.csv') or file.endswith('.TBL') or file.endswith('.txt') or file.endswith('.prj'):
                # Change was made around 09/2017 to remove ASCII raster header from gw_basins_geogrid.txt ACII raster. Thus, ASCIIToRaster is no longer usable
                shutil.copy2(infile, out_folder)
                printMessages(arcpy, ['  File Created: {0}'.format(file)])
                del file, infile
                continue
    del dirpath, dirnames, filenames

    # Remove each file from the temporary extraction directory
    for infile in dellist:
        os.remove(infile)
    printMessages(arcpy, ['Extraction of WRF routing grids completed.'])
    return out_sfolder

# The following lambda function will read a 2D grid and provide the corner values
corner_coords = lambda x: list((x[0, 0], x[-1, 0], x[-1, -1], x[0, -1]))

# Function to determine a set of corner lats and lons from a single available value
def find_latlon(arcpy, inpoint, grid_obj):
    '''
    Take a grid cell center point and produce a list of all the corners
    coordinates for that point.
    '''
    # Setup the lat/lon point CRS
    insrs = arcpy.SpatialReference(pointSR)

    # Create coordinate transform for this point and transform
    pointGeometry = arcpy.PointGeometry(inpoint, insrs)
    projpoint = pointGeometry.projectAs(grid_obj.proj)                                # Convert to latitude/longitude on the sphere

    # Gather all corner coordinates for this point in projected space
    points = grid_obj.getgrid(arcpy, projpoint.firstPoint.X, projpoint.firstPoint.Y)

    # Conver back to original coordinate system coordinates
    points_proj = []
    for num, pt in enumerate(points):
        pt2 = points[num]
        pt2.projectAs(insrs)
        points_proj.append(pt2)
    del num, pt, pt2, inpoint, points, pointGeometry
    return points_proj                                                          # Return list of out-points

def add_CRS_var(rootgrp, sr, map_pro, CoordSysVarName, grid_mapping, PE_string, GeoTransformStr=None):
    '''
    10/13/2017 (KMS):
        This function was added to generalize the creating of a CF-compliant
        coordinate reference system variable. This was modularized in order to
        create CRS variables for both gridded and point time-series CF-netCDF
        files.
    '''

    # Scalar projection variable - http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/StandardCoordinateTransforms.html
    proj_var = rootgrp.createVariable(CoordSysVarName, 'S1')                    # (Scalar Char variable)
    proj_var.transform_name = grid_mapping                                      # grid_mapping. grid_mapping_name is an alias for this
    proj_var.grid_mapping_name = grid_mapping                                   # for CF compatibility
    proj_var.esri_pe_string = PE_string                                         # For ArcGIS. Not required if esri_pe_string exists in the 2D variable attributes
    proj_var.spatial_ref = PE_string                                            # For GDAl
    proj_var.long_name = "CRS definition"                                       # Added 10/13/2017 by KMS to match GDAL format
    proj_var.longitude_of_prime_meridian = 0.0                                  # Added 10/13/2017 by KMS to match GDAL format
    if GeoTransformStr is not None:
        proj_var.GeoTransform = GeoTransformStr                                 # For GDAl - GeoTransform array

    # Projection specific parameters - http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/StandardCoordinateTransforms.html
    if map_pro == 1:
        # Lambert Conformal Conic

        # Required transform variables
        proj_var._CoordinateAxes = 'y x'                                            # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        proj_var._CoordinateTransformType = "Projection"
        proj_var.standard_parallel = sr.standardParallel1, sr.standardParallel2     # Double
        proj_var.longitude_of_central_meridian = float(sr.centralMeridianInDegrees) # Double. Necessary in combination with longitude_of_prime_meridian?
        proj_var.latitude_of_projection_origin = float(sr.latitudeOfOrigin)         # Double

        # Optional tansform variable attributes
        proj_var.false_easting = float(sr.falseEasting)                         # Double  Always in the units of the x and y projection coordinates
        proj_var.false_northing = float(sr.falseNorthing)                       # Double  Always in the units of the x and y projection coordinates
        proj_var.earth_radius = sphere_radius                                   # OPTIONAL. Parameter not read by Esri. Default CF sphere: 6371.229 km.
        proj_var.semi_major_axis = sphere_radius                                # Added 10/13/2017 by KMS to match GDAL format
        proj_var.inverse_flattening = float(0)                                  # Added 10/13/2017 by KMS to match GDAL format: Double - optional Lambert Conformal Conic parameter

    elif map_pro == 2:
        # Polar Stereographic

        # Required transform variables
        proj_var._CoordinateAxes = 'y x'                                            # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        proj_var._CoordinateTransformType = "Projection"
        proj_var.longitude_of_projection_origin = float(sr.longitudeOfOrigin)   # Double - proj_var.straight_vertical_longitude_from_pole = ''
        proj_var.latitude_of_projection_origin = float(sr.latitudeOfOrigin)     # Double
        proj_var.scale_factor_at_projection_origin = float(sr.scaleFactor)      # Double

        # Optional tansform variable attributes
        proj_var.false_easting = float(sr.falseEasting)                         # Double  Always in the units of the x and y projection coordinates
        proj_var.false_northing = float(sr.falseNorthing)                       # Double  Always in the units of the x and y projection coordinates
        proj_var.earth_radius = sphere_radius                                   # OPTIONAL. Parameter not read by Esri. Default CF sphere: 6371.229 km.
        proj_var.semi_major_axis = sphere_radius                                # Added 10/13/2017 by KMS to match GDAL format
        proj_var.inverse_flattening = float(0)                                  # Added 10/13/2017 by KMS to match GDAL format: Double - optional Lambert Conformal Conic parameter

    elif map_pro == 3:
        # Mercator

        # Required transform variables
        proj_var._CoordinateAxes = 'y x'                                            # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        proj_var._CoordinateTransformType = "Projection"
        proj_var.longitude_of_projection_origin = float(sr.longitudeOfOrigin)   # Double
        proj_var.latitude_of_projection_origin = float(sr.latitudeOfOrigin)     # Double
        proj_var.standard_parallel = float(sr.standardParallel1)                # Double
        proj_var.earth_radius = sphere_radius                                   # OPTIONAL. Parameter not read by Esri. Default CF sphere: 6371.229 km.
        proj_var.semi_major_axis = sphere_radius                                # Added 10/13/2017 by KMS to match GDAL format
        proj_var.inverse_flattening = float(0)                                  # Added 10/13/2017 by KMS to match GDAL format: Double - optional Lambert Conformal Conic parameter

    elif map_pro == 6:
        # Cylindrical Equidistant or rotated pole

        #http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#appendix-grid-mappings
        # Required transform variables
        proj_var.grid_mapping_name = "latitude_longitude"                      # or "rotated_latitude_longitude"

        # No extra parameters needed for latitude_longitude

        # If rotated pole:
        #printMessages(arcpy, ['        Cylindrical Equidistant projection not supported.'])
        #raise SystemExit

    # Added 10/13/2017 by KMS to accomodate alternate datums
    elif map_pro == 0:
        proj_var._CoordinateAxes = 'lat lon'
        proj_var.semi_major_axis = sr.semiMajorAxis
        proj_var.semi_minor_axis = sr.semiMinorAxis
        if sr.flattening != 0:
            proj_var.inverse_flattening = float(float(1)/sr.flattening)         # This avoids a division by 0 error
        else:
            proj_var.inverse_flattening = float(0)

    # Global attributes related to CF-netCDF
    rootgrp.Conventions = CFConv                                                # Maybe 1.0 is enough?
    return rootgrp

def create_CF_NetCDF(arcpy, grid_obj, rootgrp, addLatLon=False, notes='', addVars=[], latArr=None, lonArr=None):
    """This function will create the netCDF file with CF conventions for the grid
    description. Valid output formats are 'GEOGRID', 'ROUTING_GRID', and 'POINT'.
    The output NetCDF will have the XMAP/YMAP created for the x and y variables
    and the LATITUDE and LONGITUDE variables populated from the XLAT_M and XLONG_M
    variables in the GEOGRID file or in the case of the routing grid, populated
    using the getxy function."""

    tic1 = time.time()
    printMessages(arcpy, ['Creating CF-netCDF File.'])

    # Gather projection information from input raster projection
    sr = grid_obj.proj
    dim1size = grid_obj.ncols
    dim2size = grid_obj.nrows
    srType = sr.type
    PE_string = grid_obj.WKT                           # Replace ' with " so Esri can read the PE String properly when running NetCDFtoRaster
    printMessages(arcpy, ['    Esri PE String: {0}'.format(PE_string)])

    # Find name for the grid mapping
    if CF_projdict.get(grid_obj.map_pro) is not None:
        grid_mapping = CF_projdict[grid_obj.map_pro]
        printMessages(arcpy, ['    Map Projection of input raster : {0}'.format(grid_mapping)])
    else:
        #grid_mapping = sr.name
        grid_mapping = 'crs'                                                    # Added 10/13/2017 by KMS to generalize the coordinate system variable names
        printMessages(arcpy, ['    Map Projection of input raster (not a WRF projection): {0}'.format(grid_mapping)])

    # Create Dimensions
    dim_y = rootgrp.createDimension('y', dim2size)
    dim_x = rootgrp.createDimension('x', dim1size)
    printMessages(arcpy, ['    Dimensions created after {0: 3.2f} seconds.'.format(time.time()-tic1)])

    # Create coordinate variables
    var_y = rootgrp.createVariable('y', 'f8', 'y')                              # (64-bit floating point)
    var_x = rootgrp.createVariable('x', 'f8', 'x')                              # (64-bit floating point)

    # Must handle difference between ProjectionCoordinateSystem and LatLonCoordinateSystem
    if srType == 'Geographic':
        if crsVarname:
            CoordSysVarName = crsVar
        else:
            CoordSysVarName = "LatLonCoordinateSystem"

        # Set variable attributes
        #var_y.standard_name = ''
        #var_x.standard_name = ''
        var_y.long_name = "latitude coordinate"
        var_x.long_name = "longitude coordinate"
        var_y.units = "degrees_north"
        var_x.units = "degrees_east"
        var_y._CoordinateAxisType = "Lat"
        var_x._CoordinateAxisType = "Lon"

    elif srType == 'Projected':
        if crsVarname:
            CoordSysVarName = crsVar
        else:
            CoordSysVarName = "ProjectionCoordinateSystem"
        #proj_units = sr.linearUnitName.lower()                                  # sr.projectionName wouldn't work for a GEOGCS
        proj_units = 'm'                                                        # Change made 11/3/2016 by request of NWC

        # Set variable attributes
        var_y.standard_name = 'projection_y_coordinate'
        var_x.standard_name = 'projection_x_coordinate'
        var_y.long_name = 'y coordinate of projection'
        var_x.long_name = 'x coordinate of projection'
        var_y.units = proj_units                                                # was 'meter', now 'm'
        var_x.units = proj_units                                                # was 'meter', now 'm'
        var_y._CoordinateAxisType = "GeoY"                                      # Use GeoX and GeoY for projected coordinate systems only
        var_x._CoordinateAxisType = "GeoX"                                      # Use GeoX and GeoY for projected coordinate systems only
        var_y.resolution = float(grid_obj.DY)                                   # Added 11/3/2016 by request of NWC
        var_x.resolution = float(grid_obj.DX)                                   # Added 11/3/2016 by request of NWC

        # Build coordinate reference system variable
        rootgrp = add_CRS_var(rootgrp, sr, grid_obj.map_pro, CoordSysVarName, grid_mapping, PE_string, grid_obj.GeoTransformStr())

    # For prefilling additional variables and attributes on the same 2D grid, given as a list [[<varname>, <vardtype>, <long_name>],]
    for varinfo in addVars:
        ncvar = rootgrp.createVariable(varinfo[0], varinfo[1], ('y', 'x'))
        ncvar.esri_pe_string = PE_string
        ncvar.grid_mapping = CoordSysVarName
        #ncvar.long_name = varinfo[2]
        #ncvar.units = varinfo[3]

    # Get x and y variables for the netCDF file
    xmaparr, ymaparr = grid_obj.getxy()
    var_y[:] = ymaparr[:, 0]                                                    # Assumes even spacing in y across domain
    var_x[:] = xmaparr[0, :]                                                    # Assumes even spacing in x across domain
    printMessages(arcpy, ['    Coordinate variables and variable attributes set after {0: 3.2f} seconds.'.format(time.time()-tic1)])

    if addLatLon:

        printMessages(arcpy, ['    Proceeding to add LATITUDE and LONGITUDE variables after {0: 3.2f} seconds.'.format(time.time()-tic1)])

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

        # Populate netCDF variables using input numpy arrays
        lat_WRF[:] = latArr
        lon_WRF[:] = lonArr

        printMessages(arcpy, ['    Variables populated after {0: 3.2f} seconds.'.format(time.time()-tic1)])
        printMessages(arcpy, ['    Process completed without error.'])
        printMessages(arcpy, ['    LATITUDE and LONGITUDE variables and variable attributes set after {0: 3.2f} seconds.'.format(time.time()-tic1)])

    # Global attributes
    rootgrp.GDAL_DataType = 'Generic'
    rootgrp.Source_Software = 'WRF-Hydro GIS Pre-processor %s' %PpVersion
    rootgrp.proj4 = grid_obj.proj4                                              # Added 3/16/2018 (KMS) to avoid a warning in WRF-Hydro output
    rootgrp.history = 'Created %s' %time.ctime()
    rootgrp.processing_notes = notes
    printMessages(arcpy, ['    netCDF global attributes set after {0: 3.2f} seconds.'.format(time.time()-tic1)])

    # Return the netCDF file to the calling script
    return rootgrp

def create_high_res_topogaphy(arcpy, in_raster, hgt_m_raster, cellsize, sr2, projdir):
    """
    The second step creates a high resolution topography raster using a hydrologically-
    corrected elevation dataset (currently either HydroSHEDS or NHDPlusv2).

    Changes made 01/26/2018 allow methods to change depending on the ArcGIS version
    available. A bug in ArcGIS 10.4 and 10.5 (Esri BUG-000096495) prevents a Custom
    Geotransformation from being used. Thus, the script determines the ArcGIS version
    and if it is 10.4 or 10.5, will attempt to fake a NULL transformation. However,
    some slight errors in the spatial location of topographic elements has been noted
    for these versions.

    3/28/2018: ArcGIS 10.6 allows custom transformations, but arcpy does not allow
    the custom geotransformation to be specified by name.
    """

    # Get ArcGIS version information and checkout Spatial Analyst extension
    ArcVersion = arcpy.GetInstallInfo()['Version']                              # Get the ArcGIS version that is being used
    ArcProduct = arcpy.GetInstallInfo()['ProductName']
    if ArcVersion.count('.') == 2:
        ArcVersionF = float(ArcVersion.rpartition('.')[0])                      # ArcGIS major release version as a float
    else:
        ArcVersionF = float(ArcVersion)                                         # ArcGIS major release version as a float

    # Second part of the process
    printMessages(arcpy, ['Step 2 initiated...'])

    # Output raster
    mosprj = os.path.join(projdir, 'mosaicprj')

    #Get the extent information from raster object
    arcpy.MakeRasterLayer_management(hgt_m_raster, 'hgt_m_Layer')
    arcpy.env.snapRaster = 'hgt_m_Layer'
    descData = arcpy.Describe('hgt_m_Layer')
    extent = descData.Extent
    boundaryPolygon = extent.polygon                                            # Added 2016-02-11 to simplify code
    extent1 = extent                                                            # Added 2016-02-11 to simplify code

    # Make sure hgt_m is an integer multiple of supplied output resolution
    cellsize1 = descData.children[0].meanCellHeight
    cellsize2 = (cellsize1/cellsize)                                            # Target resolution is the Geogrid resolution divided by the regridding factor
    printMessages(arcpy, ['    The GEOGRID File resolution is {0}sm'.format(str(cellsize1))])
    printMessages(arcpy, ['    The High-resolution dataset will be {0}m'.format(str(cellsize2))])

    # Create a projected boundary polygon of the model domain with which to clip the in_raster
    sr3 = arcpy.Describe(in_raster).spatialReference                            # Obtain the SRS object for the input high-resolution DEM
    if ArcProduct == 'ArcGISPro':
        arcpy.CreateCustomGeoTransformation_management(geoTransfmName, sr2, sr3, customGeoTransfm)
        printMessages(arcpy, ['    Tranformation: {0}'.format(geoTransfmName)])
        projpoly = boundaryPolygon.projectAs(sr3)                               # Reproject the boundary polygon from the WRF domain to the input raster CRS using custom geotransformation
    elif ArcVersionF > 10.3 and ArcVersionF < 10.6:
        # Create two new geographic coordinate systems; one for projecting the model boundary polygon geometry,
        #   and one for projecting the input high resolution DEM (in_raster). In each one, the datum from the
        #   other dataset will be substituted, such that we can perform a 'Null'-like transformation.
        printMessages(arcpy, ['    WARNING: ArcGIS version {0} detected.  Work-around for ESRI BUG-000096495 being implemented. Check output to ensure proper spatial referencing of topographic features in output data.'.format(ArcVersion)])
        modelSR_GCS = sr2.GCS.exportToString().split(';')[0]                    # Find the GCS sub-string for the model SRS
        DEMSR_GCS = sr3.GCS.exportToString().split(';')[0]                      # Find the GCS sub-string for the input high-resolution DEM SRS
        modelSR2 = sr2.exportToString().replace(modelSR_GCS, DEMSR_GCS)         # Replace the GEOGCS portion to 'fake' a known spheroidal datum
        DEMSR2 = sr3.exportToString().replace(DEMSR_GCS, modelSR_GCS)           # Replace the GEOGCS portion to 'fake' a known spherical datum
        sr4 = arcpy.SpatialReference()                                          # Initiate the new SRS
        sr4.loadFromString(modelSR2)                                            # Load SRS from altered string
        sr5 = arcpy.SpatialReference()                                          # Initiate the new SRS
        sr5.loadFromString(DEMSR2)                                              # Load SRS from altered string
        del modelSR_GCS, DEMSR_GCS, modelSR2, DEMSR2
        projpoly = boundaryPolygon.projectAs(sr5)                               # Reproject the boundary polygon from the WRF domain to the input raster CRS
    elif ArcVersionF <= 10.3:
        arcpy.CreateCustomGeoTransformation_management(geoTransfmName, sr2, sr3, customGeoTransfm)
        printMessages(arcpy, ['    Tranformation: {0}'.format(geoTransfmName)])
        projpoly = boundaryPolygon.projectAs(sr3, geoTransfmName)               # Reproject the boundary polygon from the WRF domain to the input raster CRS using custom geotransformation
    elif ArcVersionF >= 10.6:
        # Custom geotransformation appears not to be a valid input to the .projectAs geometry tool in ArcGIS 10.6
        # However, if you create a custom geotransformation beforehand, then it will be respected bye the projectAs function
        arcpy.CreateCustomGeoTransformation_management(geoTransfmName, sr2, sr3, customGeoTransfm)
        projpoly = boundaryPolygon.projectAs(sr3) # Reproject the boundary polygon from the WRF domain to the input raster CRS
    polyextent = projpoly.extent
    del projpoly, boundaryPolygon

    # Create raster layer from input raster or mosaic dataset
    MosaicLayer = "MosaicLayer"
    arcpy.MakeRasterLayer_management(in_raster, MosaicLayer, "#", polyextent)
    printMessages(arcpy, ['    MakeRasterLayer process completed without error.'])
    printMessages(arcpy, ['    The coarse grid has {0} rows and {1} columns.'.format(arcpy.Describe(hgt_m_raster).height, arcpy.Describe(hgt_m_raster).width)])
    printMessages(arcpy, ['    The input elevation grid (before reprojection) has {0} rows and {1} columns.'.format(arcpy.Describe(MosaicLayer).height, arcpy.Describe(MosaicLayer).width)])

    # Set environments to force creation of high-res raster to have exact extent and cellsize needed
    arcpy.env.extent = extent1                                                  # using extent directly doesn't work.
    arcpy.env.outputCoordinateSystem = sr2
    arcpy.env.cellSize = cellsize2
    arcpy.env.snapRaster = hgt_m_raster                                         # Redundant?

    # Alter the method of reprojection depending on the ArcGIS version
    printMessages(arcpy, ['    Projecting input elevation data to WRF coordinate system.'])
    if ArcVersionF > 10.3 and ArcVersionF < 10.6:
        '''If the ArcGIS Version is 10.3. or 10.3.1, then a Null Custom Geotransformation
        may be used to appropriately transform data between sphere and spheroid without
        applying any translation.'''
        printMessages(arcpy, ['    ArcGIS version > 10.3 & < 10.6 found ({0}). No Custom Geotransformation used.'.format(ArcVersion)])
        arcpy.DefineProjection_management(MosaicLayer, sr5)                     # Reset the projection to the model SRS
        arcpy.ProjectRaster_management(MosaicLayer, mosprj, sr4, ElevResampleMethod, cellsize2)
        arcpy.DefineProjection_management(MosaicLayer, sr3)                     # Put it back to the way it was
    elif ArcVersionF <= 10.3 or ArcVersionF >= 10.6 or ArcProduct == 'ArcGISPro':
        printMessages(arcpy, ['    ArcGIS version {0} found. Using Custom Geotransformation ({1})'.format(ArcVersion, geoTransfmName)])
        arcpy.ProjectRaster_management(MosaicLayer, mosprj, sr2, ElevResampleMethod, cellsize2, geoTransfmName)

    printMessages(arcpy, ['    Finished projecting input elevation data to WRF coordinate system.'])
    printMessages(arcpy, ['    The fine grid (before ExtractByMask) has {0} rows and {1} columns.'.format(arcpy.Describe(mosprj).height, arcpy.Describe(mosprj).width)])

    # Extract By Mask
    arcpy.env.cellSize = cellsize2
    mosprj2 = ExtractByMask(mosprj, Con(IsNull('hgt_m_Layer') == 1, 0, 0))          # New method 10/29/2018: Use Con to eliminate NoData cells if there are any.
    #arcpy.Delete_management(mosprj)
    mosprj2.save(mosprj)                                                        # Save this new extracted raster to the same name as before

    # Check that the number of rows and columns are correct
    printMessages(arcpy, ['    Fine Grid has {0} rows and {1} columns.'.format(arcpy.Describe(mosprj2).height, arcpy.Describe(mosprj2).width)])

    # Clean up
    arcpy.Delete_management("MosaicLayer")
    del MosaicLayer

    # Finish
    printMessages(arcpy, ['    Step 2 completed without error.'])
    return mosprj

def getxy(arcpy, inraster, projdir):
    """
    This function will use the affine transformation (GeoTransform) to produce an
    array of X and Y 1D arrays. Note that the GDAL affine transformation provides
    the grid cell coordinates from the upper left corner. This is typical in GIS
    applications. However, WRF uses a south_north ordering, where the arrays are
    written from the bottom to the top. Thus, in order to flip the y array, select
    flipY = True (default).

    5/31/2019:
        This function was altered to reduce dependence on the arcgisscripting module
        single output map algebra function for $$XMAP and $$YMAP found in the getxy
        function, which is deprecated.
    """
    printMessages(arcpy, ['    Starting Process: Converting raster to XMap/YMap'])

    # Setup output rasters
    xmap = os.path.join(projdir, 'xmap2')
    ymap = os.path.join(projdir, 'ymap2')

    # Use the raster to get the boundary and grid information
    descData = arcpy.Describe(inraster)
    extent = descData.Extent
    DX = descData.meanCellWidth
    DY = descData.meanCellHeight

    # Build i,j arrays
    j = numpy.arange(descData.height) + float(0.5)                              # Add 0.5 to estimate coordinate of grid cell centers
    i = numpy.arange(descData.width) + float(0.5)                               # Add 0.5 to estimate coordinate of grid cell centers

    # col, row to x, y   From https://www.perrygeo.com/python-affine-transforms.html
    x = (i * DX) + extent.XMin
    y = (j * -DY) + extent.YMax

    # Create 2D arrays from 1D
    x2 = numpy.repeat(x[numpy.newaxis, :], y.shape, 0)
    y2 = numpy.repeat(y[:, numpy.newaxis], x.shape, 1)

    # Convert back to rasters
    xmap_arr = arcpy.NumPyArrayToRaster(x2, extent.lowerLeft, DX, DY, NoDataVal)
    ymap_arr = arcpy.NumPyArrayToRaster(y2, extent.lowerLeft, DX, DY, NoDataVal)
    xmap_arr.save(xmap)
    ymap_arr.save(ymap)

    # Clean up
    del descData, extent, xmap_arr, ymap_arr, i, j, x, y, x2, y2
    printMessages(arcpy, ['    Conversion of input raster to XMap/YMap completed without error.'])
    return xmap, ymap

def create_lat_lon_rasters(arcpy, projdir, mosprj, wkt_text):
    """The third function in the process is to create the latitude and longitude
    rasters that are necessary for running wrf-hydro.  The latitude and longitude
    that are given by the cell values in the resulting output grids are the latitude
    and longitude of the cell center."""

    # Set buffer distance so that data is still present when converting between Sphere and WGS84
    # Added 2016-02-11 for conversion to WGS84 instead of sphere
    #buffdist = 0.2        # Distance in degrees that ~22000m (maximum sphere-spheroid error) represents (really 0.199)

    # Initiate logging
    printMessages(arcpy, ['    Latitude and Longitude grid generation initiated...'])

    # Create Constant Raster in-memory so we don't have to send the hgt_m layer
    arcpy.env.outputCoordinateSystem = mosprj
    arcpy.env.snapRaster = mosprj
    arcpy.env.extent = mosprj
    arcpy.env.cellSize = mosprj
    OutRas = CreateConstantRaster(1)                                            # Create constant raster with value of 1

    # Set coordinate reference systems
    sr1 = arcpy.SpatialReference()
    sr1.loadFromString(wkt_text)
    sr2 = arcpy.Describe(mosprj).spatialReference

    # Project constant raster to geocentric coordinates system
    projraster = os.path.join(projdir, 'mosprj2')
    arcpy.ResetEnvironments()
    arcpy.env.outputCoordinateSystem = sr1                                      #arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(104128)           # EMEP sphere (same as WRF sphere)
    arcpy.CopyRaster_management(OutRas, projraster)

    # Create xmap/ymap grids
    xmap, ymap = getxy(arcpy, projraster, projdir)
    printMessages(arcpy, ['        XMAP and YMAP rasters created.'])

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
    arcpy.Delete_management(OutRas)

    printMessages(arcpy, ['        Latitude and Longitude NetCDF Files created.'])
    printMessages(arcpy, ['        Latitude and Longitude grid generation completed without error.'])
    return xout2, yout2, xmap, ymap

def build_RouteLink(arcpy, RoutingNC, order, From_To, NodeElev, Arc_To_From, Arc_From_To, NodesLL, NodesXY, Lengths, Straglers, StrOrder, sr, gageDict=None):
    '''
    8/10/2017: This function is designed to build the routiing parameter netCDF file.
                Ideally, this will be the only place the produces the file, and
                all functions wishing to write the file will reference this function.
    '''

    # To create a netCDF parameter file
    rootgrp = netCDF4.Dataset(RoutingNC, 'w', format=outNCType)

    # Create dimensions and set other attribute information
    #dim1 = 'linkDim'
    dim1 = 'feature_id'
    dim2 = 'IDLength'
    dim = rootgrp.createDimension(dim1, len(order))
    gage_id = rootgrp.createDimension(dim2, 15)

    # Create fixed-length variables
    ids = rootgrp.createVariable('link', 'i4', (dim1))                          # Variable (32-bit signed integer)
    froms = rootgrp.createVariable('from', 'i4', (dim1))                        # Variable (32-bit signed integer)
    tos = rootgrp.createVariable('to', 'i4', (dim1))                            # Variable (32-bit signed integer)
    slons = rootgrp.createVariable('lon', 'f4', (dim1))                         # Variable (32-bit floating point)
    slats = rootgrp.createVariable('lat', 'f4', (dim1))                         # Variable (32-bit floating point)
    selevs = rootgrp.createVariable('alt', 'f4', (dim1))                        # Variable (32-bit floating point)
    orders = rootgrp.createVariable('order', 'i4',(dim1))                       # Variable (32-bit signed integer)
    Qis = rootgrp.createVariable('Qi', 'f4', (dim1))                            # Variable (32-bit floating point)
    MusKs = rootgrp.createVariable('MusK', 'f4', (dim1))                        # Variable (32-bit floating point)
    MusXs = rootgrp.createVariable('MusX', 'f4', (dim1))                        # Variable (32-bit floating point)
    Lengthsnc = rootgrp.createVariable('Length', 'f4', (dim1))                  # Variable (32-bit floating point)
    ns = rootgrp.createVariable('n', 'f4', (dim1))                              # Variable (32-bit floating point)
    Sos = rootgrp.createVariable('So', 'f4', (dim1))                            # Variable (32-bit floating point)
    ChSlps = rootgrp.createVariable('ChSlp', 'f4', (dim1))                      # Variable (32-bit floating point)
    BtmWdths = rootgrp.createVariable('BtmWdth', 'f4', (dim1))                  # Variable (32-bit floating point)
    Times = rootgrp.createVariable('time', 'f4')                                # Scalar Variable (32-bit floating point)
    geo_x = rootgrp.createVariable('x', 'f4', (dim1))                           # Variable (32-bit floating point)
    geo_y = rootgrp.createVariable('y', 'f4', (dim1))                           # Variable (32-bit floating point)
    Kcs = rootgrp.createVariable('Kchan', 'i2', (dim1))                         # Variable (16-bit signed integer)
    Gages = rootgrp.createVariable('gages', 'S1', (dim1, dim2))                 # Variable (string type character) Added 07/27/2015 - 15 character strings
    LakeDis = rootgrp.createVariable('NHDWaterbodyComID', 'i4', (dim1))         # Variable (32-bit signed integer)

    # Add CF-compliant coordinate system variable
    if pointCF:
        sr = arcpy.SpatialReference(pointSR)                                    # Build a spatial reference object
        PE_string = sr.exportToString().replace("'", '"')                       # Replace ' with " so Esri can read the PE String properly when running NetCDFtoRaster
        grid_mapping = crsVar
        rootgrp = add_CRS_var(rootgrp, sr, 0, grid_mapping, 'latitude_longitude', PE_string)

    # Set variable descriptions
    ids.long_name = 'Link ID'
    froms.long_name = 'From Link ID'
    tos.long_name = 'To Link ID'
    slons.long_name = 'longitude of the start node'
    slats.long_name = 'latitude of the start node'
    selevs.long_name = 'Elevation in meters at start node'
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
    LakeDis.long_name = 'ID of the lake element that intersects this flowline'
    Gages.long_name = 'Gage ID'

    # Variable attributes for CF compliance
    slons.units = 'degrees_east'                                                # For compliance
    slats.units = 'degrees_north'                                               # For compliance
    slons.standard_name = 'longitude'                                           # For compliance
    slats.standard_name = 'latitude'                                            # For compliance
    Times.standard_name = 'time'                                                # For compliance
    Times.long_name = 'time of measurement'                                     # For compliance
    Times.units = 'days since 2000-01-01 00:00:00'                              # For compliance
    selevs.standard_name = "height"                                             # For compliance
    selevs.units = "m"                                                          # For compliance
    selevs.positive = "up"                                                      # For compliance
    selevs.axis = "Z"                                                           # For compliance
    ids.cf_role = "timeseries_id"                                               # For compliance
    geo_x.standard_name = "projection_x_coordinate"
    geo_y.standard_name = "projection_y_coordinate"
    geo_x.units = "m"
    geo_y.units = "m"
    Kcs.units = "mm h-2"
    slons.standard_name = 'longitude'                                           # For compliance with NCO
    slats.standard_name = 'latitude'                                            # For compliance with NCO

    # Apply grid_mapping and coordinates attributes to all variables
    for varname, ncVar in rootgrp.variables.items():
        if dim1 in ncVar.dimensions and varname not in ['alt', 'lat', 'lon', 'x', 'y']:
            ncVar.setncattr('coordinates', 'lat lon')                           # For CF-compliance
            if pointCF:
                ncVar.setncattr('grid_mapping', grid_mapping)                       # For CF-compliance
        del ncVar, varname

    # Fill in global attributes
    rootgrp.featureType = 'timeSeries'                                          # For compliance
    rootgrp.history = 'Created %s' %time.ctime()

    printMessages(arcpy, ['        Starting to fill in routing table NC file.'])
    ids[:] = numpy.array(order)                                                             # Fill in id field information
    fromnodes = [From_To[arcid][0] for arcid in order]                                      # The FROM node from the streams shapefile (used later as a key))
    tonodes = [From_To[arcid][1] for arcid in order]                                        # The TO node from the streams shapefile (used later as a key)
    drops = [int(NodeElev[fromnode] or 0)-int(NodeElev[tonode] or 0) for fromnode, tonode in zip(fromnodes, tonodes)]	# Fix issues related to None in NodeElev
    #drops = [NodeElev[fromnode]-NodeElev[tonode] for fromnode, tonode in zip(fromnodes, tonodes)]
    drops = [x if x > 0 else 0 for x in drops]                                                # Replace negative values with 0

    # Set variable value arrays
    fromlist = [Arc_To_From[arcid] for arcid in order]                                      # List containes 'None' values, must convert to numpy.nan
    tolist = [Arc_From_To[arcid] for arcid in order]                                        # List containes 'None' values, must convert to numpy.nan

    # Change None values to 0.  Could alternatively use numpy.nan
    froms[:] = numpy.array([0 if not x else x for x in fromlist])                    # Note that the From in this case is the ARCID of any of the immediately upstream contributing segments
    tos[:] = numpy.array([0 if not x else x for x in tolist])                        # Note that the To in this case is the ARCID of the immediately downstream segment

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
    Sos_ = numpy.round(numpy.array(drops).astype(float)/Lengthsnc[:], 3)        # Must convert list to float to result in floats
    numpy.place(Sos_, Sos_ == 0, [0.005])                                         # Set minimum slope to be 0.005
    Sos[:] = Sos_[:]

    # Set default arrays
    Qis[:] = Qi
    MusKs[:] = MusK
    MusXs[:] = MusX
    ns[:] = n
    ChSlps[:] = ChSlp
    BtmWdths[:] = BtmWdth
    Times[:] = 0
    Kcs[:] = Kc

    # Added 10/10/2017 by KMS to include user-supplied gages in reach-based routing files
    if gageDict is not None:
        Gages[:, :] = numpy.asarray([tuple(str(gageDict[arcid]).rjust(15)) if arcid in gageDict else tuple('               ') for arcid in order])
    else:
        Gages[:, :] = numpy.asarray([tuple('               ') for arcid in order])    # asarray converts from tuple to array
    del gageDict

    # Close file
    rootgrp.close()
    printMessages(arcpy, ['        Done writing NC file to disk.'])
    printMessages(arcpy, ['    Routing table created without error.'])

def assign_lake_IDs(arcpy, in_lakes, lakeIDfield=None):
    '''
    This function will either assign a new, unique lake ID field to an existing
    lake feature layer, or return an existing (provided) lake ID field if it
    exists in the input feature class.
    '''
    tic1 = time.time()

    if lakeIDfield is None:
        # Renumber lakes 1-n (addfield, calculate field)
        printMessages(arcpy, ['    Adding auto-incremented lake ID field (1...n)'])
        lakeID = defaultLakeID
        arcpy.AddField_management(in_lakes, lakeID, "LONG")
        expression = 'autoIncrement()'
        code_block = """rec = 0\ndef autoIncrement():\n    global rec\n    pStart = 1\n    pInterval = 1\n    if (rec == 0):\n        rec = pStart\n    else:\n        rec = rec + pInterval\n    return rec"""
        arcpy.CalculateField_management(in_lakes, lakeID, expression, "PYTHON_9.3", code_block)
    else:
        printMessages(arcpy, ['    Using provided lake ID field: {0}'.format(lakeIDfield)])
        lakeID = lakeIDfield                                                    # Use existing field specified by 'lakeIDfield' parameter
    printMessages(arcpy, ['  Checked lake input shapefile for unique ID in {0: 3.2f} seconds.'.format(time.time()-tic1)])
    return lakeID

def reaches_with_lakes(arcpy, FL, WB, outDir, ToSeg, sorted_Flowlinearr, FLID, in_raster):
    '''
    8/7/2017: This function will add lakes to an NCAR Reach-based routing configuration.

    1) This function is designed to be called from within the function that is
    generating the reach-based routing files. This is a design decision to allow
    the lakes to be tested against known incompatiblities with the network in WRF-Hydro
    and to be included in the reach-based routing files.
    '''
    tic1 = time.time()

    # Get information about the input domain and set environments
    descData = arcpy.Describe(in_raster)
    extent = descData.Extent                                                    # Overkill?
    arcpy.env.snapRaster = in_raster
    cellsize = descData.children[0].meanCellHeight
    arcpy.env.extent = extent
    arcpy.env.cellSize = cellsize
    arcpy.env.outputCoordinateSystem = descData.spatialReference

    # Use extent of channelgrid raster to add a feature layer of lake polygons
    outshp = os.path.join(outDir, TempLakeFIle)
    arcpy.CopyFeatures_management(WB, outshp)
    lakeID = assign_lake_IDs(arcpy, outshp, lakeIDfield=defaultLakeID)

    # Set some initial globals for the Lake Preprocessing steps
    WBFields = (FLID, lakeID)                                                  # If FC is NHDPlus-derived
    LK_Routing.FLID = FLID                                                      # Try replacing a global value in other script
    LK_Routing.save_Lake_Link_Type_arr = True                                  # Tell the script not to save the lake link type array since we will save it in this script
    LK_Routing.LakeAssoc = lakeID

    #arcpy.AddField_management(FL, lakeID, "SHORT")
    dtypes = numpy.dtype([(FLID, 'i4'), (LK_Routing.hydroSeq, 'i4')])           # Create a numpy dtype object
    FLarr = numpy.empty(len(sorted_Flowlinearr), dtype=dtypes)
    FLarr[FLID] = sorted_Flowlinearr

    # Reverse this range since the sorted array is sorted from uptream to downstream.
    FLarr[LK_Routing.hydroSeq] = numpy.arange(len(sorted_Flowlinearr))[::-1]

    # Runt he main lake pre-processor
    WaterbodyDict, Lake_Link_Type_arr, Old_New_LakeComID = LK_Routing.main(arcpy,
        outDir, FL, outshp, ToSeg, FLarr, WBFields, Subset_arr=None, NJ=True)
    print('    size of lake link type array: %s' %Lake_Link_Type_arr.shape[0])
    print('    size of waterbody dict: %s' %len(WaterbodyDict))
    del Old_New_LakeComID

    printMessages(arcpy, ['    Lakes added to reach-based routing files after {0: 3.2f} seconds.'.format(time.time()-tic1)])
    return WaterbodyDict, Lake_Link_Type_arr

def Routing_Table(arcpy, projdir, sr2, channelgrid, fdir, Elev, Strahler, gages=None, Lakes=None):
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
            while stack:
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

    printMessages(arcpy, ['    Routing table will be created...'])

    # Get grid information from channelgrid
    descdata = arcpy.Describe(channelgrid)
    extent = descdata.extent
    cellsizeY = descdata.meanCellHeight
    cellsizeX = descdata.meanCellWidth
    del descdata

    # Set output coordinate system environment
    sr1 = arcpy.SpatialReference()
    sr1.loadFromString(wkt_text)                                                # Load default sphere lat/lon CRS from global attribute wkt_text
    arcpy.env.outputCoordinateSystem = sr2

    # Output files
    outStreams = os.path.join(projdir, StreamSHP)
    RoutingNC = os.path.join(projdir, RT_nc)
    OutFC = os.path.join("in_memory", "Nodes")
    OutFC2 = os.path.join("in_memory", "NodeElev")
    OutFC3 = os.path.join("in_memory", "NodeOrder")
    outRaster = os.path.join("in_memory", "LINKID")

    # Build Stream Features shapefile
    StreamToFeature(channelgrid, fdir, outStreams, "NO_SIMPLIFY")
    printMessages(arcpy, ['        Stream to features step complete.'])

    # Set environments based on channelgrid
    arcpy.env.snapRaster = channelgrid
    arcpy.env.extent = channelgrid

    # Create a raster based on the feature IDs that were created in StreamToFeature
    arcpy.FeatureToRaster_conversion(outStreams, 'ARCID', outRaster, channelgrid)                   # Must do this to get "ARCID" field into the raster
    maxValue = arcpy.SearchCursor(outStreams, "", "", "", 'ARCID' + " D").next().getValue('ARCID')  # Gather highest "ARCID" value from field of segment IDs
    maxRasterValue = arcpy.GetRasterProperties_management(outRaster, "MAXIMUM")                     # Gather maximum "ARCID" value from raster
    if int(maxRasterValue[0]) > maxValue:
        printMessages(arcpy, ['        Setting linkid values of {0} to Null.'.format(maxRasterValue[0])])
        whereClause = "VALUE = %s" %int(maxRasterValue[0])
        outRaster = SetNull(outRaster, outRaster, whereClause)                  # This should eliminate the discrepency between the numbers of features in outStreams and outRaster
    outRaster = Con(IsNull(outRaster) == 1, NoDataVal, outRaster)               # Set Null values to -9999 in LINKID raster

    # Add gage points from input forecast points file
    frxst_linkID = {}                                                           # Create blank dictionary so that it exists and can be deleted later
    if gages is not None:
        printMessages(arcpy, ['        Adding forecast points:LINKID association.'])

        # Input forecast points raster must be forecast point IDs and NoData only
        out_frxst_linkIDs = os.path.join('in_memory', 'frxst_linkIDs')

        # Sample the LINKID value for each forecast point. Result is a table
        Sample(outRaster, gages, out_frxst_linkIDs, "NEAREST")
        frxst_linkID = {int(row[-1]):int(row[1]) for row in arcpy.da.SearchCursor(out_frxst_linkIDs, '*')}  # Dictionary of LINKID:forecast point for all forecast points

        # Clean up
        arcpy.Delete_management(gages)
        arcpy.Delete_management(out_frxst_linkIDs)
        printMessages(arcpy, ['        Found {0} forecast point:LINKID associations.'.format(len(frxst_linkID))])
        del gages, out_frxst_linkIDs

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

            if not feat.isMultipart:                                            # Make sure that each part is a single part feature
                firstpoint = feat.firstPoint                                    # First point feature geometry
                lastpoint = feat.lastPoint                                      # Last point feature geometry
                Lengths[ID] = feat.length                                       # Store length of the stream segment

                # Gather the X and Y of the top and bottom ends
                for i in firstpoint, lastpoint:                                  # Now put this geometry into new feature class
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
                printMessages(arcpy, ['This is a multipart line feature and cannot be handled.'])
    #rows.reset()                                                                # Not sure why this is necessary
    del row, rows, point, ID, feat, firstpoint, lastpoint, pointGeometry, projpoint, projpoint1, i, sr2
    printMessages(arcpy, ['        Done reading streams layer.'])

    # Make a point feature class out of the nodes
    IC = arcpy.da.InsertCursor(OutFC, ['SHAPE@', 'NODE'])
    NodesXY = Nodes                                                             # Copy the Nodes dictionary before it gets modified
    for node in list(Nodes.keys()):                                                   # Now we have to adjust the points that fall outside of the raster edge

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
    printMessages(arcpy, ['        Done building Nodes layer with adjustments.'])

    # Get the elevation values for the nodes feature class
    ExtractValuesToPoints(OutFC, Elev, OutFC2, "NONE", "VALUE_ONLY")
    printMessages(arcpy, ['        Done extracting elevations to points.'])
    NodeElev = {row[0]: row[1] for row in arcpy.da.SearchCursor(OutFC2, ['NODE', 'RASTERVALU'])}
    printMessages(arcpy, ['        Done reading node elevations.'])
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
    printMessages(arcpy, ['        Done reading Strahler stream orders.'])
    arcpy.Delete_management(OutFC)                                              # Clean up
    arcpy.Delete_management(OutFC3)                                             # Clean up

    # Add stream order into the streams shapefile
    arcpy.AddField_management(outStreams, "Order_", "SHORT")                    # Add field for "Order_"
    arcpy.AddField_management(outStreams, "index", "LONG")                      # Add field for "index" that gives the topologically sorted order of streams in the out nc file
    arcpy.AddField_management(outStreams, "GageID", "TEXT", "#", "#", 15, "#", "NULLABLE")  # Add field for gages
    with arcpy.da.UpdateCursor(outStreams, ['ARCID', 'Order_', 'GageID']) as rows:# Start UpdateCursor to add the stream order information
        for row in rows:
            row[1] = StrOrder[From_To[row[0]][0]]                               # This gets the stream order of the upstream node
            if row[0] in frxst_linkID:
                row[2] = frxst_linkID[row[0]]
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
    printMessages(arcpy, ['        Time elapsed for sorting: {0:3.2f}s'.format(time.time()-tic2)])

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

    # Make call to Lake module
    if Lakes is not None:
        ToSeg = {arcid:From_To[arcid][1] for arcid in order}
        reaches_with_lakes(arcpy, outStreams, Lakes, projdir, ToSeg, order, "ARCID", channelgrid)

    # Call function to build the netCDF parameter table
    build_RouteLink(arcpy, RoutingNC, order, From_To, NodeElev, Arc_To_From, Arc_From_To, NodesLL, NodesXY, Lengths, Straglers, StrOrder, sr1, gageDict=frxst_linkID)
    del frxst_linkID, sr1
    return outRaster

def build_LAKEPARM(arcpy, LakeNC, min_elevs, areas, max_elevs, OrificEs, cen_lats, cen_lons, WeirE_vals):
    '''
    8/10/2017: This function is designed to build the lake parameter netCDF file.
                Ideally, this will be the only place the produces the file, and
                all functions wishing to write the file will reference this function.
    '''
    min_elev_keys = list(min_elevs.keys())                                      # 5/31/2019: Supporting Python3


    # Create Lake parameter file
    printMessages(arcpy, ['    Starting to create lake parameter table.'])
    printMessages(arcpy, ['        Lakes Table: %s Lakes' %len(list(areas.keys()))])

    # Create NetCDF output table
    rootgrp = netCDF4.Dataset(LakeNC, 'w', format=outNCType)

    # Create dimensions and set other attribute information
    #dim1 = 'nlakes'
    dim1 = 'feature_id'
    dim = rootgrp.createDimension(dim1, len(min_elevs))

    # Create coordinate variables
    ids = rootgrp.createVariable('lake_id', 'i4', (dim1))                       # Variable (32-bit signed integer)
    ids[:] = numpy.array(min_elev_keys)                                         # Variable (32-bit signed integer)

    # Create fixed-length variables
    LkAreas = rootgrp.createVariable('LkArea', 'f8', (dim1))                    # Variable (64-bit floating point)
    LkMxEs = rootgrp.createVariable('LkMxE', 'f8', (dim1))                      # Variable (64-bit floating point)
    WeirCs = rootgrp.createVariable('WeirC', 'f8', (dim1))                      # Variable (64-bit floating point)
    WeirLs = rootgrp.createVariable('WeirL', 'f8', (dim1))                      # Variable (64-bit floating point)
    OrificeCs = rootgrp.createVariable('OrificeC', 'f8', (dim1))                # Variable (64-bit floating point)
    OrificeAs = rootgrp.createVariable('OrificeA', 'f8', (dim1))                # Variable (64-bit floating point)
    OrificeEs = rootgrp.createVariable('OrificeE', 'f8', (dim1))                # Variable (64-bit floating point)
    lats = rootgrp.createVariable('lat', 'f4', (dim1))                          # Variable (32-bit floating point)
    longs = rootgrp.createVariable('lon', 'f4', (dim1))                         # Variable (32-bit floating point)
    Times = rootgrp.createVariable('time', 'f8', (dim1))                        # Variable (64-bit floating point)
    WeirEs = rootgrp.createVariable('WeirE', 'f8', (dim1))                      # Variable (64-bit floating point)
    AscendOrder = rootgrp.createVariable('ascendingIndex', 'i4', (dim1))        # Variable (32-bit signed integer)
    ifd = rootgrp.createVariable('ifd', 'f4', (dim1))                           # Variable (32-bit floating point)

    # Add CF-compliant coordinate system variable
    if pointCF:
        sr = arcpy.SpatialReference(pointSR)                                    # Build a spatial reference object
        PE_string = sr.exportToString().replace("'", '"')                       # Replace ' with " so Esri can read the PE String properly when running NetCDFtoRaster
        grid_mapping = crsVar
        rootgrp = add_CRS_var(rootgrp, sr, 0, grid_mapping, 'latitude_longitude', PE_string)

    # Set variable descriptions
    ids.long_name = 'Lake ID'
    LkAreas.long_name = 'Lake area (sq. km)'
    LkMxEs.long_name = 'Maximum lake elevation (m ASL)'
    WeirCs.long_name = 'Weir coefficient'
    WeirLs.long_name = 'Weir length (m)'
    OrificeCs.long_name = 'Orifice coefficient'
    OrificeAs.long_name = 'Orifice cross-sectional area (sq. m)'
    OrificeEs.long_name = 'Orifice elevation (m ASL)'
    WeirEs.long_name = 'Weir elevation (m ASL)'
    lats.long_name = 'latitude of the lake centroid'
    longs.long_name = 'longitude of the lake centroid'
    AscendOrder.long_name = 'Index to use for sorting IDs (ascending)'
    ifd.long_name = 'Initial fraction water depth'
    longs.units = 'degrees_east'                                                # For compliance
    lats.units = 'degrees_north'                                                # For compliance
    longs.standard_name = 'longitude'                                           # For compliance
    lats.standard_name = 'latitude'                                             # For compliance
    Times.standard_name = 'time'                                                # For compliance
    Times.long_name = 'time of measurement'                                     # For compliance
    Times.units = 'days since 2000-01-01 00:00:00'                              # For compliance. Reference time arbitrary
    WeirEs.units = 'm'
    ids.cf_role = "timeseries_id"                                               # For compliance

    # Apply grid_mapping and coordinates attributes to all variables
    for varname, ncVar in rootgrp.variables.items():
        if dim1 in ncVar.dimensions and varname not in ['alt', 'lat', 'lon', 'x', 'y']:
            ncVar.setncattr('coordinates', 'lat lon')                           # For CF-compliance
            if pointCF:
                ncVar.setncattr('grid_mapping', grid_mapping)                       # For CF-compliance
        del ncVar, varname

    # Fill in global attributes
    rootgrp.featureType = 'timeSeries'                                          # For compliance
    rootgrp.history = 'Created %s' %time.ctime()

    printMessages(arcpy, ['        Starting to fill in lake parameter table NC file.'])
    AscendOrder[:] = numpy.argsort(ids[:])                                  # Use argsort to give the ascending sort order for IDs. Added by KMS 4/4/2017
    LkAreas[:] = numpy.array([float(areas[lkid])/float(1000000) for lkid in min_elev_keys])  # Divide by 1M for kilometers^2
    LkMxEs[:] = numpy.array([max_elevs[lkid] for lkid in min_elev_keys])
    WeirCs[:] = WeirC
    WeirLs[:] = WeirL
    OrificeCs[:] = OrificeC
    OrificeAs[:] = OrificA
    Times[:] = 0
    OrificeEs[:] = numpy.array([OrificEs[lkid] for lkid in min_elev_keys])   # Orifice elevation is 1/3 between 'min' and max lake elevation.
    lats[:] = numpy.array([cen_lats[lkid] for lkid in min_elev_keys])
    longs[:] = numpy.array([cen_lons[lkid] for lkid in min_elev_keys])
    WeirEs[:] = numpy.array([WeirE_vals[lkid] for lkid in min_elev_keys])    # WierH is 0.9 of the distance between the low elevation and max lake elevation
    ifd[:] = ifd_Val

    # Close file
    rootgrp.close()
    printMessages(arcpy, ['        Done writing {0} table to disk.'.format(LK_nc)])
    return

def build_LAKEPARM_ascii(LakeTBL, min_elevs, areas, max_elevs, OrificEs, cen_lats, cen_lons, WeirH_vals):

    '''
    Function to build LAKEPARM.TBL ascii-format file.
    Not currently called from any other function.
    '''

    # Create .TBL output table (ASCII)
    with open(LakeTBL, 'wb') as fp:
        a = csv.writer(fp, dialect='excel-tab', quoting=csv.QUOTE_NONE)
        #a.writerow(['lake', 'LkArea', 'LkMxH', 'WeirC', 'WeirL', 'OrificeC', 'OrificeA', 'OrificeE', 'lat', 'long', 'elevation', 'WeirH']) #
        for lkid in list(min_elevs.keys()):
            lkarea = float(areas[lkid])/float(1000000)                          # Divide by 1M for kilometers^2
            lkmaxelev = max_elevs[lkid]
            OrificeE = OrificEs[lkid]                                           # Orifice Elevation is 1/3 between 'min' and max lake elevation.
            cen_lat = cen_lats[lkid]
            cen_lon = cen_lons[lkid]
            WeirH = WeirH_vals[lkid]
            a.writerow([lkid, lkarea, lkmaxelev, WeirC, WeirL, OrificeC, OrificA, OrificeE, cen_lat, cen_lon, ifd_Val, WeirH])   #COMID?

def add_reservoirs(arcpy, channelgrid, in_lakes, flac, projdir, fill2, cellsize, sr2, lakeIDfield=None):
    """
    This function is intended to add reservoirs into the model grid stack, such
    that the channelgrid and lake grids are modified to accomodate reservoirs and
    lakes.

    This version does not attempt to subset the lakes by a size threshold, nor
    does it filter based on FTYPE.

    2/23/2018:
        Change made to how AREA paramter is calculated in LAKEPARM. Previously, it
        was based on the gridded lake area. Now it is based on the AREASQKM field
        in the input shapefile. This change was made because in NWM, the lakes
        are represented as objects, and are not resolved on a grid.

    """

    tic1 = time.time()                                                          # Set timer
    printMessages(arcpy, ['  Adding reservoirs to routing stack.'])

    # Get information about the input domain and set environments
    arcpy.env.cellSize = cellsize
    arcpy.env.extent = channelgrid
    arcpy.env.outputCoordinateSystem = sr2
    arcpy.env.snapRaster = channelgrid

    # Use extent of channelgrid raster to add a feature layer of lake polygons
    # If reach-based routing (with lakes) is selected, this shapefile may have already been created
    outshp = os.path.join(projdir, TempLakeFIle)
    if not os.path.exists(outshp):
        arcpy.CopyFeatures_management(in_lakes, outshp)
        lakeID = assign_lake_IDs(arcpy, outshp, lakeIDfield=lakeIDfield)
    elif lakeIDfield is None:
        lakeID = defaultLakeID

    # Create new area field
    Field1 = 'AREASQKM'
    if Field1 not in [field.name for field in arcpy.ListFields(outshp)]:
        arcpy.AddField_management(outshp, Field1, "FLOAT")
        arcpy.CalculateField_management(outshp, Field1, '!shape.area@squarekilometers!', "PYTHON_9.3")
    arcpy.MakeFeatureLayer_management(outshp, "Lakeslyr")

    # Gather areas from AREASQKM field (2/23/2018 altered in order to provide non-gridded areas)
    areas = {row[0]: row[1]*1000000 for row in arcpy.da.SearchCursor("Lakeslyr", [lakeID, Field1])}     # Convert to square meters
    lakeIDList = list(areas.keys())                                             # 2/23/2018: Added to find which lakes go missing after resolving on the grid

    # Create a raster from the lake polygons that matches the channelgrid layer
    outRastername = os.path.join(projdir, "Lakesras")
    outfeatures = os.path.join(projdir, LakesShp)
    arcpy.CopyFeatures_management("Lakeslyr", outfeatures)
    arcpy.PolygonToRaster_conversion(outfeatures, lakeID, outRastername, "MAXIMUM_AREA")    # This tool requires ArcGIS for Desktop Advanced OR Spatial Analyst Extension
    #arcpy.FeatureToRaster_conversion(outfeatures, lakeID, outRastername)       # This tool requires only ArcGIS for Desktop Basic, but does not allow a priority field

    # Hack to convert Lakesras to 16bit integer
    outRaster1 = (channelgrid * 0) + Raster(outRastername)                      # Create 16bit ratser for Lakesras out of channelgrid
    outRaster = Con(IsNull(outRaster1) == 1, NoDataVal, outRaster1)             # Convert Null or NoData to -9999

    # Pull flow accumulation values masked to lake polygons
    zonstat = ZonalStatistics(outRastername, "Value", flac, "MAXIMUM")          # Get maximum flow accumulation value for each lake
    lakeacc = Con(outRastername, flac)                                          # Flow accumulation over lakes only
    TestCon = Con(lakeacc == zonstat, outRastername, NoDataVal)                 # Bottom of lake channelgrid pixel = lake number
    NewChannelgrid = Con(IsNull(TestCon) == 1, channelgrid, TestCon)            # This is the new lake-routed channelgrid layer

    # Now march down a set number of pixels to get minimum lake elevation
    tolerance = int(float(arcpy.GetRasterProperties_management(channelgrid, 'CELLSIZEX').getOutput(0)) * LK_walker)
    SnapPour = SnapPourPoint(SetNull(TestCon, TestCon, "VALUE = %s" %NoDataVal), flac, tolerance)   # Snap lake outlet pixel to FlowAccumulation with tolerance
    del flac
    outTable2 = r'in_memory/zonalstat2'
    Sample(fill2, SnapPour, outTable2, "NEAREST")
    min_elevs = {row[1]: row[-1] for row in arcpy.da.SearchCursor(outTable2, "*")}

    # Convert lakes raster to polygons and gather size & elevations
    printMessages(arcpy, ['    Gathering lake parameter information.'])
    outTable = r'in_memory/zonalstat'
    zontable = ZonalStatisticsAsTable(outRastername, 'VALUE', fill2, outTable, "DATA", "MIN_MAX_MEAN")
    max_elevs = {row[0]: row[1] for row in arcpy.da.SearchCursor(outTable, ['VALUE', 'MAX'])}                   # Searchcursor on zonal stats table

    # 2/23/2018: Find the missing lakes and sample elevation at their true centroid.
    min_elev_keys = list(min_elevs.keys())                                      # 5/31/2019: Supporting Python3
    printMessages(arcpy, ['    Lakes in minimum elevation dict: {0}'.format(min_elev_keys)])
    MissingLks = [item for item in lakeIDList if item not in min_elev_keys]  # 2/23/2018: Find lakes that were not resolved on the grid
    shapes = {}
    if len(MissingLks) > 0:
        printMessages(arcpy, ['    Found {0} lakes that could not be resolved on the grid: {1}\n      Sampling elevation from the centroid of these features.'.format(len(MissingLks), str(MissingLks))])
        arcpy.SelectLayerByAttribute_management("Lakeslyr", "NEW_SELECTION", '"%s" IN (%s)' %(lakeID, str(MissingLks)[1:-1]))  # Select the missing lakes from the input shapefile
        FID_to_ID = {row[0]: row[1] for row in arcpy.da.SearchCursor("Lakeslyr", ['OID@', lakeID])}
        centroids = os.path.join('in_memory', 'lake_centroids')                     # Input lake centroid points
        out_centroids = os.path.join('in_memory', 'lake_elevs')                     # Input lake centroid points
        arcpy.FeatureToPoint_management("Lakeslyr", centroids, "INSIDE")            # Convert polygons to points inside each polygon
        Sample(fill2, centroids, out_centroids, "NEAREST", 'ORIG_FID')              # Sample the elevation value for each lake centroid point. Result is a table
        centroidElev = {FID_to_ID[row[1]]: row[-1] for row in arcpy.da.SearchCursor(out_centroids, '*')}   # Dictionary of LakeID:Centroid Elevatoin point for all forecast points
        shapes.update({row[0]: row[1] for row in arcpy.da.SearchCursor(centroids, [lakeID, 'SHAPE@XY'])})		# Get the XY values to add to attributes later
        max_elevs.update(centroidElev)                                              # Add single elevation value as max elevation
        min_elevs.update(centroidElev)                                              # Make these lakes the minimum depth
        arcpy.Delete_management(out_centroids)
        arcpy.Delete_management(centroids)
        del centroidElev, MissingLks

    # Give a minimum active lake depth to all lakes with no elevation variation
    elevRange = {key:max_elevs[key]-val for key, val in min_elevs.items()}      # Get lake depths
    noDepthLks = {key:val for key, val in elevRange.items() if val < minDepth}  # Make a dictionary of these lakes
    if len(noDepthLks) > 0:
        printMessages(arcpy, ['    Found {0} lakes with no elevation range. Providing minimum depth of {1}m for these lakes.'.format(len(noDepthLks), minDepth)])
        min_elevs.update({key:max_elevs[key]-minDepth for key, val in noDepthLks.items() if val == 0}) # Give these lakes a minimum depth
        noDepthFile = os.path.join(projdir, 'Lakes_with_minimum_depth.csv')
        with open(noDepthFile, 'wb') as f:
            w = csv.writer(f)
            w.writerows(list(noDepthLks.items()))
            del noDepthFile
    del elevRange, noDepthLks

    # Calculate the Orifice and Wier heights
    OrificEs = {x:(min_elevs[x] + ((max_elevs[x] - min_elevs[x])/3)) for x in min_elev_keys}             # Orific elevation is 1/3 between the low elevation and max lake elevation
    WeirE_vals = {x:(min_elevs[x] + ((max_elevs[x] - min_elevs[x]) * 0.9)) for x in min_elev_keys}       # WierH is 0.9 of the distance between the low elevation and max lake elevation

    #  Gather centroid lat/lons
    out_lake_raster = os.path.join(projdir, "out_lake_raster.shp")
    out_lake_raster_dis = os.path.join(projdir, "out_lake_raster_dissolve.shp")
    arcpy.RasterToPolygon_conversion(outRastername, out_lake_raster, "NO_SIMPLIFY", "VALUE")
    arcpy.Dissolve_management(out_lake_raster, out_lake_raster_dis, "GRIDCODE", "", "MULTI_PART")               # Dissolve to eliminate multipart features
    arcpy.Delete_management(out_lake_raster)                                                                    # Added 9/4/2015

    # Create a point geometry object from gathered lake centroid points
    printMessages(arcpy, ['    Starting to gather lake centroid information.'])
    sr1 = arcpy.SpatialReference()                                              # Project lake points to whatever coordinate system is specified by wkt_text in globals
    sr1.loadFromString(wkt_text)                                                # Load the Sphere datum CRS using WKT
    point = arcpy.Point()
    cen_lats = {}
    cen_lons = {}
    shapes.update({row[0]: row[1] for row in arcpy.da.SearchCursor(out_lake_raster_dis, ['GRIDCODE', 'SHAPE@XY'])})
    arcpy.Delete_management(out_lake_raster_dis)                                                                # Added 9/4/2015
    for shape in shapes:
        point.X = shapes[shape][0]
        point.Y = shapes[shape][1]
        pointGeometry = arcpy.PointGeometry(point, sr2)
        projpoint = pointGeometry.projectAs(sr1)                                # Optionally add transformation method:
        cen_lats[shape] = projpoint.firstPoint.Y
        cen_lons[shape] = projpoint.firstPoint.X
    printMessages(arcpy, ['    Done gathering lake centroid information.'])

    # Create Lake parameter file
    printMessages(arcpy, ['    Starting to create lake parameter table.'])
    printMessages(arcpy, ['        Lakes Table: {0} lakes'.format(len(list(areas.keys())))])

    # Call function to build lake parameter netCDF file
    if 'nc' in out_LKtype:
        LakeNC = os.path.join(projdir, LK_nc)
        build_LAKEPARM(arcpy, LakeNC, min_elevs, areas, max_elevs, OrificEs, cen_lats, cen_lons, WeirE_vals)
    if 'ascii' in out_LKtype:
        LakeTBL = os.path.join(projdir, LK_tbl)
        build_LAKEPARM_ascii(LakeTBL, min_elevs, areas, max_elevs, OrificEs, cen_lats, cen_lons, WeirE_vals)
        printMessages(arcpy, ['        Done writing LAKEPARM.TBL table to disk.'])

    # Process: Output Channelgrid
    channelgrid_arr = arcpy.RasterToNumPyArray(NewChannelgrid)
    outRaster_arr = arcpy.RasterToNumPyArray(outRaster)

    # Clean up by deletion
    printMessages(arcpy, ['    Lake parameter table created without error in {0:3.2f} seconds.'.format(time.time()-tic1)])

    # Clean up
    del sr1, point, outRaster1
    arcpy.Delete_management("Lakeslyr")
    arcpy.Delete_management(outshp)
    #arcpy.Delete_management(outfeatures)
    arcpy.Delete_management('in_memory')
    arcpy.Delete_management(outRastername)
    del projdir, fill2, cellsize, sr2, in_lakes, lakeIDList, outRaster, channelgrid
    return arcpy, channelgrid_arr, outRaster_arr

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
    printMessages(arcpy, ['Step 4(a) initiated...'])
    mosprj2 = os.path.join(projdir, 'mosaicprj2')

    # Set environments
    arcpy.env.outputCoordinateSystem = sr2
    descData = arcpy.Describe(in_raster)
    arcpy.env.snapRaster = in_raster
    cellsize = descData.children[0].meanCellHeight
    printMessages(arcpy, ['    Cellsize: {0}.'.format(cellsize)])
    arcpy.env.cellSize = in_raster
    printMessages(arcpy, ['    Environments set.'])

    # Adjust height to raise all land areas by depending on input z units
    if inunits == 'm':
        above = 0.5
    elif inunits == 'cm':                                                       # trap for fixing cm to m conversion
        above = 50

    # Create landmask on the fine grid
    printMessages(arcpy, ['    Creating resampled LANDMASK layer from GEOGRID file.'])
    LANDMASK2 = os.path.join(projdir, "LANDMASK")
    arcpy.Resample_management(LANDMASK, LANDMASK2, cellsize, "NEAREST")
    printMessages(arcpy, ['    Finished creating resampled LANDMASK layer from GEOGRID file.'])

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
    printMessages(arcpy, ['    Topography corrected to match coarse grid LANDMASK.'])
    return mosprj2

def build_GWBUCKPARM(arcpy, out_dir, cat_areas, cat_comids, tbl_type='.nc'):
    '''
    5/17/2017: This function will build the groundwater bucket parameter table.
               Currently, two output formats are available: netCDF or .TBL
    '''
    tic1 = time.time()
    Community = True                                                            # Switch to provide Community WRF-Hydro GWBUCKPARM outputs

    if tbl_type in ['.nc', '.nc and .TBL']:

        # Produce output in NetCDF format (binary and much faster to produce)
        out_file = os.path.join(out_dir, GW_nc)                                     # Groundwater bucket parameter table path and filename
        rootgrp = netCDF4.Dataset(out_file, 'w', format=outNCType)

        # Create dimensions and set other attribute information
        #dim1 = 'BasinDim'
        dim1 = 'feature_id'
        dim = rootgrp.createDimension(dim1, len(cat_comids))

        # Create fixed-length variables
        Basins = rootgrp.createVariable('Basin', 'i4', (dim1))                  # Variable (32-bit signed integer)
        coeffs = rootgrp.createVariable('Coeff', 'f4', (dim1))                  # Variable (32-bit floating point)
        Expons = rootgrp.createVariable('Expon', 'f4', (dim1))                  # Variable (32-bit floating point)
        Zmaxs = rootgrp.createVariable('Zmax', 'f4', (dim1))                    # Variable (32-bit floating point)
        Zinits = rootgrp.createVariable('Zinit', 'f4', (dim1))                  # Variable (32-bit floating point)
        Area_sqkms = rootgrp.createVariable('Area_sqkm', 'f4', (dim1))          # Variable (32-bit floating point)
        ComIDs = rootgrp.createVariable('ComID', 'i4', (dim1))                  # Variable (32-bit signed integer)

        # Set variable descriptions
        Basins.long_name = 'Basin monotonic ID (1...n)'
        coeffs.long_name = 'Coefficient'
        Expons.long_name = 'Exponent'
        Zmaxs.long_name = 'Bucket height'
        Zinits.long_name = 'Initial height of water in bucket'
        Area_sqkms.long_name = 'Basin area in square kilometers'
        if Community:
            ComIDs.long_name = 'Catchment Gridcode'
        else:
            ComIDs.long_name = 'NHDCatchment FEATUREID (NHDFlowline ComID)'     # For NWM
        Zmaxs.units = 'mm'
        Zinits.units = 'mm'
        Area_sqkms.units = 'km2'

        # Fill in global attributes
        rootgrp.featureType = 'point'                                           # For compliance
        rootgrp.history = 'Created %s' %time.ctime()

        # Fill in variables
        Basins[:] = cat_comids                                                  #Basins[:] = numpy.arange(1,cat_comids.shape[0]+1)
        coeffs[:] = coeff
        Expons[:] = expon
        Zmaxs[:] = zmax
        Zinits[:] = zinit
        Area_sqkms[:] = numpy.array(cat_areas)
        ComIDs[:] = numpy.array(cat_comids)

        # Close file
        rootgrp.close()
        printMessages(arcpy, ['    Created output bucket parameter table (.nc): {0}.'.format(out_file)])
        del rootgrp, dim1, dim, Basins, coeffs, Expons, Zmaxs, Zinits, Area_sqkms, ComIDs, out_file

    # Build output files
    if tbl_type in ['.TBL', '.nc and .TBL']:
        out_file = os.path.join(out_dir, GW_TBL)                                # Groundwater bucket parameter table path and filename
        out_file = out_file.replace('.nc', '.TBL')
        fp = open(out_file, 'w')                                               # Build Bucket parameter table
        if Community:
            fp.write('Basin,Coeff,Expon,Zmax,Zinit\n')                          # Header for Community WRF-Hydro
        else:
            fp.write('Basin,Coeff.,Expon.,Zmax,Zinit,Area_sqkm,ComID\n')        # Header for NWM WRF-Hydro
        for ID, AREA in zip(cat_comids, cat_areas):
            if Community:
                newline = '{0:>8},{1:>6.4f},{2:>6.3f},{3:>6.2f},{4:>7.4f}\n'.format(ID, coeff, expon, zmax, zinit)
            else:
                newline = '{0:>8},{1:>6.4f},{2:>6.3f},{3:>6.2f},{4:>7.4f},{5:>11.3f},{0:>8}\n'.format(ID, coeff, expon, zmax, zinit, AREA)
            fp.write(newline)
        fp.close()
        printMessages(arcpy, ['  Created output bucket parameter table (.TBL): {0}.'.format(out_file)])
        del newline, fp, ID, AREA, out_file

    # Print statements and return
    printMessages(arcpy, ['  Finished building groundwater bucket parameter table in {0:3.2f} seconds.'.format(time.time()-tic1)])
    del arcpy, tic1, Community, tbl_type, cat_areas, cat_comids

def build_GWBASINS_nc(arcpy, GW_BUCKS2, out_dir, grid_obj):
    '''
    5/17/2017: This function will build the 2-dimensional groundwater bucket
    grid in netCDF format.

    NOTES: Redundant to provide GW_BUCKS2 array and GWBasns_arr, but the array is
    used for building the netCDF file and the raster is used for building the ASCII raseter.
    '''
    tic1 = time.time()
    out_file = os.path.join(out_dir, GWGRID_nc)
    varList2D = [['BASIN', 'i4', 'Basin ID corresponding to GWBUCKPARM table values']]

    # Create spatial metadata file for GEOGRID/LDASOUT grids
    rootgrp = netCDF4.Dataset(out_file, 'w', format=outNCType)
    rootgrp = create_CF_NetCDF(arcpy, grid_obj, rootgrp, addLatLon=False, notes='', addVars=varList2D)

    # Array size check
    GWBasns_arr = arcpy.RasterToNumPyArray(GW_BUCKS2, '#', '#', '#', nodata_to_value=NoDataVal) # This is done because of an occasional dimension mismatch
    printMessages(arcpy, ['    NC dimensions: {0}, {1}'.format(len(rootgrp.dimensions['y']), len(rootgrp.dimensions['x']))])
    printMessages(arcpy, ['    GWBUCKS array dimensions: {0}, {1}'.format(GWBasns_arr.shape[0], GWBasns_arr.shape[1])])

    # Add groundwater buckets to the file (flip UP-DOWN?)
    ncvar = rootgrp.variables[varList2D[0][0]]
    ncvar[:] = GWBasns_arr[:]                                                   # Populate variable with groundwater bucket array
    rootgrp.close()
    printMessages(arcpy, ['    Process: {0} completed without error'.format(out_file)])
    printMessages(arcpy, ['    Finished building groundwater grid file in {0:3.2f} seconds'.format(time.time()-tic1)])
    del arcpy, GW_BUCKS2, out_dir, tic1, out_file, varList2D, ncvar, rootgrp, GWBasns_arr

def build_GWBASINS_ascii(arcpy, GW_BUCKS2, out_dir, cellsize, GW_ASCII=GW_ASCII):
    '''
    5/17/2017: This function will build the 2-dimensional groundwater bucket
    grid in netCDF format.

    NOTES: Redundant to provide GW_BUCKS2 array and GWBasns_arr, but the array is
    used for building the netCDF file and the raster is used for building the ASCII raseter.
    '''
    tic1 = time.time()

    # Output to ASCII file
    arcpy.env.cellSize = cellsize                                               # Set cellsize to coarse grid
    out_file = os.path.join(out_dir, GW_ASCII)
    arcpy.RasterToASCII_conversion(GW_BUCKS2, out_file)
    printMessages(arcpy, ['    Process: {0} completed without error'.format(GW_ASCII)])

    # Remove header (first 6 rows) from ASCII Raster. NOTE: Opening file in r+ mode and using .writelines() produces strange behavior and an invalid ASCII raster for WRF-Hydro.
    f = open(out_file, 'r')                                                     # Open the input file in read mode
    lines = f.readlines()                                                       # Read entire file
    f.close()                                                                   # Close the file for reading
    os.remove(out_file)                                                         # Delete the ascii file completely
    f = open(out_file, 'w')
    for line in lines[6:]:
        f.write(line)                                                           # Write the non-header lines back to the file
    f.close()                                                                   # Close the file
    del f, line, lines                                                          # Clean up

    printMessages(arcpy, ['    Finished building groundwater grid file in {0:3.2f} seconds'.format(time.time()-tic1)])
    del arcpy, GW_BUCKS2, out_dir, tic1, out_file, GW_ASCII

def build_GW_Basin_Raster(arcpy, in_nc, projdir, in_method, grid_obj, in_Polys=None):
    '''
    10/10/2017:
    This function was added to build the groundwater basins raster using a variety
    of methods. The result is a raster on the fine-grid which can be used to create
    groundwater bucket parameter tables in 1D and 2D for input to WRF-Hydro.
    '''
    # Get ArcGIS version information and checkout Spatial Analyst extension
    ArcVersion = arcpy.GetInstallInfo()['Version']                              # Get the ArcGIS version that is being used
    if ArcVersion.count('.') == 2:
        ArcVersionF = float(ArcVersion.rpartition('.')[0])                      # ArcGIS major release version as a float
    else:
        ArcVersionF = float(ArcVersion)                                         # ArcGIS major release version as a float

    tic1 = time.time()
    printMessages(arcpy, ['Beginning to build 2D groundwater basin inputs'])
    printMessages(arcpy, ['  Building groundwater inputs using {0}'.format(in_method)])

    # Set environments
    arcpy.env.overwriteOutput = True
    arcpy.env.workspace = projdir
    arcpy.env.scratchWorkspace = projdir
    sr = grid_obj.proj

    # Open input FullDom file
    rootgrp1 = netCDF4.Dataset(in_nc, 'r')                                      # Read-only on FullDom file

    try:
        # Determine which method will be used to generate groundwater bucket grid
        if in_method == 'FullDom basn_msk variable':

            # Create a raster layer from the netCDF
            #GWBasns = arcpy.NumPyArrayToRaster(numpy.array(rootgrp1.variables['basn_msk'][:]), point, grid_obj.DX, grid_obj.DY)
            GWBasns = grid_obj.numpy_to_Raster(arcpy, numpy.array(rootgrp1.variables['basn_msk'][:]))
            arcpy.CalculateStatistics_management(GWBasns)
            #arcpy.DefineProjection_management(GWBasns, sr)

            # Gather projection and set output coordinate system
            descData = arcpy.Describe(GWBasns)
            #sr = descData.spatialReference
            arcpy.env.outputCoordinateSystem = sr

        elif in_method == 'FullDom LINKID local basins':

            if 'LINKID' not in rootgrp1.variables:
                printMessages(arcpy, ['    Generating LINKID grid from CHANNELGRID and FLOWDIRECTION'])

                # In-memory files
                outStreams_ = r'in_memory\Streams'
                outRaster = r'in_memory\LINKID'

                # Read Fulldom file variable to raster layer
                for ncvarname in ['CHANNELGRID', 'FLOWDIRECTION']:
                    printMessages(arcpy, ['    Creating layer from netCDF variable {0}'.format(ncvarname)])
                   # nc_raster = arcpy.NumPyArrayToRaster(numpy.array(rootgrp1.variables[ncvarname][:]), point, grid_obj.DX, grid_obj.DY)
                    nc_raster = grid_obj.numpy_to_Raster(arcpy, numpy.array(rootgrp1.variables[ncvarname][:]))
                    #arcpy.DefineProjection_management(nc_raster, sr)
                    nc_raster.save(os.path.join('in_memory', ncvarname))
                    arcpy.CalculateStatistics_management(nc_raster)
                    arcpy.env.snapRaster = nc_raster
                    if ncvarname == 'CHANNELGRID':
                        chgrid = SetNull(nc_raster, '1', 'VALUE < 0')
                    elif ncvarname == 'FLOWDIRECTION':
                        fdir = Int(nc_raster)

                # Gather projection and set output coordinate system
                descData = arcpy.Describe(chgrid)
                sr = descData.spatialReference
                arcpy.env.outputCoordinateSystem = sr
                arcpy.env.snapRaster = fdir
                arcpy.env.extent = fdir

                # Convert to stream features, then back to raster
                StreamToFeature(chgrid, fdir, outStreams_, "NO_SIMPLIFY")            # Stream to feature
                if ArcVersionF >= 10.5:
                    arcidField = 'arcid'
                else:
                    arcidField = 'ARCID'
                arcpy.FeatureToRaster_conversion(outStreams_, arcidField, outRaster)    # Must do this to get "ARCID" field into the raster

                # Alter raster to handle some spurious nodata values (?)
                maxValue = arcpy.SearchCursor(outStreams_, "", "", "", arcidField + " D").next().getValue(arcidField)  # Gather highest "ARCID" value from field of segment IDs
                maxRasterValue = arcpy.GetRasterProperties_management(outRaster, "MAXIMUM")                     # Gather maximum "ARCID" value from raster
                if int(maxRasterValue[0]) > maxValue:
                    print('    Setting linkid values of %s to Null.' %maxRasterValue[0])
                    whereClause = "VALUE = %s" %int(maxRasterValue[0])
                    outRaster = SetNull(outRaster, outRaster, whereClause)          # This should eliminate the discrepency between the numbers of features in outStreams and outRaster
                outRaster = Con(IsNull(outRaster) == 1, NoDataVal, outRaster)         # Set Null values to -9999 in LINKID raster

                printMessages(arcpy, ['    Creating StreamLink grid'])
                strm = SetNull(outRaster, outRaster, "VALUE = %s" %NoDataVal)
                strm.save(os.path.join('in_memory', 'LINKID'))
                arcpy.Delete_management(outStreams_)
                arcpy.Delete_management(chgrid)
                #arcpy.Delete_management(outRaster)
                del chgrid, maxValue, maxRasterValue, outStreams_    # outRaster

            else:
                printMessages(arcpy, ['    LINKID exists in FullDom file.'])
                for ncvarname in ['LINKID', 'FLOWDIRECTION']:
                    #nc_raster = arcpy.NumPyArrayToRaster(numpy.array(rootgrp1.variables[ncvarname][:]), point, grid_obj.DX, grid_obj.DY)
                    nc_raster = grid_obj.numpy_to_Raster(arcpy, numpy.array(rootgrp1.variables[ncvarname][:]))
                    #arcpy.DefineProjection_management(nc_raster, sr)
                    nc_raster.save(os.path.join('in_memory', ncvarname))
                    arcpy.CalculateStatistics_management(nc_raster)
                    arcpy.env.snapRaster = nc_raster
                    if ncvarname == 'LINKID':
                        strm = SetNull(nc_raster, nc_raster, 'VALUE = %s' %NoDataVal)
                        strm.save(os.path.join('in_memory', 'strm'))
                    elif ncvarname == 'FLOWDIRECTION':
                        fdir = Int(nc_raster)

                # Gather projection and set output coordinate system
                descData = arcpy.Describe(strm)
                sr = descData.spatialReference
                arcpy.env.outputCoordinateSystem = sr

            # Create contributing watershed grid
            GWBasns = Watershed(fdir, strm, 'VALUE')
            arcpy.Delete_management(strm)
            arcpy.Delete_management(fdir)
            arcpy.Delete_management(nc_raster)
            del strm, fdir, nc_raster, ncvarname
            printMessages(arcpy, ['    Stream to features step complete.'])

        elif in_method == 'Polygon Shapefile or Feature Class':
            printMessages(arcpy, ['    Polygon Shapefile input: {0}'.format(in_Polys)])
            #nc_raster = arcpy.NumPyArrayToRaster(numpy.array(rootgrp1.variables['TOPOGRAPHY'][:]), point, grid_obj.DX, grid_obj.DY)
            nc_raster = grid_obj.numpy_to_Raster(arcpy, numpy.array(rootgrp1.variables['TOPOGRAPHY'][:]))
            #arcpy.DefineProjection_management(nc_raster, sr)
            arcpy.CalculateStatistics_management(nc_raster)

           # Gather projection and set output coordinate system
            descData = arcpy.Describe(nc_raster)
            sr = descData.spatialReference
            arcpy.env.outputCoordinateSystem = sr
            arcpy.env.snapRaster = nc_raster
            arcpy.env.extent = nc_raster
            arcpy.env.cellSize = grid_obj.DX                                               # Set cellsize to fine grid

            # Resolve basins on the fine grid
            descData = arcpy.Describe(in_Polys)
            GWBasnsFile = os.path.join('in_memory', 'poly_basins')
            arcpy.PolygonToRaster_conversion(in_Polys, descData.OIDFieldName, GWBasnsFile, "MAXIMUM_AREA", "", grid_obj.DX)
            GWBasns = arcpy.Raster(GWBasnsFile)                                         # Create raster object from raster layer
            arcpy.Delete_management(nc_raster)
            del nc_raster

        GWBasns_arr = arcpy.RasterToNumPyArray(GWBasns)                             # Create array from raster
        rootgrp1.close()
        printMessages(arcpy, ['Finished building groundwater basin grids in %3.2f seconds' %(time.time()-tic1)])
        del sr, descData

    except:
        rootgrp1.close()
        printMessages(arcpy, ['Finished building groundwater basin grids in %3.2f seconds' %(time.time()-tic1)])
        raise SystemExit

    del arcpy, rootgrp1, tic1
    return GWBasns, GWBasns_arr

def build_GW_buckets(arcpy, out_dir, GWBasns, GWBasns_arr, coarse_grid, tbl_type='.nc', Grid=True):
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

    tic1 = time.time()
    printMessages(arcpy, ['Beginning to build groundwater inputs'])

    # Read basin information from the array
    UniqueVals = numpy.unique(GWBasns_arr[:])                                   # Array to store the basin ID values in the fine-grid groundwater basins
    UniqueVals = UniqueVals[UniqueVals >= 0]                                    # Remove NoData, removes potential noData values (-2147483647, -9999)
    printMessages(arcpy, ['    Found {0} basins in the watershed grid'.format(UniqueVals.shape)])
    del UniqueVals

    GW_BUCKS = os.path.join('in_memory', 'GWBUCKS_orig')                        # Groundwater bucket raster on the coarse grid before re-numbering 1...n
    coarse_grid.project_to_model_grid(arcpy, GWBasns, GW_BUCKS, resampling="NEAREST")   # Resample from fine grid to coarse grid resolution
    arcpy.Delete_management(GWBasns)                                            # Delete original fine-grid groundwater basin raster
    del GWBasns, GWBasns_arr

    # Re-assign basin IDs to 1...n because sometimes the basins get lost when converting to coarse grid
    GWBasns_arr2 = arcpy.RasterToNumPyArray(GW_BUCKS, nodata_to_value=NoDataVal)# Array to store the basin ID values in the resampled coarse-grid groundwater basins
    UniqueVals2 = numpy.unique(GWBasns_arr2[:])                                 # Get the unique values, including nodata
    arcpy.Delete_management(GW_BUCKS)                                           # Delete the resampled-to-coarse-grid groundwater basin raster
    printMessages(arcpy, ['    Found {0} basins (potentially including nodata values) in the file after resampling to the coarse grid.'.format(UniqueVals2.shape[0])])

    '''Because we resampled to the coarse grid, we lost some basins. Thus, we need to
    re-assign basin ID values to conform to the required 1...n groundwater basin
    ID assignment scheme.'''
    # Fast replace loop from https://stackoverflow.com/questions/3403973/fast-replacement-of-values-in-a-numpy-array
    # This method ensures that any nodata values are issued a 0 index in sort_idx
    sort_idx = numpy.argsort(UniqueVals2)                                       # Index of each unique value, 0-based
    idx = numpy.searchsorted(UniqueVals2, GWBasns_arr2, sorter=sort_idx)        # 2D array of index values against GWBasns_arr2
    del GWBasns_arr2, sort_idx                                                  # Free up memory

    # This method requires the nodata value to be issued a 0 index in sort_idx and idx
    to_values = numpy.arange(UniqueVals2.size)                                  # 0..n values to be substituted, 0 in place of NoDataVal
    GWBasns_arr3 = to_values[idx]                                               # Same as to_values[sort_idx][idx]
    if numpy.where(UniqueVals2 == NoDataVal)[0].shape[0] > 0:
        new_ndv = int(to_values[numpy.where(UniqueVals2 == NoDataVal)[0]][0])         # Obtain the newly-assigned nodatavalue
    else:
        new_ndv = NoDataVal
        GWBasns_arr3 += 1                                                       # Add one so that the basin IDs will be 1...n rather than 0...n when there are no nodata values in the grid
    del UniqueVals2

    # Build rasters and arrays to create the NC or ASCII outputs
    GW_BUCKS = coarse_grid.numpy_to_Raster(arcpy, GWBasns_arr3, value_to_nodata=new_ndv)
    GW_BUCKS2 = os.path.join(out_dir, "GWBUCKS")                                # Raster file location
    GW_BUCKS.save(GW_BUCKS2)                                                    # Save raster
    del GWBasns_arr3, idx, to_values, new_ndv

    # If requested, create 2D gridded bucket parameter table
    if Grid:
        if 'nc' in out_2Dtype:
            build_GWBASINS_nc(arcpy, GW_BUCKS2, out_dir, coarse_grid)
        if 'ascii' in out_2Dtype:
            build_GWBASINS_ascii(arcpy, GW_BUCKS2, out_dir, coarse_grid.DX)

    # Alternate method to obtain IDs - read directly from raster attribute table
    printMessages(arcpy, ['    Calculating size and ID parameters for basin polygons.'])
    Raster_arr = arcpy.da.TableToNumPyArray(GW_BUCKS2, '*')
    cat_comids = Raster_arr['VALUE'].tolist()
    cat_areas = [float((item*(coarse_grid.DX**2))/1000000) for item in Raster_arr['COUNT'].tolist()]  # Assumes cellsize is in units of meters
    arcpy.Delete_management(GW_BUCKS)
    arcpy.Delete_management(GW_BUCKS2)
    del GW_BUCKS, GW_BUCKS2, Raster_arr

    # Build the groundwater bucket parameter table in netCDF format
    build_GWBUCKPARM(arcpy, out_dir, cat_areas, cat_comids, tbl_type)

    # Clean up and return
    del cat_comids, cat_areas
    arcpy.Delete_management('in_memory')
    printMessages(arcpy, ['Finished building groundwater parameter files in {0:3.2f} seconds'.format(time.time()-tic1)])
    del arcpy, tic1, tbl_type, out_dir, Grid, coarse_grid   # Attempt to get the script to release the Fulldom_hires.nc file
    return

def sa_functions(arcpy, rootgrp, bsn_msk, mosprj, ovroughrtfac_val, retdeprtfac_val, projdir, in_csv, threshold, routing, WKT='', in_lakes=None, lakeIDfield=None):
    """The last major function in the processing chain is to perform the spatial
    analyst functions to hydrologically process the input raster datasets."""

    arcpy.env.overwriteOutput = True

    # Fourth part of the process
    printMessages(arcpy, ['Step 4 initiated...'])

    # Set Basin mask attribute to boolean from ArcGIS text
    if bsn_msk:
        printMessages(arcpy, ['    Channelgrid will be masked to basins.'])
    else:
        printMessages(arcpy, ['    Channelgrid will not be masked to basins.'])

    if routing:
        printMessages(arcpy, ['    Reach-based routing files will be created.'])
    else:
        printMessages(arcpy, ['    Reach-based routing files will not be created.'])

    # Set environments
    arcpy.MakeRasterLayer_management(mosprj, 'mosaicprj')
    descData = arcpy.Describe('mosaicprj')
    sr2 = descData.spatialReference
    cellsize = descData.meanCellHeight
    arcpy.env.outputCoordinateSystem = sr2
    arcpy.env.workspace = projdir
    arcpy.env.scratchWorkspace = projdir

    # Process: Fill DEM
    fill = Fill(mosprj, z_limit)                                                # Sink-filling using z-limit (global)
    fill1 = Float(fill)
    #del fill

    # Process: Flow Direction
    fdir = FlowDirection(fill1)
    fdir_var = rootgrp.variables['FLOWDIRECTION']
    fdir_arr = arcpy.RasterToNumPyArray(fdir)
    fdir_var[:] = fdir_arr
    printMessages(arcpy, ['    Process: FLOWDIRECTION written to output netCDF.'])
    del fdir_arr

    # Process: Flow Accumulation (intermediate
    flac = FlowAccumulation(fdir, '#', 'FLOAT')
    flac_var = rootgrp.variables['FLOWACC']
    flac_arr = arcpy.RasterToNumPyArray(flac)
    flac_var[:] = flac_arr
    printMessages(arcpy, ['    Process: FLOWACC written to output netCDF.'])
    del flac_arr

    # Set NoData Value for topography
    fill2 = Con(IsNull(fill1) == 1, NoDataVal, fill1)
    del fill1
    nodataval = flac.noDataValue
    if nodataval is None:
        nodataval = NoDataVal                                                   # Fix to help this portion run in ArcMap
    arcpy.SetRasterProperties_management(fill2, "ELEVATION", "#", "#", [[1, nodataval]])     # Must set NoData away from -Infinityf
    fill2_var = rootgrp.variables['TOPOGRAPHY']
    fill2_arr = arcpy.RasterToNumPyArray(fill2)
    fill2_var[:] = fill2_arr
    printMessages(arcpy, ['    Process: TOPOGRAPHY written to output netCDF.'])
    del fill2_arr

    # Create stream channel raster according to threshold
    strm = SetNull(flac, '1', 'VALUE < %s' % threshold)
    channelgrid = Con(IsNull(strm) == 0, 0, NoDataVal)

    # Create initial constant raster of -9999
    constraster = CreateConstantRaster(NoDataVal, "INTEGER")
    constraster_arr = arcpy.RasterToNumPyArray(constraster)
    arcpy.Delete_management(constraster)

    # Create initial constant raster of value retdeprtfac_val
    inraster2 = CreateConstantRaster(retdeprtfac_val, "FLOAT")
    inraster2_var = rootgrp.variables['RETDEPRTFAC']
    inraster2_arr = arcpy.RasterToNumPyArray(inraster2)
    inraster2_var[:] = inraster2_arr
    printMessages(arcpy, ['    Process: RETDEPRTFAC written to output netCDF.'])
    arcpy.Delete_management(inraster2)
    del retdeprtfac_val, inraster2, inraster2_arr

    # Create initial constant raster of ovroughrtfac_val
    inraster3 = CreateConstantRaster(ovroughrtfac_val, "FLOAT")
    inraster3_var = rootgrp.variables['OVROUGHRTFAC']
    inraster3_arr = arcpy.RasterToNumPyArray(inraster3)
    inraster3_var[:] = inraster3_arr
    printMessages(arcpy, ['    Process: OVROUGHRTFAC written to output netCDF.'])
    arcpy.Delete_management(inraster3)
    del ovroughrtfac_val, inraster3, inraster3_arr

    # Create initial constant raster of LKSATFAC - added 3/30/2017
    inraster4 = CreateConstantRaster(lksatfac_val, "FLOAT")
    inraster4_var = rootgrp.variables['LKSATFAC']
    inraster4_arr = arcpy.RasterToNumPyArray(inraster4)
    inraster4_var[:] = inraster4_arr
    printMessages(arcpy, ['    Process: LKSATFAC written to output netCDF.'])
    arcpy.Delete_management(inraster4)
    del inraster4, inraster4_arr

    # Process: Stream Order
    order = StreamOrder(strm, fdir)         # Default = "STRAHLER"
    order2 = Con(IsNull(order) == 1, NoDataVal, order)
    order2_var = rootgrp.variables['STREAMORDER']
    order2_arr = arcpy.RasterToNumPyArray(order2)
    order2_var[:] = order2_arr
    printMessages(arcpy, ['    Process: STREAMORDER written to output netCDF.'])
    del order2_arr

    # Find out if forecast points are chosen, then set mask for them
    if in_csv is not None:
        # Make feature layer from CSV
        printMessages(arcpy, ['    Forecast points provided and basins being delineated.'])
        frxst_layer = 'frxst_layer'                                             # XY Event Layer name
        frsxt_FC = os.path.join('in_memory', 'frxst_FC')                        # In-memory feature class
        sr1 = arcpy.SpatialReference(4326)                                      # GCS_WGS_1984
        arcpy.MakeXYEventLayer_management(in_csv, 'LON', 'LAT', frxst_layer, sr1)
        arcpy.CopyFeatures_management(frxst_layer, frsxt_FC)                    # To support ArcGIS 10.6, an XY event layer cannot be used in SnapRaster
        tolerance = int(float(arcpy.GetRasterProperties_management('mosaicprj', 'CELLSIZEX').getOutput(0)) * walker)
        tolerance1 = int(float(arcpy.GetRasterProperties_management('mosaicprj', 'CELLSIZEX').getOutput(0)))
        frxst_raster = SnapPourPoint(frsxt_FC, flac, tolerance1, 'FID')
        frxst_raster2 = Con(IsNull(frxst_raster) == 0, frxst_raster, NoDataVal)
        frxst_raster2_var = rootgrp.variables['frxst_pts']
        frxst_raster2_arr = arcpy.RasterToNumPyArray(frxst_raster2)
        frxst_raster2_var[:] = frxst_raster2_arr
        printMessages(arcpy, ['    Process: frxst_pts written to output netCDF.'])
        del frxst_raster2_arr

        SnapPour = SnapPourPoint(frsxt_FC, flac, tolerance)
        arcpy.Delete_management(frxst_layer)
        arcpy.Delete_management(frsxt_FC)
        arcpy.Delete_management(frxst_raster2)
        del frxst_raster2

        # Delineate above points
        outWatershed = Watershed(fdir, SnapPour, 'VALUE')
        outWatershed2 = Con(IsNull(outWatershed) == 0, outWatershed, NoDataVal)
        outWatershed2_arr = arcpy.RasterToNumPyArray(outWatershed2)

        # Default groundwater method changed 10/5/2017 to LINKID local basin method
        gw_basns_var = rootgrp.variables['basn_msk']
        gw_basns_var[:] = outWatershed2_arr
        printMessages(arcpy, ['    Process: basn_msk written to output netCDF.'])
        del outWatershed2_arr

        # Set mask for future raster output
        if bsn_msk:
            channelgrid2 = Con(outWatershed2 >= 0, IsNull(strm), Con(IsNull(strm), 1, -1))  # Converts channelgrid values inside basins to 0, outside to -1
            channelgrid = Con(channelgrid2 == 1, NoDataVal, channelgrid2)
            del channelgrid2
        del outWatershed2

    else:
        # Handle the change if input forecast points are not provided
        constraster_var = rootgrp.variables['frxst_pts']
        constraster_var[:] = constraster_arr
        printMessages(arcpy, ['    Process: frxst_pts was empty. Constant value raster created.'])

        constraster_var2 = rootgrp.variables['basn_msk']
        constraster_var2[:] = constraster_arr
        printMessages(arcpy, ['    Process: gw_basns was empty. Constant value raster created.'])

    # Moved 10/9/2017 by KMS to allow masking routing files (LINKID, Route_Link, etc.) to forecast points if requested
    if routing:
        rasterExp = "Value = -9999"                                             # Default: all channels will be represented in reach-based routing file
        if in_csv is None:                                                      # Added 10/10/2017 by KMS to include forecast points in reach-based routing file
            frxst_raster = None                                                 # Default is no forecast points for reach-based routing file
        elif maskRL:
            rasterExp = "Value < 0"                                             # Only channels within the masked basin will be in reach-based routing file
        strm2 = SetNull(channelgrid, channelgrid, rasterExp)                    # Alter channelgrid such that -9999 and -1 to be NoData
        linkid = Routing_Table(arcpy, projdir, sr2, strm2, fdir, fill2, order2, gages=frxst_raster, Lakes=in_lakes)
        linkid_var = rootgrp.variables['LINKID']
        linkid_arr = arcpy.RasterToNumPyArray(linkid)
        linkid_var[:] = linkid_arr
        printMessages(arcpy, ['    Process: LINKID written to output netCDF.'])
        del linkid_arr, strm2

    # Clean up
    arcpy.Delete_management(order)
    arcpy.Delete_management(order2)
    del order, order2

    if in_csv is not None:
        arcpy.Delete_management(frxst_raster)                                   # Clean up
        del frxst_raster

    if in_lakes is None:
        # Process: Output LAKEGRID
        constraster_var3 = rootgrp.variables['LAKEGRID']
        constraster_var3[:] = constraster_arr
        printMessages(arcpy, ['    Process: LAKEGRID written to output netCDF.'])

        # Process: Output Channelgrid
        channelgrid_var = rootgrp.variables['CHANNELGRID']
        channelgrid_arr = arcpy.RasterToNumPyArray(channelgrid)
        channelgrid_var[:] = channelgrid_arr
        printMessages(arcpy, ['    Process: CHANNELGRID written to output netCDF.'])

    else:
        # Alter Channelgrid for reservoirs
        fill2.save(os.path.join(projdir, 'fill2.tif'))                          # Added 6/4/2019 to keep ArcGIS 10.6.1 from crashing later on in the lake script
        arcpy, channelgrid_arr, outRaster_arr = add_reservoirs(arcpy, channelgrid, in_lakes, flac, projdir, fill2, cellsize, sr2, lakeIDfield)

        # Process: Output LAKEGRID
        outRaster_var = rootgrp.variables['LAKEGRID']
        outRaster_var[:] = outRaster_arr
        printMessages(arcpy, ['    Process: LAKEGRID written to output netCDF.'])

        # Process: Output Channelgrid
        channelgrid_var = rootgrp.variables['CHANNELGRID']
        channelgrid_var[:] = channelgrid_arr
        printMessages(arcpy, ['    Process: CHANNELGRID written to output netCDF.'])

    # Clean up
    arcpy.Delete_management(fill)                                               # This cannot seem to be done until this portion of the function
    arcpy.Delete_management(flac)
    arcpy.Delete_management(fill2)
    arcpy.Delete_management(mosprj)
    arcpy.Delete_management('mosaicprj')
    arcpy.Delete_management('in_memory')
    del constraster, fdir, flac, strm, channelgrid, fill2, channelgrid_arr, constraster_arr, fill
    printMessages(arcpy, ['    Step 4 completed without error.'])
    return rootgrp

# --- End Functions --- #

if __name__ == '__main__':
    '''Protect this script against import from another script.'''
    pass