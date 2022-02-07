#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
conda activate geo_env

RUN GDAL_RASTERIZE FIRST. (Unless the required .tif file already exists.)

This code came from Phil Broderick strategy he described in 2022/2/03 email.

The shapefile (and generated raster files) are in:
/raid8/pkalmus/data/coral/data/location/14_001_WCMC008_CoralReefs2018_v4_1/01_Data

To create the .tif file BEFORE RUNNING THIS PYTHON CODE:
gdal_rasterize WCMC008_CoralReef2018_Py_v4_1.shp output1km.tif -te -180.005 -34.305 179.995 32.525 -tr 0.01 0.01 -burn 1
(125479, 2)

all touched: (DECIDED NOT TO DO THIS)
gdal_rasterize WCMC008_CoralReef2018_Py_v4_1.shp output1at.tif -tr 0.01 0.01 -burn 1 -at
(421207, 2)

Also, note that the gdal_rasterize step could probably also be done within this .py script
using osgeo.gdal.RasterizeLayer, which also has options like 'ALL_TOUCHED=TRUE'

Here is what ogrinfo says about the shapefile set:
(geo_env) [pkalmus@weather2 01_Data]$ ogrinfo -so WCMC008_CoralReef2018_Py_v4_1.shp WCMC008_CoralReef2018_Py_v4_1
INFO: Open of `WCMC008_CoralReef2018_Py_v4_1.shp'
      using driver `ESRI Shapefile' successful.

Layer name: WCMC008_CoralReef2018_Py_v4_1
Metadata:
  DBF_DATE_LAST_UPDATE=2021-03-30
Geometry: Polygon
Feature Count: 17504
Extent: (-179.999935, -34.298230) - (179.999936, 32.514818)
Layer SRS WKT:
GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0,
        AUTHORITY["EPSG","8901"]],
    UNIT["degree",0.0174532925199433,
        AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4326"]]

@author: pkalmus
"""

from osgeo import gdal
import numpy as np
import xarray as xr
# import coral_plotter as plotter
import pdb

basedir = '/home/pkalmus/projects/coral2/'
olddir = '/home/pkalmus/projects/coral/'
datadir = '/raid8/pkalmus/data/coral/data/location/14_001_WCMC008_CoralReefs2018_v4_1/01_Data/'
lat_extent = 35.5
lat_extent = [-34.3, 32.525]

raster = datadir+'/output1km.tif'

ds = gdal.Open(raster) # osgeo.gdal.Dataset

trans = ds.GetGeoTransform() # https://gdal.org/api/gdaldataset_cpp.html#_CPPv4N11GDALDataset15GetGeoTransformEPd
dat = ds.ReadAsArray() # 2D numpy array, (6682, 36001)

locs = np.where(dat == 1)
output = np.zeros((len(locs[0]),2))
for i in range(len(locs[0])): 
    y = locs[0][i]
    x = locs[1][i]
    output[i,0] = trans[0] + (x+0.5) * trans[1]
    output[i,1] = trans[3] + (y+0.5) * trans[5]

print(output.shape)

# dump output in csv file
csv_name = datadir+'locationsWCMCv41.csv'
np.savetxt(csv_name, output, delimiter=',', fmt='%1.3f')
print(csv_name)

# start making the new dataset with the reef_mask
mur_file = olddir + 'data/mur/201912-JPL-L4-SSTfnd-MUR_monthly-GLOB-fv04.2.nc'
mur_ds  = xr.open_dataset(mur_file)
# drop unneeded variables
mur_ds = mur_ds.drop('monthly_mean_sst')
mur_ds = mur_ds.drop('monthly_std_sst')
mur_ds = mur_ds.sel(lat=slice(lat_extent[0], lat_extent[1]))

wcmc_ds = xr.DataArray(dat, name='reef_mask', coords=[mur_ds.coords['lat'], mur_ds.coords['lon']], dims=['lat', 'lon'])  

# this was the old way, see coral/src/python/reef_locations.py
# # make reef_mask, initialized to zeros.
# null_data = np.empty((mur_ds.dims['lat'], mur_ds.dims['lon'])) # already initialized with zeros 
# reef_mask_da = xr.DataArray(null_data, coords=[mur_ds.coords['lat'], mur_ds.coords['lon']], dims=['lat', 'lon'])  
# # instead of loopoing through the 0.01 resolution data, which would take about 24 hours, a nearest neighbor interp
# interpolated = wcmc_ds.interp_like(reef_mask_da, method='nearest')

# dump output in netcdf file
netcdf_name = datadir+'locationsWCMCv41.nc'
#encoding cuts file size from 2 GB to 245 MB. encoding lat and lon to float32 prevents extra decimal garbage.
wcmc_ds.to_netcdf(netcdf_name, encoding={'reef_mask': {'dtype': 'int8', '_FillValue': -9999}, 'lat': {'dtype': 'float32', '_FillValue': -9999}, 'lon': {'dtype': 'float32', '_FillValue': -9999} } )  
print(netcdf_name)

#plotter.scatter(output[:,0], output[:,1], 'test.png', c=np.ones(output.shape[0]), marker_relsize=0.1, figsize=(15,5), cbar_fraction=0.008, cbar_orientation='vertical',fontsize=19)

pdb.set_trace()