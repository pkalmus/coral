#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Must run as GEO

Take the "vector" format (one-dimensional index) output files and put them
on a lat-lon grid.

todo: could save values as ints for compression.
todo: could maybe get rioxarray to work!

Created on Wed Mar  3 14:08:06 2021

@author: pkalmus
"""

#import rioxarray
import xarray as xr
import numpy as np
import pdb 
import sys
import subprocess

# read in user params
import importlib
paramfile = sys.argv[1]   
params = importlib.import_module(paramfile)
projectdir = params.projectdir
runname = params.runname
basedir = params.basedir
pythondir = params.pythondir
listfilename = params.listfilename
scenarios = params.scenarios
do_plots = params.do_plots
return_years = params.return_years
DHW_threshold = params.DHW_threshold
bias_correct = params.bias_correct
do_weights = params.do_weights
do_obs = params.do_obs
finerun = params.finerun

finedir = basedir+'output/%s/fine/%s/' % (runname, finerun) 
rasterdir = finedir+'raster/'

sysStr = "mkdir -p %s" % rasterdir
subprocess.call(sysStr, shell=True)

# make a little test file
#lon: -1 to 1; lat: -1 to 1. 0.01 increments
lat = np.arange(-35, 35, 0.01)
lon = np.arange(0, 360, 0.01)
longrid, latgrid = np.meshgrid(lon,lat)



for return_year in return_years:
    for scenario in scenarios:
        if do_weights:
            experiment = '%s_%iyrs_DHW%1.1f' % (scenario, return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
        else:
            experiment = 'FLAT_%s_%iyrs_DHW%1.1f' % (scenario, return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod            

        print(experiment )
        finefile = finedir+'finemur_%s_%s.nc' % (runname, experiment)
        rasterfile = rasterdir+'raster_%s_%s.nc' % (runname, experiment)
        tiff_file = rasterdir+'raster_%s_%s.tif' % (runname, experiment)
        
        # rasterfile = rasterdir+'test_sort.nc'
        # tiff_file = rasterdir+'test_sort.tif'
        
        xds = xr.open_dataset(finefile)
        lats = xds.lat.values
        lons = xds.lon.values
        vals = xds.departure_year.values
        
        # convert from "index" to regular x, y grid on 0.01 lat lon grid.
        # make new ds with x, y dimensions (all null data)
        lat = np.arange(-35, 35, 0.01)
        lon = np.arange(0, 360, 0.01)
        longrid, latgrid = np.meshgrid(lon,lat)
        
        null_data = np.empty(longrid.shape)
        null_data[:] = np.nan
        
        #dep_da = xr.DataArray(null_data, dims=['y', 'x'], coords=dict(lon=(['y', 'x'], longrid), lat=(['y', 'x'], latgrid) ))  
        #dep_da = xr.DataArray(null_data,  dims=['lat', 'lon'], coords=dict(lon=(['lon'], lon), lat=(['lat'], lat) ))  
        dep_da = xr.DataArray(null_data,  dims=['y', 'x'], coords=dict(lon=(['x'], lon), lat=(['y'], lat) )) 

        #loop through old ds and populate new one
        for (ind, ignore) in enumerate(xds.index.values):
            mylat = lats[ind]
            mylon = lons[ind]
            
            xind = int(np.rint(mylon*100)) # had trouble with numerical noise leading to "striping"
            yind = int(np.rint((mylat+35)*100))
            dep_da.loc[dict(x=xind, y=yind)] = vals[ind]
                    
            if np.mod(ind, 10000)==0:
                print(str(ind) + '/' + str(len(xds.index.values)))
        
        # note: I checked the netCDF files, both 0 to 360 and -180 to 180 versions and both seemed fine.
        # switch lon coord to -180 to 180
        dep_da = dep_da.assign_coords(lon=(((dep_da.lon + 180) % 360) - 180)).sortby('lon') # this is checked and correct.
        
        dep_ds = xr.Dataset({'departure_year':dep_da})
        dep_ds.to_netcdf(rasterfile) 
        print(rasterfile)        
        
        # here is the GDAL command I used
        # gdal_translate test2_370_3yrs_DHW4.8.nc test2_370.tif -a_srs EPSG:4326 -a_ullr 0 35 360 -35
        sysStr = "gdal_translate %s %s -a_srs EPSG:4326 -a_ullr -180 35 180 -35" % (rasterfile, tiff_file)
        #sysStr = "gdal_translate %s %s -a_srs EPSG:4326 -a_ullr 0 35 360 -35" % (rasterfile, tiff_file)
        
        subprocess.call(sysStr, shell=True)
        print(tiff_file)
