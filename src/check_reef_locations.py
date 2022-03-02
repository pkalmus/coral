#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
conda activate geo_env

Read in 
    
@author: pkalmus
"""

import numpy as np
import xarray as xr
import coral_plotter as plotter
import matplotlib.cm
import pdb
import pandas as pd

basedir = '/home/pkalmus/projects/coral2/'
datadir = '/raid8/pkalmus/data/coral/data/location/14_001_WCMC008_CoralReefs2018_v4_1/01_Data/'
netcdf_name = datadir+'locationsWCMCv41_mur.nc'
figdir = basedir+'src/fig/'

netcdf_name = datadir+'locationsWCMCv41_gdal.nc'
csv_name = datadir+'locationsWCMCv41.csv'

# plotting whole netCDF file is very slow
ds = xr.open_dataset(netcdf_name)
stack = ds.reef_mask.stack(index=([...])) # this syntax stacks (ravels) all dimensions onto index.
stack = stack.where(stack==1, drop=True)
plotter.scatter(stack.lon.values, stack.lat.values, figdir+'check_reef_locations_netcdf.png', c=np.ones(stack.lat.values.shape[0]), marker_relsize=0.1, figsize=(15,5), cbar_fraction=0.008, cbar_orientation='vertical',fontsize=19)

pdb.set_trace()

ds = ds.where(ds.reef_mask==1)
dszoom = ds.sel(lat=slice(-20, 20))
longrid, latgrid = np.meshgrid(dszoom.lon.values, dszoom.lat.values)
cmap = matplotlib.cm.jet
cmap.set_bad(color='w')
plotter.image(longrid, latgrid, dszoom.reef_mask.values, figdir+'check_reef_locations_zoom.png', cmap=cmap, figsize=(15,5))

pdb.set_trace()

# plot the CSV file
df = pd.read_csv(csv_name)
plotter.scatter(df.values[:,0], df.values[:,1], figdir+'check_reef_locations_scatter.png', c=np.ones(df.values.shape[0]), marker_relsize=0.1, figsize=(15,5), cbar_fraction=0.008, cbar_orientation='vertical',fontsize=19)
