 #!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Get coral locations for fine MUR grid as defined.

Based on coral_locations.py

Read in sample MUR ds.
Make a new ds with lat, lon, time, and mask.
For each element where mask is ocean or coast, check for a reef location.
Make a new variable  called "reef_mask" where 1 means it has a reef, 0 if not.

Make a plot from this file to test. (Read it in and use the xarray plotter).

Created on Wed Sep  4 12:29:02 2019

@author: pkalmus
"""

import matplotlib.pyplot as plt
import xarray as xr
import pdb
from sklearn.utils.extmath import cartesian

#from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as path
import coral_plotter as plotter


basedir = '/0work/projects/coral/'
basedir = '/home/pkalmus/projects/coral/'
do_plot = True

lat_extent = 35.5
#lat_extent = 0.5
#climatology_start_year = 1985 #earliest is 1985
#climatology_end_year = 2008 #latest is 2012

#gridspace = 0.5
#geolimit=([-35.5,35.5], [0,360])
gridspace = 0.06 # need to use the clim_ds grid space
geolimit=([-35.5,35.5], [0.5,355.5])
#geolimit=([-40,-2], [100,170]) # Australia

#def get_mmm_clim(point):
#    '''
#    Tried commenting this function out and just returning an integer; it's not the slow part.
#    '''
#    mylon = point[0]
#    mylat = point[1]
#    myclim = clim_ds.sel(lat=slice(mylat-gridspace/2, mylat+gridspace/2), lon=slice(mylon-gridspace/2, mylon+gridspace/2))
#    if np.nanmean(myclim.reef_mask.values) > 0:
#        myclim = myclim.where(myclim.reef_mask==1, drop=True) 
#        clim_vect=[]
#        mmm_vect =[]
#        for thislat in myclim.lat.values:
#            for thislon in myclim.lon.values:
#                thisind = myclim.sel(lat=thislat, lon=thislon).mmm_month.values - 1 # goes from 0 to 12, but I think 0 indicates fill val??
#                if ~np.isnan(thisind):
#                    mmm = myclim.sel(lat=thislat, lon=thislon).clim_monthly.values[int(thisind)]
#                    mmm_vect.append(mmm)
#                    clim = np.mean(myclim.sel(lat=thislat, lon=thislon).clim_monthly.values) # take the mean of all 12 months
#                    clim_vect.append(clim)
#        return np.mean(mmm_vect), np.mean(clim_vect)
#    else:
#        return np.nan, np.nan

def has_reef(mylon, mylat):
    #myreefs = clim_ds.sel(lat=slice(mylat-gridspace/2.0, mylat+gridspace/2), lon=slice(mylon-gridspace/2.0, mylon+gridspace/2))
    myreefs = clim_ds.sel(lat=slice(mylat-gridspace/2.0, mylat+gridspace/2.0), lon=slice(mylon-gridspace/2.0, mylon+gridspace/2.0)) # can't use "tolerance" with slice.
    return np.nanmean(myreefs.reef_mask.values) > 0.0

# read in climatology. want MMM for each location. 
clim_file = basedir + 'data/climatology/noaa_crw_thermal_history_climatology.nc'
clim_ds  = xr.open_dataset(clim_file) # lon: -180 to 180. lat: 34.14 to -35.27 (so in reverse order... be careful.)
#clim_ds = clim_ds.assign_coords(lon=((clim_ds.lon + 360) % 360)).sortby('lon')
clim_ds = clim_ds.assign_coords(lat=clim_ds.lat).sortby('lat') # fixed
#clim_ds = clim_ds.sel(lat=slice(-1*lat_extent, lat_extent))

# start making the new dataset with the reef_mask
mur_file = basedir + 'data/mur/201912-JPL-L4-SSTfnd-MUR_monthly-GLOB-fv04.2.nc'
mur_ds  = xr.open_dataset(mur_file)
# drop unneeded variables
mur_ds = mur_ds.drop('monthly_mean_sst')
mur_ds = mur_ds.drop('monthly_std_sst')
mur_ds = mur_ds.sel(lat=slice(-1*lat_extent, lat_extent))

# make reef_mask, initialized to zeros.
null_data = np.empty((mur_ds.dims['lat'], mur_ds.dims['lon'])) # already initialized with zeros 
reef_mask_da = xr.DataArray(null_data, coords=[mur_ds.coords['lat'], mur_ds.coords['lon']], dims=['lat', 'lon'])  


# instead of loopoing through the 0.01 resolution data, which would take about 24 hours, a nearest neighbor interp
interpolated = clim_ds.reef_mask.interp_like(reef_mask_da, method='nearest')

# save netcdf
netcdf_name = 'coral_locations_mur.nc'
#encoding cuts file size from 2 GB to 245 MB. encoding lat and lon to float32 prevents extra decimal garbage.
interpolated.to_netcdf(netcdf_name, encoding={'reef_mask': {'dtype': 'int8', '_FillValue': -9999}, 'lat': {'dtype': 'float32', '_FillValue': -9999}, 'lon': {'dtype': 'float32', '_FillValue': -9999} } )  
print(netcdf_name)

# diagnostic plots. note the 0.01 resolution plot takes a few minutes.
print('sanity check plots...')
plotter.plot_data_array(interpolated.sel(lat=slice(-40,2), lon=slice(100,170) ), 'reef_grid_mur.png')
plotter.plot_data_array(clim_ds.reef_mask.sel(lat=slice(-40,2), lon=slice(100,170) ), 'reef_grid_clim.png')




