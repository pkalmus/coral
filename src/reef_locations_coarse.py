#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Get locations AND climatology for coarse grid as defined.
Also list of nearest-neighbor pixels to the coral locations, that doesn't include coral locations.

Climatolgy used in Hooidonk2016: 1982-2008. I have it from 1985-2012. Let's go from 1985-2008.

We have two ways to get coral locations:
    1. The NOAA climatology file, which is about 5km resolution.
    2. The WCM shapefile
    
Cartopy is for Python 3, but does not seem ready for prime time yet.
It's nice because you can make plots directly from xarray data objects (see example below, "reef_mask.plot...")
A good way to limit the map in Cartopy is to limit the dataset and then just plot that dataset.
But there is no easy way to do a land mask. Grid lines, labels, colorbars seem to still be an uphill struggle.


Quality check:
reef_grid_1x1.txt and reef_clim_1x1.txt are ordered the same: by increasing latitude    
These currently match the way the model regridding order.
When comparing these .txt files with e.g. indexed models, always do an explicit check (e.g. in coarse_project.py)

Some documentation:
https://xarray.pydata.org/en/stable/generated/xarray.plot.pcolormesh.html

from pdb: dir(reef_mask.where(reef_mask==1).plot())

https://matplotlib.org/api/pyplot_api.htmlf


Created on Wed Sep  4 12:29:02 2019

@author: pkalmus
"""

import xarray as xr
import pdb
from sklearn.utils.extmath import cartesian
import numpy as np
import coral_plotter as plotter

basedir = '/0work/projects/coral/'
basedir = '/home/pkalmus/projects/coral/'
finelocation_file = '/raid8/pkalmus/data/coral/data/location/14_001_WCMC008_CoralReefs2018_v4_1/01_Data/locationsWCMCv41.nc'
coarselocation_dir = '/raid8/pkalmus/data/coral/data/location/coarse/'

do_plot = False # be sure to use geo_env if you set to True!

lat_extent = 35

gridspace = 1.0
geolimit=([-35.5,35.5], [0.5,355.5])
#geolimit=([-5.5,5.5], [0.5,355.5])

# def get_mmm_clim_coral_only(point):
#     '''
#     This will need some updating, from the coral1 version.
#     To use mur_clim_ds.
#     '''
#     mylon = point[0]
#     mylat = point[1]
#     myclim = clim_ds.sel(lat=slice(mylat-gridspace/2.0, mylat+gridspace/2.0), lon=slice(mylon-gridspace/2.0, mylon+gridspace/2.0))

#     if np.nanmean(myclim.reef_mask.values) > 0:
#         myclim = myclim.where(myclim.reef_mask==1, drop=True) 
#         clim_vect=[]
#         mmm_vect =[]
#         for thislat in myclim.lat.values:
#             for thislon in myclim.lon.values:
#                 thisind = myclim.sel(lat=thislat, lon=thislon).mmm_month.values - 1 # goes from 0 to 12, but I think 0 indicates fill val??
#                 if ~np.isnan(thisind):
#                     mmm = myclim.sel(lat=thislat, lon=thislon).clim_monthly.values[int(thisind)]
#                     if mmm == 0.0:
#                         return np.nan, np.nan
#                     mmm_vect.append(mmm)
#                     clim = np.mean(myclim.sel(lat=thislat, lon=thislon).clim_monthly.values) # take the mean of all 12 months
#                     clim_vect.append(clim)
#         return np.mean(mmm_vect), np.mean(clim_vect)
#     else:
#         return np.nan, np.nan


def has_reef(point):
    '''
    Count up all the fine pionts in the fine locations ds, and see if there's a reef
    Want to include "edges"
    '''
    mylon = point[0]
    mylat = point[1]
    myreefs = fine_reef_ds.sel(lat=slice(mylat-gridspace/2, mylat+gridspace/2), lon=slice(mylon-gridspace/2, mylon+gridspace/2))
    return np.nanmean(myreefs.reef_mask.values) > 0.0

fine_reef_ds = xr.open_dataset(finelocation_file)

# climatology_start_year = 1985 #earliest is 1985
# climatology_end_year = 2008 #latest is 2012
# clim_file = basedir + 'data/climatology/noaa_crw_thermal_history_climatology.nc'
# clim_ds  = xr.open_dataset(clim_file) # lon: -180 to 180. lat: 34.14 to -35.27 (so in reverse order... be careful.)
# clim_ds = clim_ds.assign_coords(lon=((clim_ds.lon + 360) % 360)).sortby('lon')
# clim_ds = clim_ds.assign_coords(lat=clim_ds.lat).sortby('lat') # fixed
# clim_ds = clim_ds.sel(years=slice(climatology_start_year, climatology_end_year))

##########################################
# determine prediction grid with reefs
##########################################
# set up prediction grid. will this work for MUR?
gridlon = np.arange(geolimit[1][0], geolimit[1][1]+gridspace, gridspace) 
gridlat = np.arange(geolimit[0][0], geolimit[0][-1]+gridspace, gridspace)
pre_grid = np.fliplr(cartesian((gridlat, gridlon))) # cartesian orders by its first argument. this gives order the same as the indexed models in regrid_reef.py

# collect the points that have reefs
post_grid = None
for point in pre_grid:
    if has_reef(point):
        print(point)
        if post_grid is None:
            post_grid = np.array([point])
        else:
            post_grid = np.append(post_grid, [point], axis=0)    
# write out reef locations
lat = post_grid[:,1]
lon = post_grid[:,0]
if do_plot:
    plotter.scatter(lon, lat, coarselocation_dir+'reef_location_coarse.png')
reef_file = coarselocation_dir+'reef_location_coarse.csv'
np.savetxt(reef_file, post_grid, fmt='%1.2f')
print(reef_file)

# collect the points that are adjacent to reef points
# it's adjacent if the abs. min. difference of the nearest point for both lat and lon is 1 or 0, but not both 0.
nn_grid = None
for point in pre_grid:
    absdiff = np.abs(post_grid-point)
    diffsum = np.sum(absdiff, axis=1)
    ind = np.argmin(diffsum)

    if (absdiff[ind,0]<=gridspace) and (absdiff[ind,1]<=gridspace) and (diffsum[ind]!=0):
        print(point)
        if nn_grid is None:
            nn_grid = np.array([point])
        else:
            nn_grid = np.append(nn_grid, [point], axis=0)    
# write out nn locations
lat = nn_grid[:,1]
lon = nn_grid[:,0]
if do_plot:
    plotter.scatter(lon, lat, coarselocation_dir+'nn_location_coarse.png')
nn_file = coarselocation_dir+'nn_location_coarse.csv'
np.savetxt(nn_file, nn_grid, fmt='%1.2f')
print(nn_file)

pdb.set_trace()


# # read in climatology. want MMM for each location. 
# clim_file = basedir + 'data/climatology/mur_mmm.nc'
# mur_clim_ds  = xr.open_dataset(clim_file)
# mur_clim_ds = mur_clim_ds.assign_coords(lon=((mur_clim_ds.lon + 360) % 360)).sortby('lon')

# # reef_clim_coarse.csv is not used, but was in coral1 for val. we will make it anyway.
# both_grid = np.vstack([post_grid, nn_grid])
# clim_grid = None
# for point in both_grid:
#     mmm, clim = get_mmm_clim_coral_only(point)
#     print(point, mmm, clim)
#     if clim_grid is None:
#         clim_grid = np.array([[point[0], point[1], mmm, clim]])
#     else:
#         clim_grid = np.append(clim_grid, [[point[0], point[1], mmm, clim]], axis=0)
# np.savetxt('reef_clim_1x1.txt', clim_grid, fmt='%1.2f')
# if do_plot:
#     lat = both_grid[:,1]
#     lon = both_grid[:,0]
#     plotter.scatter(lon, lat, 'reef_mmm_1x1.png', c=clim_grid[:,2])
#     plotter.scatter(lon, lat, 'reef_clim_1x1.png', c=clim_grid[:,3])
# print(clim_grid.shape)

# pdb.set_trace()


