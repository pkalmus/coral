#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

conda activate base
conda activate geo_env

Read in model or obs.

For each time step:
    send data, zero bias, output grid to interpolator
    save in the output file

Needs to know about the different formats from different models
    Put them all on the same "days from" time standard
    three "types"
Start them all from the same date (the latest date whatever that is)


lat lon time
ascending lat
0,360 lon
C

todo: only use basis functions at about 1 degree
todo: test the regridding. look at 1880 maps for all 5 of them.
todo: all 4 models

todo: weights for all 4 models (maybe separate .nc files)
    
todo: output grid ocean only


todo: put them all at start of month or mid month. spatiotemporal kriging?

@author: pkalmus
"""

projectdir = '/0work/projects/'
projectdir = '/home/pkalmus/projects/'

import xarray as xr
import pandas as pd
import numpy as np
import pdb
import sys
sys.path.append(projectdir+'/pySSDF/src/')
import interpolator
from sklearn.utils.extmath import cartesian
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

basedir = projectdir+'/coral/'
Data_dir = basedir+'data/tos/'


modelind = 1 
do_plot = False
do_nc_grid = False # if true, lat,lon,tos are on regular grid; if false, index (reef locations only)
models = ['HadISST_sst', 'tos_Omon_CMIP6_CanESM5_historical_ssp585_r1i1p2f1_gn_185001-210012', 'tos_Omon_bcc-csm1-1_rcp85_000', 'tos_Omon_ACCESS1-0_rcp85_000'] # how were these made - explicitly with historical too?
years_to_add = [1870, 1850, 1850, 0]

#0 = HadISST (latitude, sst),latitude descends, lon [-180,180] _FillValue = -1.e+30f, C
#1 = CMIP5 tos time lat lon, lat ascends,       lon [0,360] _FillValue = 1.e+20f, K
#2 = CMIP6 tos time j i latitude(i,j) longitude(i,j). i[0, 360] j[0,290]. lon [0,360] _FillValue = 1.e+20f, C
# i is x/longitude index (not a 1-to-1 mapping), j is y/latitude index (appears to be 1-to-1 mapping)
ds_type = [0,2,1,1] 
# tos_Omon_CanESM5_historical_r1i1p2f1_gn_185001-201412.nc tos_Omon_CanESM5_ssp585_r1i1p2f1_gn_201501-210012.nc
# tos_Omon_CanESM5_historical_ssp585_r1i1p2f1_gn_185001-210012.nc

prune_percent = .5
matching_radius = .5
#basis_file = 'lvl_1_2_3.csv'
basis_file = 'lvl_2_3_5.csv'
#basis_file = 'lvl_5.csv' # looks like crap. amplified noise?
#basis_file = 'lvl_2_3_6.csv'
# 3_5
# 5
# 5_6

year_cutoff = 1880 # throw out all years before this
year_cutoff = 1910
ntimes = 1030

geolimit=([-35,35], [0,360])
#geolimit=([-40,-2], [100,170]) # Australia

#gridspace = 0.5
#gridlat = np.arange(geolimit[0][1], geolimit[0][0]-gridspace, -1*gridspace)
#gridlon = np.arange(geolimit[1][0], geolimit[1][1]+gridspace, gridspace) 
#pre_grid = np.fliplr(cartesian((gridlat,gridlon)))

prediction_grid = np.loadtxt('reef_grid.txt')
index = np.arange(prediction_grid.shape[0])
        

#null_data = np.empty((obs_ds.dims['latitude'], obs_ds.dims['longitude']))
#null_data[:] = np.nan

model = models[modelind]
filename = Data_dir+'raw/'+model+'.nc'
if do_nc_grid:
    outfilename = Data_dir+'regrid/'+model+'_regrid.nc'
else:
    outfilename = Data_dir+'regrid/'+model+'_reef.nc'
modds_all = xr.open_dataset(filename, decode_times=False) # just keep as "days since 1870-1-1 0:0:0", don't go to datetime64
modds_all = modds_all.assign_coords(time=(modds_all.time/365. + years_to_add[modelind]))

if ds_type[modelind] == 0:
    modds_all = modds_all.assign_coords(latitude=modds_all.latitude).sortby('latitude') # put into ascending
    modds_all = modds_all.assign_coords(longitude=((modds_all.longitude + 360) % 360)).sortby('longitude') # put into 0, 360.
    modds = modds_all.sel(latitude=slice(geolimit[0][0], geolimit[0][1]), time=slice(year_cutoff,None) )
    t = modds.time.values       # "days since 1870-1-1 0:0:0"
    lats = modds.latitude.values
    lons = modds.longitude.values
    tos_all = modds.sst.values # time, latitude, longitude
    # get tos, lat, and lon in the shape required by interpolator
    latsravel = np.tile(lats, (len(lons), 1)).transpose().ravel()
    lonsravel = np.tile(lons, (len(lats), 1)).ravel() 
elif ds_type[modelind] == 1:
    modds = modds_all.sel(lat=slice(geolimit[0][0], geolimit[0][1]), time=slice(year_cutoff,None) )
    t = modds.time.values       # "days since XX
    lats = modds.lat.values
    lons = modds.lon.values
    tos_all = modds.tos.values-273.15 # time, latitude, longitude
    # get tos, lat, and lon in the shape required by interpolator
    latsravel = np.tile(lats, (len(lons), 1)).transpose().ravel()
    lonsravel = np.tile(lons, (len(lats), 1)).ravel()  
elif ds_type[modelind] == 2:
    # i,j are the coords. 
    modds = modds_all.sel(time=slice(year_cutoff,None) )
    t = modds.time.values       # "days since XX
    lats = modds.latitude.values
    lons = modds.longitude.values
    tos_all = modds.tos.values # time, latitude, longitude
    minind = np.max(np.where(lats<geolimit[0][0])[0])
    maxind = np.min(np.where(lats>geolimit[0][1])[0])
    lats = lats[minind:maxind,:]
    lons = lons[minind:maxind,:]
    tos_all = tos_all[:,minind:maxind,:]

    
    mylons = lons[1,:]
    sortind = np.argsort(mylons)
        
    lons = lons[:,sortind]
    lats = lats[:,sortind]
    tos_all = tos_all[:,:,sortind] #lats x lons

    if geolimit[1][0] > 0 and geolimit[1][1] < 360:
        minind = np.max(np.where(lons<geolimit[1][0])[1])
        maxind = np.min(np.where(lons>geolimit[1][1])[1])
        lats = lats[:,minind:maxind]
        lons = lons[:,minind:maxind]
        tos_all = tos_all[:,:,minind:maxind]    
    
    
    # get tos, lat, and lon in the shape required by interpolator
    latsravel = lats.ravel()
    lonsravel = lons.ravel()  

  
da_list = []
for loopind, looptime in enumerate(t):
    print(looptime)
    if loopind >= ntimes:
        continue
    tos = tos_all[loopind,:,:]
    tosravel = np.squeeze(tos).ravel()               
    missingind = np.where(tosravel == -1000)[0] 
    latsravel[missingind] = np.nan
    lonsravel[missingind] = np.nan
    tosravel[missingind] = np.nan
      
    nonanind = ~np.isnan(tosravel)
    latsravel_nonan = latsravel[nonanind]
    lonsravel_nonan = lonsravel[nonanind]
    tosravel_nonan = tosravel[nonanind]
    
    sstos = np.array(list(zip(*[lonsravel_nonan,latsravel_nonan])) ) #need list() in Python 3
    sstos_plot = sstos.copy()
    sigma2 = np.ones(tosravel_nonan.shape)
       
    print('kriging...')                     
    s_output, Y_bar, error = interpolator.interpolate(sstos, tosravel_nonan, sigma2, basis_file, prediction_grid, matching_radius=matching_radius)    
    print('...done!')

    
    if do_plot:
        # this only works if regular grid! (the reshape step).
        
        # maybe undo the ravel first? or maybe cartopy can handle point / ravel data?
        latDF = np.flipud(pd.unique(s_output[:,1]))
        lonDF = pd.unique(s_output[:,0])
        
        # len(lon)*len(lat) should equal Y_bar.shape(0).
        # you cannot plot satellite-type data in this way; it is for gridded data.
        imageDF = np.reshape(Y_bar, [len(lonDF),len(latDF)], 'F') # 'F' means 0th dimension varies fastest (Fortran-like)
        imageDF = np.fliplr(imageDF)  
        imageDF = imageDF.T

        # make and save maps of raw and regridded (oceans only)... maybe just 1st time step
        # common colorbar
        alldata = np.concatenate((np.squeeze(Y_bar), tosravel_nonan))
        allmin = np.nanmin(alldata)
        allmax = np.nanmax(alldata)
        maskcmap = plt.get_cmap('jet')
    
        mapfile='5auresmoo_'+model+'_'+str(looptime)+'.png'
        #plotdata = ssdf_plotter.PlotInfo(sstos, tosravel_nonan, 'raw', issat=True)
        #ssdf_plotter.save_maps(plotdata, s_output, Y_bar, error, geolimit, mapfile, maintitle='', maskland_flag=True) 
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180.0))
        #img_extent = (-180, 180, geolimit[0][0], geolimit[0][1])
        #dep_da.plot.contourf(ax=ax, extent=img_extent, transform=ccrs.PlateCarree())
        plt.contourf(lonDF, latDF, imageDF, 60, cmap=maskcmap, vmin=allmin, vmax=allmax, transform=ccrs.PlateCarree())
        ax.coastlines()
        plt.savefig(mapfile) 
        plt.close()


        mapfile='5auregrid_'+model+'_'+str(looptime)+'.png'
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180.0))
        #img_extent = (-180, 180, geolimit[0][0], geolimit[0][1])    
        plt.scatter(s_output[:,0], s_output[:,1], s=3, c=np.squeeze(Y_bar), cmap=maskcmap, edgecolor='', vmin=allmin, vmax=allmax, transform=ccrs.PlateCarree()) #, cmap=maskcmap, vmin=allmin, vmax=allmax, edgecolor='') 
        ax.coastlines()
        plt.savefig(mapfile) 
        plt.close()


        mapfile='5auraw_'+model+'_'+str(looptime)+'.png'
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180.0))
        #img_extent = (-180, 180, geolimit[0][0], geolimit[0][1])
        lat = sstos_plot[:,1]
        lon = sstos_plot[:,0]       
        plt.scatter(lon, lat, s=11, c=tosravel_nonan, cmap=maskcmap, edgecolor='', vmin=allmin, vmax=allmax, transform=ccrs.PlateCarree()) #, cmap=maskcmap, vmin=allmin, vmax=allmax, edgecolor='') 
        ax.coastlines()
        plt.savefig(mapfile) 
        plt.close()
        pdb.set_trace()
        
    # add to datarray for this model
    # note: coords must be lists, even if containing only a single value
    if do_nc_grid:
        da = xr.DataArray(imageDF[:,:,None], coords=(latDF, lonDF, [looptime]), dims=('lat', 'lon', 'time'))
    else:
        data = np.squeeze(Y_bar)
        da = xr.DataArray(data[:,None], coords=(index, [looptime]), dims=('index', 'time')) # why don't these need square brackets?
    da_list.append(da)
            
# save a netcdf file with XARRAY, one per model, lat lon pval
da = xr.concat(da_list, dim='time')

if do_nc_grid:
    ds = xr.Dataset({'tos': da})
else:   
    lon_da = xr.DataArray(s_output[:,0], coords=[index], dims=['index']) # why do these need square brackets?
    lat_da = xr.DataArray(s_output[:,1], coords=[index], dims=['index'])
    ds = xr.Dataset({'tos':da, 'lon':lon_da, 'lat':lat_da})


ds.to_netcdf(outfilename) 
print(outfilename)
pdb.set_trace()


