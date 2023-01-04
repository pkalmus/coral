#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

    conda activate geo_env

regrid_bcdp: Get model outputs using BCDP, regrid to 1x1°, and save on the local file system for regrid_reef.py to digest.
coral_locations: get reef-containing points (don't need to run every time)

regrid_reef: Standardize time dimension, find reef an neighbor locations (using previously made reef_grid_1x1.txt and nn_grid_1x1.txt),
   interpolate to missing coastal grid points, and repackage into 1-D "ravel" dataset with "index".
weighting
coarse_project

Filter out non-reef locations and save new netcdf files in "ravel" form.
Only for 1x1 regular gridded data (such as the HadISST, and some of the CMIP6 gr models.)
Read in model or obs.


TODO: shore interpoloation.
TODO: Proper time interpolation. Doing the master_time replacement with /365.25+1850 is reasonable at least to 2100: 3.5 days.

TODO: make plots to validate the time changing, as it is very delicate. It could scramble the whole analysis.
    At the end. Pick one model, and HadISST, and make a plot, at a random location, of the raw time series (native time dims) and the reef timeseries.
TODO: make maps to validate the shore interpolation. 


Needs to know about the different formats from different models
    Put them all on the same "days from" time standard
    three "types"
Start them all from the same date (the latest date whatever that is)




@author: pkalmus
"""

import xarray as xr
import pandas as pd
import numpy as np
import pdb
import sys
import os.path
import subprocess
from sklearn.utils.extmath import cartesian
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import coral_plotter as plotter
#import nc_time_axis   conda install -c conda-forge nc-time-axis    https://github.com/SciTools/nc-time-axis


# read in user params
import importlib
paramfile = sys.argv[1]   
params = importlib.import_module(paramfile)
projectdir = params.projectdir
tosrunname = params.tosrunname
basedir = params.basedir
pythondir = params.pythondir
listfilename = params.listfilename
scenarios = params.scenarios
dryRun = params.dryRun
oneMember = params.oneMember
year_cutoff = params.year_cutoff
geolimit = params.geolimit

#useHadISST = False # if true, include HadISST in the analysis

# special import from pySSDF, needs to come after "projectdir" defined
sys.path.append(projectdir+'/pySSDF/src/')
import make_s

# define in and out dirs
indir = basedir+'data/tos/%s/raw/' % (tosrunname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
outdir = basedir+'data/tos/%s/reef/' % (tosrunname)

sysStr = "mkdir -p %s" % outdir
subprocess.call(sysStr, shell=True)

# load up reef points and neighbor points
prediction_grid = np.loadtxt('reef_grid_1x1.txt') # col1 lon,    col2 lat
nn_grid = np.loadtxt('nn_grid_1x1.txt')
both_grid = np.vstack([prediction_grid, nn_grid])
index = np.arange(both_grid.shape[0])
nreef = len(prediction_grid)

with open(listfilename) as f:
    models = f.read().splitlines()

master_time = None # we are just replacing time axis. probably safer to interpolate instead. this is effectively nearest-neighbor interp.
#timefilename = basedir+'data/tos/gr/tos_Omon_GFDL-CM4_historical_ssp585_r1i1p1f1_gr_185001-210012.nc' # "days since 1850-01-01 00:00:00" noleap
#timefilename = basedir+'data/tos/nature/raw/ACCESS-CM2_r1i1p1f1_gn_ssp126.nc' # hours since 1850-01-16 12:00:00... 0, 3012 proleptic_gregorian
timefilename = basedir+'data/tos/nature/raw/UKESM1-0-LL_r1i1p1f2_gn_ssp126.nc' # days since 1850-01-16 00:00:00.000000..., 0, 3012, 360_day
timeds = xr.open_dataset(timefilename, decode_times=True) 
timeds = timeds.sel(time=slice(None, str(2099)) ) # length 3000, goes up to December 2099
#master_time = timeds.time/(365.) + 1850 # the calendar used by GFDL-CM4 wants 365 not 365.25. This puts it in years
master_time = timeds.time # datetimeds

def regrid_and_save(latsravel, lonsravel, mytime, tos_all, outfilename, do_latitudinal_interp=False):
        '''
        check to see if each point is on prediction_grid   
        save the netcdf file
        
        the ordering of locations relative to the index must be consistent across models, HadISST, and climatology.
        '''
        tosravel = tos_all.reshape(tos_all.shape[0], -1) 
        sstos = np.array(list(zip(*[lonsravel,latsravel])) ) #need list() in Python 3
    
        # find all the reef points
        ind_plus_negative_ones = make_s.ismember(sstos,prediction_grid.copy()) # the -1 values indicate "not a member"
        ind = np.where(ind_plus_negative_ones > -1)[0]
        sstos_reef = sstos[ind]
        tos_reef = tosravel[:,ind]
        tos_reef_swap = np.swapaxes(tos_reef, 0, 1) 
        lons_reef = sstos_reef[:,0]
        lats_reef = sstos_reef[:,1]
    
        # find all the nn points
        ind_plus_negative_ones = make_s.ismember(sstos, nn_grid.copy()) # the -1 values indicate "not a member"
        ind = np.where(ind_plus_negative_ones > -1)[0]
        sstos_nn = sstos[ind]
        tos_nn = tosravel[:,ind]
        tos_nn_swap = np.swapaxes(tos_nn, 0, 1) 
        lons_nn = sstos_nn[:,0]
        lats_nn = sstos_nn[:,1]
    
        tos_all = np.vstack([tos_reef_swap, tos_nn_swap]) # 4836x2220
        lons_all = np.hstack([lons_reef, lons_nn])        # 4836
        lats_all = np.hstack([lats_reef, lats_nn])
        
        # Here are the NaNs that have reef or are nearest neighbors.
        #plotter.scatter(lons_all, lats_all, outfilename.replace('.nc', '_pixels.png'))
        nanind = np.where(np.isnan(tos_all[:,0]))[0] #156 points, all along coasts.
        print('number of reefs & nn:', len(index))
        print('number of nan:', len(nanind))
        #plotter.scatter(lons_all[nanind], lats_all[nanind], outfilename.replace('.nc', '_nanpixels.png'))        
        
        if do_latitudinal_interp:
            # I don't think I'm going to implement this. Too many things can go wrong.
            # I suspect the relationships between these pixel central temperatures could be very nonlinear.
            # Could talk about it with Mark, or maybe Ian.
            # Note that currently, requiring at least 10 non-NaN coarse pixels, reduces
            # the number of fine pixels by 16%. 828639/989932. = 0.837
            raise NotImplementedError() 
            # # To do latitudinal extrapolation, get East/West nearest neighbors.
            # # Try to vectorize the interp over the time dimension or else it may be too slow
            # for ind in nanind:
            #     mylon = lons_all[ind]
            #     mylat = lats_all[ind]
            #     mytos = tos_all[ind,:]
                
            #     lonE = mylon-1
            #     indE = np.intersect1d(np.where(lons_all==lonE)[0], np.where(lats_all==mylat)[0])
            #     tosE = np.ones(mytos.shape)*np.nan
            #     if len(indE)==1:
            #         tosE = np.squeeze(tos_all[indE,:])

            #     lonEE = mylon-2
            #     indEE = np.intersect1d(np.where(lons_all==lonEE)[0], np.where(lats_all==mylat)[0])
            #     tosEE = np.ones(mytos.shape)*np.nan
            #     if len(indEE)==1:
            #         tosEE = np.squeeze(tos_all[indEE,:])
                    
            #     lonW = mylon+1
            #     indW = np.intersect1d(np.where(lons_all==lonW)[0], np.where(lats_all==mylat)[0])
            #     tosW = np.ones(mytos.shape)*np.nan
            #     if len(indW)==1:
            #         tosW = np.squeeze(tos_all[indW,:])

            #     lonWW = mylon+2
            #     indWW = np.intersect1d(np.where(lons_all==lonWW)[0], np.where(lats_all==mylat)[0])
            #     tosWW = np.ones(mytos.shape)*np.nan
            #     if len(indWW)==1:
            #         tosWW = np.squeeze(tos_all[indWW,:])

            #     pdb.set_trace()
            
        
        myyears = np.array([item.year for item in mytime])
        mymonths = np.array([item.month for item in mytime])
        myyearfrac = myyears + (mymonths-0.5)/12.0
        
        isreef = np.zeros(len(index))
        isreef[0:nreef]=1
        lon_da = xr.DataArray(lons_all, coords=[index], dims=['index']) # why do these need square brackets?
        lat_da = xr.DataArray(lats_all, coords=[index], dims=['index'])
        year_da = xr.DataArray(myyearfrac, coords=[myyearfrac], dims=['time'])
        da = xr.DataArray(tos_all, coords=(index, myyearfrac), dims=('index', 'time')) # why don't these need square brackets?
        isreef_da = xr.DataArray(isreef, coords=[index], dims=(['index']))
        ds = xr.Dataset({'tos':da, 'isreef':isreef_da, 'year':year_da, 'lon':lon_da, 'lat':lat_da})
        
        ds.to_netcdf(outfilename) 
        print(outfilename)



# HadISST
filename = '/raid8/pkalmus/data/coral/data/tos/raw/HadISST_sst.nc'
print(filename)
outfilename = outdir+'HadISST_sst_reef.nc'
modds_all = xr.open_dataset(filename, decode_times=False)  # time:units = "days since 1870-1-1 0:0:0" ; length:1794, first element: 15.5 gregorian

print(modds_all.time.units)
print(modds_all.time.calendar)
print(len(modds_all.time))

# x, y for time test plot. 
x = modds_all.time.values
x = x / 365.0 + 1870
y = modds_all.sst.values[:,40,40]

time_copy = master_time.copy()
time_had = modds_all.time
ntimes = len(time_had.values)    
modds_all = modds_all.assign_coords(time=time_copy[240:240+1794])  

years = np.array([item.year for item in modds_all.time.values])
months = np.array([item.month for item in modds_all.time.values])
days = np.array([item.day for item in modds_all.time.values])
x2 = years + months/12. + days/365.
y2 = modds_all.sst.values[:,40,40]
plotter.plot(x, y, outfilename.replace('.nc', '_testtime.png'), x2=x2, y2=y2, xlabel='time', ylabel='tos', title=None,  xlim=[1980,1985], ylim=None)

modds_all = modds_all.assign_coords(latitude=modds_all.latitude).sortby('latitude') # put into ascending
modds_all = modds_all.assign_coords(longitude=((modds_all.longitude + 360) % 360)).sortby('longitude') # put into 0, 360.

modds = modds_all.sel(latitude=slice(geolimit[0][0], geolimit[0][1]), time=slice(str(year_cutoff),None) )
t = modds.time.values       # "days since 1870-1-1 0:0:0"
lats = modds.latitude.values
lons = modds.longitude.values
tos_all = modds.sst.values # time, latitude, longitude
latsravel = np.tile(lats, (len(lons), 1)).transpose().ravel()
lonsravel = np.tile(lons, (len(lats), 1)).ravel() 

regrid_and_save(latsravel, lonsravel, t, tos_all, outfilename) # no interp for obs


for (modelind, model) in enumerate(models):
    print('\n\n***********'+str(modelind+1) + '   ' + model)
    
    # if modelind < 16:
    #     continue
    
    # if model != 'KACE-1-0-G r2i1p1f1 gr': # needs regridding
    #     continue
    # if model != 'KACE-1-0-G r3i1p1f1 gr': # needs regridding
    #     continue
    # if model != 'INM-CM4-8 r1i1p1f1 gr1': #OK
    #     continue    
    # if model != 'CESM2-WACCM r1i1p1f1 gr': #OK
    #     continue   
    # if model != 'GFDL-ESM4 r1i1p1f1 gr': #OK 4836 reefs & nn, 493 nans
    #     continue       
    # if model != 'GFDL-ESM4 r1i1p1f1 gn': #OK 4836 reefs & nn, 1182 nans
    #     continue       
    
    if model == 'CNRM-ESM2-1 r4i1p1f2 gn': # hada non-monotonic time!
        print('DELETE MODEL FROM %s: %s' % (listfilename, model))
        continue
    
    if model == 'IITM-ESM r1i1p1f1 gn': # hours since 1900-01-16 12:00:00.000000
        print('DELETE MODEL FROM %s: %s' % (listfilename, model))
        continue
      
    if model == 'KACE-1-0-G r1i1p1f1 gr': # ??? 
        print('DELETE MODEL FROM %s: %s' % (listfilename, model))
        continue
        
    
    for scenario in scenarios:
        # e.g. CanESM5_r12i1p1f1_ssp126.nc
        modelstr = model.replace(' ', '_')
        modelstr = modelstr+'_'+scenario
        filename = indir+modelstr+'.nc'
        #print(filename)
        outfilename = outdir+modelstr+'_reef.nc'
        #print(outfilename)
        
        # if os.path.exists(outfilename):
        #     print('%s already exists, continuing.' % (outfilename))
        #     continue
                        
        modds_all = xr.open_dataset(filename, decode_times=False) 
        # these attributes are only preserved if using decode_times = False
        units = modds_all.time.units
        print(units)
        print(modds_all.time.calendar)
        print(modds_all.time.values[0]) # they all start at 0.
        print(len(modds_all.time))

# timeds: days since 1850-01-01 00:00:00, 15.5, 3012
    
# hours since 1850-01-16 12:00:00   noleap, proleptic_gregorian
# hours since 1850-01-15 13:00:00.000000
# hours since 1850-01-17 00:00:00.000000
# days since 1850-01-16 00:00:00.000000   , 360_day     

# 3000

# /home/pkalmus/projects//coral/data/tos/nature/raw/CanESM5_r1i1p1f1_gn_ssp126.nc
# hours since 1850-01-16 12:00:00.000000
# 5412

# /home/pkalmus/projects//coral/data/tos/nature/raw/CanESM5_r1i1p1f1_gn_ssp585.nc
# hours since 1850-01-16 12:00:00.000000
# 5412
    
        # x, y for time test plot. 
        x = modds_all.time.values
        if modds_all.time.calendar[0] == 'd': #days since...
            x = x + 16
            x = x / (365.0) + 1850
        else:
            x = x + (16*24.0)
            x = x / (365.0*24.0) + 1850
            
        y = modds_all.tos.values[:,40,40]
    
        modds_all = xr.open_dataset(filename, decode_times=True) 
        modds_all = modds_all.sel(time=slice(None, str(2099)) ) # length 3000, goes up to December 2099
        try:
            modds_all = modds_all.assign_coords(time=master_time.copy())   

            years = np.array([item.year for item in modds_all.time.values])
            months = np.array([item.month for item in modds_all.time.values])
            days = np.array([item.day for item in modds_all.time.values])
            x2 = years + (months-1)/12. + days/365.
            y2 = modds_all.tos.values[:,40,40]
            plotter.plot(x, y, outfilename.replace('.nc', '_testtime.png'), x2=x2, y2=y2, xlabel='time', ylabel='tos', title=None,  xlim=[1850,1855], ylim=None)

        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            #raise
            pdb.set_trace()
            continue

        modds = modds_all.sel(lat=slice(geolimit[0][0], geolimit[0][1]), time=slice(str(year_cutoff),None) )
        t = modds.time.values
        lats = modds.lat.values
        lons = modds.lon.values
        tos_all = modds.tos.values # time, latitude, longitude
        latsravel = np.tile(lats, (len(lons), 1)).transpose().ravel()
        lonsravel = np.tile(lons, (len(lats), 1)).ravel() 

        regrid_and_save(latsravel, lonsravel, t, tos_all, outfilename)
        
        
        
        
        
        
        
        
        
