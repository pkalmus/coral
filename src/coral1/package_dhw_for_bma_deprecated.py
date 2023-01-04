#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Deprecated, because weight_dhw.py now saves obs annd mod dhwmaxmean in the files.
Use package_dhw_for_bma.py

Cycle through each model and the obs, save a matrix for the BMA.

1st attempt: the 1980-2000 mean (or whatever values came from params)

@author: pkalmus
"""

#projectdir = '/0work/projects/'
projectdir = '/home/pkalmus/projects/'

import sys
import xarray as xr
import numpy as np
import pdb
import coral_plotter as plotter
import subprocess
from scipy.io import savemat

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
start_year = params.start_year
nmonths = params.nmonths
meanstart=params.meanstart
meanend=params.meanend


# define in and out dirs
indir = basedir+'data/tos/%s/reef/' % (runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
outdir = basedir+'data/tos/%s/weightem/' % (runname)

sysStr = "mkdir -p %s" % outdir
subprocess.call(sysStr, shell=True)



# HadISST_sst_reef.nc  tos_Omon_CMIP6_CanESM5_historical_ssp585_r1i1p2f1_gn_185001-210012_reef.nc

# Obs: matrix, num, [1:1024, 1:2] t in fractional years, mysst in K.
# read in the Obs matrix
# send to R
# check and see if it's the same
obs_filename = indir+'HadISST_sst_reef.nc'
obs_ds = xr.open_dataset(obs_filename, decode_times=False)

# calculate obs mmm climatology from HadISST
clim_ds = obs_ds.sel(time=slice(1960, 1980))
# take annual mean
clim_tos = clim_ds.tos.values
clim_mo = np.nanmean(clim_tos.reshape(clim_tos.shape[0],-1,12), 1)
# find mmm for each index
clim_mmm = np.max(clim_mo, 1)


# limit to times
obs_ds =  obs_ds.sel(time=slice(meanstart, meanend))

with open(listfilename) as f:
    models = f.read().splitlines()

def get_dhw(mmm, sst_time_series):
    # calculate DHW. 
    dhw_mo = sst_time_series - mmm
    dhw_mo = dhw_mo*4.34 #convert from months to weeks
    dhw_mo[dhw_mo < 0] = 0
    dhw_mo = np.insert(dhw_mo, 0, np.array([0.0,0.0]))
    dhw = dhw_mo[2:]+dhw_mo[1:-1]+dhw_mo[0:-2] # 3-month running sum
    dhw = np.insert(dhw, 0, np.array([np.nan])) # shift to right once
    dhw = dhw[0:-1] 
    return dhw


obs_means = np.empty(obs_ds.index.values.shape)
obs_means[:] = np.nan
mod_matrix = np.empty((len(models), len(obs_ds.index.values))) # 127 x inndex
mod_matrix[:] = np.nan
for loopind in obs_ds.index.values: # validation loop
    print(loopind)

    # if loopind > 2:
    #     break
    
    #mmm = lon_lat_clim[loopind,2]
    mmm = clim_mmm[loopind]
    if ~np.isnan(mmm):
        obs = obs_ds.tos.values[loopind,:]
        # check if it has tos (is ocean)
        if np.isnan(obs[0]):
            print('tos nan')
            continue

        
        for (model_ind, model) in enumerate(models):
            # if model_ind > 1:
            #     continue
            modelstr = model.replace(' ', '_')
            filename = indir+modelstr+'_ssp126_reef.nc'
            print(filename)
            mod_ds = xr.open_dataset(filename, decode_times=False)
        
            # limit to times
            mod_ds =  mod_ds.sel(time=slice(meanstart, meanend))   
            
            # replace any zeros, at least until I fix this at the source in the BCDP code.
            #print(np.where(mod_ds.tos.values==0)[0])
            mod_ds = mod_ds.where(mod_ds['tos'] != 0.) # btw this takes about half a second to do - and might add more to concat / average
            mod = mod_ds.tos.values[loopind,:]    
            
            
            obs_dhw_ts = get_dhw(mmm, obs)
            mod_dhw_ts = get_dhw(mmm, mod)
           
            # select the maximum for each year (every 12 data points)
            obsmax = np.nanmax(obs_dhw_ts.reshape(-1, 12), axis=1)
            modmax = np.nanmax(mod_dhw_ts.reshape(-1, 12), axis=1)
            
            # take the mean of the annual maxima
            obsmean = np.nanmean(obsmax)
            modmean = np.nanmean(modmax)

            obs_means[loopind] = obsmean
            mod_matrix[model_ind, loopind] = modmean

# save a version before the less than N check
mdic = {'mod_matrix': mod_matrix, 'obs_vect': obs_means}
savemat('obs_mod_dhwmax_all_%i-%i.mat' % (meanstart, meanend), mdic)

# any index with fewer than 10 models w/o nan, have to nan out. e.g. 1798 has only one. 
for loopind in obs_ds.index.values: 
    num_model_estimates = np.count_nonzero(~np.isnan(mod_matrix[:,loopind]))
    if num_model_estimates < 10:
        print('too few model estimates %i for index %i' % (num_model_estimates, loopind))
        mod_matrix[:,loopind] = np.nan


# create the .mat package for the BMA. I think it's best to redo it.

# obs_means is a vect of mean max dhw values from 1980-2000, one per index.
# mod_matrix has mean max dhw values from 1980-2000, one per index, per model. 
# any index with fewer than 10 models gets nans for that index. 
mdic = {'mod_matrix': mod_matrix, 'obs_vect': obs_means}
savemat('obs_mod_dhwmax%i-%i.mat' % (meanstart, meanend), mdic)


