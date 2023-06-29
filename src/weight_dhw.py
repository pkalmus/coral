#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Run through the model list, and calculate the DHW mean annual max metric.
Also calculates the "arithmetic weight."


Take some years, find the max DHW of each year, and take the average of those.
Score: reciprocal of abs. difference, but with a maximum constraint.

After running this, run package_dhw_for_bma.py which produces the matrices for Elias to run BMA.
After he provides BMA, run weight_bma.py to create the final weights.
Then run coarse_project.py to create the files for Ayesha.

Climatology and mean windows:
Nominal: 1960-1979 for the climatology (using HadISST) - hardcoded here
         1980-1999 for the mean of DHW annual maxima 

python weight_dhw.py params 4. 
or call start_weights_NO_SSHX.sh

Note: replaced zeros in individual model timeseries with nan. 
These zeros need to be taken out at the source (BCDP stage)




Note that the DHW climatology used for both the obs and the models, is from the obs. 
There is no “bias correction” of the models relative to observation. 
Something to think about, the “right” way to do that. 

We are NOT taking into account whether a model has a hotter or colder trend (i.e. ECS). 


@author: pkalmus
"""

#projectdir = '/0work/projects/'
projectdir = '/home/pkalmus/projects/'

import sys
import xarray as xr
import numpy as np
import pdb
import subprocess

# read in user params
import importlib
paramfile = sys.argv[1]   
params = importlib.import_module(paramfile)
projectdir = params.projectdir
runname = params.runname
basedir = params.basedir
listfilename = params.modellist_filename
meanstart=params.meanstart
meanend=params.meanend
climstartyear = params.climstartyear_weight
climendyear = params.climendyear_weight

# define in and out dirs
indir = basedir+'data/tos/coral2/%s/reef/' % (runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
outdir = indir

sysStr = "mkdir -p %s" % outdir
subprocess.call(sysStr, shell=True)


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


# HadISST_sst_reef.nc  tos_Omon_CMIP6_CanESM5_historical_ssp585_r1i1p2f1_gn_185001-210012_reef.nc

# Obs: matrix, num, [1:1024, 1:2] t in fractional years, mysst in K.
# read in the Obs matrix
# send to R
# check and see if it's the same
obs_filename = indir+'HadISST_sst_reef.nc'
obs_filename = indir+'HadISST_sst_202110_reef.nc'
obs_ds = xr.open_dataset(obs_filename, decode_times=False)

# calculate obs mmm climatology from HadISST
clim_ds = obs_ds.sel(time=slice(climstartyear, climendyear))
# take annual mean
clim_tos = clim_ds.tos.values
clim_mo = np.nanmean(clim_tos.reshape(clim_tos.shape[0],-1,12), 1)
# find mmm for each index
clim_mmm = np.max(clim_mo, 1)

obs_ds =  obs_ds.sel(time=slice(meanstart, meanend))

with open(listfilename) as f:
    models = f.read().splitlines()

for (modelind, model) in enumerate(models):
    # if modelind == 0:
    #     continue

    print(str(modelind) + '   ' + model)
    
    # e.g. UKESM1-0-LL_r8i1p1f2_gn_ssp370_reef.nc
    #modelind = int(sys.argv[2])
    
    model = models[modelind]
    modelstr = model.replace(' ', '_')
    filename = indir+modelstr+'_ssp126_reef.nc'
    print(filename)


    pval_da = xr.DataArray(coords=[obs_ds.coords['index']], dims=['index'])    
    dhwmeanmax_mod_da = xr.DataArray(coords=[obs_ds.coords['index']], dims=['index'])   
    dhwmeanmax_obs_da = xr.DataArray(coords=[obs_ds.coords['index']], dims=['index'])   
    
    weightfilename = outdir+modelstr+'_weightdhw_%i-%i_%i-%i.nc' % (climstartyear, climendyear, meanstart, meanend)
    mod_ds = xr.open_dataset(filename, decode_times=False)
    
    # limit to start, end
    mod_ds =  mod_ds.sel(time=slice(meanstart, meanend))
    
    # replace any zeros, at least until I fix this at the source in the BCDP code.
    #print(np.where(mod_ds.tos.values==0)[0])
    mod_ds = mod_ds.where(mod_ds['tos'] != 0.)   
    
    
    for loopind in obs_ds.index.values: 
        # if loopind > 10:
        #     break
    
        mmm = clim_mmm[loopind]
        if ~np.isnan(mmm):
            obs = obs_ds.tos.values[loopind,:]
            # check if it has tos (is ocean)
            if np.isnan(obs[0]):
                #print('tos nan')
                continue
            obs_dhw_ts = get_dhw(mmm, obs)
            mod_dhw_ts = get_dhw(mmm, mod_ds.tos.values[loopind,:])
           
            # select the maximum for each year (every 12 data points)
            obsmax = np.nanmax(obs_dhw_ts.reshape(-1, 12), axis=1)
            modmax = np.nanmax(mod_dhw_ts.reshape(-1, 12), axis=1)
            
            # take the mean of the annual maxima
            obsmean = np.nanmean(obsmax)
            modmean = np.nanmean(modmax)
        
            # this is to make it easy to collate the seed data for BMA in pacage_dhw_for_bma.py
            dhwmeanmax_obs_da.loc[dict(index=loopind)] = obsmean
            dhwmeanmax_mod_da.loc[dict(index=loopind)] = modmean
        
            diff = np.abs(modmean - obsmean)
            
            # if ~np.isnan(modmean):
            #     pdb.set_trace()
            
            # some locations might have NaNs due to coast
            if np.isnan(diff):
                #print('diff nan')
                continue
            
            score = 1.0/diff
            if score > 10.: # this is arbitrary for now.
                score = 10.0
            print(str(modelind) + ' ' + str(loopind)+ ' ' +str(score))    
            pval_da.loc[dict(index=loopind)] = score
        
                   
    # save a netcdf file with XARRAY, one per model, lat lon pval    
    pval_ds = xr.Dataset({'pval': pval_da, 'dhwmeanmax_obs': dhwmeanmax_obs_da, 'dhwmeanmax_mod': dhwmeanmax_mod_da})
    pval_ds.to_netcdf(weightfilename) 
    print(weightfilename) 


    #pdb.set_trace()
