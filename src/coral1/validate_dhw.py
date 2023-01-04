#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This validation runs off of the model weights themselves, no need to run coarse project first.

Note: replaced zeros in individual model timeseries with nan. 
These zeros need to be taken out at the source (BCDP stage)

Any index with fewer than 10 models is NaN, to prevent e.g. a single bad model from
creating a bad overall value that sticks out. There are more of these than I thought.

Take weighted mean.
Take flat mean.
Look at difference between obs, means 2010-2014.

Created on Wed Sep  4 14:52:12 2019

@author: pkalmus
"""
import xarray as xr
import pdb
import numpy as np
import coral_plotter as plotter
import sys

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
meanstart=params.meanstart
meanend=params.meanend
do_trend=params.do_trend
min_models=params.min_models

valstartyear = params.valstartyear
valendyear = params.valendyear
climstartyear = params.climstartyear
climendyear = params.climendyear

# define in and out dirs
#weightdir = basedir+'data/tos/%s/weightem/' % (runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
    
# filetag = '1965trend'
# filetag = 'weightbma_%i-%i' % (meanstart, meanend)

#filetag = 'weightem_%i-%i' % (meanstart, meanend)
#filetag = 'bestbma_%i-%i' % (meanstart, meanend)
#filetag = 'meanbma_%i-%i' % (meanstart, meanend)

# weightdir = basedir+'data/tos/%s/weightdhw/' % (runname) 
# filetag = 'weightdhw_%i-%i_%i-%i' % (climstartyear, climendyear, meanstart, meanend)

weightdir = basedir+'data/tos/%s/weightbma/' % (runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
filetag = 'dhwbma_%i-%i_%i-%i' % (climstartyear, climendyear, meanstart, meanend)

reefdir = basedir+'data/tos/%s/reef/' % (runname) 
outdir = weightdir

# read in locations and climatology, made using coral_locations.py
# lon_lat_clim = np.loadtxt('reef_clim_1x1.txt')
# index = np.arange(lon_lat_clim.shape[0])

with open(listfilename) as f:
    models = f.read().splitlines()

# read in obs
obsfilename = reefdir+'HadISST_sst_reef.nc'
print(obsfilename)
obs_ds = xr.open_dataset(obsfilename)
# calculate obs mmm climatology from HadISST
clim_ds = obs_ds.sel(time=slice(climstartyear, climendyear))
# take annual mean
clim_tos = clim_ds.tos.values
clim_mo = np.nanmean(clim_tos.reshape(clim_tos.shape[0],-1,12), 1)
# find mmm for each index
clim_mmm = np.max(clim_mo, 1)


myobs = obs_ds.sel(time=slice(valstartyear, valendyear))


# https://github.com/pydata/xarray/issues/422. note, alternative (deleted, did not use): #https://github.com/pydata/xarray/issues/2713
def average_da(myda, dim=None, weights=None):
    """
    weighted average for DataArrays

    Parameters
    ----------
    dim : str or sequence of str, optional
        Dimension(s) over which to apply average.
    weights : DataArray
        weights to apply. Shape must be broadcastable to shape of self.

    Returns
    -------
    reduced : DataArray
        New DataArray with average applied to its data and the indicated
        dimension(s) removed.

    """
    if weights is None:
        return myda.mean(dim)
    else:
        if not isinstance(weights, xr.DataArray):
            raise ValueError("weights must be a DataArray")

        # if NaNs are present, we need individual weights
        if myda.notnull().any():
            total_weights = weights.where(myda.notnull()).sum(dim=dim, skipna=True)
        else:
            total_weights = weights.sum(dim, skipna=True)

        return (myda * weights).sum(dim, skipna=True) / total_weights

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
        

# average all models' tos
# create new DS with a new model dim, and then average over that dim.
ds_list = []
weights_ds_list = []
for (modelind, model) in enumerate(models):
    print('\n\n***********'+str(modelind) + '   ' + model)
    
    # if modelind > 11: # not less than  min_models!!
    #     continue
    
    # if model == 'CNRM-ESM2-1 r4i1p1f2 gn': # hada non-monotonic time!
    #     continue

    modelstr = model.replace(' ', '_')
    modsspstr = modelstr+'_'+'ssp126'
    print('model, ssp: ' + modsspstr)
    filename = reefdir+modsspstr+'_reef.nc'
    
    try:
        mod_ds = xr.open_dataset(filename) #tos      (time, lat, lon). lon: 0 to 360
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        print('file not found: %s' % filename)   
    if len(mod_ds.lat) > 4000: #a couple of models in the list couldn't run through the NN process.
        
        # replace any zeros, at least until I fix this at the source in the BCDP code.
        #print(np.where(mod_ds.tos.values==0)[0])
        #pdb.set_trace()
        #mod_ds = mod_ds.where(mod_ds['tos'] != 0.) # btw this takes about half a second to do - and might add more to concat / average
        mod_ds['tos'] = mod_ds.tos.where(mod_ds.tos != 0)
        
        weightfilename = weightdir+modelstr+'_%s.nc' % (filetag)
        weight_ds = xr.open_dataset(weightfilename) #tos      (time, lat, lon). lon: 0 to 360
        weights_ds_list.append(weight_ds)  # put this here for safety, to make sure it's never double-added.              
        ds_list.append(mod_ds)        
    else:
        print('model had less than 4000 pixels...')

all_ds = xr.concat(ds_list, 'model_index') # tos: (127, 4836, 2220).  compat='override'? https://github.com/pydata/xarray/issues/2217 https://xarray.pydata.org/en/stable/generated/xarray.concat.html
for loopind in obs_ds.index.values: 
    num_model_estimates = np.count_nonzero(~np.isnan(all_ds.tos.values[:,loopind, 1]))
    if num_model_estimates < min_models:
        print('too few model estimates %i for index %i' % (num_model_estimates, loopind))
        all_ds.tos.values[:,loopind,:] = np.nan
flat_ds = all_ds.mean(dim='model_index', skipna=True, keep_attrs=True)
flat_tos_da = flat_ds.tos
flat_tos_da = flat_tos_da.sel(time=slice(valstartyear, valendyear))

weights_ds = xr.concat(weights_ds_list, 'model_index')
# any index with fewer than 10 models w/o nan, have to nan out. e.g. 1798 has only one. 
for loopind in obs_ds.index.values: 
    num_model_estimates = np.count_nonzero(~np.isnan(weights_ds.pval.values[:,loopind]))
    if num_model_estimates < min_models:
        print('too few model estimates %i for index %i' % (num_model_estimates, loopind))
        weights_ds.pval.values[:,loopind] = np.nan
weighted_tos_da = average_da(all_ds.tos, dim='model_index', weights=weights_ds.pval)
weighted_tos_da = weighted_tos_da.sel(time=slice(valstartyear, valendyear))

flatdiffs = []
weightdiffs = []
for loopind in obs_ds.index.values: # validation loop
    print(loopind)
    mmm = clim_mmm[loopind]
    if ~np.isnan(mmm):
        obs = myobs.tos.values[loopind,:]
        # check if it has tos (is ocean)
        if np.isnan(obs[0]):
            print('tos nan')
            continue
        

        obs_dhw_ts = get_dhw(mmm, obs)
        flat_dhw_ts = get_dhw(mmm, flat_tos_da.values[loopind,:])
        weight_dhw_ts = get_dhw(mmm, weighted_tos_da.values[loopind,:])

       
        # select the maximum for each year (every 12 data points)
        obsmax = np.nanmax(obs_dhw_ts.reshape(-1, 12), axis=1)
        flatmax = np.nanmax(flat_dhw_ts.reshape(-1, 12), axis=1)
        weightmax = np.nanmax(weight_dhw_ts.reshape(-1, 12), axis=1)
        
        # take the mean of the annual maxima
        obsmean = np.nanmean(obsmax)
        flatmean = np.nanmean(flatmax)
        weightmean = np.nanmean(weightmax)
        
        # there a few weightmean values that are very big.
        # e.g. index 1798
        # if weightmean > 50: 
        #     pdb.set_trace()

        flatdiff = flatmean-obsmean
        weightdiff = weightmean-obsmean

       
        # flatdiff = np.nanmean(flat_dhw_ts)-np.nanmean(obs_dhw_ts)
        # weightdiff = np.nanmean(weight_dhw_ts)-np.nanmean(obs_dhw_ts)
        
        flatdiffs.append(flatdiff)
        weightdiffs.append(weightdiff)

data = np.array(flatdiffs)
plotter.hist(data, 'degree heating weeks', outdir+'hist_flatdiffs_dhwmax_%s_%i.png' % (filetag, min_models), xlim=[-15,40], logy=True, fontsize=16)
data = np.array(weightdiffs)
plotter.hist(data, 'degree heating weeks', outdir+'hist_weightdiffs_dhwmax_%s_%i.png' % (filetag, min_models), xlim=[-15,40], logy=True, fontsize=16)


