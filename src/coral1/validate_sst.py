#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This validation runs off of the model weights themselves, no need to run coarse project first.

I think for doing DHW validation, it is better to run coarse_project first.

Take weighted mean.
Take flat mean.
Look at difference between obs, means 2010-2014.

Created on Wed Sep  4 14:52:12 2019

@author: pkalmus
"""
import xarray as xr
import pdb
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import coral_plotter as plotter
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
start_year = params.start_year
nmonths = params.nmonths
R = params.R
do_plots = params.do_plots
J0 = params.J0 # to find the proper weight files
return_years = params.return_years
DHW_threshold = params.DHW_threshold
bias_correct = params.bias_correct
do_weights = params.do_weights
do_obs = params.do_obs
meanstart=params.meanstart
meanend=params.meanend
do_trend=params.do_trend

# define in and out dirs
weightdir = basedir+'data/tos/%s/weightbma/' % (runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
weightdir = basedir+'data/tos/%s/weightem/' % (runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
    
# filetag = '1965trend'
# filetag = 'weightbma_%i-%i' % (meanstart, meanend)

filetag = 'weightem_%i-%i' % (meanstart, meanend)
#filetag = 'bestbma_%i-%i' % (meanstart, meanend)
#filetag = 'meanbma_%i-%i' % (meanstart, meanend)

reefdir = basedir+'data/tos/%s/reef/' % (runname) 
outdir = weightdir

use_best_only = False

# read in locations and climatology, made using coral_locations.py
lon_lat_clim = np.loadtxt('reef_clim_1x1.txt')
index = np.arange(lon_lat_clim.shape[0])

with open(listfilename) as f:
    models = f.read().splitlines()


# read in obs
obsfilename = reefdir+'HadISST_sst_reef.nc'
print(obsfilename)
obs_ds = xr.open_dataset(obsfilename)
myobs = obs_ds.sel(time=slice(2010, 2014))

#https://github.com/pydata/xarray/issues/2713
def weighted_mean(data_da, dim, weights):
    r"""Computes the weighted mean.

    We can only do the actual weighted mean over the dimensions that
    ``data_da`` and ``weights`` share, so for dimensions in ``dim`` that aren't
    included in ``weights`` we must take the unweighted mean.

    This functions skips NaNs, i.e. Data points that are NaN have corresponding
    NaN weights.

    Args:
        data_da (xarray.DataArray):
            Data to compute a weighted mean for.
        dim (str | list[str]):
            dimension(s) of the dataarray to reduce over
        weights (xarray.DataArray):
            a 1-D dataarray the same length as the weighted dim, with dimension
            name equal to that of the weighted dim. Must be nonnegative.
    Returns:
        (xarray.DataArray):
            The mean over the given dimension. So it will contain all
            dimensions of the input that are not in ``dim``.
    Raises:
        (IndexError):
            If ``weights.dims`` is not a subset of ``dim``.
        (ValueError):
            If ``weights`` has values that are negative or infinite.
    """
    if isinstance(dim, str):
        dim = [dim]
    else:
        dim = list(dim)

    if not set(weights.dims) <= set(dim):
        dim_err_msg = (
            "`weights.dims` must be a subset of `dim`. {} are dimensions in "
            "`weights`, but not in `dim`."
        ).format(set(weights.dims) - set(dim))
        raise IndexError(dim_err_msg)
    else:
        pass  # `weights.dims` is a subset of `dim`

    if (weights < 0).any() or xr.ufuncs.isinf(weights).any():
        negative_weight_err_msg = "Weight must be nonnegative and finite"
        raise ValueError(negative_weight_err_msg)
    else:
        pass  # `weights` are nonnegative

    weight_dims = [
        weight_dim for weight_dim in dim if weight_dim in weights.dims
    ]

    if np.isnan(data_da).any():
        expanded_weights, _ = xr.broadcast(weights, data_da)
        weights_with_nans = expanded_weights.where(~np.isnan(data_da))
    else:
        weights_with_nans = weights

    mean_da = ((data_da * weights_with_nans).sum(weight_dims, skipna=True)
               / weights_with_nans.sum(weight_dims))
    other_dims = list(set(dim) - set(weight_dims))
    return mean_da.mean(other_dims, skipna=True)

# https://github.com/pydata/xarray/issues/422
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
    

# average all models' tos
# create new DS with a new model dim, and then average over that dim.
ds_list = []
weights_ds_list = []
for (modelind, model) in enumerate(models):
    print('\n\n***********'+str(modelind) + '   ' + model)
    
    # if modelind > 2:
    #     continue
    
    if model == 'CNRM-ESM2-1 r4i1p1f2 gn': # hada non-monotonic time!
        continue

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
    if len(mod_ds.lat) > 4000: # TODO: a couple of models in the list couldn't run through the NN process.
        weightfilename = weightdir+modelstr+'_%s.nc' % (filetag)
        weight_ds = xr.open_dataset(weightfilename) #tos      (time, lat, lon). lon: 0 to 360
        weights_ds_list.append(weight_ds)  # put this here for safety, to make sure it's never double-added.              
        ds_list.append(mod_ds)        
    else:
        print('model had less than 4000 pixels...')

all_ds = xr.concat(ds_list, 'model_index') # compat='override'? https://github.com/pydata/xarray/issues/2217 https://xarray.pydata.org/en/stable/generated/xarray.concat.html
flat_ds = all_ds.mean(dim='model_index', skipna=True, keep_attrs=True)
flat_tos_da = flat_ds.tos
flat_tos_da = flat_tos_da.sel(time=slice(2010,2014))
weights_ds = xr.concat(weights_ds_list, 'model_index')
weighted_tos_da = average_da(all_ds.tos, dim='model_index', weights=weights_ds.pval)
weighted_tos_da = weighted_tos_da.sel(time=slice(2010,2014))

flatdiffs = []
weightdiffs = []
for loopind in index:    
    looplon = lon_lat_clim[loopind,0]
    looplat = lon_lat_clim[loopind,1]
    obsmean = np.mean(myobs.tos.values[loopind,:])
    flatmean = np.mean(flat_tos_da.values[loopind,:])
    weightmean = np.mean(weighted_tos_da.values[loopind,:])
    
    flatdiff = np.abs(flatmean-obsmean)
    weightdiff = np.abs(weightmean-obsmean)
    
    # flatdiff = flatmean-obsmean
    # weightdiff = weightmean-obsmean
    
    flatdiffs.append(flatdiff)
    weightdiffs.append(weightdiff)

plotter.hist(np.array(flatdiffs), 'flat', outdir+'hist_flatdiffs_abs_%s.png' % (filetag))
plotter.hist(np.array(weightdiffs), 'weight', outdir+'hist_weightdiffs_abs_%s.png' % (filetag))


