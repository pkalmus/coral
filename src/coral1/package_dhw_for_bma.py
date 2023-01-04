#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Cycle through each model and the obs, save a matrix for the BMA
Depends on running weight_dhw.py

After Elias provides BMA using hte .mat file created here, run weight_bma.py to create the final weights.
Then run coarse_project.py to create the files for Ayesha.

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
meanstart=params.meanstart
meanend=params.meanend
climstartyear = params.climstartyear
climendyear = params.climendyear
min_models=params.min_models

# define in and out dirs
indir = basedir+'data/tos/%s/reef/' % (runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
weightdir = basedir+'data/tos/%s/weightdhw/' % (runname) 
filetag = 'weightdhw_%i-%i_%i-%i' % (climstartyear, climendyear, meanstart, meanend)

# HadISST_sst_reef.nc  tos_Omon_CMIP6_CanESM5_historical_ssp585_r1i1p2f1_gn_185001-210012_reef.nc

# Obs: matrix, num, [1:1024, 1:2] t in fractional years, mysst in K.
# read in the Obs matrix
# send to R
# check and see if it's the same
# just getting this for the index
obs_filename = indir+'HadISST_sst_reef.nc'
obs_ds = xr.open_dataset(obs_filename, decode_times=False)

with open(listfilename) as f:
    models = f.read().splitlines()

# load in the weight files
weights_ds_list = []
for (modelind, model) in enumerate(models):
    print(str(modelind) + '   ' + model)
    
    # if modelind > 2:
    #     continue

    modelstr = model.replace(' ', '_')
    weightfilename = weightdir+modelstr+'_%s.nc' % (filetag)       
    weight_ds = xr.open_dataset(weightfilename) 
    weights_ds_list.append(weight_ds)            

weights_ds = xr.concat(weights_ds_list, 'model_index')
obs_means = weights_ds.dhwmeanmax_obs.values[4,:] # any loop index will do
mod_matrix = weights_ds.dhwmeanmax_mod.values

# any index with fewer than 10 models w/o nan, have to nan out. e.g. 1798 has only one. 
for loopind in obs_ds.index.values: 
    num_model_estimates = np.count_nonzero(~np.isnan(mod_matrix[:,loopind]))
    if num_model_estimates <= min_models:
        print('too few model estimates %i for index %i' % (num_model_estimates, loopind))
        mod_matrix[:,loopind] = np.nan

# create the .mat package for the BMA. I think it's best to redo it.

# obs_means is a vect of mean max dhw values from 1980-2000, one per index.
# mod_matrix has mean max dhw values from 1980-2000, one per index, per model. 
# any index with fewer than 10 models gets nans for that index. 
mdic = {'mod_matrix': mod_matrix, 'obs_vect': obs_means}
savemat('obs_mod_dhwmax_%i-%i_%i-%i.mat' % (climstartyear, climendyear, meanstart, meanend), mdic)


