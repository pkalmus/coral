#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

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

# limit to times
obs_ds =  obs_ds.sel(time=slice(meanstart, meanend))
obs_means = np.nanmean(obs_ds.tos.values, axis=1)

with open(listfilename) as f:
    models = f.read().splitlines()

mod_matrix = []
for (model_ind, model) in enumerate(models):
    # if model_ind > 1:
    #     continue
    modelstr = model.replace(' ', '_')
    filename = indir+modelstr+'_ssp126_reef.nc'
    print(filename)
    mod_ds = xr.open_dataset(filename, decode_times=False)

    # limit to times
    mod_ds =  mod_ds.sel(time=slice(meanstart, meanend))
    mod_means = np.nanmean(mod_ds.tos.values, axis=1)
    # if np.nanmin(mod_means) == 0:
    #     pdb.set_trace()
    mod_matrix.append(mod_means.tolist())

mdic = {'mod_matrix': mod_matrix, 'obs_vect': obs_means}
savemat('obs_mod_%i-%i.mat' % (meanstart, meanend), mdic)




