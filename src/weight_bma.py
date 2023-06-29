#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Read in Elias' BMA .mat file (which he makes from output of package_dhw_for_bma.py)
 and create weight .nc files.

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
import scipy.io as sio

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
# /home/pkalmus/projects/coral/data/tos/nature/weightbma/
indir = basedir+'data/tos/%s/reef/' % (runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
outdir = basedir+'data/tos/%s/weightbma/' % (runname)

sysStr = "mkdir -p %s" % outdir
subprocess.call(sysStr, shell=True)


# matfile = outdir + 'mean_posterior_BMA_allpixels.mat'
# filetag = 'weightbma_%i-%i' % (meanstart, meanend)
# mat_contents = sio.loadmat(matfile)
# weight_mat = mat_contents['mean_posterior_BMA_allpixels']

# matfile = outdir + 'BMA_best_weights.mat'
# filetag = 'bestbma_%i-%i' % (meanstart, meanend)
# mat_contents = sio.loadmat(matfile)
# weight_mat = mat_contents['BMA_best_weights']

# matfile = outdir + 'BMA_mean_weights.mat'
# filetag = 'meanbma_%i-%i' % (meanstart, meanend)
# mat_contents = sio.loadmat(matfile)
# weight_mat = mat_contents['BMA_mean_weights']

matfile = outdir + 'best_BMA_combo_DHW.mat'
filetag = 'dhwbma_%i-%i_%i-%i' % (climstartyear, climendyear, meanstart, meanend)
mat_contents = sio.loadmat(matfile)
weight_mat = mat_contents['best_BMA_combo_DHW']

# matfile = outdir + 'best_BMA_combo_DHW_LongerPast.mat'
# filetag = 'dhwbma_%i-%i_%i-%i' % (climstartyear, climendyear, meanstart, meanend)
# mat_contents = sio.loadmat(matfile)
# weight_mat = mat_contents['best_BMA_combo_DHW_LongerPast']

# matfile = outdir + 'best_BMA_combo_DHW_Production.mat'
# filetag = 'dhwbma_%i-%i_%i-%i' % (climstartyear, climendyear, meanstart, meanend)
# mat_contents = sio.loadmat(matfile)
# weight_mat = mat_contents['best_BMA_combo_DHW_Production']
 

# Obs: matrix, num, [1:1024, 1:2] t in fractional years, mysst in K.
# read in the Obs matrix
# send to R
# check and see if it's the same
obs_filename = indir+'HadISST_sst_reef.nc'
obs_ds = xr.open_dataset(obs_filename, decode_times=False)

with open(listfilename) as f:
    models = f.read().splitlines()

    
for (modelind, model) in enumerate(models):
    model = models[modelind]
    modelstr = model.replace(' ', '_')
    filename = indir+modelstr+'_ssp126_reef.nc'

    weightfilename = outdir+modelstr+'_%s.nc' % (filetag)
    print(weightfilename+str(modelind))

    null_data = np.empty((obs_ds.dims['index']))
    null_data[:] = np.nan
    pval_da = xr.DataArray(null_data, coords=[obs_ds.coords['index']], dims=['index'])  


    for loopind in obs_ds.index.values: 
        # if loopind > 10:
        #     break
        
        score = weight_mat[loopind, modelind]
        pval_da.loc[dict(index=loopind)] = score
               
                
    # save a netcdf file with XARRAY, one per model, lat lon pval    
    pval_ds = xr.Dataset({'pval': pval_da})
    pval_ds.to_netcdf(weightfilename) 
    print(weightfilename) 



