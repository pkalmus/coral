#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
python validate_coarse_sst.py params

conda activate geo_env

Not technically validation. Instead a bunch of sanity checks.

Make 12-month running means.

x make sst plot for the mean
x make sst plots for individual models
look at interannual variability from each model vs. the mean
    to measure interannual variability: 
        seasonal_decompose doesn't work b/c it folds interannual variability into "trend"
        fourier space, or bandpass.
look at the range across models at each pixel (make distribution).


@author: pkalmus
"""
import xarray as xr
import pdb
import numpy as np
import glob
import coral_plotter
import sys
import subprocess

# read in user params
import importlib
paramfile = sys.argv[1]   
params = importlib.import_module(paramfile)
projectdir = params.projectdir
runname = params.runname
tosrunname = params.tosrunname
basedir = params.basedir
pythondir = params.pythondir
listfilename = params.listfilename
scenarios = params.scenarios
do_plots = params.do_plots
bias_correct = params.bias_correct
do_obs = params.do_obs
finerun = params.finerun
includemur = params.include_mur

lat_extent = 35
include_coarse = False

# define in and out dirs
reefdir = basedir+'data/tos/%s/reef/' % (tosrunname) # the reef (individual coarse) files
indir = basedir+'data/tos/%s/coarse/' % (runname)    # the mean files
outdir = basedir+'output/%s/coarse/' % (runname)

sysStr = "mkdir -p %s" % outdir
subprocess.call(sysStr, shell=True)



ds_list = []
coarse_list = []
for scenario in scenarios:
    experiment = '%s' % (scenario) # this will be replaced by individual model if we're in single_mod
    print('++++++starting SST ' + experiment)

    coarsefile = indir+'mean_tos_%s_%s.nc' % (runname, scenario)
    print(coarsefile)
    coarse_ds = xr.open_dataset(coarsefile)
    coarse_ds = coarse_ds.where(coarse_ds.isreef==1, drop=True)
    coarse_ds = coarse_ds.where(coarse_ds.time>=2010, drop=True)
    coarse_list.append(coarse_ds)
    #pdb.set_trace()
coral_plotter.sstlist(None, coarse_list, None, 'Year', outdir+'sst_%s_%s.png' % (runname, experiment), fontsize=17, ylabel='SST ($^\circ$C)', xlim=[2010,2099],axvspan=[2010,2020],label=None)


# look at individual models for SSP126
scenario = 'ssp126'
with open(listfilename) as f:
    models = f.read().splitlines()
    
dslist = []
for (modelind, model) in enumerate(models):
    print('\n\n***********'+str(modelind+1) + '   ' + model)

    # e.g. CanESM5 r11i1p2f1 gn
    # UKESM1-0-LL_r8i1p1f2_gn_ssp126_reef.nc
    modelstr = model.replace(' ', '_')
    modsspstr = modelstr+'_'+scenario
    print('model, ssp: ' + modsspstr)
    filename = reefdir+modsspstr+'_reef.nc'
    
    try:
        mod_ds = xr.open_dataset(filename) #tos      (time, lat, lon). lon: 0 to 360  
        dslist.append(mod_ds)
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        print('file not found: %s' % filename)  

meands = xr.open_dataset(indir+'mean_tos_%s_%s.nc' % (runname, scenario))
coral_plotter.sstlistmodels(dslist, meands, 'Year', outdir+'sstmodels_%s_%s.png' % (runname, experiment), fontsize=17, ylabel='SST ($^\circ$C)', xlim=[2010,2099],axvspan=[2010,2020],label=None)
