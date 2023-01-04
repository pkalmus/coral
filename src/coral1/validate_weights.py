#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

conda activate geo_env

Make hist and map of weights, for each model.

todo: add param file and new model list.

@author: pkalmus
"""

#projectdir = '/0work/projects/'
projectdir = '/home/pkalmus/projects/'

import sys
import xarray as xr
import numpy as np
import pdb
import coral_plotter as plotter

projectdir = '/home/pkalmus/projects/'
basedir = projectdir+'/coral/'
datadir = basedir+'data/tos/' # data is a soft link to /raid8/pkalmus/data/coral/data
outdir = basedir+'output/initial/'
outdir = basedir+'output/lowfreq/'
runname = 'test'

with open('model_list_for_weights.txt') as f:
    models = f.read().splitlines()

for model in models:
    print(model)
    try:
        tosfilename = datadir+'regrid/CMIP6_1x1/'+model+'_reef.nc'
        tos_ds = xr.open_dataset(tosfilename)
        weightfilename = datadir+'regrid/CMIP6_1x1/'+model+'_weights_biascorrect.nc'
        weights_ds = xr.open_dataset(weightfilename)
        
        plotter.hist(weights_ds.pval.values, 'pval', outdir+'hist_pvals_%s_%s.png' % (model, runname), logy=True, logx= False, ylim=[1,1000], bins=20)
        plotter.scatter(tos_ds.lon.values, tos_ds.lat.values, outdir+'map_pvals_%s_%s.png' % (model, runname), c=weights_ds.pval.values, vmin=0.0, vmax=1.0)
    except:
        continue


