#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

conda activate geo_env

2022/03/01:
Ayesha found an issue with the HadISST reef file not going past 180 longitude. The other files had it too.
The problem was in reef_locations_coarse.py.

@author: pkalmus
"""

import xarray as xr
import pandas as pd
import numpy as np
import pdb
import sys
import subprocess
import coral_plotter as plotter

# read in user params
import importlib
paramfile = sys.argv[1]   
params = importlib.import_module(paramfile)
projectdir = params.projectdir
tosrunname = params.tosrunname
basedir = params.basedir
srcdir = params.srcdir
figdir = params.figdir
modellist_filename = params.modellist_filename
scenarios = ['ssp126']
dryRun = params.dryRun
oneMember = params.oneMember


# define in and out dirs
indir = basedir+'data/tos/coral2/%s/raw/' % (tosrunname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
outdir = basedir+'data/tos/coral2/%s/reef/' % (tosrunname)

reefmodellist_filename = modellist_filename.replace('.txt', '_reef.txt')
df = pd.read_csv(reefmodellist_filename, delimiter=' ', header=None)
models = df.iloc[:,0].values



for (modelind, modelstr) in enumerate(models):
    print('\n\n***********'+str(modelind+1) + '   ' + modelstr)
      
    for scenario in scenarios:
        if modelind == 0:
            ssp_modelstr = modelstr            
        else:
            # e.g. CanESM5_r12i1p1f1_ssp126.nc
            ssp_modelstr = modelstr+'_'+scenario
        outfilename = outdir+ssp_modelstr+'_reef.nc'
        print(outfilename)
        ds = xr.open_dataset(outfilename)
        plotter.scatter(ds.lon.values, ds.lat.values, figdir+'check_regrid_reef_%s.png' % (ssp_modelstr), c=ds.tos.values[:,0], marker_relsize=0.1, figsize=(15,5), cbar_fraction=0.008, cbar_orientation='vertical')