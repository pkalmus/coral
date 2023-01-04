#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

tas_from_bcdp.py
write_gmst_from_tas.py.

Read in files.
Take the two mean over ensemble.

Save netcdf files with the means.

to check:
    1 they are all in K
    2 can I use bcdp to do time homogenization?
    
Once the means are made, need to convert to temp anomaly -


Created on Mon Mar 22 18:15:43 2021

@author: pkalmus
"""


import xarray as xr
import numpy as np
import pdb
import sys
import traceback
import subprocess
import os.path
import logging
import glob
logging.basicConfig(filename='run.log',level=logging.WARNING)

#import cftime

# read in user params
import importlib
paramfile = sys.argv[1]   
params = importlib.import_module(paramfile)
runname = params.runname
projectdir = '/home/pkalmus/projects/'
basedir = projectdir+'/coral/'
listfilename = params.listfilename
scenarios = params.scenarios
dryRun = params.dryRun
oneMember = True # only get the first member from a source (group)
grids = params.grids


outdir = basedir+'data/tas/%s/' % (runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
listfilename = outdir+'models_gmst.txt'

# simplified float times
time = []
for year in np.arange(1850, 2101):
    for month in np.arange(1,13):
        time.append(year+month/100.)
time = np.floor(time) + ((time - np.floor(time))*30*100 - 15)/365. 

for scenario in scenarios:
    #/home/pkalmus/projects//coral/data/tas/nature/INM-CM5-0_r1i1p1f1_gr1_tas_ssp585.nc
    infiles = glob.glob(outdir+'*tas_%s.nc' % (scenario))
 
    
    
    ds_list = [] 
    for file in infiles:
        mod_ds = xr.open_dataset(file)
        #pdb.set_trace()
        #print(mod_ds.tas.values[0])
        #print(mod_ds.time.values[0])
        #print(len(mod_ds.time.values))
        
        if len(mod_ds.time.values) == 3012: # there are 4 oddballs
            mod_ds = mod_ds.assign_coords(time=np.copy(time))  
            ds_list.append(mod_ds)

     
    all_ds = xr.concat(ds_list, 'model_index', coords='minimal', compat='override') # don't care about 'lev.' tos: (127, 4836, 2220).  compat='override'? https://github.com/pydata/xarray/issues/2217 https://xarray.pydata.org/en/stable/generated/xarray.concat.html  
    mean_ds = all_ds.mean(dim='model_index', skipna=True, keep_attrs=True)
    mean_ds['tas60'] = mean_ds.tas.rolling(time=60).mean() # 60 month rolling mean
    # calculate anomaly relative to 1880-1900 
    baseline = mean_ds.tas.sel(time=slice(1880, 1900)).mean(dim='time').values
    mean_ds['anom'] = mean_ds.tas60 - baseline
    outfile = outdir + 'mean_%s.nc' % (scenario)
    print(outfile)
    mean_ds.to_netcdf(outfile)