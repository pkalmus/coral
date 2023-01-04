#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""



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
from scipy.ndimage.filters import uniform_filter1d

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

with open(listfilename) as f:
    models = f.read().splitlines()

if len(sys.argv) != 3:
    print('You must supply the integer index for a model in models_<runname>.txt as an argument, e.g. "weighting.py params 3".')
    sys.exit(0)
else:
    # e.g. UKESM1-0-LL_r8i1p1f2_gn_ssp370_reef.nc
    modelind = int(sys.argv[2])
    model = models[modelind]
    modelstr = model.replace(' ', '_')
    filename = indir+modelstr+'_ssp126_reef.nc'
    print(filename)

null_data = np.empty((obs_ds.dims['index']))
null_data[:] = np.nan
pval_da = xr.DataArray(null_data, coords=[obs_ds.coords['index']], dims=['index'])    

weightfilename = outdir+modelstr+'_weightem_1965trend.nc'
mod_ds = xr.open_dataset(filename, decode_times=False)


# based on 0-2000m global ocean heat content, a clear trend starts at 1990. We'll take 5-year running mean
obs_ds =  obs_ds.sel(time=slice(1960, 2015))
mod_ds =  mod_ds.sel(time=slice(1960, 2015))    
N = 60 # 5-year (60 month) running average

for loopind in obs_ds.index.values: 
    # if loopind > 10:
    #     break

    myobs_ds = obs_ds.sel(index=loopind)
    sstobs = myobs_ds.tos.values
    avobs = uniform_filter1d(sstobs, size=N)     # running average
    trendobs = np.mean(avobs[-120:-59]) - np.mean(avobs[60:121])
    
    mymod_ds = mod_ds.sel(index=loopind)
    sstmod = mymod_ds.tos.values
    avmod = uniform_filter1d(sstmod, size=N)     # running average
    trendmod = np.mean(avmod[-120:-59]) - np.mean(avmod[60:121])

    diff = np.abs(trendmod - trendobs)
    
    # some locations might have NaNs due to coast
    if np.isnan(diff):
        continue
    
    score = 1.0/diff
    print(score)    
    pval_da.loc[dict(index=loopind)] = score


    # if modelind==3 and (pval==0.002) and not wrote_low_sample:
    # #if modelind==3 and loopind==2:
    #     plotter.timeseries(obs, weightdir+model+'_'+str(loopind)+'_%1.3f_1024.png' % (pval), t2=mod, title='pval: %1.3f' % (pval))
    #     np.savetxt(weightdir+model+'_obs'+str(loopind)+'_%1.3f_1024.txt' % (pval), obs, fmt='%1.2f')
    #     np.savetxt(weightdir+model+'_mod'+str(loopind)+'_%1.3f_1024.txt' % (pval), mod, fmt='%1.2f')
    #     wrote_low_sample = True
    #     print(loopind)
    #     pdb.set_trace()
        
#     if modelind==5 and (pval > 0.5) and not wrote_high_sample:
#         plotter.timeseries(obs, weightdir+model+'_'+str(loopind)+'_%1.3f_1024.png' % (pval), t2=mod, title='pval: %1.3f' % (pval))
#         np.savetxt(weightdir+model+'_obs'+str(loopind)+'_%1.3f_1024.txt' % (pval), obs, fmt='%1.2f')
#         np.savetxt(weightdir+model+'_mod'+str(loopind)+'_%1.3f_1024.txt' % (pval), mod, fmt='%1.2f')
#         wrote_high_sample = True
#         # print(loopind)
#         # pdb.set_trace()
       
            
# save a netcdf file with XARRAY, one per model, lat lon pval    
pval_ds = xr.Dataset({'pval': pval_da})
pval_ds.to_netcdf(weightfilename) 
print(weightfilename) 

# this was for appending to the model file. but you would need to duplicate weights for each ssp.
#mod_ds.close() # before appending, you must close
#pval_ds.to_netcdf(filename, mode='a')


