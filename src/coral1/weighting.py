#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

conda activate geo_env, SSH to server without X, to avoid a million plots from the R code.
python weighting.py 4 # provide integer corresponding to model in the model list, or call start_weights_NO_SSHX.sh


regrid_bcdp
coral_locations
regrid_reef
weighting
coarse_project
validating code

run with start_weights.sh or python weighting.py i, where i is the ith model in model_list_for_weights.txt
run without -X in ssh -X, or else many plots pop up from R.

Read in reef HadISST and the reef model netcdf (created by regrid_reef.py)
Call R from python.
Calculate the weights at each location. Lowest value is 1/R. Highest value is (R-1)/R, R as in EstimateWeight.R
Store in netcdf.



******************************************************** 1942 0.002
Traceback (most recent call last):
  File "weighting.py", line 137, in <module>
    clim = lon_lat_clim[loopind,3] # 1985-2012 inclusive
IndexError: index 1946 is out of bounds for axis 0 with size 1946



todo: add to existing netcdf (according to index), not a separate file

todo: make sure coords are sorted correctly (ascending... latitude..., lon 0 to 360)
todo: maybe just add the weights to existing STANDARDIZED mod files? 

todo: validate the weights 

todo: do the grid points match?
todo: do the time points match?

todo: cmip6?

Created on Tue Oct  8 14:43:52 2019

@author: pkalmus
"""

#projectdir = '/0work/projects/'
projectdir = '/home/pkalmus/projects/'

import sys
import xarray as xr
import numpy as np
import pdb
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import coral_plotter as plotter
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
J0 = params.J0

# define in and out dirs
indir = basedir+'data/tos/%s/reef/' % (runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
outdir = basedir+'data/tos/%s/weights/' % (runname)

sysStr = "mkdir -p %s" % outdir
subprocess.call(sysStr, shell=True)

# these two lines needed to pass np array into R
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

utils = importr('utils')
base = importr('base')
print(base._libPaths())

#utils.install_packages('wavethresh, WiSEBoot, normwhn.test, tseries, forecast', repos='http://r-forge.r-project.org') # only need to load once
# SEE here for alternative: https://stackoverflow.com/questions/15419740/calling-custom-functions-from-python-using-rpy2
robjects.r('''
       source('../r/EstimateWeight.R')
''')
r_func = robjects.globalenv['EstimateWeight']

wrote_high_sample = False # get sample time series to examine weighting methodology
wrote_low_sample = False

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
    filename = indir+modelstr+'_ssp585_reef.nc'
    print(filename)

null_data = np.empty((obs_ds.dims['index']))
null_data[:] = np.nan
pval_da = xr.DataArray(null_data, coords=[obs_ds.coords['index']], dims=['index'])    

weightfilename = outdir+modelstr+'_weightsJ0%i_%i.nc' % (J0, nmonths)
mod_ds = xr.open_dataset(filename, decode_times=False)

for loopind in obs_ds.index.values: 
    # if loopind > 10:
    #     break
    
    mymod_ds = mod_ds.sel(index=loopind)
    mysst = mymod_ds.tos.values

    t = mymod_ds.year.values
    myind = np.where(t >= start_year)[0]
    myind = myind[0:nmonths]
    mod = np.array([t[myind], mysst[myind]])
    mod = mod.T
    
    # some locations might have NaNs due to coast
    if np.isnan(mod[0,1]):
        continue

    myobs_ds = obs_ds.sel(index=loopind)
    mysst = myobs_ds.tos.values
    t = myobs_ds.year.values
    myind = np.where(t >= start_year)[0]
    myind = myind[0:nmonths]
    obs = np.array([t[myind], mysst[myind]])
    obs = obs.T 
    
    # some locations might have NaNs due to coast - even where mod doesn't.
    if np.isnan(obs[0,1]):
        continue
    if np.max(mod[:,1])==0:
        continue
    
    # call the R function and get weight
    # note: if you put in "mod, mod" it gives pvals very near 1.
    pval = r_func(obs, mod)[0]
    if pval==0:
        pval = 1/R
    if pval==1:
        pval = (R-1)/R
    print('********************************************************', loopind, pval)
    pval_da.loc[dict(index=loopind)] = pval

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
 
goodvals = pval_da.values[~np.isnan(pval_da.values)]       
print(goodvals)
print(np.mean(goodvals))        
            
# save a netcdf file with XARRAY, one per model, lat lon pval    
pval_ds = xr.Dataset({'pval': pval_da})
pval_ds.to_netcdf(weightfilename) 
print(weightfilename) 

# this was for appending to the model file. but you would need to duplicate weights for each ssp.
#mod_ds.close() # before appending, you must close
#pval_ds.to_netcdf(filename, mode='a')


