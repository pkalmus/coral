#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Just take last 20 years, average, see how it matches obs. 
python weight_endmatch.py 4 # provide integer corresponding to model in the model list, or call start_weights_NO_SSHX.sh


Could try also making an ensemble where the worst matches have highest weights, 
  to  make sure there's at least some 'skill

Test on ssp126.
Check to see if there's a discontinuity between historical and ssp.

Historical is 1850-2014. For this test, let's take the mean of 2000-2010,
and then validate off of mean for 2011-2014.

Created on Tue Oct  8 14:43:52 2019

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

weightfilename = outdir+modelstr+'_weightem_%i-%i.nc' % (meanstart, meanend)
mod_ds = xr.open_dataset(filename, decode_times=False)


# limit to 2000-2009
obs_ds =  obs_ds.sel(time=slice(meanstart, meanend))
mod_ds =  mod_ds.sel(time=slice(meanstart, meanend))

for loopind in obs_ds.index.values: 
    # if loopind > 10:
    #     break

    myobs_ds = obs_ds.sel(index=loopind)
    mysst = myobs_ds.tos.values
    obsmean = np.mean(mysst)
    
    mymod_ds = mod_ds.sel(index=loopind)
    mysst = mymod_ds.tos.values
    modmean = np.mean(mysst)

    diff = np.abs(modmean - obsmean)
    
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


