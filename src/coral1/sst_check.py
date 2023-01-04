#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
python sst_plot.py params

conda activate geo_env



@author: pkalmus
"""
import xarray as xr
import pdb
import numpy as np
import glob
import coral_plotter as plotter
import sys
import subprocess

import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import cartopy.feature
import cartopy.crs as ccrs

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
do_plots = params.do_plots
return_years = params.return_years
DHW_threshold = params.DHW_threshold
bias_correct = params.bias_correct
do_weights = params.do_weights
do_obs = params.do_obs
finerun = params.finerun
includemur = params.include_mur

lat_extent = 35


# define in and out dirs
indir = basedir+'data/tos/%s/weights/' % (runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
reefdir = basedir+'data/tos/%s/reef/' % (runname) 
outdir = basedir+'data/tos/%s/coarse/' % (runname)
#finedir = basedir+'data/tos/%s/fine/%s/' % (runname, finerun) # e.g. for scp: /raid8/pkalmus/data/coral/data/tos/nature/fine/flat/

finedir = basedir+'output/%s/fine/%s/' % (runname, finerun) # e.g. for scp: /raid8/pkalmus/data/coral/data/tos/nature/fine/flat/
murdir = basedir + 'data/mur/MUR_coral-pixels_2002-2019/'

sysStr = "mkdir -p %s" % outdir
subprocess.call(sysStr, shell=True)

use_best_only = False

fig, ax = plt.subplots()


N=12
startyear=1900
endyear=2019

# read in obs
obsfilename = reefdir+'HadISST_sst_reef.nc'
print(obsfilename)
obsds = xr.open_dataset(obsfilename)
obsds = obsds.where(obsds.isreef==1, drop=True)
obsds = obsds.where(obsds.time>=startyear, drop=True)
obsds = obsds.where(obsds.time<=endyear, drop=True)   
temp_obs = obsds.tos.values
time_obs = obsds.time.values
temp_obs = np.nanmean(temp_obs, axis=0)
temp_obs = np.convolve(temp_obs, np.ones(N)/N, mode='valid') # take 12-month running average
time_obs = time_obs[0:len(temp_obs)] 
ax.plot(time_obs, temp_obs, color='k', linewidth=4, linestyle='-')


# now load the coarse SST projections, for every reef location, to compare
coarsefile = outdir+'mean_tos_nature_ssp370_1yrs_dhwbma_1975-1995_1995-2015_20032017.nc'
print(coarsefile)
coarse_ds = xr.open_dataset(coarsefile)
coarse_ds = coarse_ds.where(coarse_ds.isreef==1, drop=True)
coarse_ds = coarse_ds.where(coarse_ds.time>=startyear, drop=True)
coarse_ds = coarse_ds.where(coarse_ds.time<=endyear, drop=True) 
temp_coarse = coarse_ds.tos.values
time_coarse = coarse_ds.time.values
temp_coarse = np.nanmean(temp_coarse, axis=0)
temp_coarse = np.convolve(temp_coarse, np.ones(N)/N, mode='valid') # take 12-month running average
time_coarse = time_coarse[0:len(temp_coarse)]  
ax.plot(time_coarse, temp_coarse, color='b', linewidth=3, linestyle='-')
  
coarsefile = '/home/pkalmus/projects//coral/output/nature/fine/bma/production_from_ayesha/mean_tos_nature_ssp370_1yrs_dhwbma_1975-1995_1995-2015.nc'
print(coarsefile)
coarse_ds = xr.open_dataset(coarsefile)
coarse_ds = coarse_ds.where(coarse_ds.isreef==1, drop=True)
coarse_ds = coarse_ds.where(coarse_ds.time>=startyear, drop=True)
coarse_ds = coarse_ds.where(coarse_ds.time<=endyear, drop=True) 
temp_coarse = coarse_ds.tos.values
time_coarse = coarse_ds.time.values
# print number of nans.
aa = temp_coarse[:,0]
print(np.sum(np.isnan(aa)))
temp_coarse = np.nanmean(temp_coarse, axis=0)
temp_coarse = np.convolve(temp_coarse, np.ones(N)/N, mode='valid') # take 12-month running average
time_coarse = time_coarse[0:len(temp_coarse)]  
ax.plot(time_coarse, temp_coarse, color='c', linewidth=3, linestyle='-')


coarsefile = 'test.nc'
print(coarsefile)
coarse_ds = xr.open_dataset(coarsefile)
coarse_ds = coarse_ds.where(coarse_ds.isreef==1, drop=True)
coarse_ds = coarse_ds.where(coarse_ds.time>=startyear, drop=True)
coarse_ds = coarse_ds.where(coarse_ds.time<=endyear, drop=True) 
temp_coarse = coarse_ds.tos.values
time_coarse = coarse_ds.time.values
# print number of nans.
aa = temp_coarse[:,0]
print(np.sum(np.isnan(aa)))
temp_coarse = np.nanmean(temp_coarse, axis=0)
temp_coarse = np.convolve(temp_coarse, np.ones(N)/N, mode='valid') # take 12-month running average
time_coarse = time_coarse[0:len(temp_coarse)]  
ax.plot(time_coarse, temp_coarse, color='m', linewidth=3, linestyle=':')


coarsefile = 'test_flat.nc'
print(coarsefile)
coarse_ds = xr.open_dataset(coarsefile)
coarse_ds = coarse_ds.where(coarse_ds.isreef==1, drop=True)
coarse_ds = coarse_ds.where(coarse_ds.time>=startyear, drop=True)
coarse_ds = coarse_ds.where(coarse_ds.time<=endyear, drop=True) 
temp_coarse = coarse_ds.tos.values
time_coarse = coarse_ds.time.values
# print number of nans.
aa = temp_coarse[:,0]
print(np.sum(np.isnan(aa)))
temp_coarse = np.nanmean(temp_coarse, axis=0)
temp_coarse = np.convolve(temp_coarse, np.ones(N)/N, mode='valid') # take 12-month running average
time_coarse = time_coarse[0:len(temp_coarse)]  
ax.plot(time_coarse, temp_coarse, color='c', linewidth=3, linestyle=':')


flatfile = outdir+'mean_tos_nature_ssp370_1yrs_flat_20032017_DHW4.8.nc'
print(flatfile)
flat_ds = xr.open_dataset(flatfile)
flat_ds = flat_ds.where(flat_ds.isreef==1, drop=True)
flat_ds = flat_ds.where(flat_ds.time>=startyear, drop=True)
flat_ds = flat_ds.where(flat_ds.time<=endyear, drop=True)   
temp_flat = flat_ds.tos.values
time_flat = flat_ds.time.values
# print number of nans.
aa = temp_flat[:,0]
print(np.sum(np.isnan(aa)))
temp_flat = np.nanmean(temp_flat, axis=0)
temp_flat = np.convolve(temp_flat, np.ones(N)/N, mode='valid') # take 12-month running average
time_flat = time_flat[0:len(temp_flat)]  
ax.plot(time_flat, temp_flat, color='r', linewidth=3, linestyle='-')



# if logx:
#     plt.xscale('log')
# plt.xlim(xlim)
# plt.ylim(ylim)
    
plt.legend(loc="upper right")          
ax.set(xlabel='year', ylabel='GMSTA')

# if fontsize is not None:
#     plt.rcParams.update({'font.size': fontsize})
ax.grid()

fig.tight_layout()
# plt.rcParams.update({'font.size': fontsize})
filename = 'sst_check.png'
plt.savefig(filename)
print(filename)
plt.close()