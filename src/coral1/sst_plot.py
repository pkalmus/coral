#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
python sst_plot.py params

conda activate geo_env

take in the MUR+laGP files and explore them. 

plot pixel-mean SST vs. time, for the SSPs.

plot some time series


@author: pkalmus
"""
import xarray as xr
import pdb
import numpy as np
import glob
import coral_plotter as plotter
import sys
import subprocess
from scipy.interpolate import interp1d


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
bias_correct = params.bias_correct
do_weights = params.do_weights
do_obs = params.do_obs
finerun = params.finerun
includemur = params.include_mur

lat_extent = 35
include_coarse = False

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


obsds = None
# read in coarse obs
obsfilename = reefdir+'HadISST_sst_reef.nc'
print(obsfilename)
obsds = xr.open_dataset(obsfilename)
obsds = obsds.where(obsds.isreef==1, drop=True)

# get estimate of warming from 1988.286 to 2008.5. 
N = 120
temp_obs = obsds.tos.values
time_obs = obsds.time.values
temp_obs = np.nanmean(temp_obs, axis=0)
temp_obs = np.convolve(temp_obs, np.ones(N)/N, mode='valid') # take 12-month running average
time_obs = time_obs[0:len(temp_obs)]     

temp1988 = interp1d(time_obs, temp_obs, bounds_error=False)(1988.286)
temp2008 = interp1d(time_obs, temp_obs, bounds_error=False)(2008.5)
print('\n*************************************\nmean 1-degree reef warming from 1988.286 to 2008.5, using 10-year running mean: %1.3f\n\n' % (temp2008-temp1988))

if include_coarse:
    obsds = obsds.where(obsds.time>=2010, drop=True)
else:
    obsds = None # so plotter knows not to plot    

ds_list = []
coarse_list = []
for scenario in scenarios:
    if do_weights:
        experiment = '%s' % (scenario) # this will be replaced by individual model if we're in single_mod
        lagp_files = glob.glob(finedir+'Globe*%s_weighted*.nc' % (scenario)) # make sure you don't glob output files from THIS code
        codestr = '1yrs_dhwbma_1975-1995_1995-2015_19932007'
    else:
        experiment = '%s_%s' % (finerun, scenario) # this will be replaced by individual model if we're in single_mod
        #lagp_files = glob.glob(finedir+'Globe*%s_uniform*.nc' % (scenario)) # make sure you don't glob output files from THIS code            

        lagp_files = glob.glob(finedir+'Globe*%s*.nc' % (scenario)) # make sure you don't glob output files from THIS code            

        codestr = '1yrs_flat'
    print('++++++starting SST ' + experiment)


    if include_coarse:
    # now load the coarse SST projections, for every reef location, to compare
        coarsefile = outdir+'mean_tos_%s_%s_%s.nc' % (runname, scenario, codestr)
        print(coarsefile)
        coarse_ds = xr.open_dataset(coarsefile)
        coarse_ds = coarse_ds.where(coarse_ds.isreef==1, drop=True)
        coarse_ds = coarse_ds.where(coarse_ds.time>=2010, drop=True)
        coarse_list.append(coarse_ds)
        #pdb.set_trace()
    
    if len(lagp_files) == 0:
        raise ValueError('no lagp files in: %s' % (finedir))
        
        
    fileds_list = []
    for lagp_file in lagp_files:
        dslagp = xr.open_dataset(lagp_file)
        dslagp = dslagp.dropna('index', subset=['temp'])            
        # merge in the 2002-2019 mur
        
        if includemur:
            # MUR goes to Dec. 2019. downscale starts in Jan. 2018. cut out overlap
            dslagp = dslagp.sel(time=slice(2020,2100))
        
            # sorted.130_150.min40_0_MUR.ncdf4.nc  Globe.130_150.min40_0.ssp126_weighted.ncdf4.nc
            murfile = murdir+'sorted'+lagp_file[lagp_file.find('Globe')+5:lagp_file.find('ssp')-1]+'_MUR.ncdf4.nc'
            dsmur = xr.open_dataset(murfile)
            dsmur['index'] = np.arange(1, len(dsmur.index.values)+1)
            dslagp['index'] = np.arange(1, len(dslagp.index.values)+1)
            #pdb.set_trace()
            ds = xr.concat([dsmur, dslagp], 'time')
            
            # for some reason, concat causes lat and lon to have (time, index) instead of (index)
            ds['lon'] = dslagp.lon
            ds['lat'] = dslagp.lat
        
            # print(dslagp.temp.shape)
            # print(dsmur.temp.shape)
            # print(ds.temp.shape)
        else:
            ds = dslagp
            
        fileds_list.append(ds)
        #pdb.set_trace()
        
    lagp_ds = xr.concat(fileds_list, 'index') # index will no longer be a monotonic index
    # added this for onepermodel. probably the "less than ten" rule but did not give Ayesha a new location file.
    lagp_ds = lagp_ds.dropna('index', subset=['lat']) 
    
    
    # if finerun=='weightbc':
    #     lagp_ds = lagp_ds.where(lagp_ds.lat==slice(-2, 1), drop=True)
    
    #plotter.plot_data_array(lagp_ds.sel(time=2020.01).temp, finedir+'fine2020.png')
    #plotter.plot_data_array(lagp_ds.sel(time=2099.01).temp, finedir+'fine2099.png')
        
    lons = lagp_ds.lon.values
    lats = lagp_ds.lat.values
    #pdb.set_trace()
    sst = lagp_ds.temp.values
    sd = lagp_ds.temp_sd.values # standard deviation
    time = lagp_ds.time.values
    time = np.floor(time) + ((time - np.floor(time))*30*100 - 15)/365. 
    index = lagp_ds.index.values
    index = np.arange(len(index)) # fixed it
    lagp_ds = lagp_ds.assign_coords(index=(index))

    ds_list.append(lagp_ds)  
    #pdb.set_trace()    
plotter.sstlist(ds_list, coarse_list, obsds, 'Year', finedir+'sst_%s_%s.png' % (runname, experiment), fontsize=17, ylabel='SST ($^\circ$C)', xlim=[2010,2099],axvspan=[2010,2020],label=None)


# plot mean global surface air temp anomaly tas
dsanom_list = []
for scenario in scenarios:
    if do_weights:
        experiment = '%s' % (scenario) # this will be replaced by individual model if we're in single_mod
        lagp_files = glob.glob(finedir+'Globe*%s_weighted*.nc' % (scenario)) # make sure you don't glob output files from THIS code
    else:
        experiment = '%s_%s' % (finerun, scenario) # this will be replaced by individual model if we're in single_mod
        lagp_files = glob.glob(finedir+'Globe*%s_flat*.nc' % (scenario)) # make sure you don't glob output files from THIS code            
    print('++++++starting GMSTA (air temp)  ' + experiment)

    dsanom = xr.open_dataset('/home/pkalmus/projects//coral/data/tas/nature/mean_%s.nc' % (scenario))
    dsanom_list.append(dsanom)
plotter.gmstalist(dsanom_list, 'Year', finedir+'gmsta_%s_%s.png' % (runname, experiment), fontsize=17, ylabel='GMST Anomaly ($^\circ$C)', xlim=[2010,2099], ylim=[1,6], label=None) 
       