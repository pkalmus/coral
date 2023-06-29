#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 09:22:00 2020

@author: pkalmus
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

conda activate geo_env

This validation requires running coarse_project first.

Read in the .nc files written by coarse_project.py: obs and all the runs listed
For each run:
    Take the [mean DHW over time period] difference between the run and the obs
    Take the [max DHW over time period] difference between the run and the obs
    Make plots: a map, a histogram.


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
do_plots = params.do_plots
J0 = params.J0 # to find the proper weight files
return_years = params.return_years
DHW_threshold = params.DHW_threshold
bias_correct = params.bias_correct
do_weights = params.do_weights
do_obs = params.do_obs
startyear = params.valstartyear
endyear = params.valendyear

# define in and out dirs
outdir = basedir+'data/tos/%s/coarse/' % (runname)



obsfile = outdir+'departure_nn_%s_obs.nc' % (runname)
obs_ds = xr.open_dataset(obsfile)
obs_start_year = obs_ds.time[0]
obs_end_year = obs_ds.time[-1]

# e.g. departure_nn_nature_ssp126_1yrs_weightJ02_bc.nc
experiments = ['1yrs', '1yrs_bc', '1yrs_weightJ02', '1yrs_weightJ02_bc', '1yrs_weightJ03_bc']
experiments = ['1yrs_weightJ02_bc', '1yrs_weightJ07_bc', '1yrs_bc',] # put the weighted run first, to get the DQ20 locations for the flat mean run

pvalsum_da = None # trick to find the DQ pixels first with weights, then apply them to other experiments
for run in experiments:
    myfile = outdir+'departure_nn_%s_ssp126_%s.nc' % (runname, run) # we're only using historical parts (1850-2014) for validation so scenario doesn't matter
    mod_ds = xr.open_dataset(myfile)
    
    # only take mod_ds for years with obs
    mod_ds = mod_ds.sel(time=slice(obs_start_year, obs_end_year))
        
    # switch to the start/end time slice; plot maps of mean and max DHW
    obs_ds = obs_ds.sel(time=slice(startyear, endyear))
    mod_ds = mod_ds.sel(time=slice(startyear, endyear))
    
    # calculate mean DHW over time range for obs, mod
    # subtract data sets; it subtracts lat and lon too, so correct for this
    diff_ds = mod_ds - obs_ds
    diff_ds = diff_ds.assign_coords(lat=obs_ds.lat, lon = obs_ds.lon) 
    
    # make a map of the mean difference over the time period (average over every month, not e.g. the yearly maximums)
    plotter.scatter(diff_ds.lon.values, diff_ds.lat.values, outdir+'map_dhw_%i%i_meandiff_%s_%s.png' % (startyear, endyear, run, runname), c=diff_ds.dhw.mean(dim='time', skipna=True).values, vmin=-5, vmax=5)
    plotter.hist(diff_ds.dhw.mean(dim='time', skipna=True).values, 'DHW', outdir+'hist_dhw_%i%i_meandiff_%s_%s.png' % (startyear, endyear, run, runname), bins=40)

    # calculate max DHW over time range for obs, mod at each pixel
    # subtract data sets; it subtracts lat and lon too, so correct for this
    obs_da = obs_ds.dhw.max('time', skipna=True)
    mod_da = mod_ds.dhw.max('time', skipna=True)
    diff_da = mod_da - obs_da
    diff_da = diff_da.assign_coords(lat=obs_ds.lat, lon = obs_ds.lon) 
    
    # make a map of the mean difference over the time period (average over every month, not e.g. the yearly maximums)
    plotter.scatter(diff_da.lon.values, diff_da.lat.values, outdir+'map_dhw_%i%i_maxdiff_%s_%s.png' % (startyear, endyear, run, runname), c=diff_da.values, vmin=-5, vmax=5)
    plotter.hist(diff_da.values, 'DHW', outdir+'hist_dhw_%i%i_maxdiff_%s_%s.png' % (startyear, endyear, run, runname), bins=40)
    
    if pvalsum_da is None:
        # now only keep ones with pvalsum > 20
        pvalsum_da = mod_ds.pvalsum
    
    diff_da_dq0 = diff_da.where((pvalsum_da > 1) & (pvalsum_da < 10))
    plotter.hist(diff_da_dq0.values, 'DHW', outdir+'hist_dhw_%i%i_maxdiff_DQ1-10_%s_%s.png' % (startyear, endyear, run, runname), bins=40)
    
    
#    pdb.set_trace()
#    
#    
#    #plotter.scatter(obs_ds.lon.values, obs_ds.lat.values, outdir+'map_dhw20152018_mean_%s.png' % ('obs'), c=obs_ds.dhw.mean(dim='time', skipna=True).values, vmin=0.0, vmax=5)
#    plotter.scatter(mod_ds.lon.values, mod_ds.lat.values, outdir+'map_dhw_%i%i_mean_%s_%s.png' % (startyear, endyear, run, runname), c=mod_ds.dhw.mean(dim='time', skipna=True).values, vmin=0.0, vmax=5)
#    #plotter.scatter(obs_ds.lon.values, obs_ds.lat.values, outdir+'map_dhw20152018_max_%s.png' % ('obs'), c=obs_ds.dhw.max(dim='time', skipna=True).values, vmin=0.0, vmax=8)
#    plotter.scatter(mod_ds.lon.values, mod_ds.lat.values, outdir+'map_dhw_%i%i_max_%s_%s.png' % (startyear, endyear, run, runname), c=mod_ds.dhw.max(dim='time', skipna=True).values, vmin=0.0, vmax=8)
#    
#    #pdb.set_trace()
#    plotter.timeseries(mod_ds.time.values, mod_ds.dhw.mean(dim='index', skipna=True).values, outdir+'series_dhw_%i%i_max_%s_%s.png' % (startyear, endyear, run, runname), xlabel='Year', ylabel='Mean DHW', ylim=[0,2.25])
#
#    # plot the number of 8DHW surpassings per year
#    # could I use da.to_masked_array() here?
#    mask_8 = mod_ds.where(mod_ds.dhw >= 8.0)
#    #pdb.set_trace()
#    plotter.timeseries(mod_ds.time.values, mask_8.dhw.count(dim='index').values, outdir+'count_8dhw_%i%i_max_%s_%s.png' % (startyear, endyear, run, runname), xlabel='Year', ylabel='Number of 8 DHW events')


