#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
python fine_plot.py params
make the plots of the fine files created by fine_project


conda activate geo_env

regrid_bcdp
coral_locations
regrid_reef
weighting
coarse_project TODO: redo climatology in coral_location; why do two of the models fail to get all 4000+ NN locations?
validating code
fine_project
fine_plot

Read in fine file
Make a map of 2020 and 2100
Read in climatology (need it on same grid) 
Calculate departure
Make map


@author: pkalmus
"""
import xarray as xr
import numpy as np
import pdb
import coral_plotter as plotter
import sys

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
DHW_thresholds = params.DHW_thresholds
bias_correct = params.bias_correct
do_weights = params.do_weights
do_obs = params.do_obs
finerun = params.finerun
includemur = params.include_mur
do_add_variability = params.do_add_variability
spike_size = params.spike_size
spike_freq = params.spike_freq
include_uncert = params.include_uncert
do_slow_plots = params.do_slow_plots

lat_extent = 35
climatology_type = 'mur'

if do_obs:
    scenarios = ['obs']

# define in and out dirs
finedir = basedir+'output/%s/fine/%s/' % (runname, finerun)

for return_year in return_years:
    for DHW_threshold in DHW_thresholds:
        ds_list = []
        dsanom_list = []
        for scenario in scenarios:
            experiment = '%s_%s_%iyrs_DHW%1.1f' % (finerun, scenario, return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
            print(experiment)
            if do_add_variability:
                # this is 1m 5f
                experiment = 'SPIKE_%s_%iyrs_DHW%1.1f' % (scenario, return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
                # this is 2m 5f
                #experiment = 'SPIKE2_%s_%iyrs_DHW%1.1f' % (scenario, return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
                #experiment = 'SPIKE%im%if_%s_%iyrs_DHW%1.1f' % (spike_size, spike_freq, scenario, return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
                
            
            if includemur:
                finefile = finedir+'finemur_%s_%s.nc' % (runname, experiment)
            else:
                finefile = finedir+'fine_%s_%s.nc' % (runname, experiment)
            ds = xr.open_dataset(finefile)
            ds_list.append(ds)
            
            ###### this is for the GMST anomaly plot
            dsanom = xr.open_dataset('/home/pkalmus/projects//coral/data/tas/nature/mean_%s.nc' % (scenario))
            dsanom_list.append(dsanom)
            
            
            
            # ###### plot 1% only. first fimd the threshold, then make the subset. e.g. mask_8 = mod_ds.where(mod_ds.dhw >= 8.0)
            # departure_year_nonan = ds.departure_year.values[~np.isnan(ds.departure_year.values)]
            # thres1pct = np.sort(departure_year_nonan)[int(np.ceil(0.99 * len(departure_year_nonan)))]
            # ds_pct = ds.where(ds.departure_year >= thres1pct, drop=True)
            # print(len(ds.index))
            # print(len(departure_year_nonan))
            # print(len(ds_pct.index))
            # plotter.scatter(ds_pct.lon.values, ds_pct.lat.values, finedir+'depart1pct_%s_%s.png' % (runname, experiment), c=ds_pct.departure_year.values, vmin=2020, vmax=2080, figsize=(15,5), cbar_fraction=0.008, cbar_orientation='vertical')
    
            # ###### pick out some promising regions for 1% refugia
            # plotter.scatter(ds.lon.values, ds.lat.values, finedir+'departurezoom1_%s_%s.png' % (runname, experiment), c=ds.departure_year.values, vmin=2020, vmax=2080, cbar=False, extent=[-145,-140,-20,-15], ticksize=1, figsize=(5,5)) # (min(lon), max(lon), min(lat), max(lat))
            # plotter.scatter(ds.lon.values, ds.lat.values, finedir+'departurezoom2_%s_%s.png' % (runname, experiment), c=ds.departure_year.values, vmin=2020, vmax=2080, extent=[148.9,149.9,-21,-20], ticksize=0.2, figsize=(5,5), cbar_fraction=0.02, cbar_orientation='vertical') # (min(lon), max(lon), min(lat), max(lat))       
            # pdb.set_trace()
            
            ###### the "normal" plots
            if do_slow_plots:
                plotter.scatter(ds.lon.values, ds.lat.values, finedir+'departure_%s_%s.png' % (runname, experiment), c=ds.departure_year.values, vmin=2010, vmax=2050, marker_relsize=0.1, figsize=(15,5), cbar_fraction=0.008, cbar_orientation='vertical',fontsize=19)
            # plotter.hist(ds.departure_year.values, 'Years', finedir+'hist_departure_%s_%s.png' % (runname, experiment), xlim=[2020,2100], ylim=[1.0,1.0e6], logy=True )
            # plotter.hist(ds.departure_year.values, 'Years', finedir+'cumul_departure_%s_%s.png' % (runname, experiment), xlim=[2020,2100], cumulative=1, normed=True, label=None, histtype='step')    
            
            # plotter.scatter(ds.lon.values, ds.lat.values, finedir+'departure_high_%s_%s.png' % (runname, experiment), c=ds.departure_year_high.values, vmin=2020, vmax=2080, figsize=(15,5), cbar_fraction=0.008, cbar_orientation='vertical')
            # plotter.hist(ds.departure_year_high.values, 'Years', finedir+'hist_departure_high_%s_%s.png' % (runname, experiment), xlim=[2020,2100], ylim=[1.0,1.0e6], logy=True )
        
            # plotter.scatter(ds.lon.values, ds.lat.values, finedir+'departure_low_%s_%s.png' % (runname, experiment), c=ds.departure_year_low.values, vmin=2020, vmax=2080, figsize=(15,5), cbar_fraction=0.008, cbar_orientation='vertical')
            # plotter.hist(ds.departure_year_low.values, 'Years', finedir+'hist_departure_low_%s_%s.png' % (runname, experiment), xlim=[2020,2100], ylim=[1.0,1.0e6], logy=True )
        
            # plotter.scatter(ds.lon.values, ds.lat.values, finedir+'departure_highdiff_%s_%s.png' % (runname, experiment), c=ds.departure_year_high.values-ds.departure_year.values, vmin=0, vmax=20, figsize=(15,5), cbar_fraction=0.008, cbar_orientation='vertical')
            # plotter.hist(ds.departure_year_high.values-ds.departure_year.values, 'Years', finedir+'hist_departure_highdiff_%s_%s.png' % (runname, experiment))
        
            # plotter.scatter(ds.lon.values, ds.lat.values, finedir+'departure_lowdiff_%s_%s.png' % (runname, experiment), c=ds.departure_year_low.values-ds.departure_year.values, vmin=-20, vmax=0, figsize=(15,5), cbar_fraction=0.008, cbar_orientation='vertical')
            # plotter.hist(ds.departure_year_low.values-ds.departure_year.values, 'Years', finedir+'hist_departure_lowdiff_%s_%s.png' % (runname, experiment))
        
    
        # plotter.histlist(ds_list, 'Year', finedir+'cumul_refugia_%s_%s.png' % (runname, experiment), fontsize=19, xlim=[2021,2099], bins=10000, cumulative=-1, normed=False, label=None, histtype='step', logy=True, ylim=[1,1e6], dolines=False)
        
        
        #plotter.histlist(ds_list, 'Year', finedir+'cumul_pct_%s_%s.png' % (runname, experiment), fontsize=19, ylabel='Fraction of reefs under threshold', xlim=[2021,2099], bins=10000, cumulative=-1, normed=True, label=None, histtype='step', logy=True, ylim=[1e-5,1.0], dolines=True)

        if DHW_threshold == 8:
            title_str = 'TD%iY, 11.2 DHW$_{1988}$ / 8 DHW$_{2008}$' % (return_year)  
        else:
            title_str = 'TD%iY, 8 DHW$_{1988}$ / 4.8 DHW$_{2008}$' % (return_year)  
        title_str = '' # decided not to have these titles on the figures, put it in the caption instead.


        ## for some reason, the first plot doesn't change the font size. so do this one twice.
        # plotter.histlist(ds_list, 'Year', finedir+'cumul_pct_%s_%s.png' % (runname, experiment), fontsize=19, ylabel='Fraction of reefs remaining', title=title_str, xlim=[2010,2095], axvspan=[2010,2020], bins=10000, cumulative=-1, normed=True, label=None, histtype='step', logy=False, ylim=[0, 1.0], dolines=True, include_uncert=include_uncert)
        # plotter.histlist(ds_list, 'Year', finedir+'cumul_pct_%s_%s.png' % (runname, experiment), fontsize=19, ylabel='Fraction of reefs remaining', title=title_str, xlim=[2010,2095], axvspan=[2010,2020], bins=10000, cumulative=-1, normed=True, label=None, histtype='step', logy=False, ylim=[0, 1.0], dolines=True, include_uncert=include_uncert)

        plotter.histlist(ds_list, 'Year', finedir+'cumul_pctLOGY_%s_%s.png' % (runname, experiment), txtfilename=finedir+'year_%s_%s.txt' % (runname, experiment), fontsize=19, ylabel='Fraction of reefs remaining', title=title_str, xlim=[2010,2095], axvspan=[2010,2020], bins=10000, cumulative=-1, normed=True, label=None, histtype='step', logy=True, ylim=[1e-3, 1.0], dolines=True, include_uncert=include_uncert)
     
     
        # plotter.histlist(ds_list[1:], 'GMST Anomaly (degrees C)', finedir+'cumul_temp_%s_%s.png' % (runname, experiment), anomlist=dsanom_list, colors=['b', 'g', 'r'], fontsize=19, ylabel='Fraction of reefs remaining', title=title_str, xlim=[1.0, 2.2], bins=10000, cumulative=-1, normed=True, label=None, histtype='step', logy=False, ylim=[0, 1.0], dolines=True, include_uncert=include_uncert)
        # plotter.histlist(ds_list[1:], 'GMST Anomaly (degrees C)', finedir+'cumul_tempLOGY_%s_%s.png' % (runname, experiment), txtfilename=finedir+'temp_%s_%s.txt' % (runname, experiment), anomlist=dsanom_list, colors=['b', 'g', 'r'], fontsize=19, ylabel='Fraction of reefs remaining', title=title_str, xlim=[1.0, 2.2], bins=10000, cumulative=-1, normed=True, label=None, histtype='step', logy=True, ylim=[1e-3, 1.0], dolines=True, include_uncert=include_uncert)
            