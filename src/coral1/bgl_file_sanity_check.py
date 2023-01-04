#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

have a look at the bgl_ onepermodel run and compare to the bgl_flat run.
val_end2end gave some wonky looking histograms for onpermodel. 

@author: pkalmus
"""
import xarray as xr
import pdb
import numpy as np
import glob
import sys
import subprocess
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
downscale_type = params.downscale_type
include_downscale_uncert = params.include_downscale_uncert

DHW_threshold = DHW_thresholds[0]
if len(DHW_thresholds) > 1:
    print('NOTE: ONLY RUNNING ON FIRST DHW_THRESHOLD. loop over DHW_thresholds not implemented. change params.py to run on additional.')

lat_extent = 35
climatology_type = 'mur'

if do_obs:
    scenarios = ['obs']

# define in and out dirs
#finedir = basedir+'data/tos/%s/fine/%s/' % (runname, finerun) # e.g. for scp: /raid8/pkalmus/data/coral/data/tos/nature/fine/flat/
finedir = basedir+'output/%s/fine/%s/' % (runname, finerun) # e.g. for scp: /raid8/pkalmus/data/coral/data/tos/nature/fine/flat/
print(finedir)
murdir = basedir + 'data/mur/MUR_coral-pixels_2002-2019/'

use_best_only = False

for scenario in scenarios:
    if finerun=='bgl_flat':
        lagp_files = glob.glob(finedir+'Globe*%s_uniform*.nc' % (scenario))
    elif finerun=='bgl_onepermodel':
        lagp_files = glob.glob(finedir+'Globe*%s*.nc' % (scenario)) # make sure you don't glob output files from THIS code 

    else:
        raise(ValueError('unknown finerun: %s' % finerun))               
    print('++++++starting ' + scenario)
    print(lagp_files)
    
    if len(lagp_files) == 0:
        print('no lagp files in: %s' % (finedir))
        
    ds_list = []
    for lagp_file in lagp_files:
        dslagp = xr.open_dataset(lagp_file)
        if downscale_type=='trend':
            dslagp['temp_sd'] = (dslagp.temp.dims, np.zeros(dslagp.temp.values.shape))
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
            ds = xr.concat([dsmur, dslagp], 'time')
            #pdb.set_trace()
            
            # for some reason, concat causes lat and lon to have (time, index) instead of (index)
            ds['lon'] = dslagp.lon
            ds['lat'] = dslagp.lat
        
            print(dslagp.temp.shape)
            print(dsmur.temp.shape)
            print(ds.temp.shape)
        else:
            ds = dslagp
            
        ds_list.append(ds)
        
    lagp_ds = xr.concat(ds_list, 'index') # index will no longer be a monotonic index
    
    # added this for onepermodel. extra NaNs here, probably due to the GCM "less than ten" rule but did not give Ayesha a new location file.
    lagp_ds = lagp_ds.dropna('index', subset=['lat']) 
    
    # get 2019 only
    
    print('total min: %1.2f' % (np.min(lagp_ds.temp.values)))
    print('total max: %1.2f' % (np.max(lagp_ds.temp.values)))        

    # make a histogram
    lagp_ds = lagp_ds.sel(time=slice(2018, 2021))
    
    plotter.hist(lagp_ds.temp.values, 'Celsius', finedir+'hist_2018-2020_%s_%s.png' % (scenario, finerun), bins=100 )
    
    pdb.set_trace()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        