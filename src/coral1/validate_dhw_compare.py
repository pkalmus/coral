#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

conda activate geo_env

Read in the .nc files written by coarse_project.py: obs and model 
Over the given year range, calculate the mean difference at each pixel between model DHW and obs DHW
Plot map of these differences, and histogram





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


basedir = projectdir+'/coral/'
datadir = basedir+'data/tos/'
outdir = basedir+'output/initial/'
runname = 'test'

startyear = 1996
endyear = 2019

startyear = 2015
endyear = 2018

obsfile = outdir+'departure_obs_1.nc'
obs_ds = xr.open_dataset(obsfile)
obs_start_year = obs_ds.time[0]
obs_end_year = obs_ds.time[-1]


# compare the two runs to obs
runs = ['all_ssp585', 'all_ssp585_weighted']

myfile = outdir+'departure_%s_1.nc' % (runs[0])
flat_ds = xr.open_dataset(myfile)
myfile = outdir+'departure_%s_1.nc' % (runs[1])
weight_ds = xr.open_dataset(myfile)

# only take mod_ds for years with obs
flat_ds = flat_ds.sel(time=slice(obs_start_year, obs_end_year))
weight_ds = weight_ds.sel(time=slice(obs_start_year, obs_end_year))

# now assume that weight_ds is smaller than flat_ds, and get the common set
flat_ds = flat_ds.where(np.invert(xr.ufuncs.isnan(weight_ds)) )
obs_ds = obs_ds.where(np.invert(xr.ufuncs.isnan(weight_ds)) )

flat_dhw = np.ravel(flat_ds.dhw.values)
weight_dhw = np.ravel(weight_ds.dhw.values)
obs_dhw = np.ravel(obs_ds.dhw.values)

# make diff histograms between flat, weight and obs
nbins = 50
globalmin = 0 
globalmax = 50
#histogram is build with fixed min and max values
histflat, edges = np.histogram(flat_dhw,range=(globalmin,globalmax), bins=nbins)
histweight, edges = np.histogram(weight_dhw,range=(globalmin,globalmax), bins=nbins)
histobs, edges = np.histogram(obs_dhw,range=(globalmin,globalmax), bins=nbins)

diff_flat = histflat - histobs
plotter.bar((edges[0:-1]+edges[1:])/2.0, diff_flat, 'DHW', outdir+'hist_diff_flat.png')
diff_weight = histweight - histobs
plotter.bar((edges[0:-1]+edges[1:])/2.0, diff_weight, 'DHW', outdir+'hist_diff_weight.png')

pdb.set_trace()


#    modfile = outdir+'departure_all_ssp585_1.nc'
#    mod_ds = xr.open_dataset(modfile)


# plot histogram of DHW values for entire model set
#plotter.hist(np.ravel(obs_ds.dhw.values), 'DHW', outdir+'hist_dhw_%s.png' % ('obs'), logy=True, bins=40)
plotter.hist(np.ravel(mod_ds.dhw.values), 'DHW', outdir+'hist_dhw_%s.png' % (run), logy=True, bins=40)

# switch to the start/end time slice; plot maps of mean and max DHW
#obs_ds = obs_ds.sel(time=slice(2015,2018))
mod_ds = mod_ds.sel(time=slice(startyear, endyear))
#plotter.scatter(obs_ds.lon.values, obs_ds.lat.values, outdir+'map_dhw20152018_mean_%s.png' % ('obs'), c=obs_ds.dhw.mean(dim='time', skipna=True).values, vmin=0.0, vmax=5)
plotter.scatter(mod_ds.lon.values, mod_ds.lat.values, outdir+'map_dhw_%i%i_mean_%s.png' % (startyear, endyear, run), c=mod_ds.dhw.mean(dim='time', skipna=True).values, vmin=0.0, vmax=5)
#plotter.scatter(obs_ds.lon.values, obs_ds.lat.values, outdir+'map_dhw20152018_max_%s.png' % ('obs'), c=obs_ds.dhw.max(dim='time', skipna=True).values, vmin=0.0, vmax=8)
plotter.scatter(mod_ds.lon.values, mod_ds.lat.values, outdir+'map_dhw_%i%i_max_%s.png' % (startyear, endyear, run), c=mod_ds.dhw.max(dim='time', skipna=True).values, vmin=0.0, vmax=8)

#pdb.set_trace()
plotter.timeseries(mod_ds.time.values, mod_ds.dhw.mean(dim='index', skipna=True).values, outdir+'series_dhw_%i%i_max_%s.png' % (startyear, endyear, run), xlabel='Year', ylabel='Mean DHW', ylim=[0,2.25])

# plot the number of 8DHW surpassings per year
# could I use da.to_masked_array() here?
mask_8 = mod_ds.where(mod_ds.dhw >= 8.0)
#pdb.set_trace()
plotter.timeseries(mod_ds.time.values, mask_8.dhw.count(dim='index').values, outdir+'count_8dhw_%i%i_max_%s.png' % (startyear, endyear, run), xlabel='Year', ylabel='Number of 8 DHW events')


