#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
To run on single model, pass in the appropriate integer argument.
To run on HadISST, pass in -1.

conda activate geo_env

regrid_bcdp. create the concatenated, regridded coarse files
coral_locations
regrid_reef 
weighting
coarse_project TODO: redo climatology in coral_location; why do two of the models fail to get all 4000+ NN locations?
validating code

NOTE: the main purpose of this file was to create the mean_tos....nc files, for inputs to the LaGP
The secondary purpose is to do the sensitivity test. This uses different DHW and climatologies from HadISST
Originally this code used to use the 4 km Pathfinder climatology. 


Start with ASB
** repeat with 5-year SB.

Create an ensemble record: flat mean, weighted mean, or best.
  Use obs

Read in low-res reef locations reef_grid_1x1.txt (created by coral_locations.py, using the climatology)

For each reef location:
    Get MMM (maximum of monthly mean)
    for each model:
        Calculate DHW on a 3-month cumulative sliding window.
        Find departure over 8 DHW.
    calculate average and weighted average

save in netcdf (individual model departures, and averages)    
Make a departure time map.


12 week sliding window for accumulating DHW
https://coralreefwatch.noaa.gov/satellite/dhw.php

https://coralreefwatch.noaa.gov/satellite/methodology/methodology.php#hotspot

Nighttime temps.

get the hottest climatological month for that region - that's the baseline.
Maximum of the Monthly Mean (MMM) SST climatology

Here's how CRW does it: DHWs = 0.5 * Summation of previous 24 twice-weekly HotSpots, where HotSpots have to be at least 1.0 °C to be accumulated

Here's how we'll do it:
    for each month, compare to climatology
    monthly DHW = the difference x 4 (if hotter than climatology, else zero)
    take the monthly DHW and for each month, sum that month and the preceding two.
    then make a new time series that holds the max DHW for the whole year


?? how do we integrate with the xarray mean? probably use all_ds instead of mean_ds and loop through models?
https://docs.scipy.org/doc/numpy/reference/generated/numpy.average.html
For weights:
    read in weight file
    if all weights for a location are zero [nan or flat average], else weighted average
    proceed    


https://coralreefwatch.noaa.gov/satellite/methodology/methodology.php#clim

https://coralreefwatch.noaa.gov/satellite/thermal_history/climatology.php

https://coralreefwatch.noaa.gov/satellite/thermal_history/climatology_monthly.php


Created on Wed Sep  4 14:52:12 2019

@author: pkalmus
"""
import xarray as xr
import pdb
import numpy as np
import coral_plotter as plotter
import sys
import subprocess

# read in user params
import importlib
paramfile = sys.argv[1]   
params = importlib.import_module(paramfile)
projectdir = params.projectdir
runname = params.runname
tosrunname = params.tosrunname
basedir = params.basedir
pythondir = params.pythondir
listfilename = params.listfilename
scenarios = params.scenarios
do_plots = params.do_plots

########################################
return_years = params.return_years
########################################

DHW_threshold = params.DHW_thresholds[0]
print('Note: only using first DHW_threshold... loop over DHW_thresholds not implemented')
bias_correct = params.bias_correct
do_weights = params.do_weights
do_obs = params.do_obs

meanstart=params.meanstart
meanend=params.meanend
min_models=params.min_models
climstartyear_bma = params.climstartyear_bma
climendyear_bma = params.climendyear_bma
climstartyear_project = params.climstartyear_project
climendyear_project = params.climendyear_project




# https://github.com/pydata/xarray/issues/422
def average_da(myda, dim=None, weights=None):
    """
    weighted average for DataArrays

    Parameters
    ----------
    dim : str or sequence of str, optional
        Dimension(s) over which to apply average.
    weights : DataArray
        weights to apply. Shape must be broadcastable to shape of self.

    Returns
    -------
    reduced : DataArray
        New DataArray with average applied to its data and the indicated
        dimension(s) removed.

    """
    if weights is None:
        return myda.mean(dim)
    else:
        if not isinstance(weights, xr.DataArray):
            raise ValueError("weights must be a DataArray")

        # if NaNs are present, we need individual weights
        if myda.notnull().any():
            total_weights = weights.where(myda.notnull()).sum(dim=dim, skipna=True)
        else:
            total_weights = weights.sum(dim, skipna=True)

        return (myda * weights).sum(dim, skipna=True) / total_weights


def get_dhw(mmm, sst_time_series):
    # calculate DHW. 
    dhw_mo = sst_time_series - mmm
    dhw_mo = dhw_mo*4.34 #convert from months to weeks
    dhw_mo[dhw_mo < 0] = 0
    dhw_mo = np.insert(dhw_mo, 0, np.array([0.0,0.0]))
    dhw = dhw_mo[2:]+dhw_mo[1:-1]+dhw_mo[0:-2] # 3-month running sum
    dhw = np.insert(dhw, 0, np.array([np.nan])) # shift to right once
    dhw = dhw[0:-1] 
    return dhw





if do_obs:
    scenarios = ['obs']

reefdir = basedir+'data/tos/%s/reef/' % (tosrunname) 
weightdir = basedir+'data/tos/%s/weightbma/' % (tosrunname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/

if do_weights:
    # the specified weight file must exist.
    # to create using BMA:
    
    filetag = 'dhwbma_%i-%i_%i-%i' % (climstartyear_bma, climendyear_bma, meanstart, meanend)
else:
    filetag = 'flat'

outdir = basedir+'data/tos/%s/coarse/' % (runname)
#outdir = basedir+'data/tos/nature/coarse/subdegree/'

sysStr = "mkdir -p %s" % outdir
subprocess.call(sysStr, shell=True)

with open(listfilename) as f:
    models = f.read().splitlines()

# if you want to run on just one of these models, e.g. for validation, use an arg with the model index from model_list
if len(sys.argv) != 3:
    single_mod = None
else:
    single_mod = int(sys.argv[2])

# read in obs
obsfilename = reefdir+'HadISST_sst_reef.nc'
print(obsfilename)
obs_ds = xr.open_dataset(obsfilename)


# calculate obs mmm climatology from HadISST
clim_ds = obs_ds.sel(time=slice(climstartyear_project, climendyear_project))
# take annual mean
clim_tos = clim_ds.tos.values
clim_mo = np.nanmean(clim_tos.reshape(clim_tos.shape[0],-1,12), 1)
# find mmm for each index
clim_mmm = np.max(clim_mo, 1)

if bias_correct:
    start_year = 1918
    nmonths = 1024
    # note: this vectorized version is a million times faster than the loop version
    myobs = obs_ds.sel(time=slice(start_year, np.floor(start_year + nmonths/12.0) ))
    obs_mean = myobs.tos.mean(dim='time', skipna=True).values


    
for scenario in scenarios:
    experiment = '%s_%iyrs_%s_%i%i' % (scenario, return_years, filetag, climstartyear_project, climendyear_project) # this will be replaced by individual model if we're in single_mod
    if DHW_threshold != 8:
        experiment = experiment + '_DHW%1.1f' % (DHW_threshold)
    if bias_correct:
        experiment = experiment + '_bc'
    print('\n\n***********'+ experiment)
    
    # average all models' tos
    # create new DS with a new model dim, and then average over that dim.
    ds_list = []
    weights_ds_list = []
    
    if scenario=='obs': # we're doing obs, for validation purposes. weights have no meaning here
        filename = reefdir+'HadISST_sst_reef.nc'
        print(filename)
        mod_ds = xr.open_dataset(filename)
        ds_list.append(mod_ds)
        experiment = 'obs'
        do_weights = False
        bias_correct = False
    else:    
        used_models_file = open('models_%s_%s' % (runname, experiment), 'w')
        for (modelind, model) in enumerate(models):
            print('\n\n***********'+str(modelind+1) + '   ' + model)

            # if modelind > 11:
            #     continue
            
            # if model == 'CNRM-ESM2-1 r4i1p1f2 gn': # hada non-monotonic time!
            #     pdb.set_trace()
            #     continue
    
            if single_mod is not None:
                if modelind != single_mod:
                    continue 

            # e.g. CanESM5 r11i1p2f1 gn
            # UKESM1-0-LL_r8i1p1f2_gn_ssp126_reef.nc
            modelstr = model.replace(' ', '_')
            modsspstr = modelstr+'_'+scenario
            print('model, ssp: ' + modsspstr)
            filename = reefdir+modsspstr+'_reef.nc'
            
            try:
                mod_ds = xr.open_dataset(filename) #tos      (time, lat, lon). lon: 0 to 360   
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                print('file not found: %s' % filename)   
            if len(mod_ds.lat) > 4000: # a couple of models in the list couldn't run through the NN process.
                # add the model name to the list of actual models
                #modellistfile.write('%s %s\n' % (mysource, mymember))
                used_models_file.write('%s\n' % (model))

                # replace any zeros, at least until I fix this at the source in the BCDP code.
                #print(np.where(mod_ds.tos.values==0)[0])
                #mod_ds = mod_ds.where(mod_ds['tos'] != 0.) # btw this takes about half a second to do - and might add more to concat / average
                mod_ds['tos'] = mod_ds.tos.where(mod_ds.tos != 0)
                
                if do_weights:
                    weightfilename = weightdir+modelstr+'_%s.nc' % (filetag)
                    try:
                        weight_ds = xr.open_dataset(weightfilename) #tos      (time, lat, lon). lon: 0 to 360
                        weights_ds_list.append(weight_ds)  # put this here for safety, to make sure it's never double-added.
                    except (KeyboardInterrupt, SystemExit):
                        raise
                    except:
                        print('file not found: %s' % weightfilename)   
                        continue
                           
                        
                # Note: we could do this in regrid_reef instead?
                # Use the long obs. record to do a mean-match over the time-period of interest at each location
                if bias_correct:
                    print('starting bias correct...')                    
                    mymod = mod_ds.sel(time=slice(start_year, np.floor(start_year + nmonths/12.0) ))
                    mod_mean = mymod.tos.mean(dim='time', skipna=True).values
                    
                    diff = (mod_mean - obs_mean)
                    mod_ds.tos.values = mod_ds.tos.values - diff[:,None] # add the singleton dimension so it can be broadast across "time" dim
                ds_list.append(mod_ds)
                
            else:
                print('model had less than 4000 pixels...')

    all_ds = xr.concat(ds_list, 'model_index') # tos: (127, 4836, 2220).  compat='override'? https://github.com/pydata/xarray/issues/2217 https://xarray.pydata.org/en/stable/generated/xarray.concat.html
    for loopind in obs_ds.index.values: 
        num_model_estimates = np.count_nonzero(~np.isnan(all_ds.tos.values[:,loopind, 1]))
        if num_model_estimates < min_models:
            print('too few model estimates (<%i) %i for index %i' % (min_models, num_model_estimates, loopind))
            all_ds.tos.values[:,loopind,:] = np.nan
            
    # flat_ds = all_ds.mean(dim='model_index', skipna=True, keep_attrs=True)
    # flat_tos_da = flat_ds.tos    

    mean_ds = all_ds.mean(dim='model_index', skipna=True, keep_attrs=True)
    
    pvalsum_da = xr.full_like(mean_ds.tos, fill_value=np.nan) 
    if do_weights:
        weights_ds = xr.concat(weights_ds_list, 'model_index')

        # any index with fewer than 10 models w/o nan, have to nan out. e.g. 1798 has only one. 
        for loopind in obs_ds.index.values: 
            num_model_estimates = np.count_nonzero(~np.isnan(weights_ds.pval.values[:,loopind]))
            if num_model_estimates < min_models:
                print('too few model estimates %i for index %i' % (num_model_estimates, loopind))
                weights_ds.pval.values[:,loopind] = np.nan
        
        pvalsum_da = weights_ds.pval.sum('model_index', skipna=True)
        
        weighted_tos_da = average_da(all_ds.tos, dim='model_index', weights=weights_ds.pval)
        
        # pdb.set_trace()
        # zipped = zip(weighted_tos_da.time.values, np.nanmean(weighted_tos_da.values, axis=0))
        # np.savetxt('weighted.csv', zipped, delimiter=',')
        # pdb.set_trace()
        
      
# (Pdb) aa = weighted_tos_da.values
# (Pdb) aa.shape
# (4836, 2220)
# (Pdb) bb = mean_ds.tos.values
# (Pdb) bb.shape
# (4836, 2220)
# (Pdb) amean = np.nanmean(np.nanmean(aa))
# (Pdb) amean
# 27.881541198558043
# (Pdb) bmean = np.nanmean(np.nanmean(bb))
# (Pdb) bmean
# 27.939623111104453
        
        
        mean_ds['tos'] = weighted_tos_da
        


    
    null_data = np.empty((mean_ds.dims['index']))
    null_data[:] = np.nan
    dep_da = xr.DataArray(null_data, coords=([mean_ds.coords['index']]), dims=('index'))

    #dep_da = xr.DataArray(coords=([mean_ds.coords['index']]), dims=('index'))

    dhw_da = xr.full_like(mean_ds.tos, fill_value=np.nan) # not sure how to do w/ subset of dims, hence the null_data syntax above
    
    mmm_vect = [] # for validating / understanding
    for loopind in obs_ds.index.values:   
        departure_time = np.nan
        weight_dep = np.nan
        
        mmm = clim_mmm[loopind]
        mmm_vect.append(mmm)
        if ~np.isnan(mmm):
            mymod = mean_ds.sel(index=loopind)
            tos = mymod.tos.values # this read is lightning fast, unlike the lat/lon read from before the "index" days.
            # check if it has tos (is ocean)
            if np.isnan(tos[0]):
                print('tos nan')
                continue
    
            modlon = mymod.lon.values
            modlat = mymod.lat.values
            
            
            # check that lat and lon match, just in case. (they should!)
            # see the 115.49999999999999 issue.
#            if modlon!=looplon or modlat!=looplat:
#                pdb.set_trace()
                
            # calculate DHW. 
            dhw = get_dhw(mmm, tos)
            
            # store this time series at loopind
            dhw_da.loc[dict(index=loopind)] = dhw
            
            # we now have monthly time series of 3-monthly DHW.
            # we want a series of binaries that state "within N years around this month ("span") there's an excursion"
            N = 12 * return_years
            above_months = dhw > DHW_threshold
            # this doesn't center the span around the hot 3-month period; instead, it starts there and goes out "span" monhts.
            # that's waht we want, though.
            above_span = np.convolve(above_months, np.ones((N,))/N, mode='same') 
            above_span[above_span > 0] = 1
            # end issue: sometimes the last value is 0
            above_span = above_span[0:-1]
            # find departure date - get the index of the last zero in above_span. the last month, before departure. :-(
            zero_ind = np.where(above_span==0)[0]
                        
            if len(zero_ind)==0: #before start of timeseries
                zero_ind = 0
            last_ind = np.max(zero_ind)
            departure_time = mymod.time[last_ind].values

            print(loopind, modlon, modlat, mmm, tos[0], departure_time)
 
        # can I do this using xarray, and add to the mod_ds as I go?
        dep_da.loc[dict(index=loopind)] = departure_time # THIS IS SLOW in the loop. important to initialize this as nan outside the if ~np.nan(mmm)! 

                                
    # write .nc file
    ds = xr.Dataset({'departure_year':dep_da, 'dhw':dhw_da, 'pvalsum':pvalsum_da, 'isreef':mean_ds.isreef, 'lon':mean_ds.lon, 'lat':mean_ds.lat} )
    outfilename = outdir+'departure_nn_%s_%s.nc' % (runname, experiment)
    ds.to_netcdf(outfilename) 
    print(outfilename)

    # these are the coarse ensemble file that gets downscaled
    outfilename = outdir+'mean_tos_%s_%s.nc' % (runname, scenario) # note this does not depend on 1yrs, 5yrs or climatology
    #outfilename = 'test.nc'
    mean_ds.to_netcdf(outfilename) 
    print(outfilename)
    
    ds_reef = ds.where(ds.isreef==1, drop=True)
    
    plotter.scatter(ds_reef.lon.values, ds_reef.lat.values, outdir+'departure_%s_%s.png' % (runname, experiment), c=ds_reef.departure_year.values, vmin=2020, vmax=2080, cbar_fraction=0.025, cbar_orientation='vertical', draw_labels=False)
    plotter.scatter(ds_reef.lon.values, ds_reef.lat.values, outdir+'departurezoom_%s_%s.png' % (runname, experiment), c=ds_reef.departure_year.values, marker_relsize=50, vmin=2020, vmax=2060, extent=[-40,-38,-18,-16], figsize=(5,5), cbar_orientation='vertical')
    plotter.scatter(mean_ds.lon.values, mean_ds.lat.values, outdir+'map2100_%s_%s.png' % (runname, experiment), c=mean_ds.tos.values[:,-1], cmap='jet', vmin=20, vmax=37, cbar_fraction=0.025, cbar_orientation='vertical', draw_labels=False)
    
    ds2020 = mean_ds.sel(time=slice(2020, 2030))
    mean2020 = np.mean(ds2020.tos.values, 1)
    ds2090 = mean_ds.sel(time=slice(2090, 2100))
    mean2090 = np.mean(ds2090.tos.values, 1)    
    anomaly = mean2090-mean2020
    if scenario == 'ssp126':
        plotter.scatter(mean_ds.lon.values, mean_ds.lat.values, outdir+'map2100-2020_%s_%s.png' % (runname, experiment), c=anomaly, cmap='jet', vmin=0, vmax=1, cbar_fraction=0.025, cbar_orientation='vertical', draw_labels=False)
    else:
        plotter.scatter(mean_ds.lon.values, mean_ds.lat.values, outdir+'map2100-2020_%s_%s.png' % (runname, experiment), c=anomaly, cmap='jet', vmin=0, vmax=4, cbar_fraction=0.025, cbar_orientation='vertical', draw_labels=False)
        
    plotter.hist(ds_reef.departure_year.values, 'Years', outdir+'hist_departure_%s_%s.png' % (runname, experiment), xlim=[1915,2100] )
    
    ds_2090 = ds_reef.sel(time=slice(2090, 2100))
    dhw_max = np.nanmax(ds_2090.dhw.values, 1)
    plotter.hist(dhw_max, 'MAX DHW, 2090-2100', outdir+'hist_maxdhw_%s_%s.png' % (runname, experiment) )

    mean_ds_reef = mean_ds.where(mean_ds.isreef==1, drop=True)
    mean_ds_2090 = mean_ds_reef.sel(time=slice(2090, 2100))
    sst_max = np.nanmax(mean_ds_2090.tos.values, 1)
    plotter.hist(sst_max, 'MAX SST, 2090-2100', outdir+'hist_maxsst_%s_%s.png' % (runname, experiment) )
    
    plotter.hist(np.array(mmm_vect), 'MMM', outdir+'hist_mmm_%s_%s.png' % (runname, experiment) )
        
    if experiment!='obs':
        used_models_file.close()
        
    #pdb.set_trace()