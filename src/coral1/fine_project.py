#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
python fine_project.py params

To run on single model, pass in the appropriate integer argument.
To run on HadISST, pass in -1.

conda activate geo_env

regrid_bcdp
coral_locations
regrid_reef
weighting
coarse_project TODO: redo climatology in coral_location; why do two of the models fail to get all 4000+ NN locations?
validating code

Read in fine file
Make a map of 2020 and 2100
Read in climatology (need it on same grid) 
Calculate departure
Make map


@author: pkalmus
"""
import xarray as xr
import pdb
import numpy as np
import glob
import sys
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


# read in climatology. want MMM for each location. 
if climatology_type  ==  'mur':
    clim_file = basedir + 'data/climatology/mur_mmm.nc'
    clim_ds  = xr.open_dataset(clim_file)
    clim_ds = clim_ds.assign_coords(lon=((clim_ds.lon + 360) % 360)).sortby('lon')
else:
    climatology_start_year = 1985 #earliest possible is 1985
    climatology_end_year = 2008 #latest possible is 2012
    
    clim_file = basedir + 'data/climatology/noaa_crw_thermal_history_climatology.nc'
    clim_ds  = xr.open_dataset(clim_file) # lon: -180 to 180. lat: 34.14 to -35.27 (so in reverse order... be careful.)
    clim_ds = clim_ds.assign_coords(lon=((clim_ds.lon + 360) % 360)).sortby('lon')
    clim_ds = clim_ds.assign_coords(lat=clim_ds.lat).sortby('lat') # fixed
    clim_ds = clim_ds.sel(years=slice(climatology_start_year, climatology_end_year))

def get_mmm_clim_coral_only(point):
    '''
    Tried commenting this function out and just returning an integer; it's not the slow part.
    '''
    mylon = point[0]
    mylat = point[1]
    myclim = clim_ds.sel(lat=mylat, lon=mylon, method='nearest')
    
    if (myclim.reef_mask.values==0) or (np.abs(mylon-myclim.lon.values) > 0.03) or (np.abs(mylat-myclim.lat.values) > 0.03):
        return np.nan
    else:
        thisind = myclim.mmm_month.values - 1 # goes from 0 to 12, but I think 0 indicates fill val??
        mmm = myclim.clim_monthly.values[int(thisind)]
        return mmm

def get_mmm_mur(point):
    '''
    '''
    mylon = point[0]
    mylat = point[1]
    #pdb.set_trace()
    myclim = clim_ds.sel(lat=mylat, lon=mylon, method='nearest')
    return myclim.mmm.values
    

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

def get_dep(dhw, mmm, tos):
    # we now have monthly time series of 3-monthly DHW.
    # we want a series of binaries that state "within N years around this month ("span") there's an excursion"
    N = 12 * return_year
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
    departure_time = time[last_ind]            
        
    if np.isnan(mmm) or (np.isnan(tos[1])): # bad failure mode: nans get converted to 2100
        departure_time = np.nan
    return departure_time

for return_year in return_years:
    for scenario in scenarios:
        experiment = '%s_%s_%iyrs_DHW%1.1f' % (finerun, scenario, return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
        if do_add_variability:
            experiment = 'SPIKE%im%if_%s_%s_%iyrs_DHW%1.1f' % (spike_size, spike_freq, finerun, scenario, return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
            

        if finerun=='bgl_flat':
            lagp_files = glob.glob(finedir+'Globe*%s_uniform*.nc' % (scenario))
        elif finerun=='lagp_beo':
            lagp_files = glob.glob(finedir+'Globe*%s_weighted*.nc' % (scenario)) # make sure you don't glob output files from THIS code
        elif finerun=='trend_beo':
            lagp_files = glob.glob(finedir+'Globe*%s_uniform.Trend*.nc' % (scenario)) # make sure you don't glob output files from THIS code 
        elif finerun=='lagp_flat':
            lagp_files = glob.glob(finedir+'Globe*%s_flat*.nc' % (scenario)) # make sure you don't glob output files from THIS code 
        elif finerun=='trend_flat':
            lagp_files = glob.glob(finedir+'Globe*%s_uniform.Trend*.nc' % (scenario)) # make sure you don't glob output files from THIS code 
        elif finerun=='bgl_onepermodel':
            lagp_files = glob.glob(finedir+'Globe*%s*.nc' % (scenario)) # make sure you don't glob output files from THIS code 

        else:
            raise(ValueError('unknown finerun: %s' % finerun))               

        print('++++++starting ' + experiment)

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
            #pdb.set_trace()
            
        lagp_ds = xr.concat(ds_list, 'index') # index will no longer be a monotonic index
        
        # added this for onepermodel. extra NaNs here, probably due to the GCM "less than ten" rule but did not give Ayesha a new location file.
        lagp_ds = lagp_ds.dropna('index', subset=['lat'])  
        
        # if finerun=='weightbc':
        #     lagp_ds = lagp_ds.where(lagp_ds.lat==slice(-2, 1), drop=True)
        
        #plotter.plot_data_array(lagp_ds.sel(time=2020.01).temp, finedir+'fine2020.png')
        #plotter.plot_data_array(lagp_ds.sel(time=2099.01).temp, finedir+'fine2099.png')
            
        lons = lagp_ds.lon.values
        lats = lagp_ds.lat.values
        #pdb.set_trace()
        sst = lagp_ds.temp.values
        if include_downscale_uncert:
            sd = lagp_ds.temp_sd.values # standard deviation
        else:
            sd = np.zeros(sst.shape)
        time = lagp_ds.time.values
        time = np.floor(time) + ((time - np.floor(time))*30*100 - 15)/365. # to fix Ayesha's weird non-decimal time convention
        index = lagp_ds.index.values
        index = np.arange(len(index)) # fixed it
        lagp_ds = lagp_ds.assign_coords(index=(index))
        
        # this is for testing.
        #lagp_ds = lagp_ds.where(lagp_ds.index < 1500, drop=True) 
        
        latravel = []
        lonravel = []
        depravel = []
        dephighravel = []
        deplowravel = []
        mmmravel = []
        did2020 = False
        did2090 = False
        maxdhw2020 = [] # for val
        nanctr = 0
        
        null_data = np.empty((lagp_ds.dims['index']))
        null_data[:] = np.nan
        dep_da = xr.DataArray(null_data, coords=([lagp_ds.coords['index']]), dims=('index'))    
        dephigh_da = xr.DataArray(np.copy(null_data), coords=([lagp_ds.coords['index']]), dims=('index')) 
        deplow_da = xr.DataArray(np.copy(null_data), coords=([lagp_ds.coords['index']]), dims=('index')) 
        maxdhw2020_da = xr.DataArray(np.copy(null_data), coords=([lagp_ds.coords['index']]), dims=('index'))  
        for (x, ignore) in enumerate(lagp_ds.index.values):
            lat = lats[x]
            lon = lons[x]
                        
            
            if np.isnan(lat) or np.isnan(lon):
                print('NAN in LAT AND/OR LON, continuing...')
                nanctr +=1
                continue
            
            latravel.append(lat)
            lonravel.append(lon)
            tos = sst[x,:]
            #tos_orig = np.copy(tos)
            
            if do_add_variability:
                # every five years plus or minus a random interaval, randomly add or subtract a degree
                tos_ind = np.where(time > 2020)[0][0] # only inject into non-MUR data
                while tos_ind < len(tos) - spike_freq*12:
                    tos_ind = tos_ind + (spike_freq*12 + np.random.randint(-24,25)) # every 5*12 months, plus or minus 24 months
                    spike = spike_size*np.random.randn() 
                    halfduration = np.random.randint(1,4) # months
                    tos[tos_ind-halfduration:tos_ind+halfduration] = tos[tos_ind-halfduration:tos_ind+halfduration]+spike
            
            tos_sd = sd[x,:]
            # add in MUR uncertainty, 0.25 degrees
            tos_sd = np.sqrt(np.power(tos_sd,2) + np.power(0.25,2))
            if climatology_type =='mur':
                mmm = get_mmm_mur([lon, lat])            
            else:
                mmm = get_mmm_clim_coral_only([lon, lat])
    
            #mmm = mmm+1
            dhw = get_dhw(mmm, tos)
            departure_time = get_dep(dhw, mmm, tos)
            depravel.append(departure_time)

            # for val
            ind2020 = np.where((time > 2019.5) & (time < 2021))[0] #1.5 years to account for austral winter
            maxdhw2020val = np.nanmax(dhw[ind2020])
            maxdhw2020.append(maxdhw2020val)
    
            dhw = get_dhw(mmm, tos-tos_sd)
            departure_time = get_dep(dhw, mmm, tos-tos_sd)
            dephighravel.append(departure_time)
    
            dhw = get_dhw(mmm, tos+tos_sd)
            departure_time = get_dep(dhw, mmm, tos+tos_sd)
            deplowravel.append(departure_time)        
            
            mmmravel.append(mmm)
    
            #print(departure_time)
            # if (departure_time > 2099):
            #     if not did2090:
            #         plotter.plot(time, dhw, 'dhw2099.png')
            #         did2090 = True            
            # if (departure_time <2021):
            #     if not did2020:
            #         plotter.plot(time, dhw, 'dhw2020.png')
            #         did2020 = True
    
            if np.mod(x, 1000)==0:
                print(str(x) + '/' + str(len(index)) + '   nan lat/lon: ' + str(nanctr))
        dep_da.values = depravel # this is SUPER SLOW if done inside the loop w .loc[dict(index=x)]
        dephigh_da.values = dephighravel
        deplow_da.values = deplowravel   
        maxdhw2020_da.values = maxdhw2020
        print(len(np.where(np.isnan(depravel))[0]))
    
        # write .nc file
        ds = xr.Dataset({'departure_year':dep_da, 'departure_year_high':dephigh_da, 'departure_year_low':deplow_da, 'maxdhw2020':maxdhw2020_da,'lon':lagp_ds.lon, 'lat':lagp_ds.lat} )
        if includemur:
            outfilename = finedir+'finemur_%s_%s.nc' % (runname, experiment)
        else:
            outfilename = finedir+'fine_%s_%s.nc' % (runname, experiment)
        ds.to_netcdf(outfilename) 
        print(outfilename)    
          
        # #plotter.scatter(lonravel, latravel, finedir+'departure_%s_%s_zoom.png' % (runname, experiment), c=depravel, vmin=2020, vmax=2080, extent=[-40,-38,-18,-16], figsize=(5,5), cbar_orientation='vertical')
        # #plotter.scatter(lonravel, latravel, finedir+'departure_%s_%s_zoomzoom.png' % (runname, experiment), c=depravel, vmin=2020, vmax=2080, extent=[-39.3,-38.8,-18,-17.5], figsize=(5,5), cbar_orientation='vertical')
        # #plotter.scatter(lonravel, latravel, finedir+'departure_%s_%s_zoom2.png' % (runname, experiment), c=depravel, vmin=2020, vmax=2080, extent=[-37,-35,-11,-8], figsize=(5,5), cbar_orientation='vertical')
        # # plotter.scatter(lagp_ds.lon.values, lagp_ds.lat.values, finedir+'map2020_%s_%s.png' % (runname, experiment), c=mean2020, vmin=25, vmax=30, extent=[-40,-38,-18,-16], figsize=(5,5), cmap='jet', cbar_fraction=0.025, cbar_orientation='vertical')
        # # plotter.scatter(lagp_ds.lon.values, lagp_ds.lat.values, finedir+'map2050_%s_%s.png' % (runname, experiment), c=mean2050, vmin=25, vmax=30, extent=[-40,-38,-18,-16], figsize=(5,5), cmap='jet', cbar_fraction=0.025, cbar_orientation='vertical')
        # # plotter.scatter(lonravel, latravel, finedir+'mmm_%s_%s_zoom.png' % (runname, experiment), c=mmmravel, vmin=25, vmax=30, extent=[-40,-38,-18,-16], figsize=(5,5), cmap='jet', cbar_fraction=0.025, cbar_orientation='vertical')
    
        # plotter.scatter(lonravel, latravel, finedir+'departure_%s_%s.png' % (runname, experiment), c=depravel, vmin=2020, vmax=2080, figsize=(5,5), cbar_orientation='vertical')
        # ds2020 = lagp_ds.sel(time=slice(2020, 2021))
        # mean2020 = np.mean(ds2020.temp.values, 1)
        # ds2050 = lagp_ds.sel(time=slice(2050, 2051))
        # mean2050 = np.mean(ds2050.temp.values, 1)  
        # plotter.scatter(lagp_ds.lon.values, lagp_ds.lat.values, finedir+'map2020_%s_%s.png' % (runname, experiment), c=mean2020, vmin=25, vmax=30, figsize=(5,5), cmap='jet', cbar_fraction=0.025, cbar_orientation='vertical')
        # plotter.scatter(lagp_ds.lon.values, lagp_ds.lat.values, finedir+'map2050_%s_%s.png' % (runname, experiment), c=mean2050, vmin=25, vmax=30, figsize=(5,5), cmap='jet', cbar_fraction=0.025, cbar_orientation='vertical')
        # plotter.scatter(lonravel, latravel, finedir+'mmm_%s_%s_zoom.png' % (runname, experiment), c=mmmravel, vmin=25, vmax=30, figsize=(5,5), cmap='jet', cbar_fraction=0.025, cbar_orientation='vertical')
    
    
        #pdb.set_trace()