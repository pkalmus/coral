#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read in the projection files.
Read in the obs files.


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
do_weights = params.do_weights
do_obs = params.do_obs
finerun = params.finerun
includemur = params.include_mur
doval2015 = params.doval2015
useValFiles = params.useValFiles

lat_extent = 35
climatology_type = 'mur'

# define in and out dirs
indir = basedir+'data/tos/%s/weights/' % (runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
reefdir = basedir+'data/tos/%s/reef/' % (runname) 
outdir = basedir+'data/tos/%s/coarse/' % (runname)
#finedir = basedir+'data/tos/%s/fine/%s/' % (runname, finerun) # e.g. for scp: /raid8/pkalmus/data/coral/data/tos/nature/fine/flat/

if useValFiles:
    # requires a "val" dir in each finrun output dir
    valdir = basedir+'output/%s/fine/%s/val/' % (runname, finerun) # # /home/pkalmus/projects/coral/output/nature/fine/val/
else:
    valdir = basedir+'output/%s/fine/%s/' % (runname, finerun) #for bgl, where we didn't (yet) produce val-only runs

#murdir = basedir + 'data/mur/MUR_coral-pixels_2002-2020/'
murdir = valdir

sysStr = "mkdir -p %s" % outdir
subprocess.call(sysStr, shell=True)
use_best_only = False

#valtypes = ['trend_flat', 'bgl_flat', 'lagp', 'trend']

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


for scenario in scenarios:
    # Globe.130_150.min40_0.ssp370_flat(Validation).ncdf4.nc                                                 100%   30MB  12.9MB/s   00:02    
    # Globe.130_150.min40_0.ssp370_weighted(Validation).Trend.ncdf4.nc                                       100%   16MB  15.0MB/s   00:01    
    # Globe.130_150.min40_0.ssp370_weighted(Validation).ncdf4.nc  
    if doval2015:
        if finerun=='trend':
             lagp_files = glob.glob(valdir+'Globe*%s_weighted*Trend*2015*.nc' % (scenario)) # make sure you don't glob output files from THIS code
        elif finerun=='flat':
            lagp_files = glob.glob(valdir+'Globe*%s_flat*2015*.nc' % (scenario)) # make sure you don't glob output files from THIS code 
        elif finerun=='lagp':
            lagp_files = glob.glob(valdir+'Globe*%s_weighted(Validation)*2015*.ncdf4.nc' % (scenario)) # make sure you don't glob output files from THIS code  
        print('++++++starting 2015-2017' + finerun)            
    else:
        if finerun=='trend':
             lagp_files = glob.glob(valdir+'Globe*%s_weighted*Trend*.nc' % (scenario)) # make sure you don't glob output files from THIS code
        elif finerun=='flat':
            lagp_files = glob.glob(valdir+'Globe*%s_flat*.nc' % (scenario)) # make sure you don't glob output files from THIS code 
        elif finerun=='lagp':
            lagp_files = glob.glob(valdir+'Globe*%s_weighted(Validation).ncdf4.nc' % (scenario)) # make sure you don't glob output files from THIS code  
        elif finerun=='bgl_flat':
            lagp_files = glob.glob(valdir+'Globe*%s_uniform.ncdf4.nc' % (scenario)) # make sure you don't glob output files from THIS code  
        elif finerun=='trend_flat':
            lagp_files = glob.glob(valdir+'Globe*%s_uniform.Trend.ncdf4.nc' % (scenario)) # make sure you don't glob output files from THIS code 
        elif finerun=='bgl_onepermodel':
            lagp_files = glob.glob(valdir+'Globe*%s*.ncdf4.nc' % (scenario)) # make sure you don't glob output files from THIS code  
        else:
            print('UNKNOWN FINE RUN')
        print('++++++starting ' + finerun)
        
    
    if len(lagp_files) == 0:
        sys.exit('cold not find downscaled files in: %s' % (valdir))

    # load projections and obs    
    ds_list = []
    mur_list = []
    for lagp_file in lagp_files:
        dslagp = xr.open_dataset(lagp_file)
        #dslagp = dslagp.dropna('index', subset=['temp'])            
        # if 'temp_sd' in dslagp.keys(): # errors='ignore' apparently came in a later version
        #      dslagp = dslagp.drop('temp_sd') # some of the 2015 ds have temp_sd, some don't; ignore errors if it doesn't
        ds_list.append(dslagp)   
        
        murfile = murdir+'Globe'+lagp_file[lagp_file.find('Globe')+5:lagp_file.find('ssp')-1]+'_MUR.ncdf4.nc'
        dsmur = xr.open_dataset(murfile)
        dsmur = dsmur.dropna('index', subset=['temp']) 
        mur_list.append(dsmur)
        print(lagp_file)
        print(murfile)
    lagp_ds = xr.concat(ds_list, 'index') # index will no longer be a monotonic index        
    mur_ds = xr.concat(mur_list, 'index') # index will no longer be a monotonic index 
    if doval2015:
        lagp_ds = lagp_ds.sel(time=slice(2015, 2018)) 
        mur_ds = mur_ds.sel(time=slice(2015, 2018))        
    else:
        lagp_ds = lagp_ds.sel(time=slice(2018, 2021)) 
        mur_ds = mur_ds.sel(time=slice(2018, 2021))
            
    #plotter.plot_data_array(lagp_ds.sel(time=2020.01).temp, finedir+'fine2020.png')
    #plotter.plot_data_array(lagp_ds.sel(time=2099.01).temp, finedir+'fine2099.png')
        
    lons = lagp_ds.lon.values
    lats = lagp_ds.lat.values
    sst = lagp_ds.temp.values
    #sd = lagp_ds.temp_sd.values # standard deviation
    time = lagp_ds.time.values
    time = np.floor(time) + ((time - np.floor(time))*30*100 - 15)/365. 
    index = lagp_ds.index.values
    
    sstmur = mur_ds.temp.values
    latmur = mur_ds.lat.values
    
    print(sst.shape)
    print(sstmur.shape)
    if len(sstmur) != len(sst): 
        # because of new NaN values in dslagp, you must also do unto dsmur.
        #this protects against the datasets getting out of sync.
        sys.exit('sst and sstmur have different lengths, which will invalidate validation.')
        
    #pdb.set_trace()
    
    # this is for testing.
    #lagp_ds = lagp_ds.where(lagp_ds.index < 1500, drop=True) 
    
    latravel = []
    lonravel = []
    mmmravel = []
    lagpmeanvect = [] # for val
    murmeanvect = []
    lagpmaxvect = [] # for val
    murmaxvect = []
    
    null_data = np.empty((lagp_ds.dims['index']))
    null_data[:] = np.nan
    #maxdhw2020_da = xr.DataArray(np.copy(null_data), coords=([lagp_ds.coords['index']]), dims=('index'))  
    for (x, ignore) in enumerate(lagp_ds.index.values):
        lat = lats[x]
        if np.isnan(lat):
            continue
        lon = lons[x]
        latravel.append(lat)
        lonravel.append(lon)
        
        tos = sst[x,:]
        toslat = latmur[x]
        if toslat != lat:
            print('lat mismatch')
            pdb.set_trace()
        
        # this works if dslagp and dsmur have same length!
        tosmur = sstmur[x,:]            

        #pdb.set_trace()
        
        # add in MUR uncertainty
        # tos_sd = np.sqrt(np.power(tos_sd,2) + np.power(0.25,2))
        
        mmm = get_mmm_mur([lon, lat])            
        dhwlagp = get_dhw(mmm, tos)
        dhwmur = get_dhw(mmm, tosmur)   
        
        # for val
        #ind2020 = np.where((time > 2019.5) & (time < 2021))[0] #1.5 years to account for austral winter
        #maxdhw2020val = np.nanmax(dhw[ind2020])    
        
        # select the maximum for each year (every 12 data points)
        lagpmax = np.nanmax(dhwlagp.reshape(-1, 12), axis=1)
        murmax = np.nanmax(dhwmur.reshape(-1, 12), axis=1)
        
        # take the mean of the annual maxima
        lagpmean = np.nanmean(lagpmax)
        lagpmeanvect.append(lagpmean)
        murmean = np.nanmean(murmax)
        murmeanvect.append(murmean)
        
        # now do it for the max over entire period
        lagpmaxvect.append(np.nanmax(dhwlagp))
        murmaxvect.append(np.nanmax(dhwmur))


        if np.mod(x, 1000)==0:
            print(str(x) + '/' + str(len(index)))

        # if x > 2000:
        #     break
    
    dhwmean_diff = np.array(lagpmeanvect) - np.array(murmeanvect) 
    dhwmax_diff = np.array(lagpmaxvect) - np.array(murmaxvect) 

    if doval2015:
        plotter.hist(dhwmean_diff, 'Projected mean maxDHW - MUR mean maxDHW', outdir+'hist_meanannualmax_dhwdiff_%s_%s_%s2015-2017.png' % (runname, scenario, finerun), bins=100 )
        plotter.hist(dhwmax_diff, 'Projected maxDHW - MUR maxDHW', outdir+'hist_max_dhwdiff_%s_%s_%s2015-2017.png' % (runname, scenario, finerun), bins=100 )

    else:
        plotter.hist(dhwmean_diff, 'Projected mean maxDHW - MUR mean maxDHW', outdir+'hist_meanannualmax_dhwdiff_%s_%s_%s2.png' % (runname, scenario, finerun), bins=100 )
        plotter.hist(dhwmax_diff, 'Projected maxDHW - MUR maxDHW', outdir+'hist_max_dhwdiff_%s_%s_%s2.png' % (runname, scenario, finerun), bins=100 )
    #pdb.set_trace()

