#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

conda activate geo_env

reef_locations: get reef-containing points (don't need to run every time)
regrid_bcdp: Get model outputs using BCDP, regrid to 1x1°, time homogenize and save on the local file system for regrid_reef.py to digest.

regrid_reef: Find reef and neighbor locations (using previously made reef_grid_1x1.txt and nn_grid_1x1.txt)
             Repackage into 1-D "ravel" dataset with "index".


Filter out non-reef locations and save new netcdf files in "ravel" form.
Only for 1x1 regular gridded data (such as the HadISST, and some of the CMIP6 gr models.)
Read in model or obs.

@author: pkalmus
"""

import xarray as xr
import numpy as np
import pdb
import sys
import subprocess
import coral_plotter as plotter

# read in user params
import importlib
paramfile = sys.argv[1]   
params = importlib.import_module(paramfile)
projectdir = params.projectdir
tosrunname = params.tosrunname
basedir = params.basedir
pythondir = params.pythondir
modellist_filename = params.modellist_filename
scenarios = params.scenarios
dryRun = params.dryRun
oneMember = params.oneMember
year_cutoff = params.year_cutoff
geolimit = params.geolimit

#useHadISST = False # if true, include HadISST in the analysis

# special import from pySSDF, needs to come after "projectdir" defined
sys.path.append(projectdir+'/pySSDF/src/')
import make_s

# define in and out dirs
indir = basedir+'data/tos/%s/raw/' % (tosrunname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
outdir = basedir+'data/tos/%s/reef/' % (tosrunname)

sysStr = "mkdir -p %s" % outdir
subprocess.call(sysStr, shell=True)

# load up reef points and neighbor points
locationsdir = '/raid8/pkalmus/data/coral/data/location/14_001_WCMC008_CoralReefs2018_v4_1/01_Data/'
prediction_grid = np.loadtxt(locationsdir+'locationsWCMCv41.csv')





prediction_grid = np.loadtxt('reef_grid_1x1.txt') # col1 lon,    col2 lat
nn_grid = np.loadtxt('nn_grid_1x1.txt')



both_grid = np.vstack([prediction_grid, nn_grid])
index = np.arange(both_grid.shape[0])
nreef = len(prediction_grid)


with open(modellist_filename) as f:
    models = f.read().splitlines()

# append the obs. file to the models list
# HadISST
filename = '/raid8/pkalmus/data/coral/data/tos/raw/HadISST_sst.nc'
print(filename)
outfilename = outdir+'HadISST_sst_reef.nc'    


def regrid_and_save(latsravel, lonsravel, mytime, tos_all, outfilename, do_latitudinal_interp=False):
        '''
        check to see if each point is on prediction_grid   
        save the netcdf file
        
        the ordering of locations relative to the index must be consistent across models, HadISST, and climatology.
        '''
        tosravel = tos_all.reshape(tos_all.shape[0], -1) 
        sstos = np.array(list(zip(*[lonsravel,latsravel])) ) #need list() in Python 3
    
        # find all the reef points
        ind_plus_negative_ones = make_s.ismember(sstos,prediction_grid.copy()) # the -1 values indicate "not a member"
        ind = np.where(ind_plus_negative_ones > -1)[0]
        sstos_reef = sstos[ind]
        tos_reef = tosravel[:,ind]
        tos_reef_swap = np.swapaxes(tos_reef, 0, 1) 
        lons_reef = sstos_reef[:,0]
        lats_reef = sstos_reef[:,1]
    
        # find all the nn points
        ind_plus_negative_ones = make_s.ismember(sstos, nn_grid.copy()) # the -1 values indicate "not a member"
        ind = np.where(ind_plus_negative_ones > -1)[0]
        sstos_nn = sstos[ind]
        tos_nn = tosravel[:,ind]
        tos_nn_swap = np.swapaxes(tos_nn, 0, 1) 
        lons_nn = sstos_nn[:,0]
        lats_nn = sstos_nn[:,1]
    
        tos_all = np.vstack([tos_reef_swap, tos_nn_swap]) # 4836x2220
        lons_all = np.hstack([lons_reef, lons_nn])        # 4836
        lats_all = np.hstack([lats_reef, lats_nn])
        
        # Here are the NaNs that have reef or are nearest neighbors.
        #plotter.scatter(lons_all, lats_all, outfilename.replace('.nc', '_pixels.png'))
        nanind = np.where(np.isnan(tos_all[:,0]))[0] #156 points, all along coasts.
        print('number of reefs & nn:', len(index))
        print('number of nan:', len(nanind))
        #plotter.scatter(lons_all[nanind], lats_all[nanind], outfilename.replace('.nc', '_nanpixels.png'))        
        
        myyears = np.array([item.year for item in mytime])
        mymonths = np.array([item.month for item in mytime])
        myyearfrac = myyears + (mymonths-0.5)/12.0
        
        isreef = np.zeros(len(index))
        isreef[0:nreef]=1
        lon_da = xr.DataArray(lons_all, coords=[index], dims=['index']) # why do these need square brackets?
        lat_da = xr.DataArray(lats_all, coords=[index], dims=['index'])
        year_da = xr.DataArray(myyearfrac, coords=[myyearfrac], dims=['time'])
        da = xr.DataArray(tos_all, coords=(index, myyearfrac), dims=('index', 'time')) # why don't these need square brackets?
        isreef_da = xr.DataArray(isreef, coords=[index], dims=(['index']))
        ds = xr.Dataset({'tos':da, 'isreef':isreef_da, 'year':year_da, 'lon':lon_da, 'lat':lat_da})
        
        ds.to_netcdf(outfilename) 
        print(outfilename)



for (modelind, model) in enumerate(models):
    print('\n\n***********'+str(modelind+1) + '   ' + model)
    
    # if model == 'CNRM-ESM2-1 r4i1p1f2 gn': # hada non-monotonic time!
    #     print('DELETE MODEL FROM %s: %s' % (listfilename, model))
    #     continue
    
    # if model == 'IITM-ESM r1i1p1f1 gn': # hours since 1900-01-16 12:00:00.000000
    #     print('DELETE MODEL FROM %s: %s' % (listfilename, model))
    #     continue
      
    # if model == 'KACE-1-0-G r1i1p1f1 gr': # ??? 
    #     print('DELETE MODEL FROM %s: %s' % (listfilename, model))
    #     continue
        
    
    for scenario in scenarios:
        # e.g. CanESM5_r12i1p1f1_ssp126.nc
        modelstr = model.replace(' ', '_')
        modelstr = modelstr+'_'+scenario
        filename = indir+modelstr+'.nc'
        #print(filename)
        outfilename = outdir+modelstr+'_reef.nc'
        #print(outfilename)
        
        # if os.path.exists(outfilename):
        #     print('%s already exists, continuing.' % (outfilename))
        #     continue
                        
        modds_all = xr.open_dataset(filename, decode_times=False) 

    
        modds_all = xr.open_dataset(filename, decode_times=True) 
        modds_all = modds_all.sel(time=slice(None, '2099')) # length 3000, goes up to December 2099
        try:

            years = np.array([item.year for item in modds_all.time.values])
            months = np.array([item.month for item in modds_all.time.values])
            days = np.array([item.day for item in modds_all.time.values])
            x2 = years + (months-1)/12. + days/365.
            y2 = modds_all.tos.values[:,40,40]
            #plotter.plot(x, y, outfilename.replace('.nc', '_testtime.png'), x2=x2, y2=y2, xlabel='time', ylabel='tos', title=None,  xlim=[1850,1855], ylim=None)

        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            #raise
            pdb.set_trace()
            continue

        modds = modds_all.sel(lat=slice(geolimit[0][0], geolimit[0][1]), time=slice(str(year_cutoff),None) )
        t = modds.time.values
        lats = modds.lat.values
        lons = modds.lon.values
        tos_all = modds.tos.values # time, latitude, longitude
        latsravel = np.tile(lats, (len(lons), 1)).transpose().ravel()
        lonsravel = np.tile(lons, (len(lats), 1)).ravel() 

        regrid_and_save(latsravel, lonsravel, t, tos_all, outfilename)
        
        
        
        
        
        
        
        
        
