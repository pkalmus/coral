#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

conda activate geo_env

reef_locations: get reef-containing points (don't need to run every time)
regrid_bcdp: Get model outputs using BCDP, regrid to 1 degree, time homogenize and save on the local file system for regrid_reef.py to digest.
reef_locations_coarse: get coarse reef-containing points (don't need to run every time)

regrid_reef: Find reef and neighbor locations (using previously made reef_location_coarse.csv and nn_location_coarse.csv)
             Repackage into 1-D "ravel" .nc dataset with "index".
               This .nc dataset also specifies which coarse cells are reefs.

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
srcdir = params.srcdir
modellist_filename = params.modellist_filename
scenarios = params.scenarios
dryRun = params.dryRun
oneMember = params.oneMember

# special import from pySSDF, needs to come after "projectdir" defined
sys.path.append(projectdir+'/pySSDF/src/')
import make_s

# define in and out dirs
indir = basedir+'data/tos/coral2/%s/raw/' % (tosrunname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
outdir = basedir+'data/tos/coral2/%s/reef/' % (tosrunname)

reefmodellist_filename = modellist_filename.replace('.txt', '_reef.txt')
reefmodellistfile = open(reefmodellist_filename, 'w')

sysStr = "mkdir -p %s" % outdir
subprocess.call(sysStr, shell=True)

reef_coarse_file = '/raid8/pkalmus/data/coral/data/location/coarse/reef_location_coarse.csv'    # 1294 [lon lat]
reef_coarse_cells = np.loadtxt(reef_coarse_file)
nn_coarse_file = '/raid8/pkalmus/data/coral/data/location/coarse/nn_location_coarse.csv'        # 1790 [lon lat]
nn_coarse_cells = np.loadtxt(nn_coarse_file)

both_grid = np.vstack([reef_coarse_cells, nn_coarse_cells])
index = np.arange(both_grid.shape[0])
nreef = len(reef_coarse_cells)

with open(modellist_filename) as f:
    models = f.read().splitlines()
  
    
# UNNECESSARY UGLINESS
# OK, this is to fix an issue where some of the models use datetime64. we want them all to use cftime.
# Alex probably needs to fix bcdp homogenize.
timefilename = indir+'BCC-CSM2-MR_r1i1p1f1_gn_ssp126.nc' # days since 1850-01-16 00:00:00.000000..., 0, 3012, 360_day
timeds = xr.open_dataset(timefilename, decode_times=True) 
timeds = timeds.sel(time=slice(None, '2099') ) # length 3000, goes up to December 2099
master_time = timeds.time # datetimeds
# UNNECESSARY UGLINESS


def regrid_and_save(latsravel, lonsravel, mytime, tos_all, outfilename):
    '''
    check to see if each point is on reef_coarse_cells   
    save the netcdf file
    
    the ordering of locations relative to the index must be consistent across models, HadISST, and climatology.
    '''
    tosravel = tos_all.reshape(tos_all.shape[0], -1) 
    sstos = np.array(list(zip(*[lonsravel,latsravel])) ) #need list() in Python 3

    # find all the reef points
    ind_plus_negative_ones = make_s.ismember(sstos,reef_coarse_cells.copy()) # the -1 values indicate "not a member"
    ind = np.where(ind_plus_negative_ones > -1)[0]
    sstos_reef = sstos[ind]
    tos_reef = tosravel[:,ind]
    tos_reef_swap = np.swapaxes(tos_reef, 0, 1) 
    lons_reef = sstos_reef[:,0]
    lats_reef = sstos_reef[:,1]

    # find all the nn points
    ind_plus_negative_ones = make_s.ismember(sstos, nn_coarse_cells.copy()) # the -1 values indicate "not a member"
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
    nreefs = len(index)
    nnans = len(nanind)
    print('number of reefs & nn:', nreefs)
    print('number of nan:', nnans)
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
    return nreefs, nnans

# HadISST
filename = '/raid8/pkalmus/data/coral/data/tos/coral2/main/obs/HadISST_sst_202110.nc'
print(filename)
outfilename = outdir+'HadISST_sst_202110_reef.nc'  
modds = xr.open_dataset(filename, decode_times=True) 
modds = modds.drop('time_bnds')
modds = modds.assign_coords(latitude=modds.latitude).sortby('latitude') # put into ascending
modds = modds.assign_coords(longitude=((modds.longitude + 360) % 360)).sortby('longitude') # put into 0, 360.
modds = modds.sel(time=slice('1970', '2020'))
modds = modds.assign_coords(time=master_time.copy().sel(time=slice('1970', '2020'))) 
print(modds)
t = modds.time.values
lats = modds.latitude.values
lons = modds.longitude.values
tos_all = modds.sst.values # time, latitude, longitude
latsravel = np.tile(lats, (len(lons), 1)).transpose().ravel()
lonsravel = np.tile(lons, (len(lats), 1)).ravel() 
nreefs, nnans = regrid_and_save(latsravel, lonsravel, t, tos_all, outfilename)
# write line in file
reefmodellistfile.write('%s %i %i\n' % ('HadISST_sst_202110', nreefs, nnans))
reefmodellistfile.flush()

for (modelind, model) in enumerate(models):
    print('\n\n***********'+str(modelind+1) + '   ' + model)
      
    # if model != 'IITM-ESM r1i1p1f1 gn': 
    #     continue  
    
    # if modelind < 30:
    #     continue
    
    for scenario in scenarios:
        # e.g. CanESM5_r12i1p1f1_ssp126.nc
        modelstr = model.replace(' ', '_')
        ssp_modelstr = modelstr+'_'+scenario
        filename = indir+ssp_modelstr+'.nc'
        outfilename = outdir+ssp_modelstr+'_reef.nc'

        modds = xr.open_dataset(filename, decode_times=True) 
        modds = modds.sel(time=slice('1970', '2099'))
        #print(modds)        
        modds = modds.assign_coords(time=master_time.copy()) 
        t = modds.time.values
        lats = modds.lat.values
        lons = modds.lon.values
        tos_all = modds.tos.values # time, latitude, longitude
        latsravel = np.tile(lats, (len(lons), 1)).transpose().ravel()
        lonsravel = np.tile(lons, (len(lats), 1)).ravel() 
        nreefs, nnans = regrid_and_save(latsravel, lonsravel, t, tos_all, outfilename)
        
    # write line in file
    reefmodellistfile.write('%s %i %i\n' % (modelstr, nreefs, nnans))
    reefmodellistfile.flush()

reefmodellistfile.close()
print(reefmodellistfile)  
        
        



