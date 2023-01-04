#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

conda activate geo_env

Get tas global mean model outputs using BCDP, regrid to 1x1°, and save on the local file system.

Must use same models for each SSP. All 4 or bust.

todo: MRI-ESM2-0 585 goes to 2300 and uses cftime.datetime objects for time, instead of datetime64.
So before concat with da_hist, would need to figure out how to convert to datetime64 or else get xr.concat to do this.
Otherwise, when attempting to write the netcdf, get ValueError: unable to infer dtype on variable 'time'; xarray cannot serialize arbitrary Python objects
    

@author: pkalmus
"""

#projectdir = '/0work/projects/'
projectdir = '/home/pkalmus/projects/'

import xarray as xr
import pandas as pd
import numpy as np
import pdb
import sys
sys.path.append(projectdir+'/pySSDF/src/')
import make_s
from sklearn.utils.extmath import cartesian
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import coral_plotter as plotter
import bcdp
import pprint
import intake
import traceback
import logging
logging.basicConfig(filename='run.log',level=logging.WARNING)

#import cftime

basedir = projectdir+'/coral/'
datadir = '/raid8/pkalmus/data/tas/'
pythondir = basedir+'src/python/'

# this is also bcdp.constants.DEFAULT_INTAKE_ESM_CAT; hardcode here in case it changes unexpectedly in bcdp
DEFAULT_INTAKE_ESM_CAT = 'https://raw.githubusercontent.com/NCAR/intake-esm-datastore/master/catalogs/pangeo-cmip6.json'


col = intake.open_esm_datastore(DEFAULT_INTAKE_ESM_CAT)
uni_dict = col.unique(["source_id"])
pprint.pprint(uni_dict, compact=True)

#query = dict(variable_id=["tos"], member_id=['r1i1p1f1'], experiment_id=["historical", "ssp126", "ssp245", "ssp370", "ssp585"]) # 12

# grid_label: gn, gr, gr1. table_id: Oday 245 Omon 526. member_id: 68!!
query = dict(variable_id=["tas"], table_id=['Amon'], grid_label=["gn"], experiment_id=["historical", "ssp126", "ssp245", "ssp370", "ssp585"]) # 18
#query = dict(table_id=['Amon'])
#query = dict(variable_id=["tos"], experiment_id=["historical", "ssp126", "ssp245", "ssp585"]) # 20
#query = dict(variable_id=["tos"], experiment_id=["historical", "ssp126", "ssp370", "ssp585"]) # 19
#query = dict(variable_id=["tos"], experiment_id=["historical", "ssp370", "ssp585"]) # 19
#query = dict(variable_id=["tos"], experiment_id=["historical", "ssp585"]) # 23

col_subset = col.search(require_all_on=["source_id"], **query)

# myset = col_subset.df.groupby("source_id")[["experiment_id"]].nunique()
# df = col_subset.df.groupby("source_id")[["experiment_id"]]
# print(len(myset))
# print(myset)

# myset = col_subset.df.groupby("variable_id")[["table_id"]].nunique()
# print(len(myset))
# print(myset)



sources = col_subset.unique()["source_id"]["values"]
print(sources)
#['ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CAMS-CSM1-0', 'CESM2', 'CESM2-WACCM', 'CNRM-CM6-1', 
#'CNRM-ESM2-1', 'CanESM5', 'CanESM5-CanOE', 'EC-Earth3-Veg', 'IPSL-CM6A-LR', 'MCM-UA-1-0', 'MIROC-ES2L', 
#'MIROC6', 'MPI-ESM1-2-HR', 
#'MRI-ESM2-0', 'UKESM1-0-LL']

# these have gr, but also all 5 experiments
# CESM2
# CESM2-WACCM
# CNRM-CM6-1
# MRI-ESM2-0


bounds = bcdp.Bounds(lon_bnds=(-179.5,180), lat_bnds=(-89.5,90), time_bnds=(1910,1920)) # note there is an "arange" bndry convention. I don't think time_bnds matters (for this at least).
grid_ds = bcdp.utils.grid_from_res((1.0,1.0), bounds)
# one dimensional lat and lon
grid_ds = grid_ds.assign_coords(lat=grid_ds.y)
grid_ds = grid_ds.assign_coords(lon=grid_ds.x)

scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
scenarios = ['ssp585']
#scenarios = ['ssp370']

# iterate over sources. get one member per source.
# iterate over SSPs
modellistfile = open(pythondir+'tas_models.txt', 'w')
for mysource in sources:
    print('starting %s' % (mysource))
    #pdb.set_trace()
    
    # if mysource=='MRI-ESM2-0':
    #     continue
    
    # find first member_id. choose the member_id with the most experiments. 
    query = dict(source_id=mysource, experiment_id=["historical", "ssp126", "ssp245", "ssp370", "ssp585"], variable_id='tas', table_id='Amon', grid_label='gn')
    col_subset = col.search(**query)
    members = col_subset.unique()["member_id"]["values"]
    print(members)
 
    myset = col_subset.df.groupby("member_id")[["experiment_id"]].nunique() #Pandas DataFrame
    print(myset.sort_values(by='experiment_id', ascending=False))
    mymember = myset.sort_values(by='experiment_id', ascending=False).head(1).index.values[0]
    print(mymember)

    # get historical
    query = dict(source_id=mysource, member_id=mymember, experiment_id='historical', variable_id='tas', table_id='Amon', grid_label='gn')
    ens = bcdp.load_intake_esm(query, catfile=DEFAULT_INTAKE_ESM_CAT)
    da_hist = ens.regrid(output_grid=grid_ds, backend='esmf', method='bilinear', ignore_degenerate=True).first
    
    for scenario in scenarios:
        print('starting %s %s' % (mysource, scenario))
        
        try:
            # get one model, stitch together with ssp
            query = dict(source_id=mysource, member_id=mymember, experiment_id=scenario, variable_id='tas', table_id='Amon', grid_label='gn')
            ens = bcdp.load_intake_esm(query, catfile=DEFAULT_INTAKE_ESM_CAT)
            #ds = ens.regrid(output_grid=grid_ds, backend='scipy', method='linear').first
            da = ens.regrid(output_grid=grid_ds, backend='esmf', method='bilinear', ignore_degenerate=True).first
            
            # cut short any that go to 2300, such as MRI-ESM2-0 585 and CanESM5
            #da = da.sel(time=slice(None, cftime.DatetimeProlepticGregorian(2101, 1, 1)))
            
            da_concat = xr.concat([da_hist, da], 'time')
            da_concat.name = 'tas'
            da_concat = da_concat.assign_coords(lon=((da_concat.lon + 360) % 360)).sortby('lon') # put into 0, 360.
            da_concat = da_concat.rename({'x': 'lon','y': 'lat'}) # must reassign to da_concat
            
            # take global mean
            da_mean = da_concat.mean(('lat', 'lon'))
            
            #pdb.set_trace()
            
            # save netcdf
            outfilename = datadir+mysource+'_'+scenario+'.nc'
            da_mean.to_netcdf(outfilename) 
            print(outfilename) 
        except Exception as e:
            print(e.message, e.args)
            logging.error('regrid_bcdp: source %s member %s error %s' % (mysource, mymember, e.message)) 
            traceback.print_exc()
            continue             
            
    # write line in file
    modellistfile.write('%s %s\n' % (mysource, mymember))

modellistfile.close()
print('done.')
#ds = ens.regrid(output_grid=grid_ds, backend='scipy', method='linear').normalize_times(assume_gregorian=True).bundle('CMIP6').first

# CAMS-CSM1-0 has an issue. *** ValueError: conflicting sizes for dimension 'time': length 3000 on 'tos' and length 3012 on 'time'
# master_time has length 3012. this dude has 3000.

