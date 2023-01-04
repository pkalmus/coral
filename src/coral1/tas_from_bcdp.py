#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 18:45:29 2021


conda activate geo_env

Based on regrid_bcdp.py

This is for plotting against GMST anomaly in Figure 2.

Get model outputs (tas, one from each model group) using BCDP.

Regrid (because we will need to adjust for area difference with latitude)


Will take model mean in a different file (write_gmst_from_tas.py)


@author: pkalmus
"""


import xarray as xr
import numpy as np
import pdb
import sys
# sys.path.append(projectdir+'/pySSDF/src/')
# import make_s
import bcdp
import intake
import traceback
import subprocess
import os.path
import logging
logging.basicConfig(filename='run.log',level=logging.WARNING)

#import cftime

# read in user params
import importlib
paramfile = sys.argv[1]   
params = importlib.import_module(paramfile)
runname = params.runname
projectdir = '/home/pkalmus/projects/'
basedir = projectdir+'/coral/'
listfilename = params.listfilename
scenarios = params.scenarios
dryRun = params.dryRun
oneMember = True # only get the first member from a source (group)
grids = params.grids

outdir = basedir+'data/tas/%s/' % (runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
listfilename = 'models_gmst.txt'

experiments = scenarios.copy()
experiments.insert(0, 'historical') 

sysStr = "mkdir -p %s" % (outdir)
subprocess.call(sysStr, shell=True)

# this is also bcdp.constants.DEFAULT_INTAKE_ESM_CAT; hardcode here in case it changes unexpectedly in bcdp
DEFAULT_INTAKE_ESM_CAT = 'https://raw.githubusercontent.com/NCAR/intake-esm-datastore/master/catalogs/pangeo-cmip6.json'


col = intake.open_esm_datastore(DEFAULT_INTAKE_ESM_CAT)
#uni_dict = col.unique(["source_id"])
#pprint.pprint(uni_dict, compact=True)

modellistfile = open(listfilename, 'w')
for scenario in scenarios:

    for grid in grids:
        # grid_label: gn, gr, gr1. table_id: Oday 245 Omon 526. member_id: 68!!
        query = dict(variable_id=["tas"], table_id=['Amon'], grid_label=[grid], experiment_id=experiments) # 18
        col_subset = col.search(require_all_on=["source_id"], **query)
        
        #myset = col_subset.df.groupby("source_id")[["experiment_id"]].nunique()
        #df = col_subset.df.groupby("source_id")[["experiment_id"]]
        #print(len(myset))
        #print(myset)
        
        sources = col_subset.unique()["source_id"]["values"]
        print(sources)
        print(len(sources))
        
        
        # these have gr, but also all 5 experiments
        # CESM2
        # CESM2-WACCM
        # CNRM-CM6-1
        # MRI-ESM2-0
        
        bounds = bcdp.Bounds(lon_bnds=(-179.5,180), lat_bnds=(-89.5,90), time_bnds=(1910,1920)) # note there is an "arange" bndry convention
        grid_ds = bcdp.utils.grid_from_res((1.0,1.0), bounds)
            
        
        # one dimensional lat and lon
        grid_ds = grid_ds.assign_coords(lat=grid_ds.y)
        grid_ds = grid_ds.assign_coords(lon=grid_ds.x)
        
        # iterate over sources. get one member per source.
        # iterate over SSPs
        nsource = 0
        nmember = 0
        for mysource in sources:
            print('----------------------------------------------------------')
            print('starting %s %s' % (mysource, grid))
            print('----------------------------------------------------------')
            
            if mysource=='MRI-ESM2-0':
                print('***************** problem source id %s' % (mysource))
                continue
        
            if mysource=='AWI-CM-1-1-MR':
                print('***************** problem source id %s' % (mysource))
                continue
    
            #Note: KACE-1-0-G has 3 gr members, none of which are OK; one throws an error below, and two would need regridding due to non-0.5-degree latitude.
           
            
            # get all member_id that have all experiments
            query = dict(source_id=mysource, experiment_id=experiments, variable_id='tas', table_id='Amon', grid_label=grid)
            col_subset = col.search(require_all_on=["member_id"], **query) # require_all_on means each member_id needs to have all of the above experiments
            members = col_subset.unique()["member_id"]["values"]
         
            myset = col_subset.df.groupby("member_id")[["experiment_id"]].nunique() #Pandas DataFrame
            print(myset.sort_values(by='experiment_id', ascending=False))
            
            if oneMember:
                # find first member_id
                mymember = myset.sort_values(by='experiment_id', ascending=False).head(1).index.values[0]
                members = [mymember]
            print(members)
        
            for mymember in members:
                outfilename = outdir+mysource+'_'+mymember+'_'+grid+'_tas_'+scenario+'.nc'
                if os.path.exists(outfilename):
                    print('%s already exists, continuing.' % (outfilename))
                    continue
                
                print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
                print('starting %s %s %s %s' % (mysource, mymember, scenario, grid))
                print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

                # get historical
                query = dict(source_id=mysource, member_id=mymember, experiment_id='historical', variable_id='tas', table_id='Amon', grid_label=grid)
                try:
                    ens = bcdp.load_intake_esm(query, catfile=DEFAULT_INTAKE_ESM_CAT)
                except Exception as e:
                    #raise
                    logging.error('regrid_bcdp: source %s member %s error %s' % (mysource, mymember, str(e))) 
                    traceback.print_exc()
                    continue
                        
                ens = ens.regrid(output_grid=grid_ds, backend='esmf', method='bilinear', ignore_degenerate=True)
                da_hist = ens.apply(lambda x: x.weighted(np.cos(np.deg2rad(x.lat))).mean(dim=('x', 'y'))).first # a bit overkill b/c ens only has one da
                
                
                try:
                    # get one model, stitch together with ssp
                    query = dict(source_id=mysource, member_id=mymember, experiment_id=scenario, variable_id='tas', table_id='Amon', grid_label=grid)
                    ens = bcdp.load_intake_esm(query, catfile=DEFAULT_INTAKE_ESM_CAT)
                    
                    ens = ens.regrid(output_grid=grid_ds, backend='esmf', method='bilinear', ignore_degenerate=True)
                    da = ens.apply(lambda x: x.weighted(np.cos(np.deg2rad(x.lat))).mean(dim=('x', 'y'))).first # a bit overkill b/c ens only has one da
                                                
                    # cut short any that go to 2300, such as MRI-ESM2-0 585 and CanESM5
                    #da = da.sel(time=slice(None, cftime.DatetimeProlepticGregorian(2101, 1, 1)))
                    
                    da_concat = xr.concat([da_hist, da], 'time')
                    da_concat.name = 'tas'
                    
                    #da_concat = da_concat.rename({'x': 'lon','y': 'lat'}) # must reassign to da_concat
                        
                    # save netcdf
                    da_concat.to_netcdf(outfilename) 
                    print(outfilename) 
                    #pdb.set_trace()
                except (KeyboardInterrupt, SystemExit):
                    raise
                except Exception as e:
                    #raise
                    logging.error('regrid_bcdp: source %s member %s error %s' % (mysource, mymember, str(e))) 
                    traceback.print_exc()
                    break
                        
                # write line in file
                modellistfile.write('%s %s %s\n' % (mysource, mymember, grid))
                modellistfile.flush()

print(listfilename)
modellistfile.close()

print('NOW DELETE GR,GR1 DUPLICATES FROM LIST!')
print('done.')
#ds = ens.regrid(output_grid=grid_ds, backend='scipy', method='linear').normalize_times(assume_gregorian=True).bundle('CMIP6').first

# CAMS-CSM1-0 has an issue. *** ValueError: conflicting sizes for dimension 'time': length 3000 on 'tos' and length 3012 on 'time'
# master_time has length 3012. this dude has 3000.

