#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

conda activate geo_env

regrid_bcdp
regrid_reef


Get model outputs using BCDP, regrid to 1x1°, and save on the local file system for regrid_reef.py to digest.

Must use same models for each SSP. All 4 or bust.

Check gr, gr1 as well. (See params file)

Note that none of the “ripf” indices can be omitted.
Example of a variant_label: if realization_index=2, initialization_index=1, physics_index=3, and forcing_index=233, then variant_label = “r2i1p3f233”.

https://github.com/intake/intake-esm
https://intake-esm.readthedocs.io/en/latest/

todo: MRI-ESM2-0 585 goes to 2300 and uses cftime.datetime objects for time, instead of datetime64.
So before concat with da_hist, would need to figure out how to convert to datetime64 or else get xr.concat to do this.
Otherwise, when attempting to write the netcdf, get ValueError: unable to infer dtype on variable 'time'; xarray cannot serialize arbitrary Python objects
Fails only for ssp585. The other 3 files are OK. 

  
todo: AWI-CM-1-1-MR
    MemoryError: Unable to allocate 5.02 TiB for an array with shape (830305, 830305) and data type float64

Note: CMIP6 historical simulation runs from 1850 to 2014.

These models appear to be only available as gr, but do not come up in the intake-ESM search. We will not use them.
    
../../data/tos/raw/CMIP6_1x1/tos_Omon_CESM2_historical_ssp126_r1i1p1f1_gr_185001-210012.nc
../../data/tos/raw/CMIP6_1x1/tos_Omon_CESM2_historical_ssp245_r1i1p1f1_gr_185001-210012.nc
../../data/tos/raw/CMIP6_1x1/tos_Omon_CESM2_historical_ssp370_r1i1p1f1_gr_185001-210012.nc
../../data/tos/raw/CMIP6_1x1/tos_Omon_CESM2_historical_ssp585_r1i1p1f1_gr_185001-210012.nc

../../data/tos/raw/CMIP6_1x1/tos_Omon_MRI-ESM2-0_historical_ssp126_r1i1p1f1_gr_185001-210012.nc
../../data/tos/raw/CMIP6_1x1/tos_Omon_MRI-ESM2-0_historical_ssp245_r1i1p1f1_gr_185001-210012.nc
../../data/tos/raw/CMIP6_1x1/tos_Omon_MRI-ESM2-0_historical_ssp370_r1i1p1f1_gr_185001-210012.nc
../../data/tos/raw/CMIP6_1x1/tos_Omon_MRI-ESM2-0_historical_ssp585_r1i1p1f1_gr_185001-210012.nc


Here's the intake-ESM search on 8/29/2020. We will use them.
(geo_env) [pkalmus@weather2 python]$ m models_gr1.txt 
INM-CM4-8 r1i1p1f1 seems OK (runs through regrid_reef)
INM-CM5-0 r1i1p1f1

(geo_env) [pkalmus@weather2 python]$ m models_gr.txt. 
CESM2-WACCM r1i1p1f1   HAS GN. Could use as a check on the regridding. would need to keep out of means.
GFDL-ESM4 r1i1p1f1     HAS GN. Could use as a check on the regridding.  
KACE-1-0-G r1i1p1f1 error getting historical run
KACE-1-0-G r2i1p1f1 note: latitude is NOT on 0.5 degree grid.
KACE-1-0-G r3i1p1f1 note: latitude is NOT on 0.5 degree grid.


So as of 8/29 there are 4 non-gn models available, two gr1 and two gr.
INM-CM4-8 r1i1p1f1
INM-CM5-0 r1i1p1f1
KACE-1-0-G r2i1p1f1
KACE-1-0-G r3i1p1f1

BE SURE TO DELETE:
CESM2-WACCM r1i1p1f1 gr
GFDL-ESM4 r1i1p1f1 gr
KACE-1-0-G r2i1p1f1 gr
KACE-1-0-G r3i1p1f1 gr

========================
3/7/2021.

mean temperature (tas)
maximum temperature (tasmax)
minimum temperature (tasmin)
precipitation (pr)
relative humidity (hurs)
cloudiness (clt)]
incoming shortwave radiation (rsds)
sea level pressure (psl)
wind speed (sfcWind))

========================

@author: pkalmus
"""


import xarray as xr
import pandas as pd
import numpy as np
import pdb
import sys
# sys.path.append(projectdir+'/pySSDF/src/')
# import make_s
from sklearn.utils.extmath import cartesian
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import coral_plotter as plotter
import bcdp
import pprint
import intake
import traceback
import subprocess
import os.path
import logging
logging.basicConfig(filename='run.log',level=logging.WARNING)

import cftime

# read in user params
import importlib
paramfile = sys.argv[1]   
params = importlib.import_module(paramfile)
projectdir = params.projectdir
runname = params.runname
basedir = params.basedir
pythondir = params.pythondir
daylistfilename = params.daylistfilename
scenarios = params.scenarios
dryRun = params.dryRun
oneMember = params.oneMember
grids = params.grids
cmip_var = params.cmip_var

if cmip_var=='tos':
    table = '0mon'
else:
    table = 'Aday'

outdir = basedir+'data/%s/%s/raw/' % (cmip_var, runname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/

experiments = scenarios.copy()
experiments.insert(0, 'historical') 

sysStr = "mkdir -p %s" % outdir
subprocess.call(sysStr, shell=True)

# this is also bcdp.constants.DEFAULT_INTAKE_ESM_CAT; hardcode here in case it changes unexpectedly in bcdp
DEFAULT_INTAKE_ESM_CAT = 'https://raw.githubusercontent.com/NCAR/intake-esm-datastore/master/catalogs/pangeo-cmip6.json'


col = intake.open_esm_datastore(DEFAULT_INTAKE_ESM_CAT)
#uni_dict = col.unique(["source_id"])
#pprint.pprint(uni_dict, compact=True)


modellistfile = open(daylistfilename, 'w')
for grid in grids:
    # grid_label: gn, gr, gr1. table_id: Oday 245 Omon 526. member_id: 68!!
    query = dict(variable_id=[cmip_var], table_id=[table], grid_label=[grid], experiment_id=experiments) # 18
    col_subset = col.search(require_all_on=["source_id"], **query)
    
    #myset = col_subset.df.groupby("source_id")[["experiment_id"]].nunique()
    #df = col_subset.df.groupby("source_id")[["experiment_id"]]
    #print(len(myset))
    #print(myset)
    
    sources = col_subset.unique()["source_id"]["values"]
    print(sources)
    print(len(sources))
    
    pdb.set_trace()
    
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

        if mysource!='CESM2':
            print('***************** source id %s is not CESM2' % (mysource))
            continue
    
        #Note: KACE-1-0-G has 3 gr members, none of which are OK; one throws an error below, and two would need regridding due to non-0.5-degree latitude.
       
        
        # get all member_id that have all experiments
        query = dict(source_id=mysource, experiment_id=experiments, variable_id=cmip_var, table_id=table, grid_label=grid)
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
            print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            print('starting %s %s %s' % (mysource, mymember, grid))
            print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            if not dryRun:
                # see if all the scenario files have created and if so continue
                # doing this here so we can allow getting "historical" only once, below, to save time
                allScenariosMade = True
                for scenario in scenarios:
                    outfilename = outdir+mysource+'_'+mymember+'_'+grid+'_'+scenario+'.nc'
                    if not os.path.exists(outfilename):
                        allScenariosMade = False
                if allScenariosMade:
                    print('all scenario files already here for %s %s %s' % (mysource, mymember, grid))
                    modellistfile.write('%s %s %s\n' % (mysource, mymember, grid))
                    nmember+=1
                    continue
    
                # get historical
                query = dict(source_id=mysource, member_id=mymember, experiment_id='historical', variable_id='tos', table_id='Omon', grid_label=grid)
                try:
                    ens = bcdp.load_intake_esm(query, catfile=DEFAULT_INTAKE_ESM_CAT)
                except Exception as e:
                    #raise
                    logging.error('regrid_bcdp: source %s member %s error %s' % (mysource, mymember, str(e))) 
                    traceback.print_exc()
                    continue
                        
                if grid=='gn':
                    da_hist = ens.regrid(output_grid=grid_ds, backend='esmf', method='bilinear', ignore_degenerate=True).first
                else:
                    da_hist = ens.first

                # only need from 1984
                da_hist = da_hist.sel(time=slice(cftime.DatetimeProlepticGregorian(1984, 1, 1), None ))
    
                for scenario in scenarios:
                        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
                        print('starting %s %s %s %s' % (mysource, mymember, scenario, grid))
                        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

                        outfilename = outdir+mysource+'_'+mymember+'_'+grid+'_'+scenario+'_DAY.nc'

                        if os.path.exists(outfilename):
                            print('%s already exists, continuing.' % (outfilename))
                            continue
                        try:
                            # get one model, stitch together with ssp
                            query = dict(source_id=mysource, member_id=mymember, experiment_id=scenario, variable_id=cmip_var, table_id=table, grid_label=grid)
                            ens = bcdp.load_intake_esm(query, catfile=DEFAULT_INTAKE_ESM_CAT)
                            
                            if grid=='gn':
                                #ds = ens.regrid(output_grid=grid_ds, backend='scipy', method='linear').first
                                da = ens.regrid(output_grid=grid_ds, backend='esmf', method='bilinear', ignore_degenerate=True).first
                            else:
                                da = ens.first
                                
                            # cut short any that go to 2300, such as MRI-ESM2-0 585 and CanESM5
                            #da = da.sel(time=slice(None, cftime.DatetimeProlepticGregorian(2101, 1, 1)))
                            
                            da_concat = xr.concat([da_hist, da], 'time')
                            da_concat.name = cmip_var
                            if grid=='gn':  
                                da_concat = da_concat.assign_coords(lon=((da_concat.lon + 360) % 360)).sortby('lon') # put into 0, 360.
                                da_concat = da_concat.rename({'x': 'lon','y': 'lat'}) # must reassign to da_concat
                            else:
                                da_concat = da_concat.drop(labels=['lon', 'lat'])
                                da_concat = da_concat.rename({'x': 'lon','y': 'lat'}) # must reassign to da_concat
                                da_concat = da_concat.assign_coords(lon=((da_concat.lon + 360) % 360)).sortby('lon') # put into 0, 360.

                                
                            # save netcdf
                            da_concat.to_netcdf(outfilename) 
                            print(outfilename) 
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
            else:
                modellistfile.write('%s %s %s\n' % (mysource, mymember, grid))
            nmember+=1
        nsource+=1
modellistfile.close()

print('NOW DELETE GR,GR1 DUPLICATES FROM LIST!')
print('done.')
#ds = ens.regrid(output_grid=grid_ds, backend='scipy', method='linear').normalize_times(assume_gregorian=True).bundle('CMIP6').first

# CAMS-CSM1-0 has an issue. *** ValueError: conflicting sizes for dimension 'time': length 3000 on 'tos' and length 3012 on 'time'
# master_time has length 3012. this dude has 3000.

