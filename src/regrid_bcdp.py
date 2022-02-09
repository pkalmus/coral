 #!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

run on weather2
conda activate geo_env

regrid_bcdp
regrid_reef

If the regridded concatenated file already exists, continue.

Get model outputs using BCDP, regrid to 1x1°, and save on the local file system for regrid_reef.py to digest.
Must use same models for each SSP. All 4 or bust.
Check gr, gr1 as well. (See params file)

Note that none of the “ripf” indices can be omitted.
Example of a variant_label: if realization_index=2, initialization_index=1, physics_index=3, and forcing_index=233, then variant_label = “r2i1p3f233”.

https://github.com/intake/intake-esm
https://intake-esm.readthedocs.io/en/latest/


SPECIAL CASESS: NOTE: I did this with Xarray 0.16.2. There is a newer version (0.21.1, https://xarray.pydata.org/en/stable/whats-new.html)
but it was too hard to figure out how to get it installed, as conda is still "broken". I think the newer xarray can serialize DatetimeProlpticGregorian.

ACCESS-CM2_r1i1p1f1_gn_ssp585:
MRI-ESM2-0_r1i1p1f1_gn_ssp585:
  Goes to 2300 and uses cftime.DatetimeProlepticGregorian. 
      fix: use a special slice, save time dimension from 370 to workaround the xarray serialization limitation when saving the netcdf.
  
todo: AWI-CM-1-1-MR
    MemoryError: Unable to allocate 5.02 TiB for an array with shape (830305, 830305) and data type float64

Note: CMIP6 historical simulation runs from 1850 to 2014.    

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
import numpy as np
import pdb
import sys
import bcdp
import intake
import cftime
import traceback
import subprocess
import os.path
import logging
logging.basicConfig(filename='run.log',level=logging.WARNING)

# read in user params
import importlib
paramfile = sys.argv[1]   
params = importlib.import_module(paramfile)
projectdir = params.projectdir
tosrunname = params.tosrunname
basedir = params.basedir
modellist_filename = params.modellist_filename
scenarios = params.scenarios
dryRun = params.dryRun
oneMember = params.oneMember
grids = params.grids
cmip_var = params.cmip_var
regrid_resolution = params.regrid_resolution

if cmip_var=='tos':
    table = 'Omon'
else:
    table = 'Amon'

do_continue = True

outdir = basedir+'data/%s/coral2/%s/raw/' % (cmip_var, tosrunname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/

experiments = scenarios.copy()
experiments.insert(0, 'historical') 

sysStr = "mkdir -p %s" % outdir
subprocess.call(sysStr, shell=True)

# this is also bcdp.constants.DEFAULT_INTAKE_ESM_CAT; hardcode here in case it changes unexpectedly in bcdp
DEFAULT_INTAKE_ESM_CAT = 'https://raw.githubusercontent.com/NCAR/intake-esm-datastore/master/catalogs/pangeo-cmip6.json'
col = intake.open_esm_datastore(DEFAULT_INTAKE_ESM_CAT)

modellistfile = open(modellist_filename, 'w')
modelresfile = open('model_resolutions.txt', 'w')
time_ACCESS_CM2_r1i1p1f1_gn_ssp585 = None # workaround, see script-level comments
time_MRI_ESM2_0_r1i1p1f1_gn_ssp370 = None # workaround, see script-level comments
for grid in grids:
    # grid_label: gn, gr, gr1. table_id: Oday 245 Omon 526. member_id: 68!!
    query = dict(variable_id=[cmip_var], table_id=[table], grid_label=[grid], experiment_id=experiments) # 18
    col_subset = col.search(require_all_on=["source_id"], **query)
        
    sources = col_subset.unique()["source_id"]["values"]
    print(sources)
    print(len(sources))
    
    bounds = bcdp.Bounds(lon_bnds=(-179.5,180), lat_bnds=(-89.5,90), time_bnds=(1910,1920)) # note there is an "arange" bndry convention
    grid_ds = bcdp.utils.grid_from_res((regrid_resolution, regrid_resolution), bounds)
    # one dimensional lat and lon
    grid_ds = grid_ds.assign_coords(lat=grid_ds.y)
    grid_ds = grid_ds.assign_coords(lon=grid_ds.x)
    
    # iterate over sources. get one member per source.
    # iterate over SSPs
    nsource = 0
    nmember = 0
    for mysource in sources:
        source_grid = mysource+'_'+grid

        # if mysource!='BCC-CSM2-MR' and mysource!='INM-CM5-0':
        #     continue

        print('----------------------------------------------------------')
        print('starting %s' % (source_grid))
        print('----------------------------------------------------------')
 
        if source_grid == 'CESM2-WACCM_gr': 
            print('GR version of GN model: %s. continuing' % (source_grid))
            continue
        if source_grid == 'GFDL-ESM4_gr':
            print('GR version of GN model: %s. continuing' % (source_grid))
            continue
        if source_grid == 'MIROC-ES2L_gr1': 
            print('GR version of GN model: %s. continuing' % (source_grid))
            continue
        if source_grid == 'MRI-ESM2-0_gr': 
            print('GR version of GN model: %s. continuing' % (source_grid))
            continue

        if mysource=='AWI-CM-1-1-MR': # 25 km model, too much memory issue (2022/02/06). try/catch handles
            print('***************%s, *** MemoryError: Unable to allocate 5.02 TiB for an array with shape (830305, 830305) and data type float64' % (mysource))
            continue
       
        if mysource=='IITM-ESM': # SSP370 is missing 2099 (ends in 2088) (2022/02/08)
            print('***************%s, SSP370 is missing 2099 (ends in 2088)' % (mysource))
            continue        
        
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
            noError = True # for some reason the noError isn't working... don't count on it and fix when necessary. for now, skipping known problem groups above.
            if not dryRun:
                # see if all the scenario files have created and if so continue. this facilitates fetching "historical" only once, below
                allScenariosMade = True
                for scenario in scenarios:
                    outfilename = outdir+mysource+'_'+mymember+'_'+grid+'_'+scenario+'.nc'
                    print(outfilename)
                    if not os.path.exists(outfilename):
                        allScenariosMade = False
                if allScenariosMade:
                    print('all scenario files already here for %s %s %s' % (mysource, mymember, grid))
                    modellistfile.write('%s %s %s\n' % (mysource, mymember, grid))
                    nmember+=1
                    if do_continue:
                        continue
    
                # get historical
                query = dict(source_id=mysource, member_id=mymember, experiment_id='historical', variable_id='tos', table_id='Omon', grid_label=grid)
                try:
                    ens = bcdp.load_intake_esm(query, catfile=DEFAULT_INTAKE_ESM_CAT)
                except (Exception) as e:
                    logging.error('regrid_bcdp: source %s member %s error %s' % (mysource, mymember, str(e))) 
                    traceback.print_exc()
                    noError = False
                    if do_continue:
                        continue
                mymod = ens.first
                
                nx = len(mymod.x.values)
                londiff = np.abs(np.diff(mymod.lon.values[40,:]))
                londiff = londiff[londiff < 5]
                meanlondiff1 = np.mean(londiff)
                londiff = np.abs(np.diff(mymod.lon.values[70,:]))
                londiff = londiff[londiff < 5]
                meanlondiff2 = np.mean(londiff)
                print('\n\n****************%s %s %s %1.2f %1.2f %i \n\n' % (mysource, mymember, grid, meanlondiff1, meanlondiff2, nx))
                modelresfile.write('%s %s %s %1.2f %1.2f %i \n' % (mysource, mymember, grid, meanlondiff1, meanlondiff2, nx))
                modelresfile.flush()
                                
                da_hist = ens.regrid(output_grid=grid_ds, backend='esmf', method='bilinear', ignore_degenerate=True).normalize_times().first
                
                # do the SSP stitching, regridding, and homogenization
                for scenario in scenarios:
                        identifier = mysource+'_'+mymember+'_'+grid+'_'+scenario
                        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
                        print('starting %s ' % (identifier))
                        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
                        outfilename = outdir+identifier+'.nc'
                        if os.path.exists(outfilename):
                            print('%s already exists, continuing.' % (outfilename))
                            if do_continue:
                                continue
                        try:
                            # get one model, stitch together with ssp
                            query = dict(source_id=mysource, member_id=mymember, experiment_id=scenario, variable_id=cmip_var, table_id=table, grid_label=grid)
                            ens = bcdp.load_intake_esm(query, catfile=DEFAULT_INTAKE_ESM_CAT)

                            da = ens.regrid(output_grid=grid_ds, backend='esmf', method='bilinear', ignore_degenerate=True).normalize_times().first
                                
                            da_concat = xr.concat([da_hist, da], 'time')
                            da_concat.name = cmip_var
                            print(np.max(da_concat.time.values))

                            # from 1970 to 2100. cut short any that go to 2300, such as MRI-ESM2-0 585 and CanESM5
                            # the ocean heat signal really started in 1970                
                            if identifier=='ACCESS-CM2_r1i1p1f1_gn_ssp370':
                                da_concat = da_concat.sel(time=slice('1970', '2099'))
                                time_ACCESS_CM2_r1i1p1f1_gn_ssp370 = da_concat.time 
                                print('@@@@@@@@@SPECIAL CASE: storing time from ACCESS-CM2_r1i1p1f1_gn_ssp370 for use with ssp585.')
                            elif identifier=='ACCESS-CM2_r1i1p1f1_gn_ssp585': # goes to 2300 and uses this funky calendar
                                da_concat = da_concat.sel(time=slice(cftime.DatetimeProlepticGregorian(1970,1,1), cftime.DatetimeProlepticGregorian(2099,12,15)))
                                da_concat = da_concat.assign_coords(time=time_ACCESS_CM2_r1i1p1f1_gn_ssp370.copy()) 
                                print('@@@@@@@@@SPECIAL CASE: ACCESS-CM2_r1i1p1f1_gn_ssp585.')
                            elif identifier=='MRI-ESM2-0_r1i1p1f1_gn_ssp370':
                                da_concat = da_concat.sel(time=slice('1970', '2099'))
                                time_MRI_ESM2_0_r1i1p1f1_gn_ssp370 = da_concat.time 
                                print('@@@@@@@@@SPECIAL CASE: storing time from MRI-ESM2-0_r1i1p1f1_gn_ssp370 for use with ssp585.')
                            elif identifier=='MRI-ESM2-0_r1i1p1f1_gn_ssp585': # goes to 2300 and uses this funky calendar
                                da_concat = da_concat.sel(time=slice(cftime.DatetimeProlepticGregorian(1970,1,1), cftime.DatetimeProlepticGregorian(2099,12,15)))
                                da_concat = da_concat.assign_coords(time=time_MRI_ESM2_0_r1i1p1f1_gn_ssp370.copy()) 
                                print('@@@@@@@@@SPECIAL CASE: MRI-ESM2-0_r1i1p1f1_gn_ssp585.')                   
                            else:
                                da_concat = da_concat.sel(time=slice('1970', '2099'))

                            da_concat = da_concat.assign_coords(lon=((da_concat.lon + 360) % 360)).sortby('lon') # put into 0, 360.
                            da_concat = da_concat.rename({'x': 'lon','y': 'lat'}) # must reassign to da_concat

                            # save netcdf
                            da_concat.to_netcdf(outfilename) 
                            print(outfilename) 
                        except (KeyboardInterrupt, SystemExit):
                            raise
                        except Exception as e:
                            #raise
                            logging.error('regrid_bcdp: source %s member %s error %s' % (mysource, mymember, str(e))) 
                            traceback.print_exc()
                            noError = False
                            break
                if noError:
                    # write line in file
                    modellistfile.write('%s %s %s\n' % (mysource, mymember, grid))
                    modellistfile.flush()
            else:
                modellistfile.write('%s %s %s\n' % (mysource, mymember, grid))
            nmember+=1
        nsource+=1
modellistfile.close()
print(modellist_filename)
print('IF DRYRUN, members with errors (the errors have not been actualized) will still be written to the list.')
print('regrid_bcdp.py should have checks to make sure unusable models are not included in the model list. But be careful...')
print('done.')


