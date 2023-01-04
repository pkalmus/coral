#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
python repackage_mur_2019.py params filenum
filenum = [0 to 8]

then run fine_project.py with include_mur=True on.

Ayesha gave the coral pixel MUR files out of order, without NaN, and without std.
This just puts them in same format as the laGP files.


Created on Fri Mar 19 18:50:03 2021

@author: pkalmus
"""

import xarray as xr
import pdb
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import coral_plotter as plotter
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
return_years = params.return_years
DHW_threshold = params.DHW_threshold
bias_correct = params.bias_correct
do_weights = params.do_weights
do_obs = params.do_obs
finerun = params.finerun

filenum = int(sys.argv[2])

finedir = basedir+'output/%s/fine/%s/' % (runname, finerun) # e.g. for scp: /raid8/pkalmus/data/coral/data/tos/nature/fine/flat/
murdir = basedir + 'data/mur/MUR_coral-pixels_2002-2019/'

scenario = 'ssp370' # this is just for template
# scenario = 'ssp126' # has NaNs
lagp_files = np.sort(glob.glob(finedir+'Globe*%s*.nc' % (scenario))) # make sure you don't glob output files from THIS code
mur_files = np.sort(glob.glob(murdir+'Globe*.nc'))

#lagp_file = finedir+'%s.Globe.310_340.min30_5_ncdf4.nc' % (scenario)

for (ind, lagp_file) in enumerate(lagp_files):
    if ind == filenum:
        error_ctr = 0
        print(lagp_file)
        newfile = mur_files[ind].replace('Globe', 'sorted')
        print(newfile)
        ds = xr.open_dataset(lagp_file)
        ds = ds.dropna('index', subset=['temp'])
        dsmur = xr.open_dataset(mur_files[ind])
        dsmur = dsmur.dropna('index', subset=['temp'])
    
    
        # make empty containers with the index length of dsmur, and fill
        dsnew = xr.full_like(dsmur, fill_value=np.nan)
        temp_sdnew = xr.full_like(dsnew.temp, fill_value=np.nan)
        dsnew = dsnew.assign(temp_sd=temp_sdnew)
        
        print(ds.temp.shape)
        print(dsmur.temp.shape)
        #pdb.set_trace()
        
        for i in np.arange(len(ds.index.values)):
            mylon=ds.lon.values[i]
            mylat=ds.lat.values[i]
            myds=dsmur.where((dsmur.lat==mylat) & (dsmur.lon==mylon), drop=True)
            #myds=dsmur.sel(index=((dsmur.lon==mylon) & (dsmur.lat==mylat)))
            
            # tol = 0.1
            # myds=dsmur.where((dsmur.lat<mylat+tol) & (dsmur.lat>mylat-tol) & (dsmur.lon<mylon+tol) & (dsmur.lon>mylon-tol), drop=True)
            
            try:
                dsnew.temp.values[i,:]=myds.temp.values
                dsnew.temp_sd.values[i,:]=0.0 # need this to not have NaNs, which lead to no high/low dep years
            except:
                error_ctr += 1
                print('no lat lon match. errors so far: %i' % error_ctr)
                pdb.set_trace()
            dsnew.lon.values[i]=mylon
            dsnew.lat.values[i]=mylat
            if np.mod(i, 10)==0:
                print(str(i)+'/'+str(len(ds.index.values)))
            # if i > 10:
            #     break
        dsnew.to_netcdf(newfile)
        print(error_ctr)
        #pdb.set_trace()
        