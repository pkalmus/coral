#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Load in a file and make some plots to check it.

Created on Wed Jul  8 20:33:36 2020

@author: pkalmus
"""

import xarray as xr
import pdb
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import coral_plotter as plotter
import sys


filename = '/home/pkalmus/projects//coral/output/initial/mean_tos_all_ssp126_test.nc'
#filename = '/home/pkalmus/projects//coral/data/tos/regrid/CMIP6_1x1/ACCESS-CM2_ssp126_reef.nc'
#filename = '/home/pkalmus/projects//coral/data/tos/regrid/CMIP6_1x1/ACCESS-CM2_ssp585_weights_biascorrect.nc'

ds = xr.open_dataset(filename)
pdb.set_trace()
# ds = ds.sel(isreef=1)
plotter.scatter(ds.lon.values, ds.lat.values, 'test.png', c=ds.isreef.values, vmin=0, vmax=1)
