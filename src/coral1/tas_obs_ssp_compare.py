#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  3 22:46:15 2021

@author: pkalmus
"""

import numpy as np
import pandas as pd
import xarray as xr
import pdb

import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import cartopy.feature
import cartopy.crs as ccrs

scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp585'] 
colors = ['k', 'b', 'g', 'r']

# load hadcrut5
hadfile = '/raid8/pkalmus/data/coral/data/tas/obs/HadCRUT5/HadCRUT5.0Analysis_gl.txt'
had = np.loadtxt(hadfile, usecols=np.arange(0,13))
had = pd.DataFrame(had[::2])

years = had.values[:,0]
years = np.repeat(years,12)
monthfrac = np.arange(0,12)/12. + 1/24.
monthfrac = np.tile(monthfrac, years.shape[0]/12)
years = years + monthfrac

vals = had.values[:,1:]
vals = vals.ravel()


fig, ax = plt.subplots()
ax.plot(years, vals, 'k', linewidth='3')

# load SSPs
for ind, scenario in enumerate(scenarios):
    sspfile = '/raid8/pkalmus/data/coral/data/tas/nature/mean_%s.nc' % (scenario)
    ds = xr.open_dataset(sspfile)
    #pdb.set_trace()
    ax.plot(ds.time.values, ds.anom.values, colors[ind])

# make plot.
ax.set(xlabel='Years', ylabel='Anomaly (C)')
ax.grid()
plt.xlim(2010,2050)
plt.ylim(0,2)

plotfile = '/raid8/pkalmus/data/coral/data/tas/obs/HadCRUT5/tas_obs_ssp.png'
plt.savefig(plotfile)
plt.close()
print(plotfile)

pdb.set_trace()