#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 09:50:05 2021

Requires creating a bunch of coarse projection experiments at different climatological baselines
DHW levels.

@author: pkalmus
"""


import xarray as xr
import pdb
import numpy as np
import coral_plotter as plotter
from scipy.interpolate import interp1d

scenario = 'ssp585'
# scenario = 'ssp370'
scenario = 'ssp245'

climstartyear_project = [1963, 1973, 1983, 1993, 2003]
climendyear_project = [1977, 1987, 1997, 2007, 2017]
meanvect = []
yearvect = []
# load in all the 15-year experiments
for climind in np.arange(len(climstartyear_project)):
    climstart = climstartyear_project[climind]
    climend = climendyear_project[climind]
    climmid = np.mean([climstart, climend])    
    myfile = '/home/pkalmus/projects//coral/data/tos/nature/coarse/departure_nn_nature_%s_1yrs_dhwbma_1975-1995_1995-2015_%i%i.nc' % (scenario, climstart, climend)
    ds = xr.open_dataset(myfile)    
    mymean = np.nanmean(ds.departure_year.values)    
    meanvect.append(mymean)
    yearvect.append(climmid)
# do a linear fit, and plot. dep_year = m1*climatology_year+b1
cc_clim = plotter.plot(yearvect, meanvect, 'sensitivity_clim_%s.png' % (scenario), xlabel='climatology mid year', ylabel='mean departure %s' % (scenario), points_only=True, xlim=[1965, 2015], linfit=True)

meanvect = []
dhws = [5,6,7,8,9]
# load in all the DHW experiments
for dhw in dhws:  
    if dhw == 8:
        myfile = '/home/pkalmus/projects//coral/data/tos/nature/coarse/departure_nn_nature_%s_1yrs_dhwbma_1975-1995_1995-2015_20032017.nc' % (scenario)       
    else:
        myfile = '/home/pkalmus/projects//coral/data/tos/nature/coarse/departure_nn_nature_%s_1yrs_dhwbma_1975-1995_1995-2015_20032017_DHW%1.1f.nc' % (scenario, dhw)
    ds = xr.open_dataset(myfile)    
    mymean = np.nanmean(ds.departure_year.values)    
    meanvect.append(mymean)
    
print(meanvect)
# do a linear fit, and plot.  dep_year = m2*DHW+b2   m2*DHW_diff+b2 = m1*clim_year_diff+b1 
cc_dhw = plotter.plot(dhws, meanvect, 'sensitivity_dhw_%s.png' % (scenario), xlabel='DHW threshold', ylabel='mean departure %s' % (scenario), points_only=True, ylim=[np.floor(np.min(meanvect)), np.ceil(np.max(meanvect))], xlim=[4,10], linfit=True)


# figure out how to go from MUR clim (mid: 2008.5) to CRW clim (mid: 1988.286)
xmur = 2008.5
ymur = cc_clim[1]*xmur + cc_clim[0]
xcrw_vect = [1988.286, 1998.5]

for xcrw in xcrw_vect:
    ycrw = cc_clim[1]*xcrw + cc_clim[0]
    adjusty = ymur-ycrw # this is the adjustment in departure years
    
    y8dhw = cc_dhw[1]*8.0 + cc_dhw[0]
    ynew = y8dhw - adjusty
    dhwnew = (ynew - cc_dhw[0]) / cc_dhw[1]
    
    print(dhwnew)
    # 4.76913919066 from 585; 
    # 4.82266931172 from 370. 
    # 4.81153389713 from 245
    #this uncertainty translates to 0.14 years of mean departure uncertainty (at 370)
    # or : 4.82266931172 - 4.76913919066 = 0.053530121060000546; at 245: 0.053530121060000546 * cc_dhw[1] = 0.203 years
    # mean: 4.8011141331700005 = 4.80.
    
    print(0.053530121060000546 * cc_dhw[1] )
    # 585: 0.11 years
    # 370: 0.14 years
    # 245: 0.20 years
    
    # do the above in one line
    dhw8_2008 = (( (cc_dhw[1]*8.0 + cc_dhw[0]) - adjusty) - cc_dhw[0]) / cc_dhw[1]
    print('this is what 8DHW in %1.3f would be with 2008.5 climatology (MUR): %1.3f' % (xcrw, dhw8_2008))
    
    # make a plot of the relationship
    dhw_1988 = np.linspace(6, 16, 100)
    delta_y = cc_dhw[1]*dhw_1988 + cc_dhw[0]
    ynew =  delta_y - adjusty
    dhw_2008 = (ynew - cc_dhw[0]) / cc_dhw[1]
    plotter.plot(dhw_1988, dhw_2008, 'dhw%i_to_dhw2008_%s.png' % (xcrw,scenario), xlabel='DHW %1.3f' % xcrw, ylabel='DHW 2008.5', points_only=True)
    
    # now you have the relationship and you can interp
    dhw8_2008_in_1988 = interp1d(dhw_2008, dhw_1988, bounds_error=False)(8.0)
    print(dhw8_2008_in_1988) # 11.2
    pdb.set_trace()