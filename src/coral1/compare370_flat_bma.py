#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Read in flat 370 file, bma 370 file.

Make correlation plot.
Find out how many of the 1% pixels are shared by both. (1218 / 8400, 14.5%) 

Created on Fri Mar 26 21:33:22 2021

@author: pkalmus
"""

import xarray as xr
import numpy as np
import pdb
import coral_plotter as plotter
import sys

flatfile = '/home/pkalmus/projects/coral/output/nature/fine/flat/fine_nature_FLAT_ssp370_1yrs_DHW4.8.nc'
dsflat = xr.open_dataset(flatfile)
bmafile = '/home/pkalmus/projects/coral/output/nature/fine/bma/fine_nature_ssp370_1yrs_DHW4.8.nc'
dsbma = xr.open_dataset(bmafile)
dsbma = dsbma.dropna('index', subset=['departure_year'])


# this is kind of slow but not awful. 
dsbma_sort = dsbma.sortby('departure_year', ascending=False)
dsflat_sort = dsflat.sortby('departure_year', ascending=False)

yearbma = dsbma_sort.departure_year.values[0:8400]
lonbma = dsbma_sort.lon.values[0:8400]
latbma = dsbma_sort.lat.values[0:8400]
indbma = lonbma*1000.+latbma
yearflat = dsflat_sort.departure_year.values[0:8400]
lonflat = dsflat_sort.lon.values[0:8400]
latflat = dsflat_sort.lat.values[0:8400]
indflat = lonflat*1000+latflat

matchctr = 0
for ind in indbma:
    for ind2 in indflat:
        if ind==ind2:
            print('match!')
            matchctr+=1

print(matchctr)



# departure_year_nonan = dsflat.departure_year.values[~np.isnan(dsflat.departure_year.values)]
# thres1pct = np.sort(departure_year_nonan)[int(np.ceil(0.99 * len(departure_year_nonan)))]
# dsflat_pct = dsflat.where(dsflat.departure_year >= thres1pct, drop=True)

# departure_year_nonan = dsbma.departure_year.values[~np.isnan(dsbma.departure_year.values)]
# thres1pct = np.sort(departure_year_nonan)[int(np.ceil(0.99 * len(departure_year_nonan)))]
# dsbma_pct = dsbma.where(dsbma.departure_year >= thres1pct, drop=True)

pdb.set_trace()


plotter.corrPlot(dsbma_pct.departure_year.values, dsflat_pct.departure_year.values,'BMA','Unweighted', '/home/pkalmus/projects/coral/output/nature/fine/bma/bma_flat_corr.png' ,maxxy=None,dofit=False,nbins=0,ylim=None,dology=False,label=None,label2=None)
