#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Load in two time series and experiment with the pval from R code.

Created on Wed Jul  8 20:33:36 2020

@author: pkalmus
"""

projectdir = '/home/pkalmus/projects/'

import sys
import xarray as xr
import pandas as pd
import numpy as np
import pdb
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import coral_plotter as plotter
import subprocess

# these two lines needed to pass np array into R
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

utils = importr('utils')
base = importr('base')
print(base._libPaths())

#utils.install_packages('wavethresh, WiSEBoot, normwhn.test, tseries, forecast', repos='http://r-forge.r-project.org') # only need to load once
# SEE here for alternative: https://stackoverflow.com/questions/15419740/calling-custom-functions-from-python-using-rpy2
robjects.r('''
       source('../r/EstimateWeight.R')
''')
r_func = robjects.globalenv['EstimateWeight']

R = 500.0 # 
len = 512

# mod = np.loadtxt('/home/pkalmus/projects//coral/data/tos/weights/nodetrend/CESM2_ssp585_mod4_0.002_1024.txt')
# obs = np.loadtxt('/home/pkalmus/projects//coral/data/tos/weights/nodetrend/CESM2_ssp585_obs4_0.002_1024.txt')



mod = np.loadtxt('/home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_mod840.520.txt')
obs = np.loadtxt('/home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_obs840.520.txt')

# mod = np.loadtxt('/home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_mod2_0.002_1024.txt')
# obs = np.loadtxt('/home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_obs2_0.002_1024.txt')

# mod = np.loadtxt('/home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_mod20.020_1024.txt')
# obs = np.loadtxt('/home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_obs20.020_1024.txt')


# mod = np.loadtxt('/home/pkalmus/projects//coral/data/tos/weights/nn/CNRM-CM6-1_ssp585_mod209_0.740_1024.txt')
# obs = np.loadtxt('/home/pkalmus/projects//coral/data/tos/weights/nn/CNRM-CM6-1_ssp585_obs209_0.740_1024.txt')



#mod[:,1] = obs[:,1]-30.0


# mod[:,1] = np.ones(len)+np.sin(np.arange(len))+np.arange(len)
# obs[:,1] = np.ones(len)+np.sin(np.arange(len))+np.arange(len)*0.9995



# call the R function and get weight
# note: if you put in "mod, mod" it gives pvals very near 1.
pval = r_func(obs, mod)[0]
if pval==0:
    pval = 1/R
if pval==1:
    pval = (R-1)/R
print('******************************************************** %1.8f' % pval)

plotter.timeseries(obs, './test_timeseries_CESM2_ssp585_mod840.520.png', t2=mod, title='pval: %1.3f' % (pval))
# plotter.timeseries(obs, './test_timeseries_CESM2_ssp585_obs2_0.002_1024.png', t2=mod, title='pval: %1.3f' % (pval))
# plotter.timeseries(obs, './test_timeseries_CESM2_ssp585_obs20.020_1024.png', t2=mod, title='pval: %1.3f' % (pval))
# plotter.timeseries(obs, './test_timeseries_CNRM-CM6-1_ssp585_mod209_0.740_1024.png', t2=mod, title='pval: %1.3f' % (pval))

# (geo_env) [pkalmus@weather2 python]$ ls /home/pkalmus/projects//coral/data/tos/weights/nn/*.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_mod2_0.002_1024.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_mod20.002_1024.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_mod20.002.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_mod2_0.020_1024.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_mod20.020_1024.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_mod65_0.002_1024.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_mod650.002_1024.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_mod650.540.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_mod840.520.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_obs2_0.002_1024.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_obs20.002_1024.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_obs20.002.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_obs2_0.020_1024.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_obs20.020_1024.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_obs65_0.002_1024.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_obs650.002_1024.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_obs650.540.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CESM2_ssp585_obs840.520.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CNRM-CM6-1_ssp585_mod209_0.740_1024.txt
# /home/pkalmus/projects//coral/data/tos/weights/nn/CNRM-CM6-1_ssp585_obs209_0.740_1024.txt


