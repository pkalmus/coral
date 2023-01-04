#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
python make_spike_diff_table.py params

read in txt file created by fine_plot.py (via coral_plotter.py) and do a
diff with the 3 spike experiments and make a table for supplement.


@author: pkalmus
"""
import xarray as xr
import numpy as np
import pdb
import coral_plotter as plotter
import sys
import pandas as pd

# read in user params
import importlib
paramfile = sys.argv[1]   
params = importlib.import_module(paramfile)
projectdir = params.projectdir
runname = params.runname
basedir = params.basedir
pythondir = params.pythondir
listfilename = params.listfilename
scenarios = params.scenarios
do_plots = params.do_plots
return_years = params.return_years
bias_correct = params.bias_correct
do_weights = params.do_weights
do_obs = params.do_obs
finerun = params.finerun
spike_sizes = [1,2,1]
spike_freqs = [5,5,10]
DHW_threshold = 4.8

do_uncert = True

# define in and out dirs
finedir = basedir+'output/%s/fine/%s/' % (runname, finerun)

outstr = '\\begin{table}[!h]\n'
outstr = outstr+'\centering\n'
outstr = outstr+'\\begin{tabular}{|c|c|c|c||c|c|c||c|c|c|}\n'
outstr = outstr+'	\hline\n'

outstr = outstr+'	& \multicolumn{3}{|c||}{TD1Y} & \multicolumn{3}{|c||}{TD3Y} & \multicolumn{3}{|c|}{TD5Y} \\\\ \n'
outstr = outstr+'	\hline\n'
outstr = outstr+'	SSP & 30\% & 10\% & 1\% & 30\% & 10\% & 1\% & 30\% & 10\% & 1\% \\\\ \n'
outstr = outstr+'	\hline\n'
    
for spike_ind in np.arange(0,3):
    spike_size = spike_sizes[spike_ind]
    spike_freq = spike_freqs[spike_ind]

    outstr = outstr+'	\hline\n'
    outstr = outstr+'	\multicolumn{10}{|c|}{Difference (years), injections (mag. %i at %i years) minus baseline }\\\\ \n' % (spike_size, spike_freq)
    outstr = outstr+'	\hline\n'

    # do it for years
    for (ind, scenario) in enumerate(scenarios): # reload each file 4 times, once per scenario
        outstr = outstr+ '%s ' % (scenario[3:])
        for return_year in return_years:
            experiment = 'SPIKE%im%if_%s_%iyrs_DHW%1.1f' % (spike_size, spike_freq, 'ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
            txtfilename=finedir+'year_%s_%s.txt' % (runname, experiment)
            print(txtfilename)    
            df = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'pct2020', 'yr30pct', 'yr10pct', 'yr1pct'])

            experiment_nospike = '%s_%iyrs_DHW%1.1f' % ('ssp585', return_year, DHW_threshold)
            txtfilename=finedir+'year_%s_%s.txt' % (runname, experiment_nospike)
            print(txtfilename)
            df_nospike = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'pct2020', 'yr30pct', 'yr10pct', 'yr1pct'])     

            bigdash='\\vcenter{\\hbox{\\rule{10pt}{0.5pt}}}'
            # bigdash='-'
            # uncertdash='-'
            val = df.yr30pct.values[ind] - df_nospike.yr30pct.values[ind] 
            if np.isnan(val):    
                valstr = bigdash
            else:
                valstr = '%i' % (np.round(val))
            outstr = outstr+ '& %s' % valstr

    
            val = df.yr10pct.values[ind] - df_nospike.yr10pct.values[ind]
            if np.isnan(val):    
                valstr = bigdash
            else:
                valstr = '%i' % (np.round(val))
            outstr = outstr+ '& %s' % valstr

    
            val = df.yr1pct.values[ind] - df_nospike.yr1pct.values[ind]
            if np.isnan(val):    
                valstr = bigdash
            else:
                valstr = '%i' % (np.round(val))
            outstr = outstr+ '& %s' % valstr            

                    
        outstr = outstr+' \\\\ \n'
        outstr = outstr+'	\hline\n'
     
    
    # do it for temps
    outstr = outstr+'	\multicolumn{10}{|c|}{Difference (\\degreec), injections (mag. %i at %i years) minus baseline }\\\\ \n' % (spike_size, spike_freq)
    outstr = outstr+'	\hline\n'
    for (ind, scenario) in enumerate(scenarios): # reload each file 4 times, once per scenario
        if ind==0: # no GMST for ssp126
            continue
        ind-=1
        outstr = outstr+ '%s ' % (scenario[3:])
        for return_year in return_years:
            experiment = 'SPIKE%im%if_%s_%iyrs_DHW%1.1f' % (spike_size, spike_freq, 'ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
            txtfilename=finedir+'temp_%s_%s.txt' % (runname, experiment)
            print(txtfilename)
            df = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'pct2020', 'yr30pct', 'yr10pct', 'yr1pct'])
    
            experiment_nospike = '%s_%iyrs_DHW%1.1f' % ('ssp585', return_year, DHW_threshold)
            txtfilename=finedir+'temp_%s_%s.txt' % (runname, experiment_nospike)
            print(txtfilename)
            df_nospike = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'pct2020', 'yr30pct', 'yr10pct', 'yr1pct'])    
    
            val = df.yr30pct.values[ind] - df_nospike.yr30pct.values[ind]
            if np.abs(val) < 0.05:
                val = np.abs(val) # to take care of -0.0's in the display
            valstr = '%1.1f' % val
            outstr = outstr+ '& %s' % valstr
    
            val = df.yr10pct.values[ind] - df_nospike.yr10pct.values[ind]
            if np.abs(val) < 0.05:
                val = np.abs(val) # to take care of -0.0's in the display            valstr = '%1.1f' % val
            outstr = outstr+ '& %s' % valstr

            val = df.yr1pct.values[ind] - df_nospike.yr1pct.values[ind]
            if np.abs(val) < 0.05:
                val = np.abs(val) # to take care of -0.0's in the display            valstr = '%1.1f' % val    
            outstr = outstr+ '& %s' % valstr
            
        outstr = outstr+' \\\\ \n'
        outstr = outstr+'	\hline\n'

       
outstr = outstr+'\end{tabular}\n'
outstr = outstr+'\caption{ Similar to Table\\,\\ref{table:results} but showing differences between experiments with random pulses injected into the downscaled SST time series, and the baseline results from Table\\,\\ref{table:results}. Pulses are added with a spacing of every five or ten years as indicated but with plus or minus 24 months chosen from an integer uniform distribution, have a duration of two, four, or six months also chosen randomly from an integer uniform distribution, and have a magnitude chosen randomly from one or two times the standard normal distribution as indicated. Dashes indicate that the milestone is not reached prior to 2100 in either the experiment or baseline.}\n'
outstr = outstr+'\label{table:resultsSpikeDiff}\n'
outstr = outstr+'\end{table}\n'
print(outstr)