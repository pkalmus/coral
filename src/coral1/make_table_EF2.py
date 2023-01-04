#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
python make_table_EF.py params

RUN fine_plot.py FIRST!!!
read in txt file created by fine_plot.py (via coral_plotter.py) and create the table for the paper.

Based on make_table EF.py. This is for the revision.
It takes 3 climatological baselines, and two return timescales.


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
includemur = params.include_mur
do_add_variability = params.do_add_variability
spike_size = params.spike_size
spike_freq = params.spike_freq
DHW_thresholds = params.DHW_thresholds
include_uncert = params.include_uncert_table



# define in and out dirs
finedir = basedir+'output/%s/fine/%s/' % (runname, finerun)

for return_year in return_years:
    # do it for projected years
    outstr = '\\begin{table}[ht]\n'
    outstr = outstr+'\centering\n'
    outstr = outstr+'\caption{Projected years and GMSTAs after which fewer than the stated percentage of 1 km$^2$ reef locations remain below the thermal thresholds, for a return timescale of %i years (TD%iY)}\n' % (return_year, return_year)
    
    
    outstr = outstr+'\\begin{tabular}{cccccccccc}\n'
    outstr = outstr+'	\hline\n'
    
    outstr = outstr+'	& \multicolumn{3}{c}{8 DHW$_{2008}$ } & \multicolumn{3}{c}{8 DHW$_{1998}$ } & \multicolumn{3}{c}{8 DHW$_{1988}$ } \\\\ \n'
    outstr = outstr+'	\hline\n'
    outstr = outstr+'	 & 30\% & 10\% & 1\% & 30\% & 10\% & 1\% & 30\% & 10\% & 1\% \\\\ \n'
    outstr = outstr+'	\hline\hline\n'
    outstr = outstr+'	\multicolumn{10}{c}{Year in twenty-first century }\\\\ \n'
    outstr = outstr+'	\hline\n'


    for (ind, scenario) in enumerate(scenarios): # reload each file 4 times, once per scenario
        outstr = outstr+ 'SSP%s ' % (scenario[3:])
        for DHW_threshold in DHW_thresholds:
            experiment = '%s_%s_%iyrs_DHW%1.1f' % (finerun, 'ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
        

            if do_add_variability:
                experiment = 'SPIKE_%s_%iyrs_DHW%1.1f' % ('ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod        
                # experiment = 'SPIKE2_%s_%iyrs_DHW%1.1f' % ('ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod        
                # experiment = 'SPIKE%im%if_%s_%iyrs_DHW%1.1f' % (spike_size, spike_freq, 'ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
    
            txtfilename=finedir+'year_%s_%s.txt' % (runname, experiment)
            print(txtfilename)    
            df = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'pct2020', 'yr30pct', 'yr10pct', 'yr1pct'])
    
            if include_uncert:
                txtfilename=finedir+'year_%s_%s_high.txt' % (runname, experiment)
                print(txtfilename)    
                dfhigh = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'pct2020', 'yr30pct', 'yr10pct', 'yr1pct'])
                txtfilename=finedir+'year_%s_%s_low.txt' % (runname, experiment)
                print(txtfilename)    
                dflow = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'pct2020', 'yr30pct', 'yr10pct', 'yr1pct'])
                
            bigdash='\\vcenter{\\hbox{\\rule{10pt}{0.5pt}}}'
            uncertdash='\\mathrm{\\vcenter{\\hbox{\\rule{8pt}{0.5pt}}}}'
            # bigdash='-'
            # uncertdash='-'
            val = df.yr30pct.values[ind]
            if np.isnan(val):    
                valstr = bigdash
            else:
                valstr = '%i' % (np.floor(val-2000)) # we take floor to report "year in which something is projected to happen"
            outstr = outstr+ '& %s' % valstr
            if include_uncert:
                val = dfhigh.yr30pct.values[ind]
                if np.isnan(val):    
                    valstrhigh = uncertdash
                else:
                    valstrhigh = '%i' % (np.floor(val-2000))
                val = dflow.yr30pct.values[ind]
                if np.isnan(val):    
                    valstrlow = uncertdash
                else:
                    valstrlow = '%i' % (np.floor(val-2000))
                outstr = outstr+ '$^{%s}_{%s}$ ' % (valstrhigh, valstrlow)
    
            val = df.yr10pct.values[ind]
            if np.isnan(val):    
                valstr = bigdash
            else:
                valstr = '%i' % (np.floor(val-2000))
            outstr = outstr+ '& %s' % valstr
            if include_uncert:
                val = dfhigh.yr10pct.values[ind]
                if np.isnan(val):    
                    valstrhigh = uncertdash
                else:
                    valstrhigh = '%i' % (np.floor(val-2000))
                val = dflow.yr10pct.values[ind]
                if np.isnan(val):    
                    valstrlow = uncertdash
                else:
                    valstrlow = '%i' % (np.floor(val-2000))
                outstr = outstr+ '$^{%s}_{%s} $' % (valstrhigh, valstrlow)
    
            val = df.yr1pct.values[ind]
            if np.isnan(val):    
                valstr = bigdash
            else:
                valstr = '%i' % (np.floor(val-2000))
            outstr = outstr+ '& %s' % valstr            
            if include_uncert:
                val = dfhigh.yr1pct.values[ind]
                if np.isnan(val):    
                    valstrhigh = uncertdash
                else:
                    valstrhigh = '%i' % (np.floor(val-2000))
                val = dflow.yr1pct.values[ind]
                if np.isnan(val):    
                    valstrlow = uncertdash
                else:
                    valstrlow = '%i' % (np.floor(val-2000)) 
                outstr = outstr+ '$^{%s}_{%s} $' % (valstrhigh, valstrlow)
            #outstr = outstr+ '& %i & %i & %i ' % (np.round(df.yr30pct.values[ind]), np.round(df.yr10pct.values[ind]), np.round(df.yr1pct.values[ind]) )
        outstr = outstr+' \\\\ \n'
    outstr = outstr+'	\hline\n'
 

    # do it for temps
    outstr = outstr+'	\multicolumn{10}{c}{Global mean surface temperature anomaly (\\degreec) }\\\\ \n'
    outstr = outstr+'	\hline\n'
    for (ind, scenario) in enumerate(scenarios): # reload each file 4 times, once per scenario
        if ind==0: # no GMST for ssp126
            continue
        ind-=1
        outstr = outstr+ 'SSP%s ' % (scenario[3:])
        for DHW_threshold in DHW_thresholds:
            experiment = '%s_%s_%iyrs_DHW%1.1f' % (finerun, 'ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
    
            if do_add_variability:
                experiment = 'SPIKE_%s_%iyrs_DHW%1.1f' % ('ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod        
                # experiment = 'SPIKE2_%s_%iyrs_DHW%1.1f' % ('ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod        
                # experiment = 'SPIKE%im%if_%s_%iyrs_DHW%1.1f' % (spike_size, spike_freq, 'ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
    
            
            txtfilename=finedir+'temp_%s_%s.txt' % (runname, experiment)
            print(txtfilename)
            df = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'pct2020', 'yr30pct', 'yr10pct', 'yr1pct'])
    
            if include_uncert:
                txtfilename=finedir+'temp_%s_%s_high.txt' % (runname, experiment)
                print(txtfilename)    
                dfhigh = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'pct2020', 'yr30pct', 'yr10pct', 'yr1pct'])
                txtfilename=finedir+'temp_%s_%s_low.txt' % (runname, experiment)
                print(txtfilename)    
                dflow = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'pct2020', 'yr30pct', 'yr10pct', 'yr1pct'])
    
            val = df.yr30pct.values[ind]
            valstr = '%1.1f' % val
            outstr = outstr+ '& %s' % valstr
            if include_uncert:
                val = dfhigh.yr30pct.values[ind]
                valstrhigh = '%1.1f' % val
                val = dflow.yr30pct.values[ind]
                valstrlow = '%1.1f' % val
                outstr = outstr+ '$^{%s}_{%s}$ ' % (valstrhigh, valstrlow)
    
            val = df.yr10pct.values[ind]
            valstr = '%1.1f' % val
            outstr = outstr+ '& %s' % valstr
            if include_uncert:
                val = dfhigh.yr10pct.values[ind]
                valstrhigh = '%1.1f' % val
                val = dflow.yr10pct.values[ind]
                valstrlow = '%1.1f' % val
                outstr = outstr+ '$^{%s}_{%s}$ ' % (valstrhigh, valstrlow)
    
            val = df.yr1pct.values[ind]
            valstr = '%1.1f' % val
            outstr = outstr+ '& %s' % valstr
            if include_uncert:
                val = dfhigh.yr1pct.values[ind]
                valstrhigh = '%1.1f' % val
                val = dflow.yr1pct.values[ind]
                valstrlow = '%1.1f' % val
                outstr = outstr+ '$^{%s}_{%s}$ ' % (valstrhigh, valstrlow)        
            #outstr = outstr+ '& %i & %i & %i ' % (np.round(df.yr30pct.values[ind]), np.round(df.yr10pct.values[ind]), np.round(df.yr1pct.values[ind]) )
        outstr = outstr+' \\\\ \n'
    outstr = outstr+'	\hline\n'

    outstr = outstr+'\end{tabular}\n'
    outstr = outstr+'\label{table:results}\n'
    outstr = outstr+'\end{table}\n'
    print(outstr)  




for return_year in return_years:
    # 2nd table using GMSTA thresholds
    # switch header
    outstr = '\\begin{table}[ht]\n'
    outstr = outstr+'\centering\n'
    outstr = outstr+'\caption{Percentages and numbers of reef locations remaining below the stated thresholds, for a return timescale of %i years (TD%iY)}\n' % (return_year, return_year)
    
    outstr = outstr+'\\begin{tabular}{cccccccccc}\n'
    outstr = outstr+'	\hline\n'
    
    outstr = outstr+'	& \multicolumn{3}{c}{8 DHW$_{2008}$ } & \multicolumn{3}{c}{8 DHW$_{1998}$ } & \multicolumn{3}{c}{8 DHW$_{1988}$ } \\\\ \n'
    outstr = outstr+'	\hline\n'
    outstr = outstr+'	 & 1.5\\degreec & 1.7\\degreec & 2.0\\degreec & 1.5\\degreec & 1.7\\degreec & 2.0\\degreec  & 1.5\\degreec & 1.7\\degreec & 2.0\\degreec  \\\\ \n'
    outstr = outstr+'	\hline\hline\n'    
    # do it for percentages above 1.5, 1.7, 2.0
    outstr = outstr+'	\multicolumn{10}{c}{Percent 1 km$^2$ reef locations remaining below threshold }\\\\ \n'
    outstr = outstr+'	\hline\n'

    for (ind, scenario) in enumerate(scenarios): # reload each file 4 times, once per scenario
        if ind==0: # no GMST for ssp126
            continue
        ind-=1
        outstr = outstr+ 'SSP%s ' % (scenario[3:])
        for DHW_threshold in DHW_thresholds:
            experiment = '%s_%s_%iyrs_DHW%1.1f' % (finerun, 'ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
    
            if do_add_variability:
                experiment = 'SPIKE_%s_%iyrs_DHW%1.1f' % ('ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod        
                # experiment = 'SPIKE2_%s_%iyrs_DHW%1.1f' % ('ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod        
                # experiment = 'SPIKE%im%if_%s_%iyrs_DHW%1.1f' % (spike_size, spike_freq, 'ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
    
            
            txtfilename=finedir+'temp_%s_%s_1p51p72p0.txt' % (runname, experiment)
            print(txtfilename)
            df = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'pct1p5', 'pct1p7', 'pct2p0'])
    
            if include_uncert:
                txtfilename=finedir+'temp_%s_%s_1p51p72p0high.txt' % (runname, experiment)
                print(txtfilename)    
                dfhigh = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'pct1p5', 'pct1p7', 'pct2p0'])
                txtfilename=finedir+'temp_%s_%s_1p51p72p0low.txt' % (runname, experiment)
                print(txtfilename)    
                dflow = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'pct1p5', 'pct1p7', 'pct2p0'])
    
            val = df.pct1p5.values[ind]*100.
            valstr = '%1.0f' % val
            #valstr = '{number:.{digits}f}'.format(number=val, digits=0)
            outstr = outstr+ '& %s\%% ' % valstr
            if include_uncert:
                val = dfhigh.pct1p5.values[ind]*100.
                valstrhigh = '%1.0f' % val
                val = dflow.pct1p5.values[ind]*100.
                valstrlow = '%1.0f' % val
                outstr = outstr+ '$^{%s}_{%s}$ ' % (valstrhigh, valstrlow)
    
            val = df.pct1p7.values[ind]*100.
            valstr = '%1.0f' % val
            outstr = outstr+ '& %s\%% ' % valstr
            if include_uncert:
                val = dfhigh.pct1p7.values[ind]*100.
                valstrhigh = '%1.0f' % val
                val = dflow.pct1p7.values[ind]*100.
                valstrlow = '%1.0f' % val
                outstr = outstr+ '$^{%s}_{%s}$ ' % (valstrhigh, valstrlow)
    
            val = df.pct2p0.values[ind]*100.
            valstr = '%1.0f' % val
            outstr = outstr+ '& %s\%% ' % valstr
            if include_uncert:
                val = dfhigh.pct2p0.values[ind]*100.
                valstrhigh = '%1.0f' % val
                val = dflow.pct2p0.values[ind]*100.
                valstrlow = '%1.0f' % val
                outstr = outstr+ '$^{%s}_{%s}$ ' % (valstrhigh, valstrlow)        
            #outstr = outstr+ '& %i & %i & %i ' % (np.round(df.yr30pct.values[ind]), np.round(df.yr10pct.values[ind]), np.round(df.yr1pct.values[ind]) )
        outstr = outstr+' \\\\ \n'
    outstr = outstr+'	\hline\n'

    # do it for N pixels above 1.5, 1.7, 2.0
    outstr = outstr+'	\multicolumn{10}{c}{Number of 1 km$^2$ reef locations remaining below threshold, out of 773K }\\\\ \n'  # 773261
    outstr = outstr+'	\hline\n'
    for (ind, scenario) in enumerate(scenarios): # reload each file 4 times, once per scenario
        if ind==0: # no GMST for ssp126
            continue
        ind-=1
        outstr = outstr+ 'SSP%s ' % (scenario[3:])
        for DHW_threshold in DHW_thresholds:

            experiment = '%s_%s_%iyrs_DHW%1.1f' % (finerun, 'ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
    
            if do_add_variability:
                experiment = 'SPIKE_%s_%iyrs_DHW%1.1f' % ('ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod        
                # experiment = 'SPIKE2_%s_%iyrs_DHW%1.1f' % ('ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod        
                # experiment = 'SPIKE%im%if_%s_%iyrs_DHW%1.1f' % (spike_size, spike_freq, 'ssp585', return_year, DHW_threshold) # this will be replaced by individual model if we're in single_mod
    
            
            txtfilename=finedir+'temp_%s_%s_nabove.txt' % (runname, experiment)
            print(txtfilename)
            df = pd.read_csv(txtfilename, sep=' ', names=['ssp', 'nabove1p5', 'nabove1p7', 'nabove2p0'])
    
            val = df.nabove1p5.values[ind]
            valstr = '%1.0f' % val
            outstr = outstr+ '& %s' % valstr
    
            val = df.nabove1p7.values[ind]
            valstr = '%1.0f' % val
            outstr = outstr+ '& %s' % valstr
    
            val = df.nabove2p0.values[ind]
            valstr = '%1.0f' % val
            outstr = outstr+ '& %s' % valstr
         
            #outstr = outstr+ '& %i & %i & %i ' % (np.round(df.yr30pct.values[ind]), np.round(df.yr10pct.values[ind]), np.round(df.yr1pct.values[ind]) )
        outstr = outstr+' \\\\ \n'
    outstr = outstr+'	\hline\n'

     
    outstr = outstr+'\end{tabular}\n'
    #outstr = outstr+'\caption{Projected years and global mean surface temperatures after which fewer than the stated percentage of 1 km$^2$ reef locations remain below the thermal thresholds (top), and percentage and number of reefs remaining below threshold beyond the specified GMSTA (bottom). Dashes indicate that the milestone is not reached prior to 2100. Superscripts and subscripts give one standard deviation uncertainty estimates.}\n'
    outstr = outstr+'\label{table:results}\n'
    outstr = outstr+'\end{table}\n'
    print(outstr)