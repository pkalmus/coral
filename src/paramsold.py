#!/usr/bin/env /raid2/sport/people/casejl/python/anaconda2/envs/nucaps/bin/python
# -*- coding: utf-8 -*-
"""
conda activate geo_env

Holds the parameters for stages of the coral analysis.

If "Memory Error" make sure geoenv

"""

tosrunname = 'nature' # runname where the tos files are in, and reef tos files.
runname = 'nature'    # where the mean coarse files are, and fine files.
# runname = 'flat'
# runname = 'onemodel'

#projectdir = '/0work/projects/'
projectdir = '/home/pkalmus/projects/'
basedir = projectdir+'/coral/'
pythondir = basedir+'src/python/'

cmip_var = 'tas' # tos
cmip_var = 'tasmax' # tos
cmip_var = 'tos'
listfilename = pythondir+'models_%s_%s.txt' % (cmip_var, runname)
listfilename = 'models_nature.txt'
listfilename = 'models_flat.txt'

# model numbers as of 8/28/2020, for gn only
scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp585'] # 29 models

#scenarios = ['ssp126', 'ssp370', 'ssp585'] 

# scenarios = ['ssp245', 'ssp370', 'ssp585']
# scenarios = ['ssp585']
scenarios = ['ssp370']
#scenarios = ['ssp245']
#scenarios = ['ssp126']

#scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp460', 'ssp585'] # ssp460 cuts models down from 29 to 7, so it's out
# scenarios = ['ssp126', 'ssp245', 'ssp585'] # 33 models
# scenarios = [ 'ssp245', 'ssp585'] # 34 models
# scenarios = ['ssp585'] #35
# scenarios = ['ssp245'] #34
# scenarios = ['ssp126'] #34
# scenarios = ['ssp370'] #31
# scenarios = ['ssp460'] #7

#############################################################
# controls whether debugging/verification plots are made
#############################################################
do_plots = False

#############################################################
# regrid_bcdp only
#############################################################
dryRun = False # only create model lists, not .nc files
oneMember = True # only use the first member with all scenarios per model. False for coral project, true for GMST.
daylistfilename = 'models_daily.txt'
#interpolate_method = 'conservative' # 'bilinear', 'conservative', 'patch', 'nearest_s2d', 'nearest_d2s.'
interpolate_method = 'bilinear' # 'bilinear', 'conservative', 'patch', 'nearest_s2d', 'nearest_d2s.'

grids = ['gn', 'gr', 'gr1']
# grids = ['gr1', 'gr']
# grids = ['gr']
#grids = ['gn', 'gr1'] # note: as of 8/29, there are two gr1 models, and two gr models available, but the gr models also have gn versions.


#############################################################
# regrid_reef only
#############################################################
year_cutoff = 1915 # throw out all years before this
# important that lon/lat geolimit is no smaller than in coral_locations.py.
geolimit=([-35.5,35.5], [0,360]) # NOTE: I don't think geolimiting for lon is implemented yet except for #2. also should be lon, lat.
#geolimit=([-40,-2], [100,170]) # Australia



#############################################################
# weighting.py (and coarse_project)
#
# NB: J0 and R both exist in r/EstimateWeight.R. They must be changed THERE. The R here must match.
#############################################################

# for simple mean model weighting
do_trend = False # for validation, use the trend weighting in weight_trend.py

###############
# climstartyear_bma, climendyear_bma -cimatology to calculate DHW for BMA training period
# meanstart, meanend - BMA training period
###############

# for validation
meanstart=1980
meanend=2000
climstartyear_bma = 1960
climendyear_bma = 1980

# meanstart=1960
# meanend=1980
# climstartyear_bma = 1940
# climendyear_bma = 1960

# for production. NOTE: BMA is trained on mean annual DHW max values, so it does depend on climatology.
# meanstart=1995
# meanend=2015
climstartyear_weights = 1975
climendyear_weights = 1995


#################
# this is climatology for doing coarse projections
#################
# can only do the climatology sensitivity experiment on coarse data b/c the MUR data only goes back to 2002/06
# 15 year periods
# climstartyear_project = 1975
# climendyear_project = 1995

climstartyear_project = 1963
climendyear_project = 1977

climstartyear_project = 1973
climendyear_project = 1987

climstartyear_project = 1983
climendyear_project = 1997

climstartyear_project = 1993
climendyear_project = 2007

climstartyear_project = 2003
climendyear_project = 2017

#############################################################
# coarse_project  fine_project only
#############################################################
return_years = [5,10] # we'll get departure for when there's a DHW event that recurs at least every return_years years
# return_years = [5]
# return_years = [10]
DHW_thresholds = [8,6.423,4.812]
# DHW_threshold = 7
# DHW_threshold = 6
# DHW_threshold = 5
DHW_threshold = 8
# DHW_thresholds = [8]
# DHW_thresholds = [4.812]
# DHW_thresholds = [6.423] # 6.423

bias_correct = False # use the mean annual climatology to bias correct the entire tos timeseries for each reef; makes things look nearly all red.
do_obs = False # use the obs. file, to calculate DHW on it for validation. do_weights, bias_correct don't matter.
min_models = 10 # use for making model ensemble, for coarse_project and also for validate_dhw.

include_downscale_uncert = False # include temp_sd in uncert estimate. might not want to include if it does not exist (Trend) or not trusted (bgl)

do_weights = False
downscale_type = 'bgl' # 'bgl', 'trend', 'lagp'
#downscale_type = 'trend' # 'bgl', 'trend', 'lagp'

include_mur = True # if True, concat the 2002-2019 MUR pixels
do_add_variability = False # add random interannual variability as a test
spike_size = 1 # multiplier for standard normal
spike_freq = 5 # years

#############################################################
# fine_project and downstream (fine_plot, make_table)
#############################################################
special_tag = '' # for no tag, use ''. this allows special experiments
if do_weights:
    finerun = downscale_type + '_' + 'beo' + special_tag
else:
    finerun = downscale_type + '_' + 'flat' + special_tag
    finerun = downscale_type + '_' + 'onepermodel' + special_tag    
 
#############################################################
# fine_plot abnd make_table_EF only
#############################################################    
include_uncert = False # for EF paper: MUR uncert in figures
include_uncert_table = False # but no uncert in tables (best estimates only)
do_slow_plots = False # make the slow global scatter plots

#############################################################
# validate_dhw only
#############################################################
valstartyear = 2005
valendyear = 2015
doval2015 = False # this was to capture the big DHW events from that period. (can't recall why I wanted that)
useValFiles = False # if true, use fine files specifically prepared for val. requires a "val" dir in each finrun output dir
