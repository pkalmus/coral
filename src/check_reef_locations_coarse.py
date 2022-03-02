#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

conda activate geo_env

reef_locations: get reef-containing points (don't need to run every time)
regrid_bcdp: Get model outputs using BCDP, regrid to 1 degree, time homogenize and save on the local file system for regrid_reef.py to digest.
reef_locations_coarse: get coarse reef-containing points (don't need to run every time)

regrid_reef: Find reef and neighbor locations (using previously made reef_location_coarse.csv and nn_location_coarse.csv)
             Repackage into 1-D "ravel" .nc dataset with "index".
               This .nc dataset also specifies which coarse cells are reefs.

@author: pkalmus
"""

import xarray as xr
import numpy as np
import pdb
import sys
import subprocess
import coral_plotter as plotter

# read in user params
import importlib
paramfile = sys.argv[1]   
params = importlib.import_module(paramfile)
projectdir = params.projectdir
tosrunname = params.tosrunname
basedir = params.basedir
srcdir = params.srcdir
figdir = params.figdir
modellist_filename = params.modellist_filename
scenarios = params.scenarios
dryRun = params.dryRun
oneMember = params.oneMember


# define in and out dirs
indir = basedir+'data/tos/coral2/%s/raw/' % (tosrunname) # note: "data" is a symlink to /raid8/pkalmus/data/coral/data/
outdir = basedir+'data/tos/coral2/%s/reef/' % (tosrunname)


reef_coarse_file = '/raid8/pkalmus/data/coral/data/location/coarse/reef_location_coarse.csv'    # 1294 [lon lat]
reef_coarse_cells = np.loadtxt(reef_coarse_file)
nn_coarse_file = '/raid8/pkalmus/data/coral/data/location/coarse/nn_location_coarse.csv'        # 1790 [lon lat]
nn_coarse_cells = np.loadtxt(nn_coarse_file)

both_grid = np.vstack([reef_coarse_cells, nn_coarse_cells])
index = np.arange(both_grid.shape[0])
nreef = len(reef_coarse_cells)

plotter.scatter(reef_coarse_cells[:,0], reef_coarse_cells[:,1], figdir+'check_reef_locations_coarse_reef_grid.png', c=np.ones(reef_coarse_cells[:,0].shape), marker_relsize=0.1, figsize=(15,5), cbar_fraction=0.008, cbar_orientation='vertical')
plotter.scatter(nn_coarse_cells[:,0], nn_coarse_cells[:,1], figdir+'check_reef_locations_coarse_nn_grid.png', c=np.ones(nn_coarse_cells[:,0].shape), marker_relsize=0.1, figsize=(15,5), cbar_fraction=0.008, cbar_orientation='vertical')
plotter.scatter(both_grid[:,0], both_grid[:,1], figdir+'check_reef_locations_coarse_both_grid.png', c=np.ones(both_grid[:,0].shape), marker_relsize=0.1, figsize=(15,5), cbar_fraction=0.008, cbar_orientation='vertical')
