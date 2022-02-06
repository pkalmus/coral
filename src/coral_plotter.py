#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Plotting functions (with Cartopy) for coral project.

Created on Sun Feb  2 17:05:22 2020

@author: pkalmus
"""

import pdb

#from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import cartopy.feature
import cartopy.crs as ccrs
from sklearn.metrics import mean_squared_error
from statsmodels.tsa.seasonal import seasonal_decompose


# from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
# # https://stackoverflow.com/questions/49956355/adding-gridlines-using-cartopy
# import matplotlib.ticker as mticker
# from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


import cartopy.mpl.ticker as cticker


from numpy.polynomial import polynomial as P

# from scipy.interpolate import interp1d
# import tauc_p_cluster


from scipy.stats import linregress
from scipy.optimize import curve_fit
# import fittools
# import logging
# import traceback


def plot_data_array(da, filename):
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180.0))
    da.plot(ax=ax, transform=ccrs.PlateCarree(), add_colorbar=True, cmap='ocean')
    ax.coastlines(resolution='50m') #Currently can be one of “110m”, “50m”, and “10m”.
    ax.add_feature(cartopy.feature.LAND,facecolor='white') 
    
    
    #plt.savefig('reef_mask%i.pdf' % (lat_extent)) #PDF is MUCH slower than png
    #plt.savefig('reef_mask_coralsea.pdf') #PDF is MUCH slower than png
    plt.savefig(filename) #PDF is MUCH slower than png
    print(filename)
    plt.close()

def image(lat, lon, data, filename, vmin=None, vmax=None, extent=None, figsize=None, cbar_orientation='vertical', cbar_fraction=None, marker_relsize=1, cmap='jet', draw_labels=True):
    '''
    For when lon, lat are gridded.

    Parameters
    ----------
    lon : TYPE
        DESCRIPTION.
    lat : TYPE
        DESCRIPTION.
    data : TYPE
        DESCRIPTION.
    filename : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    if figsize is None:
        figsize=(15,5)
    plt.figure(figsize=figsize)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180.0))
    ax.add_feature(cartopy.feature.LAND)
    ax.coastlines(resolution='10m') #Currently can be one of “110m”, “50m”, and “10m”.    

    sc = plt.pcolormesh(lat, lon, data, cmap=cmap, edgecolor='', transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax) #, cmap=maskcmap, vmin=allmin, vmax=allmax, edgecolor='') 
    if cbar_fraction is not None:        
        plt.colorbar(sc, extend='max', fraction=cbar_fraction, aspect=30, pad=0.05, orientation=cbar_orientation) #fix the annoying way the colorbar is way too tall
    else:
        plt.colorbar(sc, extend='max', aspect=30, pad=0.05, orientation=cbar_orientation)


    gl = ax.gridlines(draw_labels=draw_labels)
    gl.ylabels_right = False
    gl.xlabels_top = False

    if (extent is not None):
        ax.set_extent(extent, crs=ccrs.PlateCarree()) # use this to zoom into small regions. NOTE: MUST COME **AFTER** gridlines and labels.

    plt.savefig(filename)
    print(filename)
    plt.close()
        
    
def scatter(lon, lat, filename, c=None, vmin=None, vmax=None, extent=None, figsize=None, cbar=True, cbar_orientation='horizontal', cbar_fraction=None, marker_relsize=1, cmap='jet_r', draw_labels=True, ticksize=None, fontsize=None):
    '''
    

    Parameters
    ----------
    lon : TYPE
        DESCRIPTION.
    lat : TYPE
        DESCRIPTION.
    filename : TYPE
        DESCRIPTION.
    c : TYPE, optional
        DESCRIPTION. The default is None.
    vmin : TYPE, optional
        DESCRIPTION. The default is None.
    vmax : TYPE, optional
        DESCRIPTION. The default is None.
    extent : TYPE, optional
        DESCRIPTION. The default is None. (min(lon), max(lon), min(lat), max(lat)) or [min(lon), max(lon), min(lat), max(lat)]
    figsize : TYPE, optional
        DESCRIPTION. The default is None.
    cbar_orientation : TYPE, optional
        DESCRIPTION. The default is 'horizontal'.
    cbar_fraction : TYPE, optional
        DESCRIPTION. The default is None.
    marker_relsize : TYPE, optional
        DESCRIPTION. The default is 1.
    cmap : TYPE, optional
        DESCRIPTION. The default is 'jet_r'.
    draw_labels : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None.

    '''
    if figsize is None:
        figsize=(15,5)
    plt.figure(figsize=figsize)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180.0))
    ax.add_feature(cartopy.feature.LAND)
    ax.coastlines(resolution='10m') #Currently can be one of “110m”, “50m”, and “10m”.    
    sc = plt.scatter(lon, lat, c=c, s=11, cmap=cmap, edgecolor='', transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax) #, cmap=maskcmap, vmin=allmin, vmax=allmax, edgecolor='') 
    if cbar:
        if cbar_fraction is not None:        
            plt.colorbar(sc, extend='max', fraction=cbar_fraction, aspect=30, pad=0.05, orientation=cbar_orientation) #fix the annoying way the colorbar is way too tall
        else:
            plt.colorbar(sc, extend='max', aspect=30, pad=0.05, orientation=cbar_orientation)
    sc.set_sizes(sc.get_sizes()*marker_relsize)

    # better labels
    if ticksize is None:
        ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
        ax.set_yticks([-30, -15, 0, 15, 30], crs=ccrs.PlateCarree())
    else:
        # must also have extent set
        if extent is None:
            print('warning, cannot set ticksize labels when extent is None; please set extent')
        else:
            ax.set_xticks(np.arange(extent[0], extent[1], ticksize), crs=ccrs.PlateCarree())
            ax.set_yticks(np.arange(extent[2], extent[3], ticksize), crs=ccrs.PlateCarree())        
    lon_formatter = cticker.LongitudeFormatter(zero_direction_label=True)
    lat_formatter = cticker.LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    if (extent is not None): 
        ax.set_extent(extent, crs=ccrs.PlateCarree()) # use this to zoom into small regions. NOTE: MUST COME **AFTER** gridlines and labels.

    if fontsize is not None:
        plt.rcParams.update({'font.size': fontsize})
    ax.grid(linewidth=1, color='gray', alpha=0.6)

    plt.savefig(filename)
    print(filename)
    plt.close()
    
def timeseries(t1, filename, t2=None, xlabel='time', ylabel='tos', title=None,  xlim=None, ylim=None, savefiles=False):
    '''
    Plot one, or two, timeseries. (quick and dirty)
    
    t1, t2 is Nx2 array with [time, val] points.
    
    '''
    fig, ax = plt.subplots()
    ax.plot(t1[:,0], t1[:,1], color='k')
    if t2 is not None:
        ax.plot(t2[:,0], t2[:,1], color='r')

    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.grid()
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.savefig('/home/pkalmus/projects/coral/src/python/'+filename)
    #pdb.set_trace()
    plt.close()
    if savefiles:
        np.savetxt(filename+'_obs.txt', t1, fmt='%1.2f')
        np.savetxt(filename+'_mod.txt', t2, fmt='%1.2f')
    print(filename)

def plot(x, y, filename, x2=None, y2=None, xlabel='time', ylabel='tos', title='',  xlim=None, ylim=None, points_only=False, linfit=False):
    '''
    Same as timeseries but takes x,y separately (and no savefiles)
    
    t1, t2 is Nx2 array with [time, val] points.
    
    '''
    cc = None
    fig, ax = plt.subplots()
    
    if points_only:
        ax.plot(x, y, 'k*')
        if x2 is not None:
            ax.plot(x2, y2, 'r*')
    else:
        ax.plot(x, y, color='k')
        if x2 is not None:
            ax.plot(x2, y2, color='r')

    # The solution is the coefficients of the polynomial p that minimizes the sum of the weighted squared errors (least squares polynomial fit)
    # https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyfit.html     
    if linfit: # only implemented for x,y; not x2, y2.
        cc = P.polyfit(x,y,1)
        ax.plot(x, cc[1]*np.array(x)+cc[0], 'b')
        ax.text(0.8, 0.1, 'y=%1.2fx+%1.1f' % (cc[1], cc[0]), horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
    
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.grid()
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.savefig(filename)
    plt.close()
    print(filename)
    return cc
    

   
def plot_all_timeseries(all_ds, mean_ds, filename, xlim=None, xlabel='time (s)', ylabel='tos'):
    '''
    '''
    nmodels = len(all_ds.model_index.values)
    print(nmodels)
    fig, ax = plt.subplots()
    
    for model_ind in np.arange(nmodels):
        mod_ds = all_ds.sel(model_index=model_ind)
        ax.plot(mod_ds.time.values, mod_ds.tos.values, color='r')
        plt.xlim(xlim)
        #plt.show()
    ax.plot(mean_ds.time.values, mean_ds.tos.values, color='k')
    #plt.show()

    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.grid()
    fig.savefig(filename)
    print(filename)
    plt.close()
    

    
def plot_departure(years, tos, mmm, departure_year, filename):
    fig, ax = plt.subplots()
    ax.plot(years, tos, color='k')
    ax.axhline(mmm, color='r')
    ax.axvline(departure_year, color='r')
    plt.xlim(departure_year-25, departure_year+25)
    plt.ylim(14,30)
    plt.savefig(filename)
    #plt.show()
    print(filename)
    plt.close()
  
def hist(data, xlabel, filename, xlim=None, ylim=None, logy=False, logx=False, bins=20, title='', label='stats', fontsize=None, cumulative=0, normed=False, histtype='bar'):
    fig, ax = plt.subplots()
    data = data[~np.isnan(data)]
    if label=='stats':
        label='mean:%1.1f\nstd:%1.1f\nrmse:%1.1f\ncount:%i' % (np.nanmean(data), np.nanstd(data), np.sqrt(mean_squared_error(data, np.zeros(data.shape))), np.count_nonzero(~np.isnan(data)) )
    ax.hist(data, log=logy, bins=bins, label=label, cumulative=cumulative, normed=normed, histtype=histtype)
    if logx:
        plt.xscale('log')
    plt.xlim(xlim)
    plt.ylim(ylim)
        
    plt.legend(loc="upper right")          
    ax.set(xlabel=xlabel, title=title )

    if fontsize is not None:
        plt.rcParams.update({'font.size': fontsize})
    ax.grid()
    plt.savefig(filename)
    #plt.show()
    print(filename)
    plt.close()

def fix_hist_step_vertical_line_at_end(ax):
    ''' helper fctn for histlist.
    https://stackoverflow.com/questions/39728723/vertical-line-at-the-end-of-a-cdf-histogram-using-matplotlib
    '''
    axpolygons = [poly for poly in ax.get_children() if isinstance(poly, matplotlib.patches.Polygon)]
    for poly in axpolygons:
        poly.set_xy(poly.get_xy()[1:])
      
def histlist(dslist, xlabel, filename, txtfilename=None, anomlist=None, ylabel=None, xlim=None, ylim=None, logy=False, logx=False, dolines=False, bins=100, title='', label=None, fontsize=None, axvspan=None, colors=None, cumulative=0, normed=False, histtype='bar', include_uncert=True):
    fig, ax = plt.subplots()
    
    if colors is None:
        colors = ['k', 'b', 'g', 'r']
        
    if len(dslist)==4: # for the txt file
        scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
    else:
        scenarios = ['ssp245', 'ssp370', 'ssp585']
    
    pct1p5str = ''    
    pct2p0str = ''
    pct15vect = []
    pct15Hvect = []
    pct15Lvect = []
    pct2vect = []
    pct2Hvect = []
    pct2Lvect = []
    for i in np.arange(len(dslist)):
        ds = dslist[i]
        data = ds.departure_year.values
        data = data[~np.isnan(data)]
        
        anom2021 = 0.0
        if anomlist is not None:
            anomds = anomlist[i]
            anom = anomds.anom.values
            N=120 # 10 years running mean to stabilize predictions
            anom = np.convolve(anom, np.ones((N,))/N, mode='same') 
            anomyear = anomds.time.values
            anomdata = interp1d(anomyear, anom, bounds_error=False)(data)
            data = anomdata
            anom2021 = interp1d(anomyear, anom, bounds_error=False)(2021)
            tempdata = np.copy(anomdata) # for use below in determining number of reef locations post 1.5 and post 2.0
            
        try:    
            y, bins, ig = ax.hist(data, linewidth=2, color=colors[i], log=logy, bins=bins, label=label, cumulative=cumulative, normed=normed, histtype=histtype)
        except:
            pdb.set_trace()
        x = (bins[:-1] + bins[1:])/2

        if include_uncert:
            # high unertainty
            data = ds.departure_year_high.values
            data = data[~np.isnan(data)]
            if anomlist is not None:
                anomds = anomlist[i]
                anom = anomds.anom.values
                anomyear = anomds.time.values
                anomdata = interp1d(anomyear, anom, bounds_error=False)(data)
                data = anomdata
            yhigh, bins, ig = ax.hist(data, linewidth=0.1, linestyle='-', color=colors[i], log=logy, bins=bins, label=label, cumulative=cumulative, normed=normed, histtype=histtype)

            # low uncertainty
            data = ds.departure_year_low.values
            data = data[~np.isnan(data)]
            
            if anomlist is not None: # convert from year to temperature anomaly
                anomds = anomlist[i]
                anom = anomds.anom.values
                anomyear = anomds.time.values 
                anomdata = interp1d(anomyear, anom, bounds_error=False)(data)
                data = anomdata
            ylow, bins, ig = ax.hist(data, linewidth=0.1, linestyle='-', color=colors[i], log=logy, bins=bins, label=label, cumulative=cumulative, normed=normed, histtype=histtype)
            ax.fill_between(x, ylow, yhigh, where=None, color=colors[i], alpha=0.4)
        else:
            yhigh = y
            ylow = y
            
        fix_hist_step_vertical_line_at_end(ax)
        if dolines: # vertical and horiz lines
            # add 1% line and estimate numbers
            ax.axhline(y=0.01, color='m', linestyle='-')
            ax.axhline(y=0.1, color='c', linestyle='-')
            # find year of intersection
            yr1pct = interp1d(y, x, bounds_error=False)(0.01)
            if ~np.isnan(yr1pct):
                ax.axvline(x=yr1pct, ymin=0, ymax=0.1, color=colors[i], linewidth=3) # this is not drawing to 0.1 on the axis; it is 0.1 fraction of the plot height.
            yr10pct = interp1d(y, x, bounds_error=False)(0.1)
            if ~np.isnan(yr10pct):
                ax.axvline(x=yr10pct, ymin=0, ymax=0.1, color=colors[i], linewidth=3)
            yr30pct = interp1d(y, x, bounds_error=False)(0.3)
            pct2021 = interp1d(x, y, bounds_error=False)(2021)
            # fraction departed at 2020; year of 30% departure; 10% ; 1%. 
            print('******\n%s %s\npct2021: %1.4f\nanom2021: %1.2f\nyr30pct: %1.2f\nyr10pct: %1.2f\nyr1pct%1.2f\n******\n\n' % (title, scenarios[i], pct2021, anom2021, yr30pct, yr10pct, yr1pct))
            
            
            if txtfilename is not None:
                # print to a text file for creating the table in make_table.py
                #print(txtfilename)
                if i == 0:
                    writefile = open(txtfilename, 'w')
                else:
                    writefile = open(txtfilename, 'a')
                writefile.write('%s %1.2f %1.2f %1.2f %1.2f\n' % (scenarios[i], pct2021, yr30pct, yr10pct, yr1pct)) 
                
                if include_uncert:
                    # make separate file for high uncert to make table
                    yr1pct = interp1d(yhigh, x, bounds_error=False)(0.01)
                    yr10pct = interp1d(yhigh, x, bounds_error=False)(0.1)
                    yr30pct = interp1d(yhigh, x, bounds_error=False)(0.3)
                    pct2021 = interp1d(x, yhigh, bounds_error=False)(2021)
                    txtfilenamehigh = txtfilename.replace('.txt', '_high.txt')
                    #print(txtfilenamehigh)
                    if i == 0:
                        writefile = open(txtfilenamehigh, 'w')
                    else:
                        writefile = open(txtfilenamehigh, 'a')
                    writefile.write('%s %1.2f %1.2f %1.2f %1.2f\n' % (scenarios[i], pct2021, yr30pct, yr10pct, yr1pct))
    
                    # make separate file for low uncert to make table
                    yr1pct = interp1d(ylow, x, bounds_error=False)(0.01)
                    yr10pct = interp1d(ylow, x, bounds_error=False)(0.1)
                    yr30pct = interp1d(ylow, x, bounds_error=False)(0.3)
                    pct2021 = interp1d(x, ylow, bounds_error=False)(2021)
                    txtfilenamelow = txtfilename.replace('.txt', '_low.txt')
                    #print(txtfilenamelow)
                    if i == 0:
                        writefile = open(txtfilenamelow, 'w')
                    else:
                        writefile = open(txtfilenamelow, 'a')
                    writefile.write('%s %1.2f %1.2f %1.2f %1.2f\n' % (scenarios[i], pct2021, yr30pct, yr10pct, yr1pct))

                # this is for more lines in table
                if anomlist is not None:
                    # find fraction left at 1.5, 1.7, 2.0
                    pct1p5 = interp1d(x, y, bounds_error=False)(1.5)     
                    pct1p7 = interp1d(x, y, bounds_error=False)(1.7) 
                    pct2p0 = interp1d(x, y, bounds_error=False)(2.0)   
                    pct1p5h = interp1d(x, yhigh, bounds_error=False)(1.5)  
                    pct1p7h = interp1d(x, yhigh, bounds_error=False)(1.7) 
                    pct2p0h = interp1d(x, yhigh, bounds_error=False)(2.0) 
                    pct1p5l = interp1d(x, ylow, bounds_error=False)(1.5)  
                    pct1p7l = interp1d(x, ylow, bounds_error=False)(1.7) 
                    pct2p0l = interp1d(x, ylow, bounds_error=False)(2.0) 
                    pct1p5str = pct1p5str + '**************%s 1.5: %1.4f (%1.4e %1.4f)\n' % (scenarios[i], pct1p5, pct1p5l, pct1p5h)
                    pct2p0str = pct2p0str + '**************%s 2.0: %1.4f (%1.4e %1.4f)\n' % (scenarios[i], pct2p0, pct2p0l, pct2p0h)
                    pct15vect.append(pct1p5)
                    pct15Hvect.append(pct1p5h)
                    pct15Lvect.append(pct1p5l)
                    pct2vect.append(pct2p0)
                    pct2Hvect.append(pct2p0h)
                    pct2Lvect.append(pct2p0l)
                    
                    # get num locations 
                    nabove1p5 = len(np.where(tempdata > 1.5)[0])
                    nabove1p7 = len(np.where(tempdata > 1.7)[0])
                    nabove2p0 = len(np.where(tempdata > 2.0)[0])
                    txtfilenamenabove = txtfilename.replace('.txt', '_nabove.txt')
                    #print(txtfilenamenabove)
                    if i == 0:
                        writefile = open(txtfilenamenabove, 'w')
                    else:
                        writefile = open(txtfilenamenabove, 'a')
                    writefile.write('%s %i %i %i\n' % (scenarios[i], nabove1p5, nabove1p7, nabove2p0))
                    
                    
                    # get the intersections with the vertical 1.5 and 2.0 anom lines
                    txtfilename1p52p0 = txtfilename.replace('.txt', '_1p51p72p0.txt')
                    #print(txtfilename1p52p0)
                    if i == 0:
                        writefile = open(txtfilename1p52p0, 'w')
                    else:
                        writefile = open(txtfilename1p52p0, 'a')
                    writefile.write('%s %1.4f %1.4f %1.4f\n' % (scenarios[i], pct1p5, pct1p7, pct2p0))

                    txtfilename1p52p0high = txtfilename.replace('.txt', '_1p51p72p0high.txt')
                    #print(txtfilename1p52p0high)
                    if i == 0:
                        writefile = open(txtfilename1p52p0high, 'w')
                    else:
                        writefile = open(txtfilename1p52p0high, 'a')
                    writefile.write('%s %1.4f %1.4f %1.4f\n' % (scenarios[i], pct1p5h, pct1p7h, pct2p0h))
                    
                    txtfilename1p52p0low = txtfilename.replace('.txt', '_1p51p72p0low.txt')
                    #print(txtfilename1p52p0low)
                    if i == 0:
                        writefile = open(txtfilename1p52p0low, 'w')
                    else:
                        writefile = open(txtfilename1p52p0low, 'a')
                    writefile.write('%s %1.4f %1.4f %1.4f\n' % (scenarios[i], pct1p5l, pct1p7l, pct2p0l))                    


    if axvspan is not None:
        ax.axvspan(axvspan[0], axvspan[1], alpha=0.3, color='gray')

    if logx:
        plt.xscale('log')
    plt.xlim(xlim)
    plt.ylim(ylim)
        
    plt.legend(loc="upper right")          
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title )

    if fontsize is not None:
        plt.rcParams.update({'font.size': fontsize})
    ax.grid()

    fig.tight_layout()
    plt.rcParams.update({'font.size': fontsize})
    plt.savefig(filename)
    #plt.show()
    print(filename)
    plt.close()
    print(pct1p5str) 
    print('%s 1.5 percent (not fraction): %1.4f (%1.4f %1.4f)\n' % (scenarios[i], np.mean(pct15vect)*100., np.mean(pct15Lvect)*100., np.mean(pct15Hvect)*100. ))
    
    print(pct2p0str) 
    print('%s 2.0 percent (not fraction): %1.4f (%1.4f %1.4f)\n' % (scenarios[i], np.mean(pct2vect)*100., np.mean(pct2Lvect)*100., np.mean(pct2Hvect)*100. ))



def sstlistmodels(modellist, meands, xlabel, filename, ylabel=None, xlim=None, ylim=None, logy=False, logx=False, dolines=False, bins=100, title='', label=None, fontsize=None, axvspan=None, colors=None, cumulative=0, normed=False, histtype='bar'):
    '''
    For one SSP, put in all the models, and the mean.
    
    Dimensions:  (index: 828639, time: 1171)
    Coordinates:
      * time     (time) float64 2.002e+03 2.002e+03 ... 2.099e+03 2.099e+03
      * index    (index) int64 0 1 2 3 4 5 ... 828634 828635 828636 828637 828638
    Data variables:
        temp     (index, time) float32 26.408 24.95 25.491 ... 32.5832 30.326
        lon      (index) float32 130.13 130.14 130.15 ... 129.64 129.65 129.66
        lat      (index) float32 -12.91 -12.91 -12.91 -12.91 ... 0.0 0.0 0.0 0.0
        temp_sd  (index, time) float32 0.0 0.0 0.0 ... 0.204394 0.814594 0.112754
    '''

    fig, ax = plt.subplots()
    
    N=12
    
    ax.set_prop_cycle('color',plt.cm.Spectral(np.linspace(0,1,len(modellist))))
    if modellist is not None:
        for i in np.arange(len(modellist)):
            ds_coarse = modellist[i]
            temp_coarse = ds_coarse.tos.values
            time_coarse = ds_coarse.time.values
            temp_coarse = np.nanmean(temp_coarse, axis=0)
            temp_coarse_smooth = np.convolve(temp_coarse, np.ones(N)/N, mode='valid') # take 12-month running average
            time_coarse_smooth = time_coarse[0:len(temp_coarse_smooth)]   
            ax.plot(time_coarse_smooth, temp_coarse_smooth, linewidth=1, linestyle='-')

            # # doesn't work b/c it folds the interannual variability into the "trend"
            # decompose_result_mult = seasonal_decompose(temp_coarse, period=12)
            # trend = decompose_result_mult.trend
            # seasonal = decompose_result_mult.seasonal
            # residual = decompose_result_mult.resid
            # myfig = decompose_result_mult.plot()
            # filename = 'seasonal_%i.png' % (i)
            # print(filename)
            # myfig.savefig(filename)

    temp_coarse = meands.tos.values
    time_coarse = meands.time.values
    temp_coarse = np.nanmean(temp_coarse, axis=0)
    temp_coarse = np.convolve(temp_coarse, np.ones(N)/N, mode='valid') # take 12-month running average
    time_coarse = time_coarse[0:len(temp_coarse)]   
    ax.plot(time_coarse, temp_coarse, 'r', linewidth=3, linestyle='-')


    if axvspan is not None:
        ax.axvspan(axvspan[0], axvspan[1], alpha=0.3, color='gray')

    if logx:
        plt.xscale('log')
    plt.xlim(xlim)
    plt.ylim(ylim)
        
    plt.legend(loc="upper left")          
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title )

    if fontsize is not None:
        plt.rcParams.update({'font.size': fontsize})
    ax.grid()

    fig.tight_layout()
    plt.rcParams.update({'font.size': fontsize})
    plt.savefig(filename)
    #plt.show()
    print(filename)
    plt.close()



def sstlist(dslist, coarselist, obsds, xlabel, filename, txtfilename=None, anomlist=None, ylabel=None, xlim=None, ylim=None, logy=False, logx=False, dolines=False, bins=100, title='', label=None, fontsize=None, axvspan=None, colors=None, cumulative=0, normed=False, histtype='bar'):
    '''
    Dimensions:  (index: 828639, time: 1171)
    Coordinates:
      * time     (time) float64 2.002e+03 2.002e+03 ... 2.099e+03 2.099e+03
      * index    (index) int64 0 1 2 3 4 5 ... 828634 828635 828636 828637 828638
    Data variables:
        temp     (index, time) float32 26.408 24.95 25.491 ... 32.5832 30.326
        lon      (index) float32 130.13 130.14 130.15 ... 129.64 129.65 129.66
        lat      (index) float32 -12.91 -12.91 -12.91 -12.91 ... 0.0 0.0 0.0 0.0
        temp_sd  (index, time) float32 0.0 0.0 0.0 ... 0.204394 0.814594 0.112754
    '''

    fig, ax = plt.subplots()
    
    N=12
            
    if colors is None:
        colors = ['k', 'b', 'g', 'r']
        
    # if len(dslist)==4: # for the txt file
    #     scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
    # else:
    #     scenarios = ['ssp245', 'ssp370', 'ssp585']
        
    if dslist is not None:
        for i in np.arange(len(dslist)):
            ds = dslist[i]
            temp = ds.temp.values
            time = ds.time.values
            time = np.floor(time) + ((time - np.floor(time))*30*100 - 15)/365. # to fix Ayesha's weird non-decimal time convention
            temp = np.nanmean(temp, axis=0)
            temp = np.convolve(temp, np.ones(N)/N, mode='valid') # take 12-month running average
            time = time[0:len(temp)]        
            ax.plot(time, temp, color=colors[i])

    if coarselist is not None:
        for i in np.arange(len(coarselist)):
            ds_coarse = coarselist[i]
            temp_coarse = ds_coarse.tos.values
            time_coarse = ds_coarse.time.values
            temp_coarse = np.nanmean(temp_coarse, axis=0)
            temp_coarse = np.convolve(temp_coarse, np.ones(N)/N, mode='valid') # take 12-month running average
            time_coarse = time_coarse[0:len(temp_coarse)]   
            ax.plot(time_coarse, temp_coarse, color=colors[i], linewidth=3, linestyle='-')

    if obsds is not None: #this is HadISST.
        #######
        # add hadisst from 2010-2020
        # Dimensions:  (index: 4836, time: 1254)
        # Coordinates:
        #   * index    (index) int64 0 1 2 3 4 5 6 ... 4829 4830 4831 4832 4833 4834 4835
        #   * time     (time) float64 1.915e+03 1.915e+03 ... 2.019e+03 2.019e+03
        # Data variables:
        #     tos      (index, time) float32 ...
        #     isreef   (index) float64 ...
        #     year     (time) float64 ...
        #     lon      (index) float32 ...
        #     lat      (index) float32 ...
        temp_obs = obsds.tos.values
        time_obs = obsds.time.values
        temp_obs = np.nanmean(temp_obs, axis=0)
        temp_obs = np.convolve(temp_obs, np.ones(N)/N, mode='valid') # take 12-month running average
        time_obs = time_obs[0:len(temp_obs)] 
        ax.plot(time_obs, temp_obs, color='k', linewidth=4, linestyle='-')

    if axvspan is not None:
        ax.axvspan(axvspan[0], axvspan[1], alpha=0.3, color='gray')

    if logx:
        plt.xscale('log')
    plt.xlim(xlim)
    plt.ylim(ylim)
        
    plt.legend(loc="upper left")          
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title )

    if fontsize is not None:
        plt.rcParams.update({'font.size': fontsize})
    ax.grid()

    fig.tight_layout()
    plt.rcParams.update({'font.size': fontsize})
    plt.savefig(filename)
    #plt.show()
    print(filename)
    plt.close()

def gmstalist(dslist, xlabel, filename, txtfilename=None, ylabel=None, xlim=None, ylim=None, logy=False, logx=False, dolines=False, bins=100, title='', label=None, fontsize=None, axvspan=None, colors=None):
    '''
Dimensions:  (time: 3012)
Coordinates:
    lev      float64 ...
  * time     (time) float64 1.85e+03 1.85e+03 1.85e+03 ... 2.101e+03 2.101e+03
Data variables:
    tas      (time) float64 ...
    tas60    (time) float64 ...
    anom     (time) float64 ...Dimensions:  (time: 3012)

    '''

    fig, ax = plt.subplots()
    
    if colors is None:
        colors = ['k', 'b', 'g', 'r']
        lables = ['SSP126', 'SSP245', 'SSP370', 'SSP585']
        
        
    for i in np.arange(len(dslist)):        
        ds = dslist[i]        
        temp = ds.anom.values
        time = ds.time.values
        #N=12
        # temp = np.nanmean(temp, axis=0)
        # temp = np.convolve(temp, np.ones(N)/N, mode='valid') # take 12-month running average
        # time = time[0:len(temp)]

        ax.plot(time, temp, color=colors[i], label=lables[i])

    if axvspan is not None:
        ax.axvspan(axvspan[0], axvspan[1], alpha=0.3, color='gray')

    if logx:
        plt.xscale('log')
    plt.xlim(xlim)
    plt.ylim(ylim)
        
    plt.legend(loc="upper left")          
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title )

    if fontsize is not None:
        plt.rcParams.update({'font.size': fontsize})
    ax.grid()

    fig.tight_layout()
    plt.rcParams.update({'font.size': fontsize})
    plt.savefig(filename)
    #plt.show()
    print(filename)
    plt.close()



def bar(x, data, xlabel, filename, xlim=None, ylim=None, logy=False, logx=False, bins=None):
    fig, ax = plt.subplots()
    ax.bar(x, data)
    if logx:
        plt.xscale('log')
    plt.xlim(xlim)
    plt.ylim(ylim)
    ax.set(xlabel=xlabel, title='mean:%1.2f sigma:%1.2f' % (np.nanmean(data), np.nanstd(data)) )
    plt.savefig(filename)
    #plt.show()
    print(filename)
    plt.close()


def corrPlot(x,y,xlabel,ylabel,savestr,y2=None,maxxy=None,dofit=False,nbins=0,
             docloudsatfit=False,cloudsatx0=None,ylim=None,dology=False,label=None,
             label2=None,write_r2=True,showpoints=True):
    ''' Simple correlation plot with rsq in title.
        If nbins>0, will bin up data according to x.
    '''
    
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    # r_value
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    if dofit and not docloudsatfit:
        #pfit = np.polyfit(x, y, 1)
        ax.plot(x, slope*x+intercept, 'b-')
    
    if showpoints:
        ax.plot(x, y, 'b.', alpha=0.2)
        if y2 is not None:
            ax.plot(x, y2, 'r.', alpha=0.2)
    
    centers = None
    bindata = None
    if nbins > 0:
        bins = np.linspace(np.nanmin(x), np.nanmax(x), nbins)
        centersorig = (bins[:-1] + bins[1:]) / 2
        bininds = np.digitize(x, bins)
        bindata = np.array([np.nanmean(y[bininds == i], axis=0) for i in range(1, len(bins))])
        binstd = np.array([np.nanstd(y[bininds == i], axis=0) for i in range(1, len(bins))])
        goodind = ~np.isnan(bindata)
        bindata = bindata[goodind]
        binstd = binstd[goodind]
        centers = centersorig[goodind]
        if y2 is not None:
            bindata2 = np.array([np.nanmean(y2[bininds == i], axis=0) for i in range(1, len(bins))])
            binstd2 = np.array([np.nanstd(y2[bininds == i], axis=0) for i in range(1, len(bins))])
            goodind2 = ~np.isnan(bindata2)
            bindata2 = bindata2[goodind2]
            binstd2 = binstd2[goodind2]
            centers2 = centersorig[goodind2]
            
        # r_value
        slope, intercept, r_value, p_value, std_err = linregress(centers, bindata)
        if dofit:
            if docloudsatfit:
#                def fitmodel(x, *p): 
#                    return p[0] + p[1]*(x-p[2])**2
                #f = lambda x, *p: p[0] + p[1]*(x-p[2])**2
                #popt, pcov = curve_fit(fitmodel, centers, bindata, [300, 1.3, 35])

                def fitmodel(x, *x0): 
                    return x0[0] * np.exp(x0[1] * x) + x0[2]
                def fitmodel_cutoff(x, *x0): 
                    #return p[0] * np.exp(p[1] * x) + p[2]  # was like this in original cloudsat_precip submission                  
                    return x0[2] + x0[0]*np.exp(x0[1] * x) * np.exp(-1*(x/90.)**12)

                #pdb.set_trace()
                binstd[binstd==0]=np.inf # otherwise bins with 1 point could pull the fit
                
                popt, pcov = curve_fit(fitmodel, centers, bindata, [0.005, 0.2, 400], sigma=binstd, absolute_sigma=True, maxfev=80000)
                ax.plot(centers, fitmodel(centers, *popt), 'r-')
                if cloudsatx0 is not None:
                    ax.plot(centers, fitmodel_cutoff(centers, *cloudsatx0), 'g-', linewidth=3)
                    ax.text(40, 10, '%1.3e* exp(%1.3e*x) + %1.2f' % (cloudsatx0[0], cloudsatx0[1], cloudsatx0[2]), color='g', size=20)
                else:
                    ax.text(40, 10, '%1.3e* exp(%1.3e*x) + %1.2f' % (popt[0], popt[1], popt[2]), color='r', size=20)
                
#                popt = np.polyfit(centers, bindata, 3)
#                print 'curve fit: ' + str(popt)
#                ax.plot(centers, np.polyval(popt, centers), 'r-')
#                ax.text(40, 4000, '%1.2f*x^3 + %1.2f*x^2 + %1.2f*x + %1.2f' % (popt[0], popt[1], popt[2], popt[3]), color='r')
            else:
                ax.plot(centers, slope*centers+intercept, 'b-')
        ax.errorbar(centers, bindata, fmt='none', yerr=binstd, ecolor='b')
        ax.plot(centers, bindata, 'b*', markersize=10, label=label)
        if y2 is not None:
            ax.errorbar(centers2, bindata2, fmt='none', yerr=binstd2, ecolor='r')
            ax.plot(centers2, bindata2, 'r*', markersize=10, label=label2)            

    if maxxy is not None:
        ax.plot(list(range(0,maxxy)), list(range(0,maxxy)), 'k-')
        ax.set_ylim([0,maxxy])
        ax.set_xlim([0,maxxy])
       
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True)
    
    
    if dofit:
        #ax.set_title(r'r$^2$=%1.2f, m=%1.1e' % (r_value**2, slope))
        ax.set_title('r=%1.2f, m=%1.1e' % (r_value, slope))
    else:
        if write_r2:
            #ax.set_title(r'r$^2$=%1.2f' % r_value**2)
            ax.set_title('r=%1.2f' % r_value)
    if docloudsatfit:
        ax.set_title('')
        ax.set_yscale('log')
    if dology:
        ax.set_yscale('log')
    ax.legend(loc='best')    
    #setFontSize(18, ax)
    plt.tight_layout()
    print(savestr)
    fig1.savefig(savestr) 
    plt.close(fig1)

    return centers, bindata

