#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 08:52:34 2022

@author: danilocoutodsouza
"""

import matplotlib.pyplot as plt
import cmocean.cm as cmo
import numpy as np 
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature, COASTLINE
from cartopy.feature import BORDERS
import pandas as pd
from matplotlib import cm
import matplotlib.colors as colors
import sys
import glob
import os
import xarray as xr

def map_features(ax):
    ax.add_feature(COASTLINE)
    ax.add_feature(BORDERS, edgecolor='#383838')
    return ax

def Brazil_states(ax):    
    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                  name='admin_1_states_provinces_lines')
    _ = ax.add_feature(states, edgecolor='#383838')
    
    cities = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                  name='populated_places')
    _ = ax.add_feature(cities)
    
def grid_labels_params(ax,i):
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5,linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    if i not in [0,3]:
        gl.left_labels = False
    gl.xlabel_style = {'size': 14, 'color': '#383838'}
    gl.ylabel_style = {'size': 14, 'color': '#383838'}
    ax.spines['geo'].set_edgecolor('#383838')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return ax

def map_variable(VariableData,figures_subdirectory,fname):
    """
    Parameters
    ----------
    VariableData : DataObj
        Object containing meteorological data from a NetCDF file
    figures_subdirectory : str
        Directory where images will be saved.
    fname : str
        Name to append to outfile.

    Returns
    -------
    Create maps for each time step for model levels closest to 1000,850,700,
    500,300 and 200 hPa.
    """
    # limits for the map (represents the maximum limits of the lagrangian box)
    westernlimit = VariableData[LonIndexer].min()
    easternlimit = VariableData[LonIndexer].max()
    southernlimit = VariableData[LatIndexer].min()
    northernlimit = VariableData[LatIndexer].max()
    # projection
    proj = ccrs.PlateCarree() 
    # create figure
    plt.close('all')
    fig = plt.figure(constrained_layout=False,figsize=(24,18))
    gs = gridspec.GridSpec(2, 3, hspace=0, wspace=0.1,
                                   left=0.05, right=0.95)
    # Find vertical levels that closely matches the desired levels for plotting
    MatchingLevels = []
    for pres in [1000,850,700,500,300,200]:
        MatchingLevels.append(min(VariableData[LevelIndexer].values, key=lambda x:abs(x-pres)))
    # loop through pressure levels
    for p,i in zip(MatchingLevels,range(len(MatchingLevels))):
        # Get current time and time strings
        # create subplot
        ax = fig.add_subplot(gs[i], projection=proj)
        ax.set_extent([westernlimit,easternlimit,southernlimit,northernlimit]) 
        # Add decorators and Brazil states
        grid_labels_params(ax,i)
        Brazil_states(ax)
        iData = VariableData.sel({LevelIndexer:p})
        # get latitude and longitude
        lon,lat = iData[LonIndexer], iData[LatIndexer]
        # get data range
        max1,min1 = float(np.amax(iData)), float(np.amin(iData))
        # (if data is an annomaly)
        if min1 > 0:
            cmap = cmo.amp
            norm = cm.colors.Normalize(vmax=max1,vmin=min1)
        # (if data is not an annomaly)
        else:
            cmap = cmo.balance
            norm = colors.TwoSlopeNorm(vmin=min1, vcenter=0, vmax=max1)
        # plot contours
        cf1 = ax.contourf(lon, lat, iData, cmap=cmap,norm=norm) 
        ax.contour(lon, lat, iData, cf1.levels,colors='#383838',
                   linewidths=0.25)
        # Title
        title = VariableData.name
        title += ' at '+str(p)+' hPa'
        ax.text(0.01,1.01,title,
                transform=ax.transAxes, fontsize=16)
        # colorbar
        cbar = plt.colorbar(cf1,fraction=0.046,pad=0.07,orientation='horizontal',
                            norm=norm)
        cbar.ax.tick_params(labelsize=10) 
        for t in cbar.ax.get_yticklabels():
             t.set_fontsize(10) 
        # decorators
        map_features(ax)
    # save file
    outfile = figures_subdirectory+'/'+fname
    plt.savefig(outfile, bbox_inches='tight')
    print(outfile+' created!')
    
if __name__ == "__main__":

    results_subdirectory = sys.argv[1]
    FigsSubDirectory = results_subdirectory+'/Figures/maps/'; os.makedirs(
        FigsSubDirectory, exist_ok=True)
    
    results_data = glob.glob(results_subdirectory+'/*.nc')[0]
    ds = xr.open_dataset(results_data, engine='netcdf4')
    
    periods = pd.read_csv('../inputs/periods',sep= ';',header=0)
    dfVars = pd.read_csv('../inputs/fvars',sep= ';',index_col=0,header=0)
    LonIndexer,LatIndexer,TimeIndexer,LevelIndexer = \
      dfVars.loc['Longitude']['Variable'],dfVars.loc['Latitude']['Variable'],\
      dfVars.loc['Time']['Variable'],dfVars.loc['Vertical Level']['Variable']
    
    for var in ['AdvHTemp','Sigma', 'Omega','dTdt','ResT']:
        
        for i in range(len(periods)):
            start,end = periods.iloc[i]['start'],periods.iloc[i]['end']
            period =  periods.iloc[i]['Period']
            dates = pd.to_datetime(ds[var][TimeIndexer])
            # Get data for selected periods
            data_period = ds[var].sel({TimeIndexer:slice(start,end)}).mean(
                dim=TimeIndexer)
            try:
                map_variable(data_period,FigsSubDirectory,var+'_'+period)
            except:
                print('issues when plotting periods, check iputs/periods')