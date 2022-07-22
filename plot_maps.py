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


def map_features(ax):
    ax.add_feature(COASTLINE)
    ax.add_feature(BORDERS, edgecolor='white')
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

def LagrangianMaps(VariableData,FigsDirectory,fname):
    """
    Parameters
    ----------
    VariableData : DataObj
        Object containing meteorological data from a NetCDF file
    FigsDirectory : str
        Directory where images will be saved.
    fname : str
        Name to append to outfile.

    Returns
    -------
    Create maps for each time step for model levels closest to 1000,850,700,
    500,300 and 200 hPa.
    """
    dfVars = pd.read_csv('./fvars',sep= ';',index_col=0,header=0)
    LonIndexer,LatIndexer,TimeIndexer,LevelIndexer = \
      dfVars.loc['Longitude']['Variable'],dfVars.loc['Latitude']['Variable'],\
      dfVars.loc['Time']['Variable'],dfVars.loc['Vertical Level']['Variable']
    trackfile = './track'
    track = pd.read_csv(trackfile,parse_dates=[0],delimiter=';',index_col='time')
    # limits for the map (represents the maximum limits of the lagrangian box)
    westernlimit, easternlimit = track['Lon'].min()-7.55,track['Lon'].max()+7.55
    southernlimit, northernlimit = track['Lat'].min()-7.55,track['Lat'].max()+7.55
    # projection
    proj = ccrs.PlateCarree() 
    # create figure
    plt.close('all')
    fig = plt.figure(constrained_layout=False,figsize=(12,10))
    gs = gridspec.GridSpec(2, 3, hspace=0.2, wspace=0.1,
                                   left=0.05, right=0.95)
    # Find vertical levels that closely matches the desired levels for plotting
    MatchingLevels = []
    for pres in [1000,850,700,500,300,200]:
        MatchingLevels.append(min(VariableData[LevelIndexer].values, key=lambda x:abs(x-pres)))
    # loop through pressure levels
    for p,i in zip(MatchingLevels,range(len(MatchingLevels))):
        # Get current time and time strings
        itime = VariableData[TimeIndexer].values
        min_lon, max_lon = track.loc[itime]['Lon']-7.5,track.loc[itime]['Lon']+7.5
        min_lat, max_lat = track.loc[itime]['Lat']-7.5,track.loc[itime]['Lat']+7.5
        # create subplot
        ax = fig.add_subplot(gs[i], projection=proj)
        ax.set_extent([westernlimit,easternlimit,southernlimit,northernlimit]) 
        # Add decorators and Brazil states
        grid_labels_params(ax,i)
        Brazil_states(ax)
        # Slice data for the desired domain and pressure level
        iData = VariableData.sel({LevelIndexer:p}).sel(
            **{LatIndexer:slice(max_lat,min_lat),
               LonIndexer: slice(min_lon,max_lon)})
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
        # plot the selected box around the contours
        ax.plot([min_lon,min_lon,max_lon,max_lon,min_lon],
                [min_lat,max_lat,max_lat,min_lat,min_lat],
                linewidth=1,c='#383838',zorder=500)
        # plot contours
        cf1 = ax.contourf(lon, lat, iData, cmap=cmap,norm=norm,extend='both') 
        ax.contour(lon, lat, iData, cf1.levels,colors='#383838',
                   linewidths=0.25)
        # get time string
        timestr = pd.to_datetime(str(iData[TimeIndexer].values))
        date = timestr.strftime('%Y-%m-%dT%H%MZ')
        ax.text(0.05,1.01,str(p)+' '+str(VariableData[LevelIndexer].units), transform=ax.transAxes, fontsize=16)
        # plot the cyclone center
        ax.scatter(track.loc[itime]['Lon'],track.loc[itime]['Lat'],
                   zorder=1000, color='#383838',linewidth=2,edgecolor='k')
        # colorbar
        cbar = plt.colorbar(cf1)
        cbar.ax.tick_params(labelsize=10) 
        # decorators
        map_features(ax)
    # save file
    outfile = FigsDirectory+'map_'+fname+'_'+str(date)
    plt.savefig(outfile)
    print(outfile+' created!')