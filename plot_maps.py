#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 08:52:34 2022

@author: danilocoutodsouza
"""

import xarray as xr
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import numpy as np 
import matplotlib.gridspec as gridspec
from matplotlib import colors
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature, COASTLINE
from cartopy.feature import BORDERS
from celluloid import Camera
import pandas as pd
import sys
from main import check_create_folder
from main import convert_lon

def grid_labels_params(ax):
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5,linestyle='--')
    gl.top_labels = False
    gl.left_labels = False
    gl.xlabel_style = {'size': 14, 'color': 'gray'}
    gl.ylabel_style = {'size': 14, 'color': 'gray'}
    ax.spines['geo'].set_edgecolor('gray')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return ax

def map_features(ax):
    ax.add_feature(COASTLINE)
    ax.add_feature(BORDERS, edgecolor='white')
    return ax

def Brazil_states(ax):    
    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                  name='admin_1_states_provinces_lines')
    _ = ax.add_feature(states, edgecolor='white')
    
    cities = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                  name='populated_places')
    _ = ax.add_feature(cities)
    
def grid_labels_params(ax):
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5,linestyle='--')
    gl.top_labels = False
    gl.left_labels = False
    gl.xlabel_style = {'size': 14, 'color': 'gray'}
    gl.ylabel_style = {'size': 14, 'color': 'gray'}
    ax.spines['geo'].set_edgecolor('gray')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return ax

def LagrangianMaps(VariableData):
    
    dfVars = pd.read_csv('./fvars',sep= ';',index_col=0,header=0)
    LonIndexer,LatIndexer,TimeIndexer,LevelIndexer = \
      dfVars.loc['Longitude']['Variable'],dfVars.loc['Latitude']['Variable'],\
      dfVars.loc['Time']['Variable'],dfVars.loc['Vertical Level']['Variable']
      
    trackfile = './track'
    track = pd.read_csv(trackfile,parse_dates=[0],delimiter=';',index_col='time')
    westernlimit, easternlimit = track['Lon'].min()-15,track['Lon'].max()+15
    southernlimit, northernlimit = track['Lat'].min()-15,track['Lat'].max()+15
     
    lon,lat = VariableData[LonIndexer], VariableData[LatIndexer]
    
    proj = ccrs.PlateCarree() 
    cmap = cmo.balance
    # levels = np.arange(1000,1036,2)
    # norm = colors.TwoSlopeNorm(vmin=1010, vcenter=1014, vmax=1030)
    lims = [westernlimit,easternlimit,southernlimit,northernlimit]
    
    plt.close('all')
    fig = plt.figure(constrained_layout=False,figsize=(10,10))
    gs = gridspec.GridSpec(2, 3, hspace=0.2, wspace=0,
                                   left=0, right=0.95)
    
    # Find vertical levels that closely matches the desired levels for plotting
    MatchingLevels = []
    for pres in [1000,850,700,500,300,200]:
        MatchingLevels.append(min(VariableData[LevelIndexer].values, key=lambda x:abs(x-pres)))
    
    for p,i in zip(MatchingLevels,range(len(MatchingLevels))):
        
        itime = VariableData[TimeIndexer].values
        min_lon, max_lon = track.loc[itime]['Lon']-7.5,track.loc[itime]['Lon']+7.5
        min_lat, max_lat = track.loc[itime]['Lat']-7.5,track.loc[itime]['Lat']+7.5
        
        ax = fig.add_subplot(gs[i], projection=proj)
        ax.set_extent(lims) 
        
        iData = VariableData.sel({LevelIndexer:p}).sel(
            **{LatIndexer:slice(max_lat,min_lat),
               LonIndexer: slice(min_lon,max_lon)})
        
        cf1 = ax.contourf(lon, lat, iData) 
        ax.contour(lon, lat, iData, cf1.levels,colors='grey', linewidths=1)
        grid_labels_params(ax)
        map_features(ax)
        Brazil_states(ax)
        timestr = pd.to_datetime(str(iData[TimeIndexer].values))
        date = timestr.strftime('%Y-%m-%dT%H%MZ')
        ax.text(0.05,1.01,str(p)+' '+str(VariableData[LevelIndexer].units), transform=ax.transAxes, fontsize=16)
        