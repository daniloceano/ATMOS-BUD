#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 16:54:18 2023

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from shapely.geometry.polygon import Polygon
from sklearn import preprocessing
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature, COASTLINE
from cartopy.feature import BORDERS
import cartopy
import cmocean.cm as cmo
import matplotlib.dates as mdates

def map_features(ax):
    ax.add_feature(COASTLINE,edgecolor='#283618',linewidth=1)
    ax.add_feature(BORDERS,edgecolor='#283618',linewidth=1)
    return ax

def Brazil_states(ax):    
    
    _ = ax.add_feature(cfeature.NaturalEarthFeature('physical',
                        'land', '50m', edgecolor='face', facecolor='#a4ab98'))
    
    states = NaturalEarthFeature(category='cultural', scale='50m', 
                                 facecolor='none',
                                  name='admin_1_states_provinces_lines')
    _ = ax.add_feature(states, edgecolor='#283618',linewidth=1)
    
    cities = NaturalEarthFeature(category='cultural', scale='50m',
                                 facecolor='none',
                                  name='populated_places')
    _ = ax.add_feature(cities, edgecolor='#283618',linewidth=1)

    
def plot_fixed_domain(min_lon, max_lon, min_lat, max_lat, outdir):

    plt.close('all')
    datacrs = ccrs.PlateCarree() # projection
    fig = plt.figure(figsize=(8, 8.5))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.83], projection=datacrs,
                  frameon=True)
    ax.set_extent([min_lon-20, max_lon+20, max_lat+20, min_lat-20], crs=datacrs)
    Brazil_states(ax)
    map_features(ax)
    
    # plot selected domain
    # create a sample polygon, `pgon`
    pgon = Polygon(((min_lon, min_lat),
            (min_lon, max_lat),
            (max_lon, max_lat),
            (max_lon, min_lat),
            (min_lon, min_lat)))
    ax.add_geometries([pgon], crs=datacrs, 
                      facecolor='None', edgecolor='#BF3D3B', linewidth = 3,
                      alpha=1, zorder = 3)
    gl = ax.gridlines(draw_labels=True,zorder=2)    
    gl.xlabel_style = {'size': 16}
    gl.ylabel_style = {'size': 16}

    plt.title('Box defined for compuations \n', fontsize = 22)
    plt.savefig(outdir+'Figures/box.png')
    print('\nCreated figure with box defined for computations at '
          +outdir+'Figures/box.png')
    
def plot_track(track, FigsDir):
        
    min_lon,max_lon,min_lat,max_lat=min(track['Lon']),max(track['Lon']),\
        min(track['Lat']),max(track['Lat'])
    
    plt.close('all')
    datacrs = ccrs.PlateCarree() # projection
    # If the figure lenght > height the fig size and legend position have
    # ot be adjusted otherwise they won't appear
    lenx, leny = max_lon - min_lon, max_lat - min_lat
    if lenx <= leny:
        fig = plt.figure(figsize=(12, 9))
    else:
        fig = plt.figure(figsize=(14, 9))
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9], projection=datacrs,
                  frameon=True)
    ax.set_extent([min_lon-10, max_lon+10, max_lat+10, min_lat-10], crs=datacrs)
    ax.coastlines(zorder = 1)
    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN,facecolor=("lightblue"))
    gl = ax.gridlines(draw_labels=True,zorder=2,linestyle='dashed',alpha=0.8,
                 color='#383838')
    gl.xlabel_style = {'size': 14, 'color': '#383838'}
    gl.ylabel_style = {'size': 14, 'color': '#383838'}
    gl.bottom_labels = None
    gl.right_labels = None
    
    for i in range(len(track)):
        
        # Model timestep
        itime = str(track.index[i])
        
        itime_track = track.index[track.index==itime]
        lon = track.loc[itime_track]['Lon'].values
        lat = track.loc[itime_track]['Lat'].values
        
        min_lon, max_lon = lon-7.5,lon+7.5
        min_lat, max_lat = lat-7.5,lat+7.5
            
    plt.plot(track['Lon'], track['Lat'],c='#383838')
    
    if 'min_zeta_850' and 'max_wind_850' in track.columns:
        normalized = preprocessing.normalize(
            track['max_wind_850'].values.reshape(1, -1))
        normalized = normalized**4 * 100000
        scatter = ax.scatter(track['Lon'].loc[track.index],
                           track['Lat'].loc[track.index],
                           zorder=100,c=track['min_zeta_850'],
                           cmap=cmo.deep_r, edgecolor='gray',
                           s=normalized, label=normalized)
        # produce a legend with a cross section of sizes from the scatter
        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
        ax.legend(handles, labels, loc="upper right",
                            title="Max wind speed at 850 hPa")
        plt.colorbar(scatter, pad=0.07, orientation='vertical', shrink=0.5,
                     label= 'Minimum vorticity')
        
    else:
        plt.scatter(track['Lon'].loc[track.index],
                           track['Lat'].loc[track.index],
                           zorder=100,edgecolors='grey')
    
    # mark beginning and end of system
    start = [track.loc[track.index[0]]['Lon'],
             track.loc[track.index[0]]['Lat']]
    end =  [track.loc[track.index[-1]]['Lon'],
             track.loc[track.index[-1]]['Lat']]
    ax.text(*start,'A',zorder=101,fontsize=24,horizontalalignment='center',
            verticalalignment='center')
    ax.text(*end,'Z',zorder=101,fontsize=24,horizontalalignment='center',
            verticalalignment='center')
    
    # plt.title('System track and boxes defined for compuations \n', fontsize = 22)
    plt.savefig(FigsDir+'track_boxes.png',bbox_inches='tight')
    print('\nCreated figure with track and boxes defined for computations: '
          +FigsDir+'track_boxes.png')
    
def plot_min_zeta_hgt(track, FigsDir):
    plt.close('all')
    fig, ax1 = plt.subplots(figsize=(15,10))
    lns1 = ax1.plot(track.index, track['min_zeta_850'],c='#554348', marker='o',
             label= '850 hPa minimum vorticity')
    ax2 = ax1.twinx()
    lns2 = ax2.plot(track.index, track['min_hgt_850'],c='#6610F2', marker='s',
             label= '850 hPa minimum geopotential height')
    # added these three lines
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc='best',prop={'size': 18})
    ax1.grid(c='gray',linewidth=0.25,linestyle='dashdot', axis='x')
    ax1.tick_params(axis='x', labelrotation=20)
    ax1.xaxis.set_tick_params(labelsize=16)
    ax1.yaxis.set_tick_params(labelsize=16)
    ax2.yaxis.set_tick_params(labelsize=16)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %HZ'))
    
    # plt.title('System track and boxes defined for compuations \n', fontsize = 22)
    plt.savefig(FigsDir+'timeseries-min_zeta_hgt.png',bbox_inches='tight')
    print('\nCreated:',FigsDir+'track_boxes.png')