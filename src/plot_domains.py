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
import os 
from select_area import plot_zeta 

def check_create_folder(DirName, verbose=True):
    """

    Check if directory exists and if not, creates it.
    
    Parameters
    ----------
    DirName : str
        directory name.

    Returns
    -------
    None. 

    """
    if not os.path.exists(DirName):
                os.makedirs(DirName)
                print(DirName+' created')
    else:
        if verbose:
            print(DirName+' directory exists')

def map_features(ax):
    ax.add_feature(COASTLINE,edgecolor='#283618',linewidth=1)
    ax.add_feature(BORDERS,edgecolor='#283618',linewidth=1)
    return ax

def Brazil_states(ax, facecolor='#a4ab98'):    
    
    _ = ax.add_feature(cfeature.NaturalEarthFeature('physical',
                        'land', '50m', edgecolor='face', facecolor=facecolor))
    
    states = NaturalEarthFeature(category='cultural', scale='50m', 
                                 facecolor='none',
                                  name='admin_1_states_provinces_lines')
    _ = ax.add_feature(states, edgecolor='#283618',linewidth=1)
    
    cities = NaturalEarthFeature(category='cultural', scale='50m',
                                 facecolor='none',
                                  name='populated_places')
    _ = ax.add_feature(cities, edgecolor='#283618',linewidth=1)

def plot_fixed_domain(limits, data850, ResultsSubDirectory, time):

    # Central coordinates
    max_lon, min_lon = limits['max_lon'], limits['min_lon']
    max_lat, min_lat = limits['max_lat'], limits['min_lat']
    central_lon = (max_lon+min_lon)/2
    central_lat = (max_lat+min_lat)/2

    # Create figure
    plt.close('all')
    fig, ax = plt.subplots(figsize=(8, 8.5), subplot_kw=dict(projection=ccrs.PlateCarree()))

    # Set map extent and features
    ax.set_extent([min_lon-20, max_lon+20, max_lat+20, min_lat-20], crs=ccrs.PlateCarree())
    map_features(ax)
    Brazil_states(ax, facecolor='None')
    
    # Plot selected domain
    # Create a sample polygon, `pgon`
    pgon = Polygon(((min_lon, min_lat),
                    (min_lon, max_lat),
                    (max_lon, max_lat),
                    (max_lon, min_lat),
                    (min_lon, min_lat)))
    ax.add_geometries([pgon], crs=ccrs.PlateCarree(), 
                      facecolor='None', edgecolor='#BF3D3B', linewidth=3,
                      alpha=1, zorder=3)

    # Add gridlines
    gl = ax.gridlines(draw_labels=True,zorder=2)    
    gl.xlabel_style = {'size': 16}
    gl.ylabel_style = {'size': 16}

    # Add title
    plt.title('Box defined for computations\n', fontsize=22)

    # Plot central point, mininum vorticity, minimum hgt and maximum wind 
    plot_zeta(ax, data850['min_zeta']['data'], data850['lat'], data850['lon'], data850['min_hgt']['data'])
    map_features(ax)
    Brazil_states(ax, facecolor='None')

    # Plot central point, mininum vorticity, minimum hgt and maximum wind
    ax.scatter(central_lon, central_lat,  marker='o', c='#31332e', s=100, zorder=4)
    ax.scatter(data850['min_zeta']['longitude'], data850['min_zeta']['latitude'],
                marker='s', c='#31332e', s=100, zorder=4, label='min zeta')
    ax.scatter(data850['min_hgt']['longitude'], data850['min_hgt']['latitude'],
                marker='x', c='#31332e', s=100, zorder=4, label='min hgt')
    ax.scatter(data850['max_wind']['longitude'], data850['max_wind']['latitude'],
                marker='^', c='#31332e', s=100, zorder=4, label='max wind')
    plt.legend(loc='upper left', frameon=True, fontsize=14, bbox_to_anchor=(1.1,1.2))

    # Save figure
    if time:
        boxes_directory = os.path.join(ResultsSubDirectory, 'Figures', 'boxes')
        check_create_folder(boxes_directory, verbose=False)
        filename = os.path.join(boxes_directory, f'box_{time}.png')
    else:
        filename = os.path.join(ResultsSubDirectory, 'Figures', 'box.png')
    plt.savefig(filename)
    print(f'\nCreated figure with box defined for computations at {filename}')    
    
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
    
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.dates as mdates

def plot_min_zeta_hgt(track, figs_dir, max_ticks=10):
    fig, ax1 = plt.subplots(figsize=(15, 10))
    line1 = ax1.plot(track.index, track['min_zeta_850'], c='#554348', marker='o',
                     label='850 hPa minimum vorticity')
    ax2 = ax1.twinx()
    line2 = ax2.plot(track.index, track['min_hgt_850'], c='#6610F2', marker='s',
                     label='850 hPa minimum geopotential height')
    
    # Combine lines and labels
    lines = line1 + line2
    labels = [line.get_label() for line in lines]
    
    # Add legend
    ax2.legend(lines, labels, loc='best', prop={'size': 18})
    
    # Customize axes
    ax1.set_title('System track and boxes defined for computations', fontsize=22)
    ax1.grid(c='gray', linewidth=0.25, linestyle='dashdot', axis='x')
    ax1.xaxis.set_major_locator(MaxNLocator(max_ticks))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %HZ'))
    ax1.tick_params(axis='x', labelrotation=20)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    
    # Save the figure
    filename = f"{figs_dir}timeseries-min_zeta_hgt.png"
    plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)
    
    print(f"Created: {filename}")


if __name__ == '__main__':
    import pandas as pd
    track = pd.read_csv('../inputs/track-test',parse_dates=[0],
                            delimiter=';',index_col='time')
    plot_min_zeta_hgt(track, './')