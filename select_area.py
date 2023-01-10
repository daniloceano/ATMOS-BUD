#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 12:33:00 2023

@author: danilocoutodsouza
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import time
from shapely.geometry.polygon import Polygon
import matplotlib.colors as colors
import cmocean.cm as cmo

nclicks = 2

# define all CRS
crs_longlat = ccrs.PlateCarree() 
crs_3857 = ccrs.epsg(3857)

# Transformation function
def coordXform(orig_crs, target_crs, x, y):
    return target_crs.transform_points( orig_crs, x, y )

def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()
    
def draw_box(ax, limits, crs):
    max_lon, min_lon = limits['max_lon'], limits['min_lon']
    max_lat, min_lat = limits['max_lat'], limits['min_lat']
    pgon = Polygon(((min_lon, min_lat),
            (min_lon, max_lat),
            (max_lon, max_lat),
            (max_lon, min_lat),
            (min_lon, min_lat)))
    ax.add_geometries([pgon], crs=crs, 
                      facecolor='None', edgecolor='#BF3D3B', linewidth = 3,
                      alpha=1, zorder = 3)

def plot_slp(ax, slp, lat, lon):
    if float(slp.min()) < float(slp.max()):
             norm = colors.TwoSlopeNorm(vmin=float(slp.min()),
                               vcenter=1014, vmax=float(slp.max()))
    else:
        norm = colors.TwoSlopeNorm(vmin=float(slp.min()),
                          vcenter=1014, vmax=1015)
    cmap = cmo.balance
    # plot contours
    cf1 = ax.contourf(lon, lat, slp, cmap=cmap,norm=norm, transform=crs_longlat) 
    ax.contour(lon, lat, slp, cf1.levels,colors='#383838'
               ,linewidths=0.25,transform=crs_longlat)
    plt.colorbar(cf1,pad=0.07,orientation='vertical', norm=norm, shrink=0.5)

def plot_wind(ax, u, v, lat, lon):
    ax.streamplot(lon.values, lat.values, u.values, v.values, 
              transform=crs_longlat)

def initial_domain(slp, lat, lon):
    plt.close('all')
    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.epsg(3857))
    fig.add_axes(ax)
    ax.coastlines()
    ax.set_global()
    plot_slp(ax, slp, lat, lon)
    plt.subplots_adjust(bottom=0, top=1.2)
    
    while True:
        pts = []
        while len(pts) < nclicks:
            tellme('Select an initial spatial domain')
            pts = np.asarray(plt.ginput(nclicks, timeout=15,
                                        mouse_stop='MouseButton.MIDDLE'))
            if len(pts) < nclicks:
                tellme('Too few points, starting over')
                time.sleep(1)  # Wait a second
                
        xmin,xmax = min([pts[0,0],pts[1,0]]),max([pts[0,0],pts[1,0]])
        ymin,ymax = min([pts[0,1],pts[1,1]]),max([pts[0,1],pts[1,1]])
        xs, ys = np.array((xmin,xmax)), np.array((ymin,ymax))
        lls = coordXform(crs_3857, crs_longlat, xs, ys)
        min_lon,min_lat = lls[0,0], lls[0,1]
        max_lon, max_lat = lls[1,0], lls[1,1]
        
        print('min_lon,min_lat: ',min_lon,min_lat)
        print('max_lon, max_lat: ',max_lon, max_lat)
        limits = {'min_lon':min_lon,'max_lon':max_lon,
                'min_lat':min_lat, 'max_lat':max_lat}
        draw_box(ax, limits, crs_longlat)
        
        tellme('Happy? Key press any keyboard key for yes, mouse click for no')
        if plt.waitforbuttonpress():
            break
        
    return limits
    

def draw_box_map(u, v, slp, lat, lon, timestr, domain_limits):
    plt.close('all')
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.epsg(3857))
    fig.add_axes(ax)
    ax.coastlines()
    ax.set_extent([domain_limits['min_lon'], domain_limits['max_lon'],
                  domain_limits['min_lat'],domain_limits['max_lat']]) 
    
    plot_slp(ax, slp, lat, lon)
    plot_wind(ax, u, v, lat, lon)
    
    while True:
        pts = []
        while len(pts) < nclicks:
            tellme('Select box up-left and bottom-right corners \n\
model time step: '+timestr)
            pts = np.asarray(plt.ginput(nclicks, timeout=15,
                                        mouse_stop='MouseButton.MIDDLE'))
            if len(pts) < nclicks:
                tellme('Too few points, starting over')
                time.sleep(1)  # Wait a second
        
        xmin,xmax = min([pts[0,0],pts[1,0]]),max([pts[0,0],pts[1,0]])
        ymin,ymax = min([pts[0,1],pts[1,1]]),max([pts[0,1],pts[1,1]])
        xs, ys = np.array((xmin,xmax)), np.array((ymin,ymax))
        lls = coordXform(crs_3857, crs_longlat, xs, ys)
        min_lon,min_lat = lls[0,0], lls[0,1]
        max_lon, max_lat = lls[1,0], lls[1,1]
        
        print('min_lon,min_lat: ',min_lon,min_lat)
        print('max_lon, max_lat: ',max_lon, max_lat)
        limits = {'min_lon':min_lon,'max_lon':max_lon,
                'min_lat':min_lat, 'max_lat':max_lat}
        draw_box(ax, limits, crs_longlat)
    
        tellme('Happy? Key press any keyboard key for yes, mouse click for no')
    
        if plt.waitforbuttonpress():
            break
        
    return limits