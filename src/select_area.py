#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 12:33:00 2023

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
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
# crs_3857 = ccrs.epsg(3857)

# Transformation function
def coordXform(orig_crs, target_crs, x, y):
    return target_crs.transform_points( orig_crs, x, y )

def tellme(s):
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
                      facecolor='None', edgecolor='#6c757d', linewidth = 3,
                      alpha=1, zorder = 3)

def plot_slp(ax, slp, lat, lon):
    norm = colors.TwoSlopeNorm(vmin=-5e-5, vcenter=0, vmax=5e-5)
    cmap = cmo.balance
    # plot contours
    cf1 = ax.contourf(lon, lat, slp, cmap=cmap,norm=norm, transform=crs_longlat) 
    plt.colorbar(cf1, pad=0.07, orientation='vertical', shrink=0.5)
    ax.contour(lon, lat, slp, cf1.levels,colors='#383838',
               linewidths=0.25,transform=crs_longlat)
    
def map_decorators(ax):
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True,zorder=2,linestyle='dashed',alpha=0.8,
                 color='#383838')
    gl.xlabel_style = {'size': 14, 'color': '#383838'}
    gl.ylabel_style = {'size': 14, 'color': '#383838'}
    gl.top_labels = None
    gl.right_labels = None
    
def plot_min_slp(ax,slp, lat, lon, limits):
    max_lon, min_lon = limits['max_lon'], limits['min_lon']
    max_lat, min_lat = limits['max_lat'], limits['min_lat']
    # Plot mininum SLP point whithin box
    islp = slp.sel({lon.dims[0]:slice(min_lon,max_lon),
                    lat.dims[0]:slice(min_lat,max_lat)})
    slp_min = islp.min()
    slp_min_loc = islp.where(islp==slp_min, drop=True).squeeze()
    # sometimes there are multiple minimuns
    if slp_min_loc.shape:
        if len(slp_min_loc.shape) >1:
            for points in slp_min_loc:
                for point in points:
                    ax.scatter(point[lon.dims[0]], point[lat.dims[0]],
                           marker='o', facecolors='none', linewidth=3,
                           edgecolor='#a5a58d',  s=200)
        else:
            for point in slp_min_loc:
                ax.scatter(point[lon.dims[0]], point[lat.dims[0]],
                       marker='o', facecolors='none', linewidth=3,
                       edgecolor='#a5a58d',  s=200)
    else:
        ax.scatter(slp_min_loc[lon.dims[0]], slp_min_loc[lat.dims[0]],
               marker='o', facecolors='none', linewidth=3,
               edgecolor='#a5a58d',  s=200)

def initial_domain(slp, lat, lon):
    plt.close('all')
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=crs_longlat)
    fig.add_axes(ax)
    ax.set_global()
    plot_slp(ax, slp, lat, lon)
    map_decorators(ax)
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
        lls = coordXform(crs_longlat, crs_longlat, xs, ys)
        min_lon,min_lat = lls[0,0], lls[0,1]
        max_lon, max_lat = lls[1,0], lls[1,1]
        
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
    ax = plt.axes(projection=crs_longlat)
    fig.add_axes(ax)
    ax.set_extent([domain_limits['min_lon'], domain_limits['max_lon'],
                  domain_limits['min_lat'],domain_limits['max_lat']]) 
    
    plot_slp(ax, slp, lat, lon)
    ax.quiver(lon.values, lat.values, u.values, v.values,
              transform=crs_longlat)
    map_decorators(ax)
    
    while True:
        pts = []
        while len(pts) < nclicks:
            tellme('Select box corners \nModel time step: '+timestr[:-13])
            pts = np.asarray(plt.ginput(nclicks, timeout=15,
                                        mouse_stop='MouseButton.MIDDLE'))
            if len(pts) < nclicks:
                tellme('Too few points, starting over')
                time.sleep(1)  # Wait a second
        
        xmin,xmax = min([pts[0,0],pts[1,0]]),max([pts[0,0],pts[1,0]])
        ymin,ymax = min([pts[0,1],pts[1,1]]),max([pts[0,1],pts[1,1]])
        xs, ys = np.array((xmin,xmax)), np.array((ymin,ymax))
        lls = coordXform(crs_longlat, crs_longlat, xs, ys)
        min_lon,min_lat = lls[0,0], lls[0,1]
        max_lon, max_lat = lls[1,0], lls[1,1]
                
        limits = {'min_lon':min_lon,'max_lon':max_lon,
                'min_lat':min_lat, 'max_lat':max_lat}
        draw_box(ax, limits, crs_longlat)
        plot_min_slp(ax, slp, lat, lon, limits)
    
        tellme('Happy? Key press any keyboard key for yes, mouse click for no')
        
    
        if plt.waitforbuttonpress():
            break
        
    return limits