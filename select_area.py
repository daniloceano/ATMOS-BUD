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


# Transformation function
def coordXform(orig_crs, target_crs, x, y):
    return target_crs.transform_points( orig_crs, x, y )

def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()

def draw_box_map(u, v, slp):
    fig = plt.figure(figsize=(15, 10))
    ax = plt.axes(projection=ccrs.epsg(3857))
    fig.add_axes(ax)
    ax.coastlines()
    ax.set_global()
    nclicks = 2
    # define all CRS
    crs_longlat = ccrs.PlateCarree() 
    crs_3857 = ccrs.epsg(3857)
    
    while True:
        pts = []
        while len(pts) < nclicks:
            tellme('Select box up-left and bottom-right corners')
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
        
        # plot selected domain
        # create a sample polygon, `pgon`
        pgon = Polygon(((min_lon, min_lat),
                (min_lon, max_lat),
                (max_lon, max_lat),
                (max_lon, min_lat),
                (min_lon, min_lat)))
        ax.add_geometries([pgon], crs=crs_longlat, 
                          facecolor='None', edgecolor='#BF3D3B', linewidth = 3,
                          alpha=1, zorder = 3)
    
        tellme('Happy? Key press any keyboard key for yes, mouse click for no')
    
        if plt.waitforbuttonpress():
            break
        
        return {'min_lon':min_lon,'max_lon':max_lon,
                'min_lat':min_lat, 'max_lat':max_lat}