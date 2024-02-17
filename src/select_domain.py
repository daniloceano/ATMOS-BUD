# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    select_domain.py                                   :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/01/06 12:33:00 by daniloceano       #+#    #+#              #
#    Updated: 2024/02/17 00:45:15 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
A module for selecting and visualizing meteorological data domains.
Designed to assist in the analysis and interpretation of atmospheric phenomena.

Author: Danilo Couto de Souza
Affiliation: Universidade de São Paulo (USP), Instituto de Astornomia, Ciências Atmosféricas e Geociências, São Paulo - Brazil
Contact: danilo.oceano@gmail.com
"""

from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import time
from shapely.geometry.polygon import Polygon
import matplotlib.colors as colors
import cmocean.cm as cmo
import matplotlib.ticker as ticker

nclicks = 2
crs_longlat = ccrs.PlateCarree() 

# Transformation function
def coordXform(orig_crs, target_crs, x, y):
    """
    Transforms coordinates from one CRS to another.

    Parameters:
    - orig_crs: The original coordinate reference system.
    - target_crs: The target coordinate reference system.
    - x, y: The coordinates to transform.

    Returns:
    - Transformed coordinates.
    """
    return target_crs.transform_points( orig_crs, x, y )

def tellme(s):
    """
    Displays a message on the plot title.

    Parameters:
    - s: The message to display.
    """
    plt.title(s, fontsize=16)
    plt.draw()

def fmt(x, pos):
    """
    Formats the colorbar labels.

    Parameters:
    - x: The value to format.
    - pos: The position (unused, but required by FuncFormatter).

    Returns:
    - Formatted label.
    """
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
    
def draw_box(ax, limits, crs):
    """
    Draws a rectangular box on the plot based on specified limits.

    Parameters:
    - ax: The matplotlib axes object.
    - limits: A dictionary with 'min_lon', 'max_lon', 'min_lat', 'max_lat' keys.
    - crs: The coordinate reference system for the plot.
    """
    max_lon, min_lon = limits['max_lon'], limits['min_lon']
    max_lat, min_lat = limits['max_lat'], limits['min_lat']
    pgon = Polygon([(min_lon, min_lat), (min_lon, max_lat), (max_lon, max_lat), (max_lon, min_lat), (min_lon, min_lat)])
    ax.add_geometries([pgon], crs=crs, facecolor='None', edgecolor='k', linewidth=3, alpha=1, zorder=3)

def plot_zeta(ax, zeta, lat, lon, hgt=None):
    """
    Plots the vorticity field and optionally the geopotential height.

    Parameters:
    - ax: The matplotlib axes object for plotting.
    - zeta: The vorticity data array.
    - lat: Latitude coordinates.
    - lon: Longitude coordinates.
    - hgt: Optional geopotential height data array.
    """
    norm = colors.TwoSlopeNorm(vmin=min(-zeta.max(), zeta.min()), vcenter=0, vmax=max(zeta.max(), -zeta.min()))
    cmap = cmo.balance
    cf1 = ax.contourf(lon, lat, zeta, cmap=cmap, norm=norm, levels=51, transform=crs_longlat) 
    plt.colorbar(cf1, pad=0.1, orientation='vertical', shrink=0.5, format=ticker.FuncFormatter(fmt))
    if hgt is not None:
        cs = ax.contour(lon, lat, hgt, levels=11, colors='#747578', linestyles='dashed', linewidths=1, transform=crs_longlat)
        ax.clabel(cs, cs.levels, inline=True, fontsize=10)
    
def map_decorators(ax):
    """
    Adds coastlines and gridlines to the map.

    Parameters:
    - ax: The matplotlib axes object for the map.
    """
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True, zorder=2, linestyle='dashed', alpha=0.7, linewidth=0.5, color='#383838')
    gl.xlabel_style = {'size': 14, 'color': '#383838'}
    gl.ylabel_style = {'size': 14, 'color': '#383838'}
    gl.top_labels = None
    gl.right_labels = None
    
def plot_min_max_zeta(ax, zeta, lat, lon, limits):
    """
    Plots the minimum or maximum zeta point within a specified domain.

    Parameters:
    - ax: The matplotlib axes object for plotting.
    - zeta: The vorticity data array.
    - lat: Latitude coordinates.
    - lon: Longitude coordinates.
    - limits: A dictionary with 'min_lon', 'max_lon', 'min_lat', 'max_lat' keys specifying the domain.
    """
    max_lon, min_lon = limits['max_lon'], limits['min_lon']
    max_lat, min_lat = limits['max_lat'], limits['min_lat']
    # Plot mininum zeta point whithin box
    izeta = zeta.sel({lon.dims[0]:slice(min_lon,max_lon),
                    lat.dims[0]:slice(min_lat,max_lat)})
    if float(min_lat) < 0:
        min_max_zeta = izeta.min()
    else:
        min_max_zeta = izeta.max()
    min_max_zeta_loc = izeta.where(izeta==min_max_zeta, drop=True).squeeze()
    # sometimes there are multiple minimuns
    if min_max_zeta_loc.shape:
        if len(min_max_zeta_loc.shape) >1:
            for points in min_max_zeta_loc:
                for point in points:
                    ax.scatter(point[lon.dims[0]], point[lat.dims[0]],
                           marker='o', facecolors='none', linewidth=3,
                           edgecolor='k',  s=200)
        else:
            for point in min_max_zeta_loc:
                ax.scatter(point[lon.dims[0]], point[lat.dims[0]],
                       marker='o', facecolors='none', linewidth=3,
                       edgecolor='k',  s=200)
    else:
        ax.scatter(min_max_zeta_loc[lon.dims[0]], min_max_zeta_loc[lat.dims[0]],
               marker='o', facecolors='none', linewidth=3,
               edgecolor='k',  s=200)

def initial_domain(zeta, lat, lon):
    """
    Interactively selects an initial spatial domain on a map.

    Parameters:
    - zeta: The vorticity data array to be plotted for domain selection.
    - lat: Latitude coordinates for the plot.
    - lon: Longitude coordinates for the plot.

    Returns:
    - limits: A dictionary containing the selected domain's 'min_lon', 'max_lon', 'min_lat', 'max_lat'.
    """
    plt.close('all')
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=crs_longlat)
    fig.add_axes(ax)
    ax.set_global()
    plot_zeta(ax, zeta, lat, lon)
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

def draw_box_map(u, v, zeta, hgt, lat, lon, timestr):
    """
    Draws a map with streamlines and allows for the interactive selection of a domain.

    Parameters:
    - u: Eastward wind component data array.
    - v: Northward wind component data array.
    - zeta: Vorticity data array.
    - hgt: Geopotential height data array.
    - lat: Latitude coordinates.
    - lon: Longitude coordinates.
    - timestr: The timestep as a string for display.

    Returns:
    - limits: A dictionary with the selected domain's limits.
    """
    plt.close('all')
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=crs_longlat)
    fig.add_axes(ax)
    
    plot_zeta(ax, zeta, lat, lon, hgt)
    ax.streamplot(lon.values, lat.values, u.values, v.values, color='#2A1D21',
              transform=crs_longlat)
    map_decorators(ax)

    # Convert to datetime object
    date_obj = datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S%z')
    # Format datetime object to desired string format (YYYY-MM-DD HH)
    formatted_date_str = date_obj.strftime('%Y-%m-%d %H')
    
    while True:
        pts = []
        while len(pts) < nclicks:
            tellme(f'Select box corners \nModel time step: {formatted_date_str}')
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
        plot_min_max_zeta(ax, zeta, lat, lon, limits)
    
        tellme('Happy? Key press any keyboard key for yes, mouse click for no')
        
    
        if plt.waitforbuttonpress():
            break
        
    return limits

def get_domain_limits(args, *variables_at_850hpa, track=None):
    """
    Determines the domain limits based on track data or user selection.

    Parameters:
    - args: Command-line arguments or options.
    - variables_at_850hpa: A tuple containing meteorological variables at 850 hPa.
    - track: Optional DataFrame containing track data for domain selection.

    Returns:
    - current_domain_limits: A dictionary with the calculated domain limits.
    """
    iu_850, iv_850, zeta, ight_850, lat, lon, itime = variables_at_850hpa

    if args.track:
        if 'width'in track.columns:
            width, length = track.loc[itime]['width'],track.loc[itime]['length']

        else:
            width, length = 15, 15
        central_lon = track.loc[itime]['Lon']
        central_lat = track.loc[itime]['Lat']
        min_lon = central_lon-(width/2)
        max_lon = central_lon+(width/2)
        min_lat = central_lat-(length/2)
        max_lat = central_lat+(length/2)
        current_domain_limits = {
            'min_lon': min_lon,
            'max_lon': max_lon,
            'min_lat': min_lat,
            'max_lat': max_lat,
            'central_lat': central_lat,
            'central_lon': central_lon
            }
        
    elif args.choose:
        # Draw maps and ask user to specify corners for specifying the box
        limits = draw_box_map(iu_850, iv_850, zeta, ight_850,
                            lat, lon, itime)
        
        # Store system position and attributes
        min_lon, max_lon = limits['min_lon'],  limits['max_lon']
        min_lat, max_lat = limits['min_lat'],  limits['max_lat']
        width, length = max_lon - min_lon, max_lat - min_lat
        central_lat = (max_lat + min_lat) / 2
        central_lon = (max_lon + min_lon) / 2
        current_domain_limits = {
            'min_lon': min_lon,
            'max_lon': max_lon,
            'min_lat': min_lat,
            'max_lat': max_lat,
            'central_lat': central_lat,
            'central_lon': central_lon
            }

    elif args.fixed:
        dfbox = pd.read_csv('inputs/box_limits',header=None, delimiter=';',index_col=0)
        min_lon = float(lon.sel({lon.name:float(dfbox.loc['min_lon'].values)}, method='nearest'))
        max_lon = float(lon.sel({lon.name:float(dfbox.loc['max_lon'].values)}, method='nearest'))
        min_lat = float(lat.sel({lat.name:float(dfbox.loc['min_lat'].values)}, method='nearest'))
        max_lat = float(lat.sel({lat.name:float(dfbox.loc['max_lat'].values)}, method='nearest'))
        current_domain_limits = {
            'min_lon': min_lon,
            'max_lon': max_lon,
            'min_lat': min_lat,
            'max_lat': max_lat,
            'central_lat': (max_lat + min_lat) / 2,
            'central_lon': (max_lon + min_lon) / 2
            }
        
    return current_domain_limits