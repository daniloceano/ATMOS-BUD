# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    visualization.py                                   :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/01/11 16:54:18 by daniloceano       #+#    #+#              #
#    Updated: 2025/05/23 11:44:47 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

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
import os 
import cmocean.cm as cmo
import datetime
import pandas as pd
import numpy as np
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import TwoSlopeNorm
from shapely.geometry.polygon import Polygon
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature, COASTLINE
from cartopy.feature import BORDERS
from src.select_domain import plot_zeta 

def map_features(ax):
    """
    Adds coastline and border features to the matplotlib Axes.

    Parameters:
    - ax (matplotlib.axes._subplots.AxesSubplot): The Axes object to add features to.
    
    Returns:
    - ax (matplotlib.axes._subplots.AxesSubplot): The modified Axes object with added features.
    """
    ax.add_feature(COASTLINE, edgecolor='#283618', linewidth=1)
    ax.add_feature(BORDERS, edgecolor='#283618', linewidth=1)
    return ax

def Brazil_states(ax, facecolor='#a4ab98'):
    """
    Adds Brazilian state boundaries and major cities to the plot.

    Parameters:
    - ax (matplotlib.axes._subplots.AxesSubplot): The Axes object to add the features to.
    - facecolor (str): The color to use for land areas. Defaults to a grayish color.
    """
    land = cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor=facecolor)
    ax.add_feature(land)
    
    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none', name='admin_1_states_provinces_lines')
    ax.add_feature(states, edgecolor='#283618', linewidth=1)
    
    cities = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none', name='populated_places')
    ax.add_feature(cities, edgecolor='#283618', linewidth=1)

def plot_fixed_domain(limits, data_plevel, args, results_subdirectory, time, app_logger):
    """
    Creates and saves a plot of the specified domain with meteorological overlays.

    Parameters:
    - limits (dict): Dictionary containing the 'min_lon', 'max_lon', 'min_lat', and 'max_lat' of the domain.
    - data_plevel (dict): Dictionary containing data at chosen pressure level for plotting.
    - args (Namespace): Namespace object containing command-line arguments.
    - results_subdirectory (str): Directory path to save the plot.
    - time (str): Timestamp or identifier for the plot filename.
    - app_logger (logging.Logger): Logger object for logging messages.

    Raises:
    - Exception: If an error occurs during plotting or file saving.
    """
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
    pgon = Polygon(((min_lon, min_lat), (min_lon, max_lat), (max_lon, max_lat), (max_lon, min_lat), (min_lon, min_lat)))
    ax.add_geometries([pgon], crs=ccrs.PlateCarree(), facecolor='None', edgecolor='#BF3D3B', linewidth=3, alpha=1, zorder=3)

    # Add gridlines
    gl = ax.gridlines(draw_labels=True,zorder=2)    
    gl.xlabel_style = {'size': 16}
    gl.ylabel_style = {'size': 16}

    # Add title
    datestr = datetime.datetime.strptime(time, '%Y%m%d%H%M')
    timestr = datestr.strftime('%Y-%m-%d %H:%MZ')
    plt.title(f'Box defined for computations at {args.level} hPa\n{timestr}\n', fontsize=22)

    # Plot central point, mininum vorticity, minimum hgt and maximum wind 
    plot_zeta(ax, data_plevel[f'{args.track_vorticity}_zeta_{args.level}']['data'], data_plevel['lat'], data_plevel['lon'], data_plevel[f'{args.track_geopotential}_hgt_{args.level}']['data'])
    map_features(ax)
    Brazil_states(ax, facecolor='None')

    # Plot central point, mininum vorticity, minimum hgt and maximum wind
    ax.scatter(central_lon, central_lat,  marker='o', c='#31332e', s=100, zorder=4)
    ax.scatter(data_plevel[f'{args.track_vorticity}_zeta_{args.level}']['longitude'], data_plevel[f'{args.track_vorticity}_zeta_{args.level}']['latitude'],
                marker='s', c='#31332e', s=100, zorder=4, label=f'{args.track_vorticity} zeta')
    ax.scatter(data_plevel[f'{args.track_geopotential}_hgt_{args.level}']['longitude'], data_plevel[f'{args.track_geopotential}_hgt_{args.level}']['latitude'],
                marker='x', c='#31332e', s=100, zorder=4, label=f'{args.track_geopotential} hgt')
    ax.scatter(data_plevel['max_wind']['longitude'], data_plevel['max_wind']['latitude'],
                marker='^', c='#31332e', s=100, zorder=4, label='max wind')
    plt.legend(loc='upper left', frameon=True, fontsize=14, bbox_to_anchor=(1.1,1.2))

    # Save figure
    if time:
        boxes_directory = os.path.join(results_subdirectory, 'Figures', 'boxes')
        os.makedirs(boxes_directory, exist_ok=True)
        filename = os.path.join(boxes_directory, f'box_{time}.png')

    else:
        filename = os.path.join(results_subdirectory, 'Figures', 'box.png')

    plt.tight_layout()
    plt.savefig(filename)
    app_logger.info(f'\nCreated figure with box defined for computations at {os.path.basename(filename)}')    
    
def plot_track(track, args, figures_directory, app_logger):
    """
    Plots the track of a weather system and saves the figure.

    Parameters:
    - track (pandas.DataFrame): DataFrame containing the system's track data including longitude ('Lon'),
      latitude ('Lat'), and optionally 'min_max_zeta_plevel' and 'max_wind_plevel' for enhanced visualization.
    - args (argparse.Namespace): Command-line arguments or parameters specifying calculation options.
    - figures_directory (str): Directory path to save the plot.
    - app_logger (logging.Logger): Logger object for logging messages.

    Raises:
    - Exception: If an error occurs during plotting or file saving.
    """
    try:
        if not os.path.exists(figures_directory):
            os.makedirs(figures_directory)

        min_lon, max_lon = track['Lon'].min(), track['Lon'].max()
        min_lat, max_lat = track['Lat'].min(), track['Lat'].max()
        
        plt.close('all')
        fig, ax = plt.subplots(figsize=(12, 9) if (max_lon - min_lon) <= (max_lat - min_lat) else (14, 9), 
                               subplot_kw={"projection": ccrs.PlateCarree()})
        
        ax.set_extent([min_lon - 10, max_lon + 10, min_lat - 10, max_lat + 10], crs=ccrs.PlateCarree())
        ax.coastlines()
        ax.add_feature(cartopy.feature.LAND)
        ax.add_feature(cartopy.feature.OCEAN, facecolor="lightblue")
        
        gl = ax.gridlines(draw_labels=True, linestyle='--', alpha=0.5, color='darkgray')
        gl.top_labels = gl.right_labels = False
        
        # Plotting the track
        ax.plot(track['Lon'], track['Lat'], color='#383838')
        
        # Enhanced visualization with normalized wind speed and vorticity
        if f'{args.track_vorticity}_zeta_{args.level}' in track and f'max_wind_{args.level}' in track:
            normalized_sizes = normalize(track[[f'max_wind_{args.level}']].values.reshape(1, -1))**4 * 100000
            scatter = ax.scatter(track['Lon'], track['Lat'], zorder=3, c=track[f'{args.track_vorticity}_zeta_{args.level}'].values, 
                                 cmap=cmo.deep_r, edgecolor='gray', s=normalized_sizes.flatten(), 
                                 label=f'Normalized max wind speed at {args.level} hPa')
            plt.colorbar(scatter, pad=0.1, orientation='vertical', shrink=0.5, label=f'{args.track_vorticity} vorticity')
        
        # Marking the start and end points
        ax.text(track.iloc[0]['Lon'], track.iloc[0]['Lat'], 'A', fontsize=12, ha='center', va='center', color='green')
        ax.text(track.iloc[-1]['Lon'], track.iloc[-1]['Lat'], 'Z', fontsize=12, ha='center', va='center', color='red')
        
        filename = os.path.join(figures_directory, 'track_boxes.png')
        plt.savefig(filename, bbox_inches='tight')
        app_logger.info(f'Created figure with track and boxes defined for computations: {filename}')
    
    except Exception as e:
        app_logger.error(f'Failed to plot track: {e}')

def plot_min_max_zeta_hgt(track_plotting, args, figs_dir, app_logger, max_ticks=10):
    """
    Creates a dual-axis time series plot for min/max vorticity and minimum geopotential height at a specific pressure level.
    
    Parameters:
    - track_plotting (pandas.DataFrame): DataFrame containing the time series data to be plotted.
    - args (argparse.Namespace): Command-line arguments containing the pressure level and min/max indicator.
    - figs_dir (str): Directory where the plot image will be saved.
    - max_ticks (int): Maximum number of ticks to display on the x-axis.
    - app_logger (logging.Logger): Logger object for logging messages.
    
    The function saves the generated plot as a PNG file within the specified directory.
    """
    try:
        os.makedirs(figs_dir, exist_ok=True)

        fig, ax1 = plt.subplots(figsize=(15, 10))
        ax1.plot(pd.to_datetime(track_plotting.index), track_plotting[f'{args.track_vorticity}_zeta_{args.level}'], color='#554348', marker='o',
                 label=fr'{args.level} hPa {args.track_vorticity} $\zeta$')
        ax2 = ax1.twinx()
        ax2.plot(pd.to_datetime(track_plotting.index), track_plotting[f'{args.track_geopotential}_hgt_{args.level}'], color='#6610F2', marker='s',
                 label=f'{args.level} hPa {args.track_geopotential} gh')
        
        # Get handles and labels for each axis
        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()

        # Combine the handles and labels
        combined_handles = handles1 + handles2
        combined_labels = labels1 + labels2

        # Create a single legend with the combined handles and labels
        ax2.legend(combined_handles, combined_labels, loc='best', prop={'size': 18})
        
        # Customize axes and title
        ax1.set_title(fr'System Track: {args.track_vorticity} $\zeta$ and {args.track_geopotential} Geopotential Height at {args.level} hPa', fontsize=22)
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
        ax1.xaxis.set_major_locator(MaxNLocator(max_ticks))
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %HZ'))
        plt.setp(ax1.xaxis.get_majorticklabels(), rotation=20)
        
        # Save the figure
        filename = os.path.join(figs_dir,f"timeseries_{args.track_vorticity}_zeta_{args.track_geopotential}_hgt_{args.level}hPa.png")
        plt.savefig(filename, bbox_inches='tight')
        plt.close(fig)
        
        app_logger.info(f"Time series plot created and saved: {filename}")
    except Exception as e:
        app_logger.error(f"Failed to create and save time series plot: {e}")

from matplotlib.colors import TwoSlopeNorm

def hovmoller_mean_zeta(Zeta, figures_subdirectory, app_logger):
    """
    Creates and saves a Hovmöller diagram of mean relative vorticity for each pressure level
    over time within the system domain. Centered at zero (white).

    Parameters:
    - Zeta (pandas.DataFrame): DataFrame with index as pressure levels (Pa) and columns as timestamps (datetime),
      containing mean vorticity (1/s) for each level and time.
    - figures_subdirectory (str): Directory path to save the Hovmöller plot.
    - app_logger (logging.Logger): Logger object for logging messages.

    Raises:
    - Exception: If an error occurs during plotting or file saving.
    """
    try:
        os.makedirs(figures_subdirectory, exist_ok=True)

        # Remove colunas completamente vazias
        Zeta = Zeta.dropna(axis=1, how='all')

        # Converter valores e calcular limites simétricos
        zeta_plot = Zeta.values * 1e5  # converte para 1e-5 s⁻¹
        vmax = np.nanmax(np.abs(zeta_plot))
        vmin = -vmax
        norm = TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax)

        # Criar figura
        fig, ax = plt.subplots(figsize=(16, 9))
        cf = ax.contourf(Zeta.columns, Zeta.index / 100., zeta_plot,
                         levels=21, cmap=cmo.balance, norm=norm, extend='both')

        ax.invert_yaxis()

        ax.set_title('Hovmöller Diagram of Mean Vorticity (System Domain)', fontsize=22)
        ax.set_xlabel('Time', fontsize=18)
        ax.set_ylabel('Pressure level (hPa)', fontsize=18)

        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %Hh'))
        ax.xaxis.set_major_locator(mdates.AutoDateLocator(maxticks=10))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha='right')

        cbar = fig.colorbar(cf, ax=ax, orientation='vertical', pad=0.02, aspect=40)
        cbar.set_label('Vorticity ($10^{-5}$ s$^{-1}$)', fontsize=16)

        filename = os.path.join(figures_subdirectory, 'hovmoller_mean_zeta.png')
        plt.savefig(filename, bbox_inches='tight')
        plt.close(fig)

        app_logger.info(f'Hovmöller diagram of mean zeta created and saved: {filename}')

    except Exception as e:
        app_logger.error(f"Failed to create Hovmöller diagram of mean zeta: {e}")