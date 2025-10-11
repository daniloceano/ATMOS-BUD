# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    utils.py                                           :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2024/02/16 18:31:30 by daniloceano       #+#    #+#              #
#    Updated: 2025/08/15 09:30:24 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import sys
import logging
import numpy as np
import pandas as pd
from metpy.calc import vorticity
from scipy.signal import savgol_filter
from src.select_domain import initial_domain

def initialize_logging(results_subdirectory, args):
    """
    Initializes the logging configuration for the application.

    Args:
        results_subdirectory (str): Directory path to save the log file.
        args (object): The argparse object containing the command line arguments.
    """

    verbose = True if args.verbose else False

    # Set root logger to higher severity level (INFO or ERROR)
    root_log_level = logging.ERROR if not verbose else logging.INFO
    logging.basicConfig(level=root_log_level, format='%(asctime)s - %(levelname)s - %(message)s')

    # Create a separate logger for the application
    app_logger = logging.getLogger('atmos_bud')
    app_log_level = logging.DEBUG if verbose else logging.INFO
    app_logger.setLevel(app_log_level)
    app_logger.propagate = False  # Prevent the logger from propagating messages to the root logger

    # Create file handler for saving logs
    log_file_name = f'log.{os.path.basename(args.infile).split(".")[0]}'
    log_file = os.path.join(results_subdirectory, log_file_name)
    file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')  # Ensure file uses UTF-8
    file_handler.setLevel(app_log_level)
    file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    app_logger.addHandler(file_handler)

    # Create a console handler for app logger with UTF-8 encoding
    console_handler = logging.StreamHandler(stream=sys.stdout)  # Use sys.stdout explicitly
    console_handler.setLevel(app_log_level)
    console_handler.setFormatter(file_formatter)
    # Set UTF-8 encoding for console output
    if sys.platform.startswith('win'):
        # On Windows, wrap the stream to ensure UTF-8 encoding
        console_handler.stream = open(sys.stdout.fileno(), mode='w', encoding='utf-8', errors='replace')
    app_logger.addHandler(console_handler)

    return app_logger

def convert_lon(input_data, longitude_indexer):
    """
    
    Convert longitudes from 0:360 range to -180:180

    Parameters
    ----------
    xr : xarray.DataArray 
        gridded data.
    longitude_indexer : str
        corrdinate indexer used for longitude.

    Returns
    -------
    xr : xarray.DataArray 
        gridded data with longitude converted to desired format.

    """
    input_data.coords[longitude_indexer] = (input_data.coords[longitude_indexer] + 180) % 360 - 180
    input_data = input_data.sortby(input_data[longitude_indexer])
    return input_data

def handle_track_file(input_data, times, longitude_indexer, latitude_indexer, app_logger):
    """
    Handles the track file by validating its time and spatial limits against the provided dataset.

    Args:
        data (xr.Dataset): A Xarray Dataset containing the data to compute the energy cycle.
        times (pd.DatetimeIndex): The time series of the dataset.
        LonIndexer (str): The name of the longitude coordinate in the dataset.
        LatIndexer (str): The name of the latitude coordinate in the dataset.
        args (argparse.Namespace): Arguments provided to the script.
        app_logger (logging.Logger): Logger for logging messages.

    Returns:
        pd.DataFrame: DataFrame containing the track information if the track file is valid.

    Raises:
        FileNotFoundError: If the track file is not found.
        ValueError: If the time or spatial limits of the track file do not match the dataset.
    """
    
    trackfile = "inputs/track"

    try:
        track = pd.read_csv(trackfile, parse_dates=[0], delimiter=';', index_col='time')
    except FileNotFoundError:
        app_logger.error(f'Track file not found: {trackfile}')
        raise
    except pd.errors.ParserError:
        app_logger.error(f'Error parsing track file: {trackfile}')
        raise

    # Remove timezone information from track index if it's tz-aware
    if track.index.tz is not None:
        track.index = track.index.tz_localize(None)

    try:
        data_lon_max, data_lon_min = float(input_data[longitude_indexer].max()), float(input_data[longitude_indexer].min())
        data_lat_max, data_lat_min = float(input_data[latitude_indexer].max()), float(input_data[latitude_indexer].min())

        app_logger.debug(f"Data spatial limits --> lon_min: {data_lon_min:.2f}, lon_max: {data_lon_max:.2f}, lat_min: {data_lat_min:.2f}, lat_max: {data_lat_max:.2f}")
        app_logger.debug(f"Track spatial limits --> lon_min: {track['Lon'].min():.2f}, lon_max: {track['Lon'].max():.2f}, lat_min: {track['Lat'].min():.2f}, lat_max: {track['Lat'].max():.2f}")

        # Check time if track time limits match with data time limits
        if track.index[0] < times.min() or track.index[-1] > times.max():
            app_logger.error("Track time limits do not match with data time limits.")
            raise ValueError("Track time limits do not match with data time limits.")

        # Check longitude limits
        if track['Lon'].max() > data_lon_max:
            app_logger.error(f"Track file longitude max limit ({track['Lon'].max():.2f}) exceeds data max longitude limit ({data_lon_max:.2f}).")
            raise ValueError(f"Track file longitude max limit ({track['Lon'].max():.2f}) exceeds data max longitude limit ({data_lon_max:.2f}).")
        if track['Lon'].min() < data_lon_min:
            app_logger.error(f"Track file longitude min limit ({track['Lon'].min():.2f}) is below data min longitude limit ({data_lon_min:.2f}).")
            raise ValueError(f"Track file longitude min limit ({track['Lon'].min():.2f}) is below data min longitude limit ({data_lon_min:.2f}).")

        # Check latitude limits
        if track['Lat'].max() > data_lat_max:
            app_logger.error(f"Track file latitude max limit ({track['Lat'].max():.2f}) exceeds data max latitude limit ({data_lat_max:.2f}).")
            raise ValueError(f"Track file latitude max limit ({track['Lat'].max():.2f}) exceeds data max latitude limit ({data_lat_max:.2f}).")
        if track['Lat'].min() < data_lat_min:
            app_logger.error(f"Track file latitude min limit ({track['Lat'].min():.2f}) is below data min latitude limit ({data_lat_min:.2f}).")
            raise ValueError(f"Track file latitude min limit ({track['Lat'].min():.2f}) is below data min latitude limit ({data_lat_min:.2f}).")

        return track
    
    except FileNotFoundError:
        app_logger.error(f"Track file {trackfile} not found.")
        raise


def find_extremum_coordinates(ds_data, lat, lon, variable, args):
    """
    Finds the indices of the extremum values for a given variable.

    Args:
    ds_data: An xarray DataArray containing the data to compute the energy cycle.
    lat: An xarray DataArray containing the latitudes of the data.
    lon: An xarray DataArray containing the longitudes of the data.
    variable: A string containing the name of the variable to find the indices for.

    Returns:
    A tuple containing the indices of the extremum values for the specified variable.
    """
    
    lat_values = lat.values
    lon_values = lon.values

    if (f'{args.track_vorticity}_zeta' in variable):
        if args.track_vorticity == "min":
            index = np.unravel_index(np.argmin(ds_data.data), ds_data.shape)

        else:
            index = np.unravel_index(np.argmax(ds_data.data), ds_data.shape)

    elif (f'{args.track_geopotential}_hgt' in variable):
        if args.track_geopotential == "min":
            index = np.unravel_index(np.argmin(ds_data.data), ds_data.shape)
        else:
            index = np.unravel_index(np.argmax(ds_data.data), ds_data.shape)

    elif variable == 'max_wind':
        index = np.unravel_index(np.argmax(ds_data.data), ds_data.shape)
    else:
        raise ValueError("Invalid variable specified.")

    lat_index = index[0]
    lon_index = index[1]

    lat_value = lat_values[lat_index]
    lon_value = lon_values[lon_index]

    return lat_value, lon_value

def slice_domain(input_data, args, namelist_df):
    """
    Slices the input dataset according to the specified domain. The domain can be defined
    based on fixed boundaries, track information, or an interactively selected area.

    Parameters:
    - input_data (xr.Dataset): The dataset containing meteorological variables.
    - args (argparse.Namespace): Parsed command-line arguments indicating the slicing method.
    - namelist_df (pd.DataFrame): DataFrame mapping variable names to their respective dataset labels.

    Returns:
    - xr.Dataset: The sliced dataset within the specified domain.
    """
    
    # Data indexers
    longitude_indexer = namelist_df.loc['Longitude']['Variable']
    latitude_indexer = namelist_df.loc['Latitude']['Variable']
    time_indexer = namelist_df.loc['Time']['Variable']
    vertical_level_indexer = namelist_df.loc['Vertical Level']['Variable']

    # Get choosen vertical level from args and convert to Pa
    plevel_Pa = int(args.level) * 100

    # Input data limits
    input_data_lon_min = input_data[longitude_indexer].min().values
    input_data_lon_max = input_data[longitude_indexer].max().values
    input_data_lat_min = input_data[latitude_indexer].min().values
    input_data_lat_max = input_data[latitude_indexer].max().values
    
    if args.fixed:
        dfbox = pd.read_csv('inputs/box_limits', header=None, delimiter=';', index_col=0)
        dfbox_min_lon = float(dfbox.loc['min_lon'].values.item())
        dfbox_max_lon = float(dfbox.loc['max_lon'].values.item())
        dfbox_min_lat = float(dfbox.loc['min_lat'].values.item())
        dfbox_max_lat = float(dfbox.loc['max_lat'].values.item())
        
        # Check if dfbox limits are within input_data limits
        if (dfbox_min_lon < input_data_lon_min or dfbox_max_lon > input_data_lon_max or
            dfbox_min_lat < input_data_lat_min or dfbox_max_lat > input_data_lat_max):
            raise ValueError("Specified domain limits in 'inputs/box_limits' are outside the range of 'input_data' limits.")

        WesternLimit = float(input_data[longitude_indexer].sel(
            {longitude_indexer: dfbox_min_lon}, method='nearest'))
        EasternLimit = float(input_data[longitude_indexer].sel(
            {longitude_indexer: dfbox_max_lon}, method='nearest'))
        SouthernLimit = float(input_data[latitude_indexer].sel(
            {latitude_indexer: dfbox_min_lat}, method='nearest'))
        NorthernLimit = float(input_data[latitude_indexer].sel(
            {latitude_indexer: dfbox_max_lat}, method='nearest'))

    elif args.track:
        trackfile = 'inputs/track'
        track = pd.read_csv(trackfile,parse_dates=[0],
                            delimiter=';',index_col='time')
        if 'width' in track.columns:
            max_width = track['width'].max()
            max_length = track['length'].max()
        else:
            max_width = 15
            max_length = 15
            
        WesternLimit = track['Lon'].min()-(max_width) - 5
        EasternLimit = track['Lon'].max()+(max_width) + 5
        SouthernLimit = track['Lat'].min()-(max_length) - 5
        NorthernLimit = track['Lat'].max()+(max_length) + 5
        
    elif args.choose:
        iu_plevel = input_data.isel({time_indexer:0}).sel({vertical_level_indexer:plevel_Pa}, method='nearest'
                        )[namelist_df.loc['Eastward Wind Component']['Variable']]
        iv_plevel = input_data.isel({time_indexer:0}).sel({vertical_level_indexer:plevel_Pa}, method='nearest'
                        )[namelist_df.loc['Northward Wind Component']['Variable']]
        zeta = vorticity(iu_plevel, iv_plevel).metpy.dequantify()
        
        lat, lon = iu_plevel[latitude_indexer], iu_plevel[longitude_indexer]
        domain_limits = initial_domain(zeta, lat, lon)
        WesternLimit = domain_limits['min_lon']
        EasternLimit = domain_limits['max_lon']
        SouthernLimit = domain_limits['min_lat']
        NorthernLimit = domain_limits['max_lat']
        
    input_data = input_data.sel(
        **{latitude_indexer:slice(SouthernLimit,NorthernLimit),
           longitude_indexer: slice(WesternLimit,EasternLimit)})
    
    return input_data

def get_domain_extreme_values(itime, args, slices_plevel, track=None):
    """
    Retrieves or calculates extreme values (minimum/maximum vorticity, minimum geopotential height,
    and maximum wind speed) within a specified domain at chosen pressure level.

    Parameters:
    - track (pd.DataFrame): Track data potentially containing extreme values.
    - itime (str or pd.Timestamp): The specific time step for which to retrieve or calculate extremes.
    - args (argparse.Namespace): Parsed command-line arguments indicating the slicing method.
    - slices_plevel (tuple): Tuple containing slices of vorticity, geopotential height, and wind speed at chosen pressure level.

    Returns:
    - tuple: Containing minimum/maximum vorticity, minimum geopotential height, and maximum wind speed.
    """
    izeta_plevel_slice, ight_plevel_slice, iwspd_plevel_slice = slices_plevel

    if args.track:
        # Check if 'min_max_zeta_plevel', 'min_max_hgt_plevel' and 'max_wind_plevel' columns exists in the track file.
        # If they exist, then retrieve and convert the value from the track file.  
        # If they do not exist, calculate them.
        if f'{args.track_vorticity}_zeta_{args.level}' in track.columns:
            min_max_zeta = float(track.loc[itime][f'{args.track_vorticity}_zeta_{args.level}'])
        else:
            if args.track_vorticity == "min":
                min_max_zeta = float(izeta_plevel_slice.min())
            else:
                min_max_zeta = float(izeta_plevel_slice.max())

        if f'{args.track_geopotential}_hgt_{args.level}' in track.columns:
            min_max_hgt = float(track.loc[itime][f'{args.track_geopotential}_hgt_{args.level}'])
        else:
            if args.track_geopotential == "min":
                min_max_hgt = float(ight_plevel_slice.min())
            else:
                min_max_hgt = float(ight_plevel_slice.max())

        if f'max_wind_{args.level}' in track.columns:
            max_wind = float(track.loc[itime][f'max_wind_{args.level}'])
        else:
            max_wind = float(iwspd_plevel_slice.max())
    
    else:
        if args.track_vorticity == "min":
            min_max_zeta = float(izeta_plevel_slice.min())

        else:
            min_max_zeta = float(izeta_plevel_slice.max())

        min_max_hgt = float(ight_plevel_slice.min())
        max_wind = float(iwspd_plevel_slice.max())

    return min_max_zeta, min_max_hgt, max_wind