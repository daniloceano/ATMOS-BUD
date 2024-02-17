# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    data_handling.py                                   :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2024/02/16 15:22:52 by daniloceano       #+#    #+#              #
#    Updated: 2024/02/16 21:13:04 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import numpy as np
import xarray as xr
import dask
from metpy.units import units
from src.utils import convert_lon
from src.utils import slice_domain

def load_data(infile, longitude_indexer, args, app_logger):
    """
    Loads data from a specified NetCDF file, handling both single files and multiple GFS files.
    
    Parameters:
        infile (str): Path to the input .nc file or a pattern matching multiple files.
        args: Parsed command-line arguments containing flags and options.
        app_logger (logging.Logger): Logger for recording messages about script progress and issues.
        
    Returns:
        xr.Dataset: The loaded dataset.
        
    Raises:
        FileNotFoundError: If the input file or files specified by infile do not exist.
        Exception: For any other issues encountered during data loading.
    """
    try:
        app_logger.info(f'Loading {infile}...')
        with dask.config.set(array={'slicing': {'split_large_chunks': True}}):
            if args.gfs:
                data = convert_lon(
                    xr.open_mfdataset(
                        infile,
                        engine='cfgrib',
                        parallel=True,
                        filter_by_keys={'typeOfLevel': 'isobaricInhPa'},
                        combine='nested',
                        concat_dim='time'
                    ),
                    longitude_indexer
                )
            else:
                data = convert_lon(xr.open_dataset(infile), longitude_indexer)

        app_logger.info(f'Loaded {infile} successfully!')
        return data

    except FileNotFoundError:
        app_logger.error(f'File not found: {infile}')
        raise
    except Exception as e:
        app_logger.error(f'Failed to load data from {infile}: {e}')
        raise

def preprocess_data(data, df_namelist, args, app_logger):
    """
    Preprocesses the loaded data by sorting, slicing, and adjusting units as necessary.
    
    Parameters:
        data (xr.Dataset): The loaded dataset to preprocess.
        df_namelist (pd.DataFrame): DataFrame containing namelist information such as variable names.
        args: Parsed command-line arguments containing flags and options.
        app_logger (logging.Logger): Logger for recording messages about script progress and issues.
        
    Returns:
        xr.Dataset: The preprocessed dataset.
        
    Raises:
        ValueError: If critical namelist variables are missing or if data preprocessing encounters an issue.
        Exception: For any other issues encountered during data preprocessing.
    """
    longitude_indexer = df_namelist.loc['Longitude']['Variable']
    latitude_indexer = df_namelist.loc['Latitude']['Variable']
    vertical_level_indexer = df_namelist.loc['Vertical Level']['Variable']

    app_logger.info('Preprocessing data...')

    # Ensure critical namelist variables are present
    if not all([longitude_indexer, latitude_indexer, vertical_level_indexer]):
        raise ValueError('Missing critical namelist variables.')

    # Sort data coordinates as data from distinc sources might have different arrangements and this might affect the results from the integrations
    app_logger.debug('Sorting data by longitude, vertical level and latitude...')
    data = data.sortby(longitude_indexer).sortby(vertical_level_indexer).sortby(latitude_indexer)

    # Slice data so the code runs faster
    app_logger.debug('Slicing data...')
    data = slice_domain(data, args, df_namelist)
    
    # Assign lat and lon as radians, for calculations
    app_logger.debug('Assigning lat and lon as radians...')
    data = data.assign_coords({"rlats": np.deg2rad(data[latitude_indexer])})
    data = data.assign_coords({"coslats": np.cos(np.deg2rad(data[latitude_indexer]))})
    data = data.assign_coords({"rlons": np.deg2rad(data[longitude_indexer])})
    
    # Force vertical levels to be in Pa
    app_logger.debug('Force vertical levels to be in Pa...')
    new_pressure = (data[vertical_level_indexer]).metpy.convert_units('Pa') * units('Pa')
    data = data.assign_coords({vertical_level_indexer: new_pressure})
    
    app_logger.info('Done.')

    return data
