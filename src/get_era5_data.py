"""
ERA5 Data Download Module

This module provides functionality for downloading ERA5 atmospheric reanalysis data 
from the Copernicus Climate Data Store (CDS). It handles authentication, data requests, 
and automatic file management for ATMOS-BUD workflows.

Author: Danilo Couto de Souza
Date: 2024
"""

import cdsapi
import os
import logging
from datetime import datetime, timedelta
from typing import List, Optional, Union


def download_era5_data(
    variables: List[str],
    pressure_levels: List[int],
    start_date: str,
    end_date: str,
    area: List[float],
    output_file: str,
    hours: Optional[List[str]] = None,
    logger: Optional[logging.Logger] = None
) -> None:
    """
    Downloads ERA5 reanalysis data from the Copernicus Climate Data Store using the modern CDSAPI.
    
    This function provides a flexible interface for downloading ERA5 atmospheric data
    with configurable variables, pressure levels, time periods, and spatial domains.
    Uses the modern CDSAPI 0.7.6+ syntax with list-based parameters.
    
    Parameters
    ----------
    variables : List[str]
        List of ERA5 variable names to download. Examples:
        - 'temperature'
        - 'u_component_of_wind' 
        - 'v_component_of_wind'
        - 'geopotential'
        - 'vorticity'
        - 'specific_humidity'
    
    pressure_levels : List[int]
        List of pressure levels in hPa. Examples: [1000, 925, 850, 700, 500, 300]
    
    start_date : str
        Start date in 'YYYY-MM-DD' format (e.g., '2023-01-01')
    
    end_date : str  
        End date in 'YYYY-MM-DD' format (e.g., '2023-01-31')
    
    area : List[float]
        Spatial domain as [North, West, South, East] in decimal degrees.
        Example: [20, -80, -60, -20] for South America region
    
    output_file : str
        Output filename for the downloaded NetCDF file
        
    hours : Optional[List[str]], default None
        List of hours in 'HH:MM' format. If None, uses ['00:00', '06:00', '12:00', '18:00']
        
    logger : Optional[logging.Logger], default None
        Logger object for progress tracking and error reporting
        
    Returns
    -------
    None
        Downloads file directly to specified output path
        
    Raises
    ------
    Exception
        If download fails due to authentication, network, or API errors
        
    Notes
    -----
    - Requires valid CDS API credentials in ~/.cdsapirc
    - Uses modern CDSAPI syntax (version 0.7.6+)
    - Automatically handles date range conversion to required format
    - Includes progress monitoring and comprehensive error handling
    
    Examples
    --------
    >>> # Basic usage for atmospheric budget analysis
    >>> variables = ['temperature', 'u_component_of_wind', 'v_component_of_wind', 'geopotential']
    >>> levels = [850, 700, 500, 300]
    >>> download_era5_data(
    ...     variables=variables,
    ...     pressure_levels=levels,
    ...     start_date='2023-01-01',
    ...     end_date='2023-01-31', 
    ...     area=[10, -80, -40, -30],
    ...     output_file='era5_january2023.nc'
    ... )
    
    >>> # Custom time selection
    >>> download_era5_data(
    ...     variables=['vorticity', 'geopotential'],
    ...     pressure_levels=[850],
    ...     start_date='2023-06-15',
    ...     end_date='2023-06-15',
    ...     area=[-10, -70, -30, -40], 
    ...     output_file='era5_single_day.nc',
    ...     hours=['00:00', '12:00']
    ... )
    """
    
    # Set up logging
    if logger is None:
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
    
    # Default hours if not specified
    if hours is None:
        hours = ["00:00", "06:00", "12:00", "18:00"]
    
    try:
        # Log the download initiation
        logger.info(f"🌍 Starting ERA5 data download for {start_date} to {end_date}")
        logger.info(f"📍 Spatial domain: {area}")
        logger.info(f"📊 Variables: {', '.join(variables)}")
        logger.info(f"🎚️ Pressure levels: {pressure_levels} hPa")
        
        # Convert dates to required format and generate date range
        start_dt = datetime.strptime(start_date, '%Y-%m-%d')
        end_dt = datetime.strptime(end_date, '%Y-%m-%d')
        
        # Generate all dates in range
        date_range = []
        current_date = start_dt
        while current_date <= end_dt:
            date_range.append(current_date.strftime('%Y-%m-%d'))
            current_date += timedelta(days=1)
        
        # Prepare the request using modern CDSAPI syntax
        request = {
            "product_type": "reanalysis",
            "variable": variables,
            "year": list(set([date.split('-')[0] for date in date_range])),
            "month": list(set([date.split('-')[1] for date in date_range])),
            "day": list(set([date.split('-')[2] for date in date_range])),
            "time": hours,
            "pressure_level": [str(level) for level in pressure_levels],
            "area": area,  # [North, West, South, East]
            "data_format": "netcdf",
            "download_format": "unarchived"
        }
        
        # Log request details
        logger.info(f"📅 Date range: {len(date_range)} days from {start_date} to {end_date}")
        logger.info(f"⏰ Time steps: {len(hours)} per day ({', '.join(hours)})")
        logger.info(f"📁 Output file: {output_file}")
        
        # Create CDS client and initiate download
        logger.info("🔗 Connecting to Copernicus Climate Data Store...")
        client = cdsapi.Client()
        
        # Retrieve and download data
        logger.info("⬇️ Initiating data retrieval...")
        result = client.retrieve("reanalysis-era5-pressure-levels", request)
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', exist_ok=True)
        
        # Download the file
        logger.info(f"💾 Downloading to {output_file}...")
        result.download(output_file)
        
        # Verify download
        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file) / (1024 * 1024)  # Size in MB
            logger.info(f"✅ Download completed successfully!")
            logger.info(f"📊 File size: {file_size:.2f} MB")
            logger.info(f"📁 Output location: {os.path.abspath(output_file)}")
        else:
            raise Exception("Download completed but file was not found")
            
    except Exception as e:
        error_msg = f"❌ Failed to download ERA5 data: {str(e)}"
        logger.error(error_msg)
        raise Exception(error_msg)


def download_era5_data_legacy():
    """
    Legacy download function for backward compatibility.
    
    Downloads specific ERA5 data for the 2005-08-08 to 2005-08-14 case study
    over the South America region. This function maintains compatibility with
    existing workflows while the main download_era5_data function provides
    more flexibility.
    
    Returns
    -------
    None
        Downloads data to 'system-20050808_ERA5.nc'
        
    Notes
    -----
    This function is maintained for backward compatibility. For new workflows,
    use the main download_era5_data function which provides more flexibility.
    """
    
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    
    # Legacy parameters for the original case study
    variables = [
        "geopotential",
        "specific_humidity", 
        "temperature",
        "u_component_of_wind",
        "v_component_of_wind",
        "vertical_velocity"
    ]
    
    pressure_levels = [
        10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250,
        300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 
        800, 825, 850, 875, 900, 925, 950, 975, 1000
    ]
    
    # Use the main download function with legacy parameters
    download_era5_data(
        variables=variables,
        pressure_levels=pressure_levels,
        start_date='2005-08-08',
        end_date='2005-08-14',
        area=[-17.5, -60, -42.5, -30],  # [North, West, South, East]
        output_file='system-20050808_ERA5.nc',
        hours=["00:00", "06:00", "12:00", "18:00"],
        logger=logger
    )


if __name__ == "__main__":
    # Run the legacy function for backward compatibility
    download_era5_data_legacy()