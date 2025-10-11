ERA5 Data Module
================

The get_era5_data module provides functionality for downloading ERA5 atmospheric reanalysis data from the Copernicus Climate Data Store (CDS). It handles authentication, data requests, and automatic file management for ATMOS-BUD workflows.

Key Features
------------

* **Modern CDSAPI Integration**: Uses the latest CDSAPI 0.7.6 with updated authentication
* **Flexible Data Selection**: Configure variables, pressure levels, and time periods
* **Automatic File Management**: Handles output file naming and organization
* **Error Handling**: Robust error handling for network and API issues
* **Progress Tracking**: Built-in progress monitoring for large downloads

Dependencies
------------

The module requires specific versions for compatibility:

* **cdsapi >= 0.7.6**: Latest version with updated API syntax
* **requests**: HTTP request handling
* **os, logging**: Standard Python libraries for file and logging operations

Authentication Setup
---------------------

Before using this module, you must set up CDS API credentials:

1. **Create CDS Account**: Register at https://cds.climate.copernicus.eu/
2. **Generate API Key**: Get your key from your CDS profile page
3. **Setup Credentials**: Create ``~/.cdsapirc`` file:

   .. code-block:: bash
   
      url: https://cds.climate.copernicus.eu/api/v2
      key: <your-uid>:<your-api-key>

Functions Overview
-------------------

.. autofunction:: src.get_era5_data.download_era5_data

   Downloads ERA5 reanalysis data from the Copernicus Climate Data Store using the modern CDSAPI.
   
   **Key Features:**
   
   * Modern CDSAPI syntax with list-based parameters
   * Flexible variable and level selection
   * Automatic temporal range handling
   * Progress monitoring and error handling
   * Standard NetCDF output format
   
   **Usage Example:**
   
   .. code-block:: python
   
      from src.get_era5_data import download_era5_data
      import logging
      
      # Setup logging
      logger = logging.getLogger(__name__)
      
      # Define download parameters
      variables = ['temperature', 'u_component_of_wind', 'v_component_of_wind']
      pressure_levels = [850, 700, 500, 300]
      
      download_era5_data(
          variables=variables,
          pressure_levels=pressure_levels,
          start_date='2023-01-01',
          end_date='2023-01-31',
          area=[20, -80, -60, -20],  # North, West, South, East
          output_file='era5_data_january2023.nc',
          logger=logger
      )

Data Request Configuration
--------------------------

Variable Selection
~~~~~~~~~~~~~~~~~

The module supports all ERA5 atmospheric variables:

* **Temperature**: temperature, potential_temperature
* **Wind Components**: u_component_of_wind, v_component_of_wind, w_component_of_wind
* **Geopotential**: geopotential, geopotential_height
* **Vorticity**: vorticity, absolute_vorticity, potential_vorticity
* **Humidity**: specific_humidity, relative_humidity
* **Surface Variables**: surface_pressure, mean_sea_level_pressure

Pressure Level Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Standard pressure levels in hPa:

.. code-block:: python

   # Common atmospheric levels
   standard_levels = [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100]
   
   # Tropospheric focus
   tropospheric_levels = [1000, 850, 700, 500, 300]
   
   # Single level analysis
   single_level = [850]  # Must be list format for CDSAPI 0.7.6

Temporal Configuration
~~~~~~~~~~~~~~~~~~~~~~

Flexible time period specification:

.. code-block:: python

   # Single day
   start_date = end_date = '2023-01-15'
   
   # Month-long period
   start_date, end_date = '2023-01-01', '2023-01-31'
   
   # Specific hours (if supported)
   hours = ['00:00', '06:00', '12:00', '18:00']

Spatial Domain Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Area specification follows [North, West, South, East] format:

.. code-block:: python

   # South America focus
   south_america = [15, -85, -60, -30]
   
   # Brazil region
   brazil_region = [10, -75, -35, -30]
   
   # Custom analysis domain
   analysis_domain = [max_lat, min_lon, min_lat, max_lon]

Integration Examples
--------------------

Complete ERA5 Download Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   """
   Complete workflow for downloading ERA5 data for atmospheric budget analysis
   """
   import logging
   from datetime import datetime, timedelta
   from src.get_era5_data import download_era5_data
   
   # Setup logging
   logging.basicConfig(level=logging.INFO)
   logger = logging.getLogger(__name__)
   
   # Define analysis period
   start_date = '2023-01-01'
   end_date = '2023-01-31'
   
   # Analysis domain (South America)
   domain = [15, -85, -60, -30]  # [North, West, South, East]
   
   # Required variables for budget analysis
   variables = [
       'temperature',
       'u_component_of_wind',
       'v_component_of_wind',
       'w_component_of_wind',
       'geopotential',
       'vorticity',
       'specific_humidity'
   ]
   
   # Atmospheric levels for analysis
   levels = [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200]
   
   # Download data
   output_file = f'era5_budget_analysis_{start_date}_{end_date}.nc'
   
   try:
       download_era5_data(
           variables=variables,
           pressure_levels=levels,
           start_date=start_date,
           end_date=end_date,
           area=domain,
           output_file=output_file,
           logger=logger
       )
       
       logger.info(f"✅ Successfully downloaded ERA5 data: {output_file}")
       
   except Exception as e:
       logger.error(f"❌ Failed to download ERA5 data: {e}")

Multi-Case Download with Error Handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   """
   Download ERA5 data for multiple analysis cases with robust error handling
   """
   import time
   from src.get_era5_data import download_era5_data
   
   # Define multiple analysis cases
   analysis_cases = [
       {
           'name': 'Summer_Case_2023',
           'start': '2023-01-15',
           'end': '2023-01-20',
           'domain': [-10, -70, -30, -40],
           'variables': ['temperature', 'geopotential', 'vorticity']
       },
       {
           'name': 'Winter_Case_2023',
           'start': '2023-07-10',
           'end': '2023-07-15', 
           'domain': [-15, -65, -35, -35],
           'variables': ['temperature', 'u_component_of_wind', 'v_component_of_wind']
       }
   ]
   
   pressure_levels = [850, 700, 500, 300]
   
   for case in analysis_cases:
       output_file = f"era5_{case['name']}.nc"
       
       try:
           logger.info(f"🌍 Downloading data for case: {case['name']}")
           
           download_era5_data(
               variables=case['variables'],
               pressure_levels=pressure_levels,
               start_date=case['start'],
               end_date=case['end'],
               area=case['domain'],
               output_file=output_file,
               logger=logger
           )
           
           logger.info(f"✅ Completed: {case['name']}")
           
           # Add delay between requests to be respectful to CDS
           time.sleep(10)
           
       except Exception as e:
           logger.error(f"❌ Failed to download {case['name']}: {e}")
           continue  # Continue with next case

Custom Variable Selection
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   """
   Custom variable selection for specific research needs
   """
   from src.get_era5_data import download_era5_data
   
   # Vorticity and wind analysis
   vorticity_variables = [
       'vorticity',
       'u_component_of_wind',
       'v_component_of_wind',
       'geopotential'
   ]
   
   # Thermodynamic analysis
   thermal_variables = [
       'temperature',
       'potential_temperature',
       'specific_humidity',
       'relative_humidity'
   ]
   
   # Download for different analysis types
   analysis_types = {
       'vorticity_analysis': vorticity_variables,
       'thermal_analysis': thermal_variables
   }
   
   for analysis_name, variables in analysis_types.items():
       output_file = f'era5_{analysis_name}_202301.nc'
       
       download_era5_data(
           variables=variables,
           pressure_levels=[850, 700, 500],
           start_date='2023-01-01',
           end_date='2023-01-31',
           area=[10, -80, -40, -30],
           output_file=output_file,
           logger=logger
       )

Technical Notes
----------------

API Limitations
~~~~~~~~~~~~~~~~

Be aware of CDS API limitations:

* **Request Size**: Large requests may be queued or rejected
* **Rate Limits**: Respect API rate limits to avoid blocking
* **Concurrent Requests**: Limit simultaneous downloads
* **Data Volume**: Monitor your CDS quota usage

Error Handling
~~~~~~~~~~~~~~~

The module handles various error conditions:

* **Authentication Errors**: Invalid or missing API credentials
* **Network Issues**: Connection timeouts and interruptions  
* **API Errors**: Server-side processing failures
* **File System Errors**: Disk space and permission issues

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

For optimal performance:

* **Request Sizing**: Balance between large requests and API limits
* **Sequential Downloads**: Process requests sequentially to avoid conflicts
* **Local Caching**: Avoid re-downloading existing data
* **Progress Monitoring**: Use logging to track download progress

Data Quality
~~~~~~~~~~~~~

ERA5 data characteristics:

* **Temporal Resolution**: Hourly data available
* **Spatial Resolution**: Approximately 31 km (0.25° × 0.25°)
* **Vertical Levels**: 37 pressure levels from 1000 to 1 hPa
* **Data Format**: NetCDF with CF conventions
* **Quality Control**: Extensive quality assurance in reanalysis

.. automodule:: src.get_era5_data
   :members:
   :undoc-members:
   :show-inheritance:
