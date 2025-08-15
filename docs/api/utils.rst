Utilities Module
================

The utilities module contains essential helper functions and common operations used throughout ATMOS-BUD for logging, coordinate transformations, data processing, and atmospheric feature detection.

Main Functions
--------------

.. automodule:: src.utils
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Core Functions
--------------

Logging and Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~

**initialize_logging()** - Comprehensive logging system setup

* Configurable verbosity levels (DEBUG, INFO, ERROR)
* Dual output: console and log file
* Timestamped entries with level identification
* Application-specific logger separation
* Log file naming based on input data

Data Preprocessing Utilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**convert_lon()** - Longitude coordinate standardization

* Converts longitude from 0°-360° to -180°-180° format
* Automatic coordinate sorting after conversion
* Maintains data integrity during transformation
* Compatible with all xarray datasets

**slice_domain()** - Spatial domain extraction

* Supports fixed, tracking, and interactive domain selection
* Automatic coordinate matching and boundary handling
* Preserves metadata and coordinate attributes
* Optimizes memory usage through selective data loading

Track File Management
~~~~~~~~~~~~~~~~~~~~~

**handle_track_file()** - Storm track data validation and processing

* Time series validation against input datasets
* Spatial boundary checking for track coverage
* Automatic reindexing for temporal alignment
* Comprehensive error handling and logging
* Support for multiple track file formats

Atmospheric Feature Detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**find_extremum_coordinates()** - Meteorological extrema identification

* Vorticity minimum/maximum detection for storm centers
* Geopotential height extrema for pressure systems
* Maximum wind speed identification
* Configurable detection criteria (min/max selection)
* Precise coordinate extraction with grid alignment

**get_domain_extreme_values()** - Domain-specific extreme value extraction

* Integration with track files for pre-computed values
* On-demand calculation for missing track data
* Multi-variable extrema processing (vorticity, height, wind)
* Pressure-level specific analysis
* Optimized for time series processing

Key Features
------------

Coordinate System Management
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Longitude Convention Handling:**
* Automatic detection of longitude format (0-360° vs -180-180°)
* Seamless conversion between conventions
* Proper dateline crossing management
* Consistent coordinate ordering

**Spatial Domain Processing:**
* Grid-aligned domain boundaries
* Nearest neighbor coordinate matching
* Memory-efficient spatial subsetting
* Coordinate metadata preservation

Logging System
~~~~~~~~~~~~~~

**Multi-Level Logging:**
* **DEBUG** - Detailed processing information
* **INFO** - General progress messages  
* **ERROR** - Critical error reporting
* **Console + File** - Dual output streams

**Configuration Options:**
* Verbose mode for detailed debugging
* Quiet mode for production runs
* Custom log file naming conventions
* Timestamp formatting for analysis tracking

Track Data Integration
~~~~~~~~~~~~~~~~~~~~~~

**File Format Support:**
* CSV files with temporal indexing
* Variable column naming flexibility
* Missing data handling and interpolation
* Automatic coordinate system detection

**Validation Features:**
* Temporal range verification
* Spatial coverage checking
* Data quality assessment
* Error reporting and logging

Storm Analysis Tools
~~~~~~~~~~~~~~~~~~~~

**Feature Detection:**
* Vorticity-based storm center identification
* Pressure system tracking via geopotential height
* Wind maximum detection for intensity analysis
* Multi-criteria extrema finding

**Track Integration:**
* Pre-computed track value utilization
* Dynamic calculation for missing data
* Multi-level analysis (different pressure levels)
* Time series consistency checking

Usage Examples
--------------

Logging Setup
~~~~~~~~~~~~~

.. code-block:: python

   from src.utils import initialize_logging
   import argparse
   
   # Configure logging
   parser = argparse.ArgumentParser()
   parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
   parser.add_argument('--infile', required=True, help='Input data file')
   args = parser.parse_args()
   
   # Initialize logging system
   logger = initialize_logging('results/', args)
   
   # Use logger throughout application
   logger.info('Starting atmospheric budget analysis')
   logger.debug('Loading configuration parameters')
   logger.error('Critical error encountered')

Coordinate Conversion
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from src.utils import convert_lon
   import xarray as xr
   
   # Load dataset with 0-360° longitude format
   data = xr.open_dataset('data_0to360.nc')
   
   # Convert to -180 to 180° format
   data_converted = convert_lon(data, 'longitude')
   
   print(f"Original range: {data.longitude.min():.1f} to {data.longitude.max():.1f}")
   print(f"Converted range: {data_converted.longitude.min():.1f} to {data_converted.longitude.max():.1f}")

Domain Slicing
~~~~~~~~~~~~~~

.. code-block:: python

   from src.utils import slice_domain
   import pandas as pd
   
   # Load namelist configuration
   namelist_df = pd.read_csv('inputs/namelist', index_col=0)
   
   # Configure domain selection
   args.fixed = True  # or args.track = True, args.choose = True
   
   # Extract spatial domain
   sliced_data = slice_domain(full_dataset, args, namelist_df)
   
   print(f"Original shape: {full_dataset.dims}")
   print(f"Sliced shape: {sliced_data.dims}")

Track File Processing
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from src.utils import handle_track_file
   import pandas as pd
   
   # Process track file
   track_data = handle_track_file(
       input_data=dataset,
       times=time_series,
       longitude_indexer='longitude',
       latitude_indexer='latitude', 
       app_logger=logger
   )
   
   # Access track information
   for timestamp in track_data.index:
       lat = track_data.loc[timestamp, 'Lat']
       lon = track_data.loc[timestamp, 'Lon']
       print(f"{timestamp}: Storm center at {lat:.1f}°N, {lon:.1f}°E")

Feature Detection
~~~~~~~~~~~~~~~~~

.. code-block:: python

   from src.utils import find_extremum_coordinates, get_domain_extreme_values
   
   # Find vorticity extremum
   lat_center, lon_center = find_extremum_coordinates(
       ds_data=vorticity_slice,
       lat=latitude,
       lon=longitude,
       variable='min_zeta',  # or 'max_zeta', 'max_wind'
       args=args
   )
   
   # Get domain extreme values
   min_zeta, min_hgt, max_wind = get_domain_extreme_values(
       itime=timestamp,
       args=args,
       slices_plevel=(vorticity_slice, height_slice, wind_slice),
       track=track_data
   )
   
   print(f"Storm center: {lat_center:.2f}°N, {lon_center:.2f}°E")
   print(f"Min vorticity: {min_zeta:.2e} s⁻¹")
   print(f"Max wind: {max_wind:.1f} m/s")

Complete Workflow Integration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from src.utils import *
   import xarray as xr
   import pandas as pd
   
   # Setup logging
   logger = initialize_logging('results/', args)
   
   # Load and preprocess data
   data = xr.open_dataset(input_file)
   data = convert_lon(data, 'longitude')
   
   # Load configuration
   namelist = pd.read_csv('inputs/namelist', index_col=0)
   
   # Process track file if needed
   if args.track:
       track = handle_track_file(data, time_series, 'longitude', 'latitude', logger)
   
   # Slice domain
   domain_data = slice_domain(data, args, namelist)
   
   # Process each time step
   for timestamp in time_series:
       # Find atmospheric features
       if args.track:
           extrema = get_domain_extreme_values(timestamp, args, data_slices, track)
       
       logger.info(f"Processed timestamp: {timestamp}")

Error Handling and Validation
------------------------------

The utilities module implements robust error handling:

**File Operations:**
* FileNotFoundError for missing track files
* Validation of file formats and contents
* Graceful handling of corrupted data

**Data Validation:**
* Coordinate system consistency checks
* Time series alignment verification
* Spatial boundary validation
* Missing data detection and reporting

**Logging Integration:**
* Comprehensive error logging with context
* Debug information for troubleshooting
* User-friendly error messages
* Stack trace preservation for development
