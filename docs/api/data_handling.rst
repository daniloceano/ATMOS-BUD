Data Handling Module
====================

The data handling module provides essential functions for loading, preprocessing, and preparing atmospheric data for budget analysis. It handles both single NetCDF files and multiple GFS files with proper coordinate system management.

Main Functions
--------------

.. automodule:: src.data_handling
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Core Functions
--------------

Data Loading
~~~~~~~~~~~~

**load_data()** - Primary function for loading atmospheric data from NetCDF files

* Supports single NetCDF files and multi-file GFS datasets
* Handles GRIB format conversion using cfgrib engine  
* Implements Dask parallel processing for large datasets
* Automatic longitude coordinate conversion
* Robust error handling and logging

Data Preprocessing  
~~~~~~~~~~~~~~~~~~

**preprocess_data()** - Comprehensive preprocessing pipeline for atmospheric data

* Unit standardization (pressure levels to Pascal)
* Coordinate sorting for consistent data arrangement
* Domain slicing for computational efficiency
* Radian coordinate assignment for mathematical calculations
* Cosine latitude weighting preparation

Key Features
------------

Multi-Format Support
~~~~~~~~~~~~~~~~~~~~

**NetCDF Files:**
* Standard atmospheric reanalysis data (ERA5, NCEP, etc.)
* Single file or concatenated multi-file datasets
* Automatic metadata preservation

**GFS GRIB Files:**
* Operational weather model data
* Multi-file time series handling
* Pressure level filtering (isobaricInhPa)
* Parallel loading with nested concatenation

Coordinate Management
~~~~~~~~~~~~~~~~~~~~~

**Longitude Conversion:**
* Automatic detection of longitude conventions (0-360° vs -180-180°)
* Standardization to consistent coordinate system
* Proper handling of dateline crossing

**Coordinate Sorting:**
* Longitude, latitude, and pressure level ordering
* Ensures consistent integration results
* Handles data from different sources uniformly

**Radian Coordinates:**
* Conversion to radians for mathematical operations
* Cosine latitude weighting for area calculations
* Proper spherical coordinate handling

Data Optimization
~~~~~~~~~~~~~~~~~

**Domain Slicing:**
* Reduces memory footprint by extracting relevant regions
* Faster processing for regional analysis
* Configurable through command-line arguments

**Unit Standardization:**
* Pressure levels converted to Pascal (Pa)
* Consistent physical units throughout analysis
* MetPy integration for unit handling

**Dask Integration:**
* Lazy loading for large datasets
* Parallel processing capabilities
* Memory-efficient chunked operations
* Large chunk splitting configuration

Error Handling
~~~~~~~~~~~~~~

Comprehensive error management:

* **FileNotFoundError** - Missing input files
* **ValueError** - Invalid namelist configurations  
* **Exception** - General processing errors
* Detailed logging for debugging and monitoring

Preprocessing Pipeline
----------------------

The ``preprocess_data()`` function implements a standardized pipeline:

1. **Validation** - Check critical namelist variables
2. **Unit Conversion** - Standardize pressure coordinates to Pa
3. **Sorting** - Order coordinates consistently  
4. **Domain Slicing** - Extract relevant spatial/temporal regions
5. **Radian Assignment** - Add mathematical coordinate systems
6. **Weighting Preparation** - Compute cosine latitude factors


Usage Examples
--------------

Basic Data Loading
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from src.data_handling import load_data, preprocess_data
   import pandas as pd
   import argparse
   import logging
   
   # Setup logging
   logger = logging.getLogger(__name__)
   
   # Load atmospheric data
   dataset = load_data(
       infile='era5_data.nc',
       longitude_indexer='longitude',
       args=args,
       app_logger=logger
   )

GFS Multi-File Loading
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Load GFS GRIB files (args.gfs = True)
   gfs_dataset = load_data(
       infile='gfs_*.grib2',
       longitude_indexer='longitude', 
       args=args,
       app_logger=logger
   )

Complete Preprocessing
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Load namelist configuration
   namelist_df = pd.read_csv('namelist.csv', index_col=0)
   
   # Preprocess the loaded data
   processed_data = preprocess_data(
       data=dataset,
       df_namelist=namelist_df,
       args=args,
       app_logger=logger
   )
   
   # Access processed coordinates
   print(f"Pressure levels: {processed_data.level.values}")
   print(f"Latitude range: {processed_data.latitude.values[[0,-1]]}")
   print(f"Radian coordinates available: {'rlats' in processed_data.coords}")

Data Pipeline Integration
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Complete data preparation workflow
   def prepare_atmospheric_data(input_file, namelist_path, args, logger):
       # Load namelist
       namelist_df = pd.read_csv(namelist_path, index_col=0)
       longitude_var = namelist_df.loc['Longitude']['Variable']
       
       # Load and preprocess data
       raw_data = load_data(input_file, longitude_var, args, logger)
       processed_data = preprocess_data(raw_data, namelist_df, args, logger)
       
       return processed_data, namelist_df

Supported Data Sources
----------------------

ATMOS-BUD can work with **any atmospheric dataset** that contains the required meteorological variables, as long as the ``inputs/namelist`` file is configured correctly to map the variable names and coordinate systems.

**Dataset Flexibility:**
* Any NetCDF or GRIB format atmospheric dataset
* Custom variable names supported through namelist configuration
* Flexible coordinate system handling (longitude, latitude, pressure, time)
* Automatic unit conversion and standardization

**Configuration Requirements:**

To use any dataset, configure the ``inputs/namelist`` file as follows:

.. code-block:: text

    ;standard_name;Variable;Units
    Air Temperature;air_temperature;T;K
    Geopotential;geopotential;Z;m**2/s**2
    Specific Humidity;specific_humidity;Q;kg/kg
    Omega Velocity;omega;W;Pa/s
    Eastward Wind Component;eastward_wind;U;m/s
    Northward Wind Component;northward_wind;V;m/s
    Longitude;;longitude
    Latitude;;latitude
    Time;;time
    Vertical Level;;level


**Required Variables:**
* Temperature
* Specific humidity
* Horizontal wind components (u, v)
* Vertical velocity (omega)
* Geopotential or geopotential height
* Coordinate arrays (longitude, latitude, pressure, time)
