Data Object Module
==================

The data object module defines the main `DataObject` class for processing meteorological data and computing atmospheric budget terms including thermodynamic, vorticity, and water budget calculations.

Main Classes
------------

.. automodule:: src.data_object
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

DataObject Class
----------------

The `DataObject` class is the core component for processing meteorological data and computing terms of atmospheric budget equations. It calculates thermodynamic and vorticity terms, with diabatic heating estimated as a residual.

Key Features
------------

* **Variable extraction and unit conversion** from input datasets
* **Thermodynamic budget calculations** including temperature advection and diabatic heating
* **Vorticity budget computations** with geostrophic and ageostrophic components
* **Water budget analysis** with moisture transport and storage terms
* **Coordinate system handling** for latitude, longitude, pressure levels, and time
* **Physical constant calculations** (Coriolis parameter, grid spacing, etc.)

Main Methods
------------

Variable Processing
~~~~~~~~~~~~~~~~~~~

* ``extract_variables()``: Extracts and processes variables from input dataset
* ``convert_units()``: Handles unit conversions using MetPy
* ``calculate_geopotential_height()``: Computes geopotential height from geopotential
* ``calculate_additional_properties()``: Calculates grid spacing and Coriolis parameter

Budget Calculations
~~~~~~~~~~~~~~~~~~~

* ``calculate_thermodynamic_terms()``: Computes temperature budget terms
* ``calculate_vorticity_terms()``: Computes vorticity budget terms  
* ``calculate_water_budget_terms()``: Computes moisture budget terms


Data Structure
--------------

The DataObject contains:

**Input Variables:**
* Temperature, humidity, wind components (u, v)
* Vertical velocity (omega), geopotential height
* Coordinate arrays (lat, lon, pressure, time)

**Calculated Terms:**
* Budget equation components for each variable
* Residual terms and diabatic heating estimates
* Grid properties (dx, dy, Coriolis parameter)

**Metadata:**
* Variable mappings from namelist
* Unit conversions and physical constants
* Processing parameters and logging

Usage Examples
--------------

.. code-block:: python

   from src.data_object import DataObject
   import xarray as xr
   import pandas as pd
   
   # Load input data and tendencies
   input_data = xr.open_dataset('era5_data.nc')
   dTdt = calculate_temperature_tendency(input_data)
   dZdt = calculate_geopotential_tendency(input_data)  
   dQdt = calculate_humidity_tendency(input_data)
   
   # Load variable mapping
   namelist_df = pd.read_csv('namelist.csv', index_col=0)
   
   # Create DataObject instance
   data_obj = DataObject(
       input_data=input_data,
       dTdt=dTdt, 
       dZdt=dZdt,
       dQdt=dQdt,
       namelist_df=namelist_df,
       args=args,
       app_logger=logger
   )
   
   # Access calculated budget terms
   horizontal_temp_adv = data_obj.HorizontalAdvT
   diabatic_heating = data_obj.ResT * Cp_d  # Q = Cp_d * ResT
   vorticity_advection = data_obj.HorizontalAdvZ
