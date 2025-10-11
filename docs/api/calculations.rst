Calculations Module
===================

The calculations module contains the main computational functions for performing atmospheric budget analysis, including area averaging, zonal averaging, and the complete budget calculation workflow.

Main Functions
--------------

.. automodule:: src.calculations
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Core Functions
--------------

Statistical Operations
~~~~~~~~~~~~~~~~~~~~~~

* **CalcZonalAverage()** - Computes zonal averages for atmospheric variables
* **CalcAreaAverage()** - Computes area averages with optional zonal averaging
* **perform_calculations()** - Main function orchestrating the complete analysis workflow

Main Workflow Function
~~~~~~~~~~~~~~~~~~~~~~

The ``perform_calculations()`` function is the central orchestrator that:

* Processes meteorological data for each time step
* Creates DataObject instances for budget term calculations  
* Handles domain selection (fixed, tracking, or interactive)
* Computes area averages for all budget terms
* Saves results in CSV and NetCDF formats
* Generates diagnostic plots and visualizations

Key Features
------------

Zonal Averaging
~~~~~~~~~~~~~~~

Computes longitudinal averages of atmospheric variables:

* Maintains all vertical levels and time steps
* Uses proper coordinate weighting
* Based on Brennan & Vincent (1980) methodology

Area Averaging
~~~~~~~~~~~~~~

Computes spatial averages over specified domains:

* Optional zonal averaging preprocessing
* Cosine latitude weighting for proper area integration
* Handles both rectangular and irregular domains

Budget Term Processing
~~~~~~~~~~~~~~~~~~~~~~

The main calculation workflow processes:

**Thermodynamic Terms:**
* ``AdvHTemp`` - Horizontal temperature advection
* ``AdvVTemp`` - Vertical temperature advection  
* ``Sigma`` - Sigma coordinate term
* ``Omega`` - Vertical velocity effects
* ``ResT`` - Thermodynamic residual (diabatic heating)

**Vorticity Terms:**
* ``AdvHZeta`` - Horizontal vorticity advection
* ``AdvVZeta`` - Vertical vorticity advection
* ``ZetaDivH`` - Vorticity stretching term
* ``fDivH`` - Coriolis stretching term
* ``Tilting`` - Tilting term
* ``vxBeta`` - Beta effect term
* ``ResZ`` - Vorticity residual

**Water Budget Terms:**
* ``dQdt`` - Moisture tendency
* ``divQ`` - Moisture flux divergence
* ``WaterBudgetResidual`` - Water budget residual

Domain Handling
~~~~~~~~~~~~~~~

Supports multiple domain selection methods:

* **Fixed domains** - User-specified rectangular regions
* **Storm tracking** - Dynamic domains following atmospheric features
* **Interactive selection** - Manual domain specification

Time Series Processing
~~~~~~~~~~~~~~~~~~~~~~

Processes complete time series with:

* Automatic time step iteration
* Consistent domain tracking across time
* Progressive result accumulation
* Memory-efficient data handling

Output Generation
~~~~~~~~~~~~~~~~~

Produces comprehensive output including:

* **CSV files** - Time series of area-averaged budget terms
* **NetCDF files** - Gridded results preserving spatial structure  
* **Diagnostic plots** - Domain maps, time series, vertical profiles
* **Track files** - Storm center coordinates and characteristics

Usage Examples
--------------

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   from src.calculations import perform_calculations, CalcAreaAverage
   
   # Complete budget analysis workflow
   perform_calculations(
       input_data=era5_dataset,
       namelist_df=variable_mapping,
       dTdt=temperature_tendency,
       dZdt=geopotential_tendency, 
       dQdt=moisture_tendency,
       args=analysis_args,
       app_logger=logger,
       results_dir, figures_dir, output_filename
   )

Area Averaging
~~~~~~~~~~~~~~

.. code-block:: python

   from src.calculations import CalcAreaAverage, CalcZonalAverage
   
   # Compute zonal average first
   zonal_avg = CalcZonalAverage(temperature_data)
   
   # Then area average
   area_avg = CalcAreaAverage(temperature_data, ZonalAverage=True)
   
   # Direct area average without zonal preprocessing
   direct_avg = CalcAreaAverage(temperature_data, ZonalAverage=False)

Mathematical Background
-----------------------

The calculations are based on established atmospheric budget methodologies:

**References:**
* Brennan, F. E., & Vincent, D. G. (1980). Zonal and Eddy Components of the Synoptic-Scale Energy Budget during Intensification of Hurricane Carmen (1974). *Monthly Weather Review*, 108(7), 954-965.

The module implements proper:

* Spherical coordinate system handling
* Area weighting for accurate spatial averaging
* Pressure coordinate vertical integration
* Conservation of physical units throughout calculations
