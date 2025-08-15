Command Line Interface
======================

The CLI module provides command-line access to ATMOS-BUD functionality.

Main Functions
--------------

.. automodule:: src.cli_interface
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Key Features
------------

* Command-line argument parsing
* Configuration management
* Batch processing capabilities
* Progress reporting
* Error handling and logging
* Interactive mode support

Command Structure
-----------------

Basic Usage
~~~~~~~~~~~

.. code-block:: bash

   atmos-bud --input data.nc --domain fixed --output results/

Available Arguments
~~~~~~~~~~~~~~~~~~~

Core Arguments:
   * ``--input``: Input NetCDF file path
   * ``--output``: Output directory path
   * ``--domain``: Domain selection method (fixed, track, choose)
   * ``--config``: Configuration file path

Domain Options:
   * ``--fixed``: Use fixed domain analysis
   * ``--track``: Enable storm tracking
   * ``--choose``: Interactive domain selection

Analysis Options:
   * ``--heat``: Calculate heat budget
   * ``--vorticity``: Calculate vorticity budget
   * ``--humidity``: Calculate humidity budget

Output Options:
   * ``--csv``: Export results to CSV format
   * ``--netcdf``: Save results as NetCDF
   * ``--plots``: Generate visualization plots

Configuration Files
-------------------

YAML Configuration
~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   input:
     file: "data/era5_data.nc"
     variables: ["temperature", "humidity", "wind"]
   
   domain:
     type: "fixed"
     bounds:
       lon_min: -60
       lon_max: -30
       lat_min: -40
       lat_max: -20
   
   analysis:
     budgets: ["heat", "vorticity", "humidity"]
     levels: [1000, 850, 500, 200]
   
   output:
     directory: "results/"
     formats: ["csv", "netcdf"]
     plots: true

Usage Examples
--------------

.. code-block:: bash

   # Basic analysis with fixed domain
   atmos-bud --input era5_data.nc --fixed --heat --output results/
   
   # Full analysis with configuration file
   atmos-bud --config analysis_config.yaml
   
   # Interactive domain selection
   atmos-bud --input data.nc --choose --all-budgets
   
   # Storm tracking analysis
   atmos-bud --input data.nc --track --center-coords "lat,lon"
