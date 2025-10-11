Visualization Module
====================

The visualization module provides comprehensive plotting capabilities for atmospheric data visualization in ATMOS-BUD. It offers functions for creating maps, plotting system tracks, generating time series plots, and creating Hovmöller diagrams for atmospheric budget analysis.

Key Features
------------

* **Map Visualization**: Create publication-quality maps with meteorological overlays
* **Track Plotting**: Visualize weather system trajectories with enhanced features
* **Time Series Analysis**: Generate dual-axis plots for vorticity and geopotential height
* **Hovmöller Diagrams**: Display time-pressure cross-sections of atmospheric variables
* **Cartographic Features**: Integrated coastlines, boundaries, and geographic references

Dependencies
------------

The module relies on several specialized libraries:

* **Cartopy**: Geographic projections and map features
* **Matplotlib**: Core plotting functionality
* **Pandas**: Data handling and time series processing
* **NumPy**: Numerical computations
* **cmocean**: Oceanographic colormaps
* **sklearn**: Data normalization utilities

Functions Overview
-------------------

Map Features and Setup
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: src.visualization.map_features

   Adds standard cartographic features to map plots including coastlines, land/ocean boundaries, 
   and geographic reference lines.

   **Usage Example:**
   
   .. code-block:: python
   
      import matplotlib.pyplot as plt
      import cartopy.crs as ccrs
      from src.visualization import map_features
      
      fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
      map_features(ax)

.. autofunction:: src.visualization.Brazil_states

   Adds Brazilian state boundaries to map visualizations with customizable styling.

   **Usage Example:**
   
   .. code-block:: python
   
      from src.visualization import Brazil_states
      
      # Add state boundaries without fill
      Brazil_states(ax, facecolor='None')

Domain and Track Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: src.visualization.plot_fixed_domain

   Creates comprehensive domain visualization with meteorological overlays showing the analysis box,
   critical points (vorticity extrema, geopotential minima, wind maxima), and atmospheric context.

   **Key Features:**
   
   * Analysis domain boundary highlighting
   * Critical point identification and labeling  
   * Meteorological field overlays
   * Automatic timestamp formatting
   * Publication-ready styling
   
   **Usage Example:**
   
   .. code-block:: python
   
      from src.visualization import plot_fixed_domain
      
      limits = {
          'min_lon': -60, 'max_lon': -30,
          'min_lat': -35, 'max_lat': -15
      }
      
      plot_fixed_domain(
          limits=limits,
          data_plevel=pressure_level_data,
          args=configuration_args,
          results_subdirectory='./results',
          time='202301011200',
          app_logger=logger
      )

.. autofunction:: src.visualization.plot_track

   Generates enhanced track visualizations with normalized wind speed markers and vorticity coloring.
   Automatically adjusts map extent and provides start/end point identification.

   **Key Features:**
   
   * Adaptive figure sizing based on track extent
   * Enhanced visualization with wind speed and vorticity
   * Start (A) and end (Z) point markers
   * Automatic extent calculation with buffer zones
   * Integrated colorbar for vorticity values
   
   **Usage Example:**
   
   .. code-block:: python
   
      import pandas as pd
      from src.visualization import plot_track
      
      # Track DataFrame with required columns
      track_data = pd.DataFrame({
          'Lon': [-45.2, -44.8, -44.1],
          'Lat': [-23.1, -22.8, -22.3],
          'min_zeta_850': [1.2e-5, 1.5e-5, 1.8e-5],
          'max_wind_850': [15.2, 18.4, 22.1]
      })
      
      plot_track(
          track=track_data,
          args=args,
          figures_directory='./figures',
          app_logger=logger
      )

Time Series Analysis
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: src.visualization.plot_min_max_zeta_hgt

   Creates dual-axis time series plots showing the evolution of vorticity extrema and geopotential 
   height minima at specified pressure levels.

   **Key Features:**
   
   * Dual-axis plotting for different variable scales
   * Automatic date formatting and tick management
   * Combined legend handling
   * Publication-quality styling
   * Configurable maximum tick count
   
   **Usage Example:**
   
   .. code-block:: python
   
      from src.visualization import plot_min_max_zeta_hgt
      
      # Time series DataFrame with datetime index
      track_plotting = pd.DataFrame({
          'min_zeta_850': vorticity_values,
          'min_hgt_850': geopotential_values
      }, index=datetime_index)
      
      plot_min_max_zeta_hgt(
          track_plotting=track_plotting,
          args=args,
          figs_dir='./figures',
          app_logger=logger,
          max_ticks=10
      )

Hovmöller Diagrams
~~~~~~~~~~~~~~~~~~~

.. autofunction:: src.visualization.hovmoller_mean_zeta

   Generates Hovmöller diagrams displaying time-pressure cross-sections of mean relative vorticity
   within the system domain, with symmetric scaling centered at zero.

   **Key Features:**
   
   * Symmetric color scaling centered at zero
   * Automatic missing data handling
   * Pressure level inversion for proper atmospheric display
   * Publication-ready formatting with proper units
   * Automatic timestamp formatting
   
   **Usage Example:**
   
   .. code-block:: python
   
      from src.visualization import hovmoller_mean_zeta
      import pandas as pd
      
      # DataFrame with pressure levels (Pa) as index, time as columns
      zeta_data = pd.DataFrame(
          data=vorticity_array,  # Shape: (pressure_levels, time_steps)
          index=pressure_levels,  # In Pa
          columns=datetime_array
      )
      
      hovmoller_mean_zeta(
          Zeta=zeta_data,
          figures_subdirectory='./figures',
          app_logger=logger
      )

Integration Examples
--------------------

Complete Visualization Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   """
   Complete visualization workflow for atmospheric budget analysis
   """
   import logging
   import pandas as pd
   from src.visualization import (
       plot_fixed_domain, plot_track, 
       plot_min_max_zeta_hgt, hovmoller_mean_zeta
   )
   
   # Setup logging
   logger = logging.getLogger(__name__)
   
   # 1. Domain visualization
   domain_limits = {
       'min_lon': -65, 'max_lon': -35,
       'min_lat': -35, 'max_lat': -10
   }
   
   plot_fixed_domain(
       limits=domain_limits,
       data_plevel=analysis_data,
       args=config_args,
       results_subdirectory='./results/case_study',
       time='202301151800',
       app_logger=logger
   )
   
   # 2. Track visualization
   track_df = pd.read_csv('system_track.csv')
   plot_track(
       track=track_df,
       args=config_args,
       figures_directory='./results/figures',
       app_logger=logger
   )
   
   # 3. Time series analysis
   timeseries_data = pd.read_csv('track_evolution.csv', index_col=0, parse_dates=True)
   plot_min_max_zeta_hgt(
       track_plotting=timeseries_data,
       args=config_args,
       figs_dir='./results/figures',
       app_logger=logger
   )
   
   # 4. Hovmöller diagram
   vorticity_field = pd.read_pickle('hovmoller_data.pkl')
   hovmoller_mean_zeta(
       Zeta=vorticity_field,
       figures_subdirectory='./results/figures',
       app_logger=logger
   )

Multi-Level Analysis Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   """
   Generate visualizations for multiple pressure levels
   """
   pressure_levels = [850, 700, 500, 300]
   
   for level in pressure_levels:
       # Update configuration for current level
       args.level = level
       
       # Generate domain plot for this level
       plot_fixed_domain(
           limits=analysis_domain,
           data_plevel=data_by_level[level],
           args=args,
           results_subdirectory=f'./results/level_{level}hPa',
           time=current_timestamp,
           app_logger=logger
       )
       
       # Generate time series for this level
       level_timeseries = track_data[[f'min_zeta_{level}', f'min_hgt_{level}']]
       plot_min_max_zeta_hgt(
           track_plotting=level_timeseries,
           args=args,
           figs_dir=f'./results/level_{level}hPa/figures',
           app_logger=logger
       )

Technical Notes
----------------

Color Schemes
~~~~~~~~~~~~~~

The module uses scientifically appropriate color schemes:

* **cmocean.cm.balance**: For symmetric data (vorticity) centered at zero
* **cmocean.cm.deep_r**: For depth-related or intensity data
* **Standard matplotlib**: For basic geographic features

Figure Management
~~~~~~~~~~~~~~~~~~

All functions include proper figure management:

* Automatic figure closing to prevent memory leaks
* Tight layout adjustment for publication quality
* Exception handling with informative error messages
* Comprehensive logging of operations and outputs

Data Requirements
~~~~~~~~~~~~~~~~~~

Functions expect specific data structures:

* **Pandas DataFrames**: For time series and track data
* **Dictionary structures**: For pressure level data
* **Proper datetime indexing**: For time-based visualizations
* **Geographic coordinates**: In decimal degrees (longitude/latitude)

Output Management
~~~~~~~~~~~~~~~~~

* Automatic directory creation as needed
* Standardized filename conventions
* High-resolution PNG output suitable for publication
* Comprehensive logging of all generated files

.. automodule:: src.visualization
   :members:
   :undoc-members:
   :show-inheritance:
