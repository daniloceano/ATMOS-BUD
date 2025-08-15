Domain Selection Module
=======================

The domain selection module provides comprehensive functions for defining and managing spatial domains for atmospheric budget analysis. It supports three main approaches: fixed domains, storm tracking, and interactive domain selection with visualization capabilities.

Main Functions
--------------

.. automodule:: src.select_domain
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Core Functions
--------------

Domain Definition
~~~~~~~~~~~~~~~~~

**get_domain_limits()** - Main function orchestrating domain selection based on analysis type

* Handles fixed, tracking, and interactive domain selection modes
* Integrates with track files for storm following analysis  
* Returns standardized domain limit dictionaries
* Calculates central coordinates and domain dimensions

Interactive Domain Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**draw_box_map()** - Interactive map-based domain selection interface

* Creates publication-quality meteorological maps
* Overlays streamlines, vorticity, and geopotential height
* Mouse-click interface for corner selection
* Real-time domain visualization and validation

**initial_domain()** - Basic interactive domain selection

* Simplified vorticity map display
* Two-click rectangular domain selection
* Coordinate transformation handling
* User confirmation interface

Visualization Functions
~~~~~~~~~~~~~~~~~~~~~~~

**plot_zeta()** - Vorticity field plotting with optional height contours
**plot_min_max_zeta()** - Storm center identification and marking
**map_decorators()** - Cartographic elements (coastlines, gridlines, labels)
**draw_box()** - Domain boundary visualization

Utility Functions  
~~~~~~~~~~~~~~~~~

**coordXform()** - Coordinate reference system transformations
**tellme()** - User interaction messages
**fmt()** - Coordinate axis formatting

Key Features
------------

Domain Selection Methods
~~~~~~~~~~~~~~~~~~~~~~~~

**Fixed Domains**
* Read predefined boundaries from ``inputs/box_limits`` file
* Consistent analysis regions across different datasets
* Nearest neighbor coordinate matching for grid alignment

**Storm Tracking**
* Dynamic domains following atmospheric features
* Center coordinates from track files with timestamps
* Configurable domain size (width/length) or default 15°×15°
* Automatic domain translation for moving systems

**Interactive Selection**
* Visual domain specification on meteorological maps
* Real-time feedback with overlaid boundaries
* Multiple variables displayed (vorticity, streamlines, heights)
* User confirmation and revision capability

Coordinate System Handling
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Transformation Support:**
* Automatic coordinate reference system detection
* Plate Carrée projection for global datasets
* Proper handling of longitude conventions (0-360° vs -180-180°)
* Grid-aligned domain boundaries

**Spatial Calculations:**
* Domain center coordinates
* Width and length computation
* Area and aspect ratio determination

Storm Center Detection
~~~~~~~~~~~~~~~~~~~~~~

**Vorticity Extrema Finding:**
* Minimum or maximum vorticity identification
* Multiple extrema handling for complex systems
* Visual marking of detected centers
* Integration with tracking algorithms

**Configurable Tracking:**
* User-specified vorticity type (min/max)
* Pressure level selection for analysis
* Geopotential height integration for verification

Interactive Interface
~~~~~~~~~~~~~~~~~~~~~

**User Experience:**
* Clear instruction prompts
* Timeout handling for automated workflows
* Keyboard/mouse input flexibility
* Visual feedback during selection

**Map Quality:**
* High-resolution coastlines and boundaries
* Professional colormap selection (cmocean)
* Customizable map extent and projection
* Streamline visualization for wind patterns

Domain Limit Dictionary Structure
---------------------------------

All domain selection methods return standardized dictionaries:

.. code-block:: python

   domain_limits = {
       'min_lon': float,      # Western boundary (degrees)
       'max_lon': float,      # Eastern boundary (degrees)  
       'min_lat': float,      # Southern boundary (degrees)
       'max_lat': float,      # Northern boundary (degrees)
       'central_lat': float,  # Domain center latitude
       'central_lon': float   # Domain center longitude
   }

Usage Examples
--------------

Fixed Domain Analysis
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from src.select_domain import get_domain_limits
   
   # Configure arguments for fixed domain
   args.fixed = True
   args.track = False
   args.choose = False
   
   # Get domain limits from box_limits file
   domain = get_domain_limits(
       args, 
       u_wind, v_wind, vorticity, geopotential,
       latitude, longitude, timestamp
   )
   
   print(f"Domain: {domain['min_lon']:.1f}°W to {domain['max_lon']:.1f}°E")
   print(f"        {domain['min_lat']:.1f}°S to {domain['max_lat']:.1f}°N")

Storm Tracking Analysis
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd
   
   # Load storm track data
   track_df = pd.read_csv('storm_track.csv', index_col=0)
   
   # Configure for tracking
   args.track = True
   args.track_vorticity = 'min'  # For cyclones
   
   # Get dynamic domain following storm
   domain = get_domain_limits(
       args,
       u_wind, v_wind, vorticity, geopotential,
       latitude, longitude, timestamp,
       track=track_df
   )

Interactive Domain Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Configure for interactive selection
   args.choose = True
   
   # Launch interactive map
   domain = get_domain_limits(
       args,
       u_wind, v_wind, vorticity, geopotential, 
       latitude, longitude, timestamp
   )
   
   # User will see map and select domain interactively
   # Returns selected boundaries

Direct Interactive Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from src.select_domain import draw_box_map, initial_domain
   
   # Basic domain selection
   simple_domain = initial_domain(vorticity_data, lat, lon)
   
   # Advanced interactive selection with full meteorological display
   detailed_domain = draw_box_map(
       u_wind, v_wind, vorticity, geopotential,
       lat, lon, timestamp_str, args
   )

Configuration Files
-------------------

Fixed Domain Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

Create ``inputs/box_limits`` file:

.. code-block:: text

   min_lon;-60.0
   max_lon;-30.0
   min_lat;-40.0
   max_lat;-20.0

Storm Track File Format
~~~~~~~~~~~~~~~~~~~~~~~

CSV file with temporal tracking data:

.. code-block:: text

   timestamp,Lat,Lon,width,length
   2020-01-01 00:00:00,-25.0,-45.0,20.0,15.0
   2020-01-01 06:00:00,-26.0,-44.0,20.0,15.0
   2020-01-01 12:00:00,-27.0,-43.0,20.0,15.0

Integration with Analysis Workflow
-----------------------------------

The domain selection module integrates seamlessly with the main analysis pipeline:

1. **Domain Definition** - Select appropriate method (fixed/track/choose)
2. **Coordinate Extraction** - Get boundary coordinates and center
3. **Data Slicing** - Extract relevant spatial subset  
4. **Budget Calculations** - Perform analysis on selected domain
5. **Results Output** - Save domain information with results

All domain selection methods ensure:
* Consistent coordinate systems
* Grid-aligned boundaries  
* Proper metadata preservation
* Integration with visualization tools
