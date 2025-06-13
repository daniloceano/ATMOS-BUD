Running the Program
===================

ATMOS-BUD is designed to analyze the heat, vorticity, and humidity budgets of cyclones, offering three distinct frameworks for setting the analysis domain. Before running the program, ensure you have configured the ``namelist`` file in the ``inputs/namelist``. This file specifies the naming conventions for the variables contained in the NetCDF (``.nc``) file. Examples of ``namelist`` files for ERA5 data are provided, allowing users to copy and adapt them for their specific needs.

**Note:** Users can specify the pressure level (in hPa) for all diagnostics and visualizations, as well as whether to track the minimum or maximum vorticity and geopotential height, by using the ``--level``, ``--track-vorticity``, and ``--track_geopotential`` command-line options. All frameworks fully support these options.

Positional Arguments
---------------------

- ``infile``: The input ``.nc`` file containing temperature, specific humidity, meridional, zonal, and vertical components of the wind, across pressure levels.

Optional Arguments
------------------

- ``-h, --help``: Display the help message and exit.
- ``-f, --fixed``: Compute energetics for a Fixed domain specified by the ``inputs/box_limits`` file.
- ``-t, --track``: Define the box using a track file specified by the ``inputs/track`` file. The track indicates the central point of the system, and an arbitrary box of 15°x15° is constructed.
- ``-c, --choose``: For each time step, the user can interactively choose the domain by clicking on the displayed map.
- ``-l LEVEL, --level LEVEL``: Pressure level (in hPa) at which to perform the analysis and domain selection (default: 850).
- ``--track-vorticity {min,max}``: Select whether to track the minimum (default) or maximum vorticity for all tracking, diagnostics, and plots.
- ``--track_geopotential {min,max}``: Select whether to track the minimum (default) or maximum geopotential height for all tracking, diagnostics, and plots.
- ``-g, --gfs``: Open multiple GFS files at once using the cfgrib engine.
- ``-o OUTNAME, --outname OUTNAME``: Choose a name for saving results (default is the same as infile).
- ``--save_nc_file``:  Wether to save the NetCDF file containing the results for each term, for the entire domain (default is True)
- ``-v, --verbose``: Show debug messages while running.

Fixed Framework
---------------

The Fixed framework utilizes a specified domain to compute the budgets, defined within the ``inputs/box_limits`` file. This method is suitable for analyzing atmospheric systems with low displacement, such as convergence zones. To use this framework, users should refer to the example ``box_limits`` files available in the ``inputs`` directory and configure their own as needed.

To run ATMOS-BUD using the Fixed framework, execute the following command in the terminal:

.. code-block:: bash

    python atmos_bud.py path-to-file/your_input_file.nc -f

Ensure that the ``namelist`` and ``box_limits`` files are correctly set up in the ``inputs`` directory before running the command. This setup allows for the precise calculation of atmospheric budgets over the defined domain, with results stored in the ``ATMOS-BUD_Results`` directory and accompanying visualizations generated for analysis.

Semi-Lagrangian Framework
-------------------------

The Semi-Lagrangian framework in ATMOS-BUD is designed for dynamic analysis of moving atmospheric systems. It utilizes a ``inputs/track`` file, which contains information about the system's position (latitude and longitude) for given time steps. By default, a 15°x15° box is constructed around the system's central position for the analysis. The ``track`` file may also include optional columns specifying the length and width of the desired box, allowing for customization of the domain size.

Track File Format
~~~~~~~~~~~~~~~~~

The track file should include the following columns at a minimum:

- Time (formatted as YYYY-MM-DD-HHMM)
- Latitude (in degrees)
- Longitude (in degrees)

Optional columns for specifying the box dimensions include:

- Length (in degrees)
- Width (in degrees)

These additional parameters allow for greater flexibility in defining the analysis domain around the atmospheric system of interest.

Running the Program with the Semi-Lagrangian Framework
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To execute ATMOS-BUD using the Semi-Lagrangian framework, use the ``-t`` flag as follows:

.. code-block:: bash

    python atmos_bud.py path-to-file/your_input_file.nc -t

This command instructs ATMOS-BUD to dynamically adjust the analysis domain based on the system's position as defined in the ``inputs/track`` file.

.. code-block:: bash

   python atmos_bud.py your_input_file.nc -t --level 700 --track-vorticity max --track_geopotential max

This command instructs ATMOS-BUD to use the Semi-Lagrangian framework, performing all diagnostics and visualizations at the 700 hPa level, and tracking the maximum vorticity and geopotential height values within the analysis domain at each time step.

Output
~~~~~~

Running the program in this framework generates a track file in the results directory (``ATMOS-BUD_Results``). This output file includes detailed tracking information for each time step, such as:

- System position (latitude and longitude)
- Box length and width
- Minimum or maximum vorticity, minimum or maximum geopotential height, and maximum wind speeds within the defined domain at the user-selected pressure level.

An example output file demonstrating this format can be found at ``inputs/track_output-example``. This file serves as a reference for understanding the structure and type of data generated by the Semi-Lagrangian framework.

Important Note on Track File Formatting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is crucial for users to be familiar with the formatting requirements of the track file to ensure accurate analysis. Examples of properly formatted track files are provided in the ``inputs/`` directory. Users are encouraged to refer to these examples when preparing their track files for analysis with the Semi-Lagrangian framework. By adhering to the correct format, users can maximize the efficiency and accuracy of their atmospheric system analyses, ensuring that the domain of interest remains centered on the system throughout the analysis period.

Interactive Framework
---------------------

The Interactive framework within ATMOS-BUD offers an engaging, hands-on approach for analyzing atmospheric systems, allowing users to dynamically choose the analysis domain at each time step.

Running the Program with the Interactive Framework
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To execute ATMOS-BUD using the Interactive framework, use the ``-c`` flag as follows:

.. code-block:: bash

    python atmos_bud.py path-to-file/your_input_file.nc -c

This command opens an interactive graphical interface, guiding users through the data subsetting and domain selection processes for each time step based on real-time visualization of atmospheric data.

Initial Data Subsetting
^^^^^^^^^^^^^^^^^^^^^^^

Upon initiating the Interactive framework, the first step involves data subsetting to define the working domain, optimizing memory usage and computational resources. A window displaying vorticity data at the user-selected pressure level (default: 850 hPa) will guide users in selecting the desired domain:

1. A graphical interface will present vorticity data at the chosen pressure level.
2. Users can subset the data to their working domain directly through this interface, aiding in the efficient use of computational resources.

Domain Selection for Each Time Step
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For every time step in the analysis, the program provides an interactive window where users can define the computational domain using their mouse. This step is crucial for tailoring the analysis to specific atmospheric conditions and phenomena:

1. The window displays key atmospheric variables at the chosen pressure level: vorticity, geopotential height, and wind streamlines. The tracked extreme (minimum or maximum vorticity/geopotential) is determined by the user-defined configuration.
2. Instructions on screen will guide users through the process of selecting the computational domain for each time step.

Output and Replicability
^^^^^^^^^^^^^^^^^^^^^^^^

Similar to the Semi-Lagrangian framework, the Interactive framework generates a track file detailing the chosen domain's parameters for each time step. This feature enhances the replicability of the analysis, allowing for future adjustments to the domain by editing the track file.

Leveraging the Interactive Framework
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Interactive framework is designed to offer researchers and students an intuitive and flexible way to engage with atmospheric data. By allowing for dynamic domain selection based on real-time data visualization, it empowers users to conduct detailed and targeted analyses of atmospheric phenomena. Familiarity with the system under study will significantly enhance the ability to choose the most appropriate domain for analysis, leading to more meaningful and accurate results.
