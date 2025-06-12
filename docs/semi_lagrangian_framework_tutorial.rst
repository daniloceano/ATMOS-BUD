Track Framework Tutorial
##########################

This section provides a comprehensive guide to using the **Track (Semi-Lagrangian) Framework** of ATMOS-BUD. The Track framework allows users to analyze the atmospheric budgets for systems that move significantly, such as cyclonic systems.

In this tutorial, we will use data from the **Reg1 cyclone**, which originated in southeastern Brazil and was documented in the article *Dias Pinto, J. R., and R. P. da Rocha (2011), The energy cycle and structural evolution of cyclones over southeastern South America in three case studies, J. Geophys. Res., 116, D14112*. The data for Reg1 comes from the NCEP reanalysis and covers the period from August 8 to August 14, 2005.


Preparing Your Environment
***************************

Before running the program, ensure that your environment is set up correctly. This includes activating the conda environment and preparing the input files.

1. Activate the Conda Environment:
----------------------------------

Activate your conda environment to ensure all dependencies are available:

.. code-block:: bash

   conda activate atmosbud

2. Prepare Input Files:
-----------------------

The ``track`` file and the ``namelist`` file must be correctly configured in the ``inputs/`` directory. These files define the system's trajectory and the settings for the NCEP reanalysis data.

- **Track file**:  
    The ``track`` file specifies the position (latitude and longitude) of the cyclone for each time step. In this tutorial, we will use the ``track_Reg1-Representative`` file for the Reg1 cyclone.

.. code-block:: bash

   cp inputs/track_Reg1-Representative inputs/track

- **Namelist file for NCEP-R2**:  
    The ``namelist`` file contains the configuration for running the program with the NCEP reanalysis data. Ensure you use the correct namelist configuration for your data.

.. code-block:: bash

   cp inputs/namelist_NCEP-R2 inputs/namelist

Once these steps are completed, you are ready to run the program with the **Track (Semi-Lagrangian) Framework**.

2. Running ATMOS-BUD with the Track Framework
=============================================

To run ATMOS-BUD using the **Track Framework**, use the following command in your terminal:

.. code-block:: bash

   python atmos_bud.py path-to-file/your_input_file.nc -t

Replace ``path-to-file/your_input_file.nc`` with the path to your NetCDF input file. For example, analyzing the Reg1 cyclone data:

.. code-block:: bash

   python atmos_bud.py samples/Reg1-Representative_NCEP-R2.nc -t

Ensure that the ``namelist`` and ``track`` files are properly configured in the ``inputs/`` directory. Results will be stored in the ``ATMOS-BUD_Results`` directory, and visualizations will be generated for further analysis.

Example Terminal Output
-----------------------

When you run the program, you will see detailed logs in the terminal. Below is an example output for the execution of the program, explaining the key components of the process:

.. code-block:: bash

    2025-06-12 16:20:09,468 - atmos_bud - INFO - Loading samples/Reg1-Representative_NCEP-R2.nc...
    2025-06-12 16:20:10,023 - atmos_bud - INFO - Loaded samples/Reg1-Representative_NCEP-R2.nc successfully!
    2025-06-12 16:20:10,023 - atmos_bud - INFO - Preprocessing data...
    2025-06-12 16:20:10,036 - atmos_bud - INFO - Done.
    2025-06-12 16:20:10,037 - atmos_bud - INFO - Computing zeta and temperature tendencies...
    2025-06-12 16:20:10,117 - atmos_bud - INFO - Done.
    2025-06-12 16:20:10,117 - atmos_bud - INFO - Directory where results will be stored: ./Results/Reg1-Representative_NCEP-R2_track
    2025-06-12 16:20:10,118 - atmos_bud - INFO - Directory where figures will be stored: ./Results/Reg1-Representative_NCEP-R2_track/Figures
    2025-06-12 16:20:10,118 - atmos_bud - INFO - Name of the output file with results: Reg1-Representative_NCEP-R2_track


Explanation of Key Terminal Outputs
-----------------------------------

- **Loading and Preprocessing**:  
  The program first loads the input data (``samples/Reg1-Representative_NCEP-R2.nc``), preprocesses it, and then begins the main computation (e.g., computing vorticity and temperature tendencies).
  
.. code-block:: bash

    2025-06-12 16:20:09,468 - atmos_bud - INFO - Loading samples/Reg1-Representative_NCEP-R2.nc...
    2025-06-12 16:20:10,023 - atmos_bud - INFO - Preprocessing data...


- **Time Step Processing**:  
    For each time step, the program calculates atmospheric diagnostics, such as vorticity (`ζ`), geopotential height, and wind speed. It will also generate figures for each time step.

.. code-block:: bash

    2025-06-12 16:20:10,130 - atmos_bud - INFO - Processing time step: 2005-08-08 00Z

After processing each time step, the program will display details like the central latitude and longitude, the domain size, and computed diagnostics:

.. code-block:: bash

    2025-06-12 16:20:10,456 - atmos_bud - INFO - central lat: -22.5, central lon: -45.0, size: 15.0 x 15.0, lon range: -52.5 to -37.5, lat range: -30.0 to -15.0
    2025-06-12 16:20:10,456 - atmos_bud - INFO - 850 hPa diagnostics --> min ζ: -1.71e-05, min geopotential height: 1563, max wind speed: 8.91


- **Generated Figures**:  
    Figures will be created for each time step, showing the defined analysis box and relevant atmospheric diagnostics.

.. code-block:: bash

    2025-06-12 16:20:11,519 - atmos_bud - INFO - Created figure with box defined for computations at box_200508080000.png


- **Completion of Outputs**:  
    After processing all time steps, the program will generate summary plots (e.g., time series, Hovmöller diagrams) and save all output files (CSV, NetCDF) in the appropriate directories.

.. code-block:: bash

    2025-06-12 16:20:37,311 - atmos_bud - INFO - Created figure with track and boxes defined for computations: ./Results/Reg1-Representative_NCEP-R2_track/Figures/track_boxes.png
    2025-06-12 16:20:37,525 - atmos_bud - INFO - Time series plot created and saved: ./Results/Reg1-Representative_NCEP-R2_track/Figures/timeseries_min_zeta_min_hgt_850hPa.png


- **Final Outputs**:  
    Once the analysis is complete, the program will display a message showing the total time taken for the execution and list all the generated output files.

.. code-block:: bash

    2025-06-12 16:20:37,749 - atmos_bud - INFO - ./Results/Reg1-Representative_NCEP-R2_track/Reg1-Representative_NCEP-R2_track.nc created
    2025-06-12 16:20:37,749 - atmos_bud - INFO - --- 28.305777072906494 seconds for running the program ---


By interpreting this output, users can confirm the successful loading of data, processing of each time step, and generation of output files for further analysis.
