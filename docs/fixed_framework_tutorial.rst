Fixed Framework Tutorial
========================

This section provides a comprehensive guide to using the **Fixed Framework** of ATMOS-BUD. The Fixed framework allows users to analyze the atmospheric budgets within a predefined, stationary domain. This is particularly useful for analyzing systems that do not move significantly, such as convergence zones and other localized phenomena.

In this tutorial, we will use data from the **Reg1 cyclone**, which originated in the southeastern Brazil and was documented in the article *Dias Pinto, J. R., and R. P. da Rocha (2011), The energy cycle and structural evolution of cyclones over southeastern South America in three case studies, J. Geophys. Res., 116, D14112*. The data for Reg1 comes from the NCEP reanalysis and covers the period from August 8 to August 14, 2005.

1. Preparing Your Environment
=============================

Before running the program, ensure that your environment is set up correctly. This includes activating the conda environment and preparing the input files.

1. **Activate the Conda Environment:**
   Ensure that your conda environment is activated before running the program. This will make sure that all required dependencies are available.

   Run the following command to activate the environment:

   .. code-block:: bash

      conda activate atmosbud

2. **Prepare Input Files:**
   You need to ensure that the `box_limits` file and the `namelist` file are correctly configured in the `inputs/` directory. These files define the domain limits for analysis and the settings for the NCEP reanalysis data.

   - **Copy the appropriate `box_limits` file:**
     The `box_limits` file specifies the region of interest for the cyclone. In this tutorial, we will use the `box_limits_Reg1` file for the Reg1 cyclone.

     .. code-block:: bash

        cp inputs/box_limits_Reg1 inputs/box_limits

   - **Copy the `namelist` file for NCEP-R2:**
     The `namelist` file contains the configuration for running the program with the NCEP reanalysis data. Ensure you use the correct namelist configuration for your data.

     .. code-block:: bash

        cp inputs/namelist_NCEP-R2 inputs/namelist

Once these steps are completed, you are ready to run the program with the Fixed Framework.


2. Running ATMOS-BUD with the Fixed Framework
===========================================

To run ATMOS-BUD using the **Fixed Framework**, use the following command in your terminal:

.. code-block:: bash

   python atmos_bud.py path-to-file/your_input_file.nc -f

- `path-to-file/your_input_file.nc`: This is the path to your input NetCDF file containing atmospheric data.
- `-f`: This flag tells ATMOS-BUD to use the **Fixed Framework** for the analysis.

For example, to analyze the cyclone data for Reg1 (as mentioned earlier), you would run the following command:

.. code-block:: bash

   python atmos_bud.py samples/Reg1-Representative_NCEP-R2.nc -f

Ensure that the `namelist` and `box_limits` files are correctly set up in the `inputs/` directory before running the command. This setup allows ATMOS-BUD to compute the atmospheric budgets over the defined domain with the results stored in the `ATMOS-BUD_Results` directory and visualizations generated for further analysis.

3. Understanding the Output Files
=================================

After running the program, several files and plots are generated. Below is an example of the output directory structure and some of the generated visualizations.

- **Figures Directory:**
  The `Figures/` folder contains plots that visualize the defined analysis box and the atmospheric conditions within that box for each time step. 

  Example of a figure showing the domain at 850 hPa for August 8, 2005:
  
  .. image:: figs/box_200508080000.png
     :alt: Box defined for computations
     :width: 500px
     :align: center

  This plot shows the domain used for the analysis of the cyclone, with the following elements:
  - **Geopotential height** (shaded) and **vorticity contours** at 850 hPa.
  - **Maximum wind speed** (triangle), **minimum vorticity** (black circle) and **minimum geopotential height** (black cross) within the domain.

  The domain can be modified using the `-l LEVEL`, `--track-vorticity`, and `--track_geopotential` flags:
  
  - `-l LEVEL`: Allows users to choose the pressure level for the analysis (default is 850 hPa).
  - `--track-vorticity {min,max}`: Tracks the minimum or maximum vorticity (default is minimum).
  - `--track_geopotential {min,max}`: Tracks the minimum or maximum geopotential height (default is minimum).

- **CSV Files:**
  In addition to the plots, ATMOS-BUD generates CSV files containing the diagnostic results. These files are organized in subdirectories by budget category: heat, moisture, and vorticity.

  Each CSV file contains the following terms:

  - **Heat Budget (`heat_terms/`)**:
    - **`dTdt`**: Temperature tendency (K/s).
    - **`Theta`**: Potential temperature (K).
    - **`AdvHTemp`**: Horizontal advection of temperature (K/s).
    - **`AdvVTemp`**: Vertical advection of temperature (K/s).
    - **`Sigma`**: Static stability term (K/Pa).
    - **`ResT`**: Residual of the thermodynamic equation (K/s).
    - **`AdiabaticHeating`**: Estimated diabatic heating (W/kg).

  - **Vorticity Budget (`vorticity_terms/`)**:
    - **`Zeta`**: Relative vorticity (1/s).
    - **`dZdt`**: Vorticity tendency (1/s).
    - **`AdvHZeta`**: Horizontal advection of vorticity (1/s²).
    - **`AdvVZeta`**: Vertical advection of vorticity (1/s²).
    - **`Beta`**: Meridional gradient of the Coriolis parameter (1/m/s).
    - **`vxBeta`**: Meridional advection of planetary vorticity (1/s²).
    - **`DivH`**: Horizontal divergence of the wind (1/s).
    - **`ZetaDivH`**: Term ζ·div(V) (1/s²).
    - **`fDivH`**: Term f·div(V) (1/s²).
    - **`Tilting`**: Tilting term (1/s²).
    - **`ResZ`**: Residual of the vorticity budget (1/s²).

  - **Moisture Budget (`moisture_terms/`)**:
    - **`dQdt`**: Specific humidity tendency (kg/m²/s).
    - **`dQdt_integrated`**: Vertically integrated `dQdt` (kg/m²/s).
    - **`divQ`**: Horizontal divergence of moisture flux (1/s).
    - **`divQ_integrated`**: Vertically integrated `div(Q)` (kg/m²/s).
    - **`WaterBudgetResidual`**: `dQdt_integrated` + `divQ_integrated` (kg/m²/s).

These CSV files allow for further analysis and visualization of the atmospheric budgets for the cyclone or system of interest.