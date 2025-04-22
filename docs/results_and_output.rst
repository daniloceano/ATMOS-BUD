# Results and Output

Each time **ATMOS-BUD** is executed, it generates a comprehensive set of results organized within a newly created directory inside the `ATMOS-BUD_Results` folder, located in the same directory as the main program. The name of each results directory is based on the input filename, appended with `_method`, where `method` refers to the analysis framework used during the run: either `fixed`, `track`, or `choose`, corresponding to the Fixed, Semi-Lagrangian, or Interactive frameworks.

## Directory Structure and Contents

Within each results directory, you will find the following files and subdirectories:

### `heat_terms/`
Contains CSV files with vertically resolved time series for each term in the heat budget equation:

- `AdvHTemp.csv`: Horizontal advection of temperature  
- `AdvVTemp.csv`: Vertical advection of temperature  
- `dTdt.csv`: Temporal derivative of temperature  
- `Sigma.csv`: Sigma term (static stability)
- `Omega.csv`: Omega term (vertical motion)  
- `ResT.csv`: Residual term (subtaction of all other terms, if multiplied by Cp_d becomes the Adiabatic heating)

### `moisture_terms/`
Contains the following diagnostic files:

- `dQdt.csv`, `dQdt_integrated.csv`: Temporal derivative of specific humidity and its vertical integration  
- `divQ.csv`, `divQ_integrated.csv`: Horizontal divergence of moisture flux and its integrated value  
- `WaterBudgetResidual.csv`, `WaterBudgetResidual_integrated.csv`: Residual of the water vapor budget equation  

### `vorticity_terms/`
Contains terms related to the vorticity budget:

- `dZdt.csv`: Temporal derivative of relative vorticity  
- `AdvHZeta.csv`, `AdvVZeta.csv`: Horizontal and vertical advection of vorticity  
- `ZetaDivH.csv`: Term from divergence acting on vorticity  
- `vxBeta.csv`: Meridional advection of planetary vorticity (v·β)  
- `Tilting.csv`: Tilting term (baroclinic)  
- `ResZ.csv`: Residual of the vorticity budget  
- `Zeta.csv`: Raw relative vorticity data over time  

### `Figures/`
Visual outputs generated during the run:

- `timeseries-min_max_zeta_hgt.png`: Time series of the min/max vorticity and minimum geopotential height at 850 hPa  
- `hovmoller_mean_zeta.png`: Hovmöller diagram of mean relative vorticity across pressure levels and time  
- `track_boxes.png`: Map showing the cyclone track and dynamically selected analysis boxes *(only for `track` and `choose` frameworks)*  

#### `Figures/boxes/`
Contains individual maps for each time step at 850 hPa:

- Shows vorticity, geopotential height, and the selected domain  
- Marks positions of min/max vorticity, minimum height, and maximum wind speed  
- Automatically selects minimum vorticity for Southern Hemisphere and maximum for Northern Hemisphere  

### Track-Related Files *(for `track` and `choose` frameworks)*

- `<input>_track_track.csv`: CSV with the complete storm track and diagnostic fields per time step  
- `<input>_track.nc`: NetCDF file with the full set of spatial variables over time  
- `track`: Auxiliary diagnostic file related to the cyclone track  

### Log File

- `log.<input>`: Records all execution steps, including detected issues and time tracking  

---

## Understanding the Output

The output files provide both **quantitative diagnostics** (CSV, NetCDF) and **qualitative insights** (visualizations) into the evolution of the atmospheric system.

- The modular structure of heat, moisture, and vorticity terms allows researchers to isolate processes, calculate residuals, and evaluate consistency between different atmospheric budgets.  
- Track files ensure replicability in semi-lagrangian and interactive experiments.  
- The log file offers full traceability of execution.  

Together, these outputs form a complete analytical package that supports both in-depth diagnostics and high-level interpretation of atmospheric systems, particularly for synoptic and mesoscale phenomena such as cyclones and monsoon depressions.
