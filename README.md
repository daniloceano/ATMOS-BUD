# Cyclone_Thermodynamics

This program solves the Quasi-Geostrophic Thermodynamic equation for a closed domain on the atmosphere, explicitly computing each term and estimating the diabatic heating term (R<sub>T</sub>) as a residual:

![image](https://user-images.githubusercontent.com/56005607/214878079-a359b897-2388-4197-bd95-3a0d0038ceda.png)

Where the first term in the right hand size is the temperature tendency, the second one is the horiztonal temperature advection and the last one is the vertcial velocity, in pressure levels, times the static stability term.

## Selecting the computatinal domain

### Fixed domain

The most basic way of using the program is by setting a fixed domain in time. For this the user needs to use the -f [flag](#Flags) and provide a text file with its maximum/minimum longitude and latitude values:

![image](https://user-images.githubusercontent.com/56005607/206709581-34ebe0a7-ff45-4bd4-86e0-8cce8dde91ea.png)

### Moving domain

#### Interactively choosing the domain each time step

By using the -c [flag](#Flags), a pop-up window will appear for the user displaying the world map and the vorticity field. The user must first select an area for slicing the global data, which improves visualization and processing time. To select the area, the user needs to click twice on the screen. The first click selects the top-left corner of the area, and the second click selects the bottom-right corner. The order of these clicks does not matter. After this, the pop-up window will display the vorticity data for the selected area only. For each time step, the user will be prompted to select a computational area by clicking on the screen. After finishing the computations, it will be created a "track" file containing the domain central latitude and longitude, its length and width, the minimum vorticity inside the domain and the maximum windspeed.

#### Using a pre-defined domain




## Inputs

### Data

The program opens either NetCDF or Grib data. The input variables are the zonal and meridional wind components components, the vertical velocity in pressure levels (omega), air temperature and geopotential height data. Alternatively, if only geopotential data is available (instead of geopotential height), the user must use the -g [flag](#Flags) and change the [fvars](#Namelist options) text file.

### Auxiliary files
