[![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/release/python-390/)

# Cyclone_Thermodynamics

This program solves the Quasi-Geostrophic Thermodynamic equation for a closed domain on the atmosphere, explicitly computing each term and estimating the diabatic heating term (R<sub>T</sub>) as a residual:

![image](https://user-images.githubusercontent.com/56005607/214878079-a359b897-2388-4197-bd95-3a0d0038ceda.png)

Where the first term in the right hand size is the temperature tendency, the second one is the horiztonal temperature advection and the last one is the vertcial velocity, in pressure levels, times the static stability term.

## Selecting the computatinal domain

### Fixed domain

The most basic way of using the program is by setting a fixed domain in time. For this the user needs to use the -f [flag](#flags) and provide a text file with its maximum/minimum longitude and latitude values:

![image](https://user-images.githubusercontent.com/56005607/206709581-34ebe0a7-ff45-4bd4-86e0-8cce8dde91ea.png)

### Moving domain

#### Interactively choosing the domain each time step

By using the -c [flag](#flags), a pop-up window will appear for the user displaying the world map and the vorticity field. The user must first select an area for slicing the global data, which improves visualization and processing time. To select the area, the user needs to click twice on the screen. The first click selects the top-left corner of the area, and the second click selects the bottom-right corner. The order of these clicks does not matter.

![image](https://user-images.githubusercontent.com/56005607/214921907-e19d0024-08dc-4475-ab65-c953e04e7859.png)


After this, the pop-up window will display the vorticity data for the selected area only. For each time step, the user will be prompted to select a computational area by clicking on the screen.

![image](https://user-images.githubusercontent.com/56005607/214922008-5b7c094f-c160-4415-a528-07cc58730827.png)

After finishing the computations, it will be created a [track](#auxiliary files)  file containing the domain central latitude and longitude, its length and width and the minimum vorticity and maximum windspeed inside the domain.

#### Using a pre-defined domain

By using the -t [flag](#flags), the program will create the computational domain from a [track](#auxiliary files)  text file from the [inputs](inputs) directory. This file should contain, for each time step, the domain central position:

![image](https://user-images.githubusercontent.com/56005607/206721056-61fa32ce-aa5d-4f16-af28-c46ac2a9bf88.png)

The resulting domain will be 15°x15° in size. Optionally, the user can add the length and width columns to the track file, for using a distinct domain size than the default.

## Inputs

### Data

The program opens either NetCDF or Grib data. The input variables are the zonal and meridional wind components components, the vertical velocity in pressure levels (omega), air temperature and geopotential height data. Alternatively, if only geopotential data is available (instead of geopotential height), the user must use the -g [flag](#flags) and change the [fvars](#namelist options) text file.

### Auxiliary files

The [inputs](inputs) directory contains auxiliary files that the program uses, depending on the method chosen for selecting the domain. The [box_limits](inputs/box_limits) file contains the limits for defining a fixed computational domain. The [track](inputs/track) file contains the domain central position for each time step and, optionally, its width and length. The [fvars](inputs/fvars) file contains the naming system used by the input data for each variable and its corresponding units:

![image](https://user-images.githubusercontent.com/56005607/210861069-1c899cc8-860a-4212-bd44-118e308db9bd.png)

The [inputs](inputs) directory contains some examples for using those files.

## Flags

- -m, --moving

Creates a fixed domain (see [Fixed domain](#fixed domain))

- -c, --choose

For each time step, interactively select the domain (see [Moving domain](#moving domain))

- -t, --track

Creates a domain that follows the system (see [Moving domain](#moving domain))

- -o, --outname

Choose a name for saving the results (optional)

- -g, --geopotential

Use this flag when instead of Geopotential Height, Geopotential data is provided. It is required that the [fVars](#auxiliary files) file is adjusted accordingly, for example:

![image](https://user-images.githubusercontent.com/56005607/210860966-713243c8-7447-4661-a33d-a988ab1055cf.png)

