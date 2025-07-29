Downloading Data
================

ATMOS-BUD requires a single NetCDF file containing the necessary meteorological data as input. The data should include geopotential, temperature, specific humidity, vertical velocity, and the u and v components of the wind. It is essential that the data be organized in a structured grid format and that the vertical coordinate be on isobaric (pressure) levels.

In the following example, we will guide you through the process of downloading data for the ECMWF's ERA5 dataset, focusing on an atmospheric system that occurred in the Southern Atlantic region near Southeast Brazil in 2005.

Method 1: Copernicus Interface
-------------------------------

To download ERA5 data from the Copernicus Climate Data Store (CDS) using their web interface, follow these steps:

1. Access the Copernicus Climate Data Store (CDS) website at `Copernicus CDS <https://cds.climate.copernicus.eu/#!/home>`_. On the top right, click to create an account or log into your existing account.

2. Visit the ERA5 dataset page at `ERA5 dataset page <https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels?tab=overview>`_ and click on "Download".

3. Select the necessary variables for your analysis:

   - **Geopotential**
   - **Temperature**
   - **Specific Humidity**
   - **Vertical Velocity**
   - **U-Component of Wind**
   - **V-Component of Wind**

4. Choose the time range of interest (year, month, day). Make sure to select the dates that correspond to the atmospheric event you wish to analyze.

5. Define the temporal resolution. Generally, 3-hourly data is suitable for most atmospheric analyses, but you can choose hourly or daily data depending on your needs.

6. Select the desired vertical levels for your analysis. Generally, the pressure levels of interest are from 1000 hPa to 100 hPa.

7. Select sub-region extraction and specify the spatial extent by selecting the geographic region that covers the atmospheric system you are studying.

8. Choose the output file format. Select 'NetCDF' as the preferred format for the output data.

9. Select download format as "Unarchived" to ensure the data is downloaded in a single NetCDF file without compression.

10. After all selections have been made, submit the data retrieval request. You will receive the dataset in the specified format once the request has been processed.

11. Once you have downloaded the dataset, it is advisable to rename the file for easier identification and organization. Use a convention that includes the system identifier and dataset name, such as ``system-identifier_dataset-name.nc``. For example, rename your downloaded file to ``system-20050808_ERA5.nc`` to reflect the system date and data source.

Method 3: NCAR's Research Data Archive
---------------------------------------

For users who prefer to access data from the National Center for Atmospheric Research (NCAR), you can download ERA5 data from their Research Data Archive. Follow these steps:

1. Visit the ERA5 dataset page on the NCAR Research Data Archive website at `NCAR RDA <https://rda.ucar.edu/datasets/d633000/dataaccess/>`_.

2. Click on the "Get a Subset" from the "ERA5 atmospheric pressure level analysis [netCDF4]" dataset

3. On "Temporal Selection", select the time range (start and end date) that corresponds to the atmospheric event you are interested in. You can specify the year, month, and day.
 - Here, we cannot select the temporal resolution, only the start and end hours. So, if the user wishes to use a distinct temporal resolution, they will need to download the data and process it using Climate Data Operators (CDO).

4. Select the necessary variables for your analysis:

   - **Geopotential**
   - **Temperature**
   - **Specific Humidity**
   - **Vertical Velocity**
   - **U Component of Wind**
   - **V Component of Wind**

5. Click on "Continue" to proceed to the next step.

6. Check whether the selected variables are correct and select the vertical levels you need for your analysis. The pressure levels of interest are generally from 1000 hPa to 100 hPa.

7. Make a spatial selection by clicking on "DRAW BOX" bellow the map and either select the desired region or enter the latitude and longitude coordinates manually. Ensure that the selected region covers the atmospheric system you are studying.

8. Click on "Submit Request" to initiate the data retrieval process.

 - Note that the system might retrieve multiple files, depending on the time range and spatial extent you selected. If you encounter multiple files, you can merge them using Climate Data Operators (CDO) or similar tools.

Method 3: Automated Download Using `get_era5_data.py`
------------------------------------------------------

For those who prefer an automated approach to downloading ERA5 data, ATMOS-BUD includes a Python script that utilizes the ECMWF's `cdsapi` library. This method streamlines the data retrieval process, allowing you to specify your data requirements directly within the script.

1. Navigate to the `src` directory where the `get_era5_data.py` script is located.

.. code-block:: bash

    cd src

2. Open the `get_era5_data.py` script in your preferred text editor. Inside the script, locate the section where input parameters are specified, such as pressure levels, dates, and spatial coverage. Modify these inputs to match your specific data requirements.

4. Once you have configured the script with your desired parameters, save the changes, and run the script from the terminal:

.. code-block:: bash

    python get_era5_data.py

The script will automatically communicate with the Copernicus Climate Data Store, request the specified data, and download it in NetCDF format to your local machine. This method is particularly useful for automating the data retrieval process and can be easily integrated into larger data processing workflows.


.. note::
   Ensure you have the `cdsapi` library installed and configured with your CDS API key before running the script. For more information on setting up `cdsapi`, visit the official `cdsapi` installation guide and user documentation: `https://cds.climate.copernicus.eu/how-to-api`

.. note::
    Regardless of the method used, it is advisable to rename the downloaded file, choosing a name that facilitates easy identification and organization of your datasets. We recommend adopting a naming convention that incorporates the system identifier and the dataset name. For example, naming your output file as ``system-20050808_ERA5.nc`` can be beneficial, where ``20050808`` denotes the date of the atmospheric event under analysis, and ``ERA5`` represents the dataset source. This practice ensures your files are organized systematically, enhancing the efficiency of data management and retrieval for future analyses.

.. note::
    If you encounter the error message "Your request is too large, please reduce your selection" when downloading data, this indicates that you are trying to download too much data in a single request. To resolve this issue, it is recommended to download data for each day separately in a loop and then merge the individual files using the Climate Data Operators (CDO) tool. You can merge multiple NetCDF files chronologically using the following command:

    .. code-block:: bash

        cdo -f nc mergetime system-*_ERA5.nc system-merged_ERA5.nc

    This command will merge all files matching the pattern ``system-*_ERA5.nc`` in chronological order into a single output file named ``system-merged_ERA5.nc``. Make sure to have CDO installed on your system before using this command. For installation instructions and more information about CDO, visit: `https://code.mpimet.mpg.de/projects/cdo <https://code.mpimet.mpg.de/projects/cdo>`_

.. note::
