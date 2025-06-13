# benchmarks/generate_benchmarks.py

import xarray as xr
import numpy as np
import os
import sys
import logging 
import pandas as pd

from metpy.units import units
from metpy.calc import vorticity

# Add src to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.calculations import CalcZonalAverage, CalcAreaAverage
from src.data_object import DataObject
from src.utils import initialize_logging
from src.cli_interface import parse_arguments
from src.data_handling import preprocess_data

def load_and_preprocess_data():
    """
    Load the sample dataset and preprocess it by selecting the necessary slices of data
    and adding the required latitude and longitude coordinates in radians.
    """
    # Load the sample dataset
    ds = xr.open_dataset("./samples/Reg1-Representative_NCEP-R2.nc")
    
    # Select the relevant data slice (specific time and level)
    data = ds['TMP_2_ISBL'].isel(initial_time0_hours=0, lv_ISBL3=-1).isel(lat_2=slice(0, 10), lon_2=slice(0, 10))

    # Add latitude and longitude coordinates in radians and their cosine for zonal averaging
    data = data.assign_coords({"rlats": np.deg2rad(data['lat_2'])})
    data = data.assign_coords({"coslats": np.cos(np.deg2rad(data['lat_2']))})
    data = data.assign_coords({"rlons": np.deg2rad(data['lon_2'])})

    return data

def generate_zonal_average_benchmark():
    """
    Generate the benchmark for the CalcZonalAverage function.
    """
    data = load_and_preprocess_data()

    # Calculate the zonal average
    result = CalcZonalAverage(data)

    # Save the result as a benchmark
    np.save('./benchmarks/benchmark_results/zonal_average.npy', result.values)
    print("Zonal Average benchmark generated.")

def generate_area_average_benchmark():
    """
    Generate the benchmark for the CalcAreaAverage function.
    """
    data = load_and_preprocess_data()

    # Calculate the area average (assumed function)
    result = CalcAreaAverage(data)

    # Save the result as a benchmark
    np.save('./benchmarks/benchmark_results/area_average.npy', result.values)
    print("Area Average benchmark generated.")

def calculate_temporal_derivatives(namelist_df, args, app_logger):
    """
    Calculate the temporal derivatives (tendencies) for temperature, vorticity, and humidity.
    """
    app_logger.info('🔄 Starting the calculation of temporal derivatives for temperature, vorticity, and humidity...')
    
    # Load the sample dataset
    ds = xr.open_dataset("./samples/Reg1-Representative_NCEP-R2.nc")

    data = preprocess_data(ds, namelist_df, args, app_logger)
    
    # Select the relevant data slice (specific time and level)
    data = ds.isel(initial_time0_hours=slice(0, 2), lat_2=slice(0, 10), lon_2=slice(0, 10)).sel(lv_ISBL3=slice(925, 1000))

    # Add latitude and longitude coordinates in radians and their cosine for zonal averaging
    data = data.assign_coords({"rlats": np.deg2rad(data['lat_2'])})
    data = data.assign_coords({"coslats": np.cos(np.deg2rad(data['lat_2']))})
    data = data.assign_coords({"rlons": np.deg2rad(data['lon_2'])})

    # Calculate the temperature tendency (dT/dt)
    dTdt = data['TMP_2_ISBL'].differentiate('initial_time0_hours', datetime_unit='s') * units('K/s')
    
    # Calculate vorticity (zeta)
    Zeta = vorticity(data['U_GRD_2_ISBL'], data['V_GRD_2_ISBL'])
    
    # Calculate the vorticity tendency (dZ/dt)
    dZdt = Zeta.differentiate('initial_time0_hours', datetime_unit='s') / units('s')
    
    # Calculate the humidity tendency (dQ/dt)
    dQdt = data['q_2_ISBL'].differentiate('initial_time0_hours', datetime_unit='s') * units('kg/kg/s')
        
    return data, dTdt, dZdt, dQdt

def generate_budget_benchmark(data, dTdt, dZdt, dQdt, args, namelist_df, app_logger):
    """
    Generate the benchmark for the terms of heat, vorticity, and water budget.
    """
    # Create DataObject instance
    obj = DataObject(data,
                     dTdt=dTdt,
                     dZdt=dZdt,
                     dQdt=dQdt,
                     namelist_df=namelist_df,
                     app_logger=app_logger,
                     args=args)

    # Save thermodynamic results
    np.save('./benchmarks/benchmark_results/advHTemp.npy', obj.AdvHTemp.values)
    np.save('./benchmarks/benchmark_results/advVTemp.npy', obj.AdvVTemp.values)
    np.save('./benchmarks/benchmark_results/sigma.npy', obj.Sigma.values)
    np.save('./benchmarks/benchmark_results/adiabaticHeating.npy', obj.AdiabaticHeating.values)
    print("Thermodynamic terms benchmark generated.")

    # Save vorticity results
    np.save('./benchmarks/benchmark_results/AdvHZeta.npy', obj.AdvHZeta.values)
    np.save('./benchmarks/benchmark_results/AdvVZeta.npy', obj.AdvVZeta.values)
    np.save('./benchmarks/benchmark_results/fDivH.npy', obj.fDivH.values)
    np.save('./benchmarks/benchmark_results/tilting.npy', obj.Tilting.values)
    print("Vorticity terms benchmark generated.")

    # Save water budget results
    np.save('./benchmarks/benchmark_results/divQ.npy', obj.divQ.values)
    np.save('./benchmarks/benchmark_results/WaterBudgetResidual.npy', obj.WaterBudgetResidual.values)
    print("Water budget terms benchmark generated.")

def generate_readme():
    """
    Generate the README file automatically with information about the dataset,
    the preprocessing steps, and general details about the benchmark functions.
    """
    data = load_and_preprocess_data()

    # Extract relevant information about the dataset
    lat_min = data['lat_2'].min().values
    lat_max = data['lat_2'].max().values
    lon_min = data['lon_2'].min().values
    lon_max = data['lon_2'].max().values

    # Write the README file
    with open('./benchmarks/README.md', 'w') as f:
        f.write("# Benchmarks for ATMOS-BUD\n\n")
        f.write("This directory contains the benchmark results for various functions used in the ATMOS-BUD project. These benchmarks are calculated using data from the NCEP-R2 dataset for the specified date, vertical level, and latitude/longitude ranges.\n\n")
        
        f.write("## How to Obtain the Data\n\n")
        f.write("The data used for generating the benchmarks is from the NCEP-R2 (National Centers for Environmental Prediction - Reanalysis 2) dataset. The specific data for the following parameters is required:\n\n")
        f.write(f"- **Initial time**: 2005-08-12T12\n")
        f.write(f"- **Vertical level**: 1000 hPa\n")
        f.write(f"- **Latitude range**: {lat_min} to {lat_max} degrees\n")
        f.write(f"- **Longitude range**: {lon_min} to {lon_max} degrees\n\n")
        
        f.write("You can access the NCEP-R2 dataset from the [NOAA ESRL website](https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.derived.html). Download the dataset for the required date and variables.\n\n")
        
        f.write("## How to Generate the Benchmarks\n\n")
        f.write("To generate the benchmark results, run the `generate_benchmarks.py` script. This script will calculate the benchmarks for various functions such as `CalcZonalAverage` and `CalcAreaAverage` using the NCEP-R2 data.\n\n")
        f.write("### Steps:\n")
        f.write("1. **Download the NCEP-R2 dataset** as mentioned above.\n")
        f.write("2. Place the dataset in the `samples/` directory.\n")
        f.write("3. Run the following command to generate the benchmarks:\n\n")
        f.write("```bash\n")
        f.write("python benchmarks/generate_benchmarks.py\n")
        f.write("```\n\n")
        
        f.write("## How to Use the Benchmarks\n\n")
        f.write("These benchmark files are used in the tests to validate that the implemented functions (like `CalcZonalAverage` and `CalcAreaAverage`) are producing correct results. During testing, the computed results are compared to these benchmarks to ensure the accuracy of the functions.\n\n")
        f.write("Run the tests with:\n\n")
        f.write("```bash\n")
        f.write("pytest\n")
        f.write("```\n")
        f.write("This will automatically compare the test results with the stored benchmarks.\n")

    print("README generated successfully.")

def main():
    """
    Generate all benchmarks and the README by calling the appropriate functions.
    """
    debug_args = ['samples/Reg1-Representative_NCEP-R2.nc', '-f']
    args = parse_arguments(debug_args)
    app_logger = initialize_logging("./", args)
    namelist_df = pd.read_csv('inputs/namelist_NCEP-R2', sep=';', index_col=0, header=0)

    generate_zonal_average_benchmark()
    generate_area_average_benchmark()

    data, dTdt, dZdt, dQdt = calculate_temporal_derivatives(namelist_df, args, app_logger)
    generate_budget_benchmark(data, dTdt, dZdt, dQdt, args, namelist_df, app_logger) 

    generate_readme()


if __name__ == '__main__':
    main()
