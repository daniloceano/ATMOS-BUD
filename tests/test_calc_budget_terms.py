# tests/test_calc_budget_terms.py

import pytest
import numpy as np
import xarray as xr
from metpy.units import units
from metpy.calc import vorticity
import pandas as pd
from src.data_object import DataObject
from src.utils import initialize_logging
from src.cli_interface import parse_arguments
from src.data_handling import preprocess_data

# Expected benchmark results for thermodynamic, vorticity, and water budget terms
EXPECTED_ADVHTEMP = np.load('./benchmarks/benchmark_results/advHTemp.npy')
EXPECTED_ADVVTEMP = np.load('./benchmarks/benchmark_results/advVTemp.npy')
EXPECTED_SIGMA = np.load('./benchmarks/benchmark_results/sigma.npy')
EXPECTED_ADIABATICHEATING = np.load('./benchmarks/benchmark_results/adiabaticHeating.npy')

EXPECTED_ADVHZETA = np.load('./benchmarks/benchmark_results/AdvHZeta.npy')
EXPECTED_ADVVZETA = np.load('./benchmarks/benchmark_results/AdvVZeta.npy')
EXPECTED_FDIVH = np.load('./benchmarks/benchmark_results/fDivH.npy')
EXPECTED_TILTING = np.load('./benchmarks/benchmark_results/tilting.npy')

EXPECTED_DIVQ = np.load('./benchmarks/benchmark_results/divQ.npy')
EXPECTED_WATERBUDGETRESIDUAL = np.load('./benchmarks/benchmark_results/WaterBudgetResidual.npy')

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

def test_thermodynamic_terms():
    """
    Test the thermodynamic terms (advHTemp, advVTemp, sigma, adiabaticHeating).
    """
    debug_args = ['samples/Reg1-Representative_NCEP-R2.nc', '-f']
    args = parse_arguments(debug_args)
    app_logger = initialize_logging("./", args)
    namelist_df = pd.read_csv('inputs/namelist_NCEP-R2', sep=';', index_col=0, header=0)

    data, dTdt, dZdt, dQdt = calculate_temporal_derivatives(namelist_df, args, app_logger)

    # Initialize logging and arguments
    debug_args = ['samples/Reg1-Representative_NCEP-R2.nc', '-f']
    args = parse_arguments(debug_args)
    app_logger = initialize_logging("./", args)
    namelist_df = pd.read_csv('inputs/namelist_NCEP-R2', sep=';', index_col=0, header=0)

    # Create DataObject instance with all required parameters
    obj = DataObject(data,
                     dTdt=dTdt,
                     dZdt=dZdt,
                     dQdt=dQdt,
                     namelist_df=namelist_df,
                     app_logger=app_logger,
                     args=args)

    # Compare the results with the benchmarks
    assert np.allclose(obj.AdvHTemp.values, EXPECTED_ADVHTEMP, atol=1e-5), \
        f"Expected {EXPECTED_ADVHTEMP} but got {obj.AdvHTemp.values}"

    assert np.allclose(obj.AdvVTemp.values, EXPECTED_ADVVTEMP, atol=1e-5), \
        f"Expected {EXPECTED_ADVVTEMP} but got {obj.AdvVTemp.values}"

    assert np.allclose(obj.Sigma.values, EXPECTED_SIGMA, atol=1e-5), \
        f"Expected {EXPECTED_SIGMA} but got {obj.Sigma.values}"

    assert np.allclose(obj.AdiabaticHeating.values, EXPECTED_ADIABATICHEATING, atol=1e-5), \
        f"Expected {EXPECTED_ADIABATICHEATING} but got {obj.AdiabaticHeating.values}"

def test_vorticity_terms():
    """
    Test the vorticity terms (AdvHZeta, AdvVZeta, fDivH, tilting).
    """
    # Load data and calculate dZ/dt
    debug_args = ['samples/Reg1-Representative_NCEP-R2.nc', '-f']
    args = parse_arguments(debug_args)
    app_logger = initialize_logging("./", args)
    namelist_df = pd.read_csv('inputs/namelist_NCEP-R2', sep=';', index_col=0, header=0)

    data, dTdt, dZdt, dQdt = calculate_temporal_derivatives(namelist_df, args, app_logger)

    # Initialize logging and arguments
    debug_args = ['samples/Reg1-Representative_NCEP-R2.nc', '-f']
    args = parse_arguments(debug_args)
    app_logger = initialize_logging("./", args)
    namelist_df = pd.read_csv('inputs/namelist_NCEP-R2', sep=';', index_col=0, header=0)

    # Create DataObject instance with all required parameters
    obj = DataObject(data,
                     dTdt=dTdt,
                     dZdt=dZdt,
                     dQdt=dQdt,
                     namelist_df=namelist_df,
                     app_logger=app_logger,
                     args=args)

    # Compare the results with the benchmarks
    assert np.allclose(obj.AdvHZeta.values, EXPECTED_ADVHZETA, atol=1e-5), \
        f"Expected {EXPECTED_ADVHZETA} but got {obj.AdvHZeta.values}"

    assert np.allclose(obj.AdvVZeta.values, EXPECTED_ADVVZETA, atol=1e-5), \
        f"Expected {EXPECTED_ADVVZETA} but got {obj.AdvVZeta.values}"

    assert np.allclose(obj.fDivH.values, EXPECTED_FDIVH, atol=1e-5), \
        f"Expected {EXPECTED_FDIVH} but got {obj.fDivH.values}"

    assert np.allclose(obj.Tilting.values, EXPECTED_TILTING, atol=1e-5), \
        f"Expected {EXPECTED_TILTING} but got {obj.Tilting.values}"

def test_water_budget_terms():
    """
    Test the water budget terms (divQ, WaterBudgetResidual).
    """
    # Load data and calculate dQ/dt
    debug_args = ['samples/Reg1-Representative_NCEP-R2.nc', '-f']
    args = parse_arguments(debug_args)
    app_logger = initialize_logging("./", args)
    namelist_df = pd.read_csv('inputs/namelist_NCEP-R2', sep=';', index_col=0, header=0)

    data, dTdt, dZdt, dQdt = calculate_temporal_derivatives(namelist_df, args, app_logger)

    # Initialize logging and arguments
    debug_args = ['samples/Reg1-Representative_NCEP-R2.nc', '-f']
    args = parse_arguments(debug_args)
    app_logger = initialize_logging("./", args)
    namelist_df = pd.read_csv('inputs/namelist_NCEP-R2', sep=';', index_col=0, header=0)

    # Create DataObject instance with all required parameters
    obj = DataObject(data,
                     dTdt=dTdt,
                     dZdt=dZdt,
                     dQdt=dQdt,
                     namelist_df=namelist_df,
                     app_logger=app_logger,
                     args=args)

    # Compare the results with the benchmarks
    assert np.allclose(obj.divQ.values, EXPECTED_DIVQ, atol=1e-5), \
        f"Expected {EXPECTED_DIVQ} but got {obj.divQ.values}"

    assert np.allclose(obj.WaterBudgetResidual.values, EXPECTED_WATERBUDGETRESIDUAL, atol=1e-5), \
        f"Expected {EXPECTED_WATERBUDGETRESIDUAL} but got {obj.WaterBudgetResidual.values}"

