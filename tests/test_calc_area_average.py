# tests/test_calc_area_average.py

import pytest
import xarray as xr
import numpy as np
from src.calculations import CalcAreaAverage

# Expected benchmark results for the area average, manually obtained and stored as .npy file
EXPECTED_BENCHMARK_RESULTS = np.load('./benchmarks/benchmark_results/area_average.npy')  # Load benchmark results from the .npy file

def test_calc_area_average():
    """
    Test the CalcAreaAverage function to ensure it produces the expected results.

    This test verifies that the CalcAreaAverage function correctly computes the area 
    average of temperature (or another variable) over the specified latitudes and longitudes, 
    as defined in the dataset. The result is compared against the pre-generated benchmark results 
    stored in 'benchmark_results/area_average.npy'.

    Steps:
    1. Load the sample dataset.
    2. Select the relevant slice of data for testing.
    3. Add latitude and longitude coordinates in radians.
    4. Calculate the area average using the CalcAreaAverage function.
    5. Compare the result with the pre-generated benchmark using np.allclose.
    
    The test passes if the result is within a small tolerance (1e-5) of the benchmark results.
    """
    
    # Load and preprocess the data
    ds = xr.open_dataset("./samples/Reg1-Representative_NCEP-R2.nc")
    data = ds['TMP_2_ISBL'].isel(initial_time0_hours=0, lv_ISBL3=-1).isel(lat_2=slice(0, 10), lon_2=slice(0, 10))

    # Add necessary latitude and longitude coordinates for the area averaging
    data = data.assign_coords({"rlats": np.deg2rad(data['lat_2'])})
    data = data.assign_coords({"coslats": np.cos(np.deg2rad(data['lat_2']))})
    data = data.assign_coords({"rlons": np.deg2rad(data['lon_2'])})

    # Calculate the area average
    result = CalcAreaAverage(data)

    # Verify that the calculated result is close to the benchmark results
    # np.allclose checks if the values are approximately equal within a given tolerance
    assert np.allclose(result.values, EXPECTED_BENCHMARK_RESULTS, atol=1e-5), \
        f"Test failed: Expected {EXPECTED_BENCHMARK_RESULTS} but got {result.values}"

