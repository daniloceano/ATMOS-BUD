import pytest
import xarray as xr
import numpy as np
from collections import namedtuple

from src.utils import find_extremum_coordinates

Args = namedtuple('Args', ['track_vorticity', 'track_geopotential'])

def create_test_data():
    # Cria um DataArray 2x2 simples para testes
    data = xr.DataArray(
        np.array([[1, 5],
                  [3, 2]]),
        dims=('lat', 'lon')
    )
    lat = xr.DataArray(np.array([10, 20]), dims='lat')
    lon = xr.DataArray(np.array([100, 110]), dims='lon')
    return data, lat, lon

def test_find_min_vorticity():
    ds_data, lat, lon = create_test_data()
    args = Args(track_vorticity='min', track_geopotential='max')
    variable = 'min_zeta'

    lat_val, lon_val = find_extremum_coordinates(ds_data, lat, lon, variable, args)

    # O mínimo do ds_data é 1 na posição (0,0)
    assert lat_val == 10
    assert lon_val == 100

def test_find_max_vorticity():
    ds_data, lat, lon = create_test_data()
    args = Args(track_vorticity='max', track_geopotential='min')
    variable = 'max_zeta'

    lat_val, lon_val = find_extremum_coordinates(ds_data, lat, lon, variable, args)

    # O máximo do ds_data é 5 na posição (0,1)
    assert lat_val == 10
    assert lon_val == 110

def test_find_min_geopotential():
    ds_data, lat, lon = create_test_data()
    args = Args(track_vorticity='max', track_geopotential='min')
    variable = 'min_hgt'

    lat_val, lon_val = find_extremum_coordinates(ds_data, lat, lon, variable, args)

    # O mínimo do ds_data é 1 na posição (0,0)
    assert lat_val == 10
    assert lon_val == 100

def test_find_max_geopotential():
    ds_data, lat, lon = create_test_data()
    args = Args(track_vorticity='min', track_geopotential='max')
    variable = 'max_hgt'

    lat_val, lon_val = find_extremum_coordinates(ds_data, lat, lon, variable, args)

    # O máximo do ds_data é 5 na posição (0,1)
    assert lat_val == 10
    assert lon_val == 110

def test_find_max_wind():
    ds_data, lat, lon = create_test_data()
    args = Args(track_vorticity='min', track_geopotential='max')
    variable = 'max_wind'

    lat_val, lon_val = find_extremum_coordinates(ds_data, lat, lon, variable, args)

    # Máximo do ds_data é 5 na posição (0,1)
    assert lat_val == 10
    assert lon_val == 110

def test_invalid_variable_raises():
    ds_data, lat, lon = create_test_data()
    args = Args(track_vorticity='min', track_geopotential='max')
    variable = 'invalid_var'

    with pytest.raises(ValueError, match="Invalid variable specified."):
        find_extremum_coordinates(ds_data, lat, lon, variable, args)
