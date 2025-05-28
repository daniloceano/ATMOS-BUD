import pytest
from unittest.mock import MagicMock, patch
import xarray as xr
import pandas as pd

from src.data_object import DataObject

def create_fake_namelist():
    data = {
        'Variable': ['lon', 'lat', 'time', 'plev', 'T', 'u', 'v', 'omega', 'geo', 'q'],
        'Units': ['degrees_east', 'degrees_north', 'hours', 'Pa', 'K', 'm/s', 'm/s', 'Pa/s', 'm^2/s^2', 'kg/kg']
    }
    index = ['Longitude', 'Latitude', 'Time', 'Vertical Level',
             'Air Temperature', 'Eastward Wind Component', 'Northward Wind Component',
             'Omega Velocity', 'Geopotential', 'Specific Humidity']
    return pd.DataFrame(data, index=index)

@patch('src.data_object.DataObject.extract_variables')
@patch('src.data_object.DataObject.calculate_thermodynamic_terms')
@patch('src.data_object.DataObject.calculate_vorticity_terms')
@patch('src.data_object.DataObject.calculate_water_budget_terms')
def test_dataobject_init(mock_water, mock_vorticity, mock_thermo, mock_extract):
    # Mocks para argumentos
    input_data = MagicMock(spec=xr.Dataset)
    dTdt = MagicMock(spec=xr.DataArray)
    dZdt = MagicMock(spec=xr.DataArray)
    dQdt = MagicMock(spec=xr.DataArray)
    namelist_df = create_fake_namelist()
    args = MagicMock()
    app_logger = MagicMock()

    # Instanciar DataObject
    obj = DataObject(input_data, dTdt, dZdt, dQdt, namelist_df, args, app_logger)

    # Verificar se métodos foram chamados
    mock_extract.assert_called_once_with(input_data, namelist_df, args)
    mock_thermo.assert_called_once_with(dTdt)
    mock_vorticity.assert_called_once_with(dZdt)
    mock_water.assert_called_once_with(dQdt)

# Poderia seguir testando cada método isoladamente, mockando dados internos.
