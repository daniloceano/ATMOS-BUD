import pytest
from unittest.mock import MagicMock, patch
import numpy as np
import pandas as pd

from src.data_handling import preprocess_data

# Função auxiliar units que a função usa (assumindo que vem do metpy.units)
from metpy.units import units

# Simulando df_namelist para teste
def create_dummy_namelist():
    data = {'Variable': ['lon', 'lat', 'plev']}
    index = ['Longitude', 'Latitude', 'Vertical Level']
    return pd.DataFrame(data, index=index)

@patch('src.data_handling.slice_domain')
def test_preprocess_data_success(mock_slice_domain):
    # Configuração dos mocks e inputs
    df_namelist = create_dummy_namelist()
    args = MagicMock()
    app_logger = MagicMock()

    # Mock para data e seus métodos encadeados
    data = MagicMock()
    # Simulando o vertical level data com metpy
    vertical_data = MagicMock()
    # metpy.convert_units retorna um objeto que pode ser multiplicado por units('Pa')
    vertical_data.metpy.convert_units.return_value = vertical_data
    vertical_data.__mul__.return_value = 'new_pressure'

    # data[vertical_level_indexer] retorna vertical_data
    data.__getitem__.return_value = vertical_data

    # assign_coords e sortby devem retornar um novo dataset (encadeável)
    data.assign_coords.return_value = data
    data.sortby.return_value = data

    # slice_domain retorna dataset modificado (mock)
    mock_slice_domain.return_value = data

    # Chama a função
    result = preprocess_data(data, df_namelist, args, app_logger)

    # Testa se ValueError não foi levantado (sucesso)
    assert result == data

    # Verifica chamadas do logger
    app_logger.info.assert_any_call('Preprocessing data...')
    app_logger.debug.assert_any_call('Force vertical levels to be in Pa...')
    app_logger.debug.assert_any_call('Sorting data by longitude, vertical level and latitude...')
    app_logger.debug.assert_any_call('Slicing data...')
    app_logger.debug.assert_any_call('Assigning lat and lon as radians...')
    app_logger.info.assert_any_call('Done.')

    # Verifica se vertical_data.metpy.convert_units foi chamado com 'Pa'
    vertical_data.metpy.convert_units.assert_called_once_with('Pa')

    # Verifica se assign_coords foi chamado para pressão e depois para lat/lon
    assert data.assign_coords.call_count >= 4

    # Verifica se sortby foi chamado para lon, plev e lat
    calls = [call.args[0] for call in data.sortby.call_args_list]
    assert 'lon' in calls
    assert 'plev' in calls
    assert 'lat' in calls

    # Verifica se slice_domain foi chamado com args corretos
    mock_slice_domain.assert_called_once_with(data, args, df_namelist)


def test_preprocess_data_missing_variable():
    # Namelist faltando uma variável crítica
    df_namelist = pd.DataFrame({'Variable': [None, 'lat', 'plev']},
                               index=['Longitude', 'Latitude', 'Vertical Level'])
    args = MagicMock()
    app_logger = MagicMock()
    data = MagicMock()

    with pytest.raises(ValueError, match='Missing critical namelist variables.'):
        preprocess_data(data, df_namelist, args, app_logger)
