import pytest
from unittest.mock import patch, MagicMock
from collections import namedtuple
import xarray as xr

from src.data_handling import load_data

Args = namedtuple('Args', ['gfs'])

@patch('src.data_handling.convert_lon')
@patch('xarray.open_dataset')
def test_load_data_normal(mock_open_dataset, mock_convert_lon):
    # Configurações
    args = Args(gfs=False)
    infile = 'file.nc'
    longitude_indexer = "longitude"
    mock_logger = MagicMock()

    # Mock do dataset retornado
    mock_ds = MagicMock(spec=xr.Dataset)
    mock_open_dataset.return_value = mock_ds
    mock_convert_lon.return_value = 'converted_data'

    # Chama a função
    result = load_data(infile, longitude_indexer, args, mock_logger)

    # Testa se open_dataset foi chamado corretamente
    mock_open_dataset.assert_called_once_with(infile)

    # Testa se convert_lon foi chamado com o dataset e longitude_indexer
    mock_convert_lon.assert_called_once_with(mock_ds, longitude_indexer)

    # Testa retorno
    assert result == 'converted_data'

    # Testa se logger chamou info duas vezes
    assert mock_logger.info.call_count == 2
    mock_logger.info.assert_any_call(f'Loading {infile}...')
    mock_logger.info.assert_any_call(f'Loaded {infile} successfully!')

@patch('src.data_handling.convert_lon')
@patch('xarray.open_mfdataset')
def test_load_data_gfs(mock_open_mfdataset, mock_convert_lon):
    args = Args(gfs=True)
    infile = 'files_pattern.nc'
    longitude_indexer = slice(0, 10)
    mock_logger = MagicMock()

    mock_ds = MagicMock(spec=xr.Dataset)
    mock_open_mfdataset.return_value = mock_ds
    mock_convert_lon.return_value = 'converted_data_gfs'

    result = load_data(infile, longitude_indexer, args, mock_logger)

    mock_open_mfdataset.assert_called_once_with(
        infile,
        engine='cfgrib',
        parallel=True,
        filter_by_keys={'typeOfLevel': 'isobaricInhPa'},
        combine='nested',
        concat_dim='time'
    )

    mock_convert_lon.assert_called_once_with(mock_ds, longitude_indexer)

    assert result == 'converted_data_gfs'

    assert mock_logger.info.call_count == 2
    mock_logger.info.assert_any_call(f'Loading {infile}...')
    mock_logger.info.assert_any_call(f'Loaded {infile} successfully!')

@patch('xarray.open_dataset')
def test_load_data_file_not_found(mock_open_dataset):
    args = Args(gfs=False)
    infile = 'missing_file.nc'
    longitude_indexer = slice(0, 10)
    mock_logger = MagicMock()

    # Simula FileNotFoundError ao abrir dataset
    mock_open_dataset.side_effect = FileNotFoundError

    with pytest.raises(FileNotFoundError):
        load_data(infile, longitude_indexer, args, mock_logger)

    # Verifica se logger registrou erro
    mock_logger.error.assert_called_once_with(f'File not found: {infile}')

@patch('xarray.open_dataset')
def test_load_data_other_exception(mock_open_dataset):
    args = Args(gfs=False)
    infile = 'file.nc'
    longitude_indexer = slice(0, 10)
    mock_logger = MagicMock()

    # Simula outra exceção
    mock_open_dataset.side_effect = RuntimeError("fail")

    with pytest.raises(RuntimeError):
        load_data(infile, longitude_indexer, args, mock_logger)

    mock_logger.error.assert_called_once()
    assert 'Failed to load data from' in mock_logger.error.call_args[0][0]
