import pytest
import pandas as pd
from unittest.mock import MagicMock, patch
from collections import namedtuple

from src.select_domain import get_domain_limits

Args = namedtuple('Args', ['track', 'choose', 'fixed'])

def create_variables_mock():
    # mocks para os 7 elementos passados em variables_at_plevel
    iu_plevel = MagicMock()
    iv_plevel = MagicMock()
    zeta = MagicMock()
    ight_plevel = MagicMock()
    lat = MagicMock()
    lon = MagicMock()
    itime = '2023-01-01 00:00'  # pode ser string qualquer para índice
    return (iu_plevel, iv_plevel, zeta, ight_plevel, lat, lon, itime)

def test_get_domain_limits_track_with_width_length():
    args = Args(track=True, choose=False, fixed=False)
    
    variables = create_variables_mock()

    # Criar DataFrame track com width e length
    track = pd.DataFrame({
        'width': [20],
        'length': [30],
        'Lon': [100],
        'Lat': [10]
    }, index=['2023-01-01 00:00'])

    result = get_domain_limits(args, *variables, track=track)
    expected = {
        'min_lon': 100 - 10,
        'max_lon': 100 + 10,
        'min_lat': 10 - 15,
        'max_lat': 10 + 15,
        'central_lat': 10,
        'central_lon': 100
    }
    assert result == expected

def test_get_domain_limits_track_without_width_length():
    args = Args(track=True, choose=False, fixed=False)
    variables = create_variables_mock()

    # DataFrame sem colunas width e length
    track = pd.DataFrame({
        'Lon': [100],
        'Lat': [10]
    }, index=['2023-01-01 00:00'])

    result = get_domain_limits(args, *variables, track=track)
    expected = {
        'min_lon': 100 - 7.5,   # width/2 = 15/2 = 7.5
        'max_lon': 100 + 7.5,
        'min_lat': 10 - 7.5,
        'max_lat': 10 + 7.5,
        'central_lat': 10,
        'central_lon': 100
    }
    assert result == expected

@patch('src.select_domain.draw_box_map')
def test_get_domain_limits_choose(mock_draw_box_map):
    args = Args(track=False, choose=True, fixed=False)
    variables = create_variables_mock()

    # Mock retorno de draw_box_map
    mock_draw_box_map.return_value = {
        'min_lon': 90,
        'max_lon': 110,
        'min_lat': 5,
        'max_lat': 15
    }

    result = get_domain_limits(args, *variables)
    expected = {
        'min_lon': 90,
        'max_lon': 110,
        'min_lat': 5,
        'max_lat': 15,
        'central_lat': 10,
        'central_lon': 100
    }
    assert result == expected
    mock_draw_box_map.assert_called_once()

@patch('pandas.read_csv')
def test_get_domain_limits_fixed(mock_read_csv):
    args = Args(track=False, choose=False, fixed=True)
    variables = create_variables_mock()

    lat = variables[4]
    lon = variables[5]

    # Mock DataFrame retornado pelo read_csv
    mock_read_csv.return_value = pd.DataFrame({
        1: ['-50.0', '-30.0', '10.0', '20.0']
    }, index=['min_lon', 'max_lon', 'min_lat', 'max_lat'])

    # Mock para lat.sel e lon.sel
    lat.sel.return_value = -40.0
    lon.sel.return_value = -40.0

    result = get_domain_limits(args, *variables)
    expected = {
        'min_lon': -40.0,
        'max_lon': -40.0,
        'min_lat': -40.0,
        'max_lat': -40.0,
        'central_lat': (-40.0 + -40.0) / 2,
        'central_lon': (-40.0 + -40.0) / 2
    }
    assert result == expected

    mock_read_csv.assert_called_once_with('inputs/box_limits', header=None, delimiter=';', index_col=0)
    lat.sel.assert_called()
    lon.sel.assert_called()
