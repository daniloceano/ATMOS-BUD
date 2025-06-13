import pytest
import pandas as pd
from collections import namedtuple
from unittest.mock import MagicMock

from src.utils import get_domain_extreme_values

Args = namedtuple('Args', ['track', 'track_vorticity', 'track_geopotential', 'level'])

class SliceMock:
    def __init__(self, min_val, max_val):
        self._min = min_val
        self._max = max_val
    
    def min(self):
        return self._min
    
    def max(self):
        return self._max

def create_slices(min_zeta=1, max_zeta=5, min_hgt=10, max_hgt=20, max_wind=30):
    izeta = SliceMock(min_zeta, max_zeta)
    ight = SliceMock(min_hgt, max_hgt)
    iwspd = SliceMock(min_val=0, max_val=max_wind)
    return izeta, ight, iwspd

def test_track_true_with_columns():
    slices = create_slices()
    itime = '2023-01-01 00:00'
    args = Args(track=True, track_vorticity='min', track_geopotential='min', level=500)
    track = pd.DataFrame({
        'min_zeta_500': [0.1],
        'min_hgt_500': [9.9],
        'max_wind_500': [29.9],
    }, index=[itime])

    zeta, hgt, wind = get_domain_extreme_values(itime, args, slices, track)
    assert zeta == 0.1
    assert hgt == 9.9
    assert wind == 29.9

def test_track_true_without_columns():
    slices = create_slices(min_zeta=1, max_zeta=5, min_hgt=10, max_hgt=20, max_wind=30)
    itime = '2023-01-01 00:00'
    args = Args(track=True, track_vorticity='max', track_geopotential='max', level=500)
    # Track sem as colunas exigidas
    track = pd.DataFrame(index=[itime])

    zeta, hgt, wind = get_domain_extreme_values(itime, args, slices, track)
    # track_vorticity = max, logo zeta deve ser max do slice (5)
    assert zeta == 5
    # track_geopotential = max, logo hgt deve ser max do slice (20)
    assert hgt == 20
    # wind sempre max do slice (30)
    assert wind == 30

def test_track_false_min_vorticity():
    slices = create_slices(min_zeta=1, max_zeta=5, min_hgt=10, max_hgt=20, max_wind=30)
    itime = '2023-01-01 00:00'
    args = Args(track=False, track_vorticity='min', track_geopotential='max', level=500)

    zeta, hgt, wind = get_domain_extreme_values(itime, args, slices)
    assert zeta == 1  # min zeta slice
    assert hgt == 10  # min hgt slice (sempre min)
    assert wind == 30 # max wind slice

def test_track_false_max_vorticity():
    slices = create_slices(min_zeta=1, max_zeta=5, min_hgt=10, max_hgt=20, max_wind=30)
    itime = '2023-01-01 00:00'
    args = Args(track=False, track_vorticity='max', track_geopotential='min', level=500)

    zeta, hgt, wind = get_domain_extreme_values(itime, args, slices)
    assert zeta == 5  # max zeta slice
    assert hgt == 10  # min hgt slice (sempre min)
    assert wind == 30 # max wind slice
