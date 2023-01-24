#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for calculations necessary for computate Lorenz Energy Cycle such 
integrations, area avreages, partial differentiation and calculate the
static stability parameter

Created by Danilo Couto de Souza
Universidade de São Paulo (USP)
Instituto de Astornomia, Ciências Atmosféricas e Geociências
São Paulo - Brazil

danilo.oceano@gmail.com

"""

import numpy as np

from metpy import units

def CalcZonalAverage(VariableData):
    """
    Computates variable zonal average of some variable, for all z
    levels and time steps.
    
    Source:
        Brennan, F. E., & Vincent, D. G. (1980).
        Zonal and Eddy Components of the Synoptic-Scale Energy Budget
        during Intensification of Hurricane Carmen (1974),
        Monthly Weather Review, 108(7), 954-965. Retrieved Jan 25, 2022, from:
        https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_0954_zaecot_2_0_co_2.xml
    
    Parameters
    ----------
    VariableData: xarray.Dataset
        arrays containing data to be integrated. Requires dimension rlons 
        (longitude in radians)
    Returns
    -------
    zonal_ave: xarray.Dataset
        Arrays of zonal avreages for all longitudes from the passed Dataset
    """
    xlength = VariableData['rlons'][-1]- VariableData['rlons'][0]
    return VariableData.integrate("rlons")/xlength

def CalcAreaAverage(VariableData, ZonalAverage=False):
    """
    Computates the Area Average of a function.
    
    The default is to computate the zonal average and then a meridional average.
    If the input data is already some sort of zonal quantity (average or not),
    simply set LonIndexer to None
    
    Source:
        Brennan, F. E., & Vincent, D. G. (1980).
        Zonal and Eddy Components of the Synoptic-Scale Energy Budget
        during Intensification of Hurricane Carmen (1974),
        Monthly Weather Review, 108(7), 954-965. Retrieved Jan 25, 2022, from:
        https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_0954_zaecot_2_0_co_2.xml
    
    Parameters
    ----------
    VariableData: xarray.Dataset
        arrays containing data to be integrated
    ZonalAverage: string (optional)
        if passed, it will first compute zonal averages
   
    Returns
    -------
    zonal_ave: xarray.Dataset
        Arrays of area avreages for all latitudes and longitudes from
        the passed Dataset
    """
    # Compute zonal average if requested
    if ZonalAverage:
        ZA = CalcZonalAverage(VariableData)
    else:
        ZA = VariableData
    # Take the area avearge
    ylength = np.sin(
        VariableData['rlats'][-1]) - np.sin(VariableData['rlats'][0])
    return (ZA*ZA["coslats"]).integrate("rlats")/ylength
