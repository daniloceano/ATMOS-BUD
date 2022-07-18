#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 10:09:03 2022

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

from EnergyContents import EnergyContents
from ConversionTerms import ConversionTerms
from BoundaryTerms import BoundaryTerms
from GenerationDissipationTerms import GenerationDissipationTerms
from BoxData import BoxData
from thermodynamics import AdiabaticHEating
from compute_terms import calc_budget_diff,calc_residuals
from metpy.units import units
import pandas as pd
import xarray as xr
import os
import numpy as np
import argparse
from metpy.constants import g

import traceback
import logging

import time

def check_create_folder(DirName):
    if not os.path.exists(DirName):
                os.makedirs(DirName)
                print(DirName+' created')
    else:
        print(DirName+' directory exists')
 
# Convert longitudes from 0:360 range to -180:180
def convert_lon(df,LonIndexer):
    df.coords[LonIndexer] = (df.coords[LonIndexer] + 180) % 360 - 180
    df = df.sortby(df[LonIndexer])
    return df

# Function for opening the data
def get_data(infile, varlist):   
    print('Variables specified by the user in: '+varlist)
    print('Attempting to read '+varlist+' file...')
    try:
        dfVars = pd.read_csv(varlist,sep= ';',index_col=0,header=0)
    except:
        raise SystemExit("Error: verify that there is a 'fvar' text file\
 located on your working directory")
    # Print variables for the user
    print('List of variables found:')
    print(dfVars)
    # Get data indexers
    LonIndexer,LatIndexer,TimeIndexer,LevelIndexer = \
      dfVars.loc['Longitude']['Variable'],dfVars.loc['Latitude']['Variable'],\
      dfVars.loc['Time']['Variable'],dfVars.loc['Vertical Level']['Variable']
    print('Ok!')
    print('Opening inpyt data...')
    try:
        full_data = convert_lon(xr.open_dataset(infile),LonIndexer)
    except:
        raise SystemExit('ERROR!!!!!\n Could not open data. Check if path is\
 correct, fvars file and file format (should be .nc)')
    print('Ok!')
    # Sort data coordinates - data from distinc sources might have different
    # arrangements, which could affect the results from the integrations
    full_data = full_data.sortby(LonIndexer).sortby(LevelIndexer,
                            ascending=False).sortby(LatIndexer,ascending=False)
    # Fill missing values with 0
    full_data = full_data.fillna(0)
    try:
        full_data = full_data.where(full_data.apply(np.isfinite)).fillna(0.0)
    except:
        full_data = full_data.fillna(0)
    # load data into memory (code optmization)
    data = full_data.load()
    # Stores data as separated variables and give them correct units
    tair = data[dfVars.loc['Air Temperature']['Variable']] \
        * units(dfVars.loc['Air Temperature']['Units']).to('K')
    omega = data[dfVars.loc['Omega Velocity']['Variable']]*\
        units(dfVars.loc['Omega Velocity']['Units']).to('Pa/s')
    u = data[dfVars.loc['Eastward Wind Component']['Variable']]*\
        units(dfVars.loc['Eastward Wind Component']['Units'])
    v = data[dfVars.loc['Northward Wind Component']['Variable']]*\
        units(dfVars.loc['Northward Wind Component']['Units'])
    if args.geopotential:
         hgt = (data[dfVars.loc['Geopotential']['Variable']] \
        * units(dfVars.loc['Geopotential']['Units'])/g).metpy.convert_units('gpm')
    else:
        hgt = data[dfVars.loc['Geopotential Height']['Variable']]\
            *units(dfVars.loc['Geopotential Height']['Units']).to('gpm')
    return LonIndexer, LatIndexer, TimeIndexer, LevelIndexer, tair, hgt,\
            omega, u, v
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "\
Program for analysing the thermodynamics of a cyclone. \n \
The program can use two distinct frameworks:\
    1) Lagragian framework. A box is definid in the box_lims' file and then the \
       energetics are computed for a fixed domain.\
    2) Eulerian framework. The domain is not fixed and follows the system using \
       the track file.\
 Both frameworks can be applied at the same time, given the required files are\
 provided. An auxilliary 'fvars' file is also needed for both frameworks. \
 It contains the specified names used for each variable.  The results \
 are stored as csv files in the 'LEC_results' directory on ../ and it also \
 creates figures for visualising the results. \
 The use of -r flag is required while the computation of friction parameters \
 is not implemented")
    parser.add_argument("infile", help = "input .nc file with temperature,\
 geopotential and meridional, zonal and vertical components of the wind,\
 in pressure levels")
    parser.add_argument("-g", "--geopotential", default = False,
    action='store_true', help = "use the geopotential data instead of\
 geopotential height. The file fvars must be adjusted for doing so.")
    parser.add_argument("-e", "--eulerian", default = False,
    action='store_true', help = "compute the energetics for a fixed domain\
 specified by the box_lims file.")
    parser.add_argument("-l", "--lagrangian", default = False,
    action='store_true', help = "compute the energetics for a fixed domain\
 specified by the box_lims file.")
    args = parser.parse_args()
    infile  = args.infile
    # box_limits = args.box_limits
    varlist = './fvars'
    # Run the program
    start_time = time.time()
    # if args.eulerian:
    #     eulerian()
    # print("--- %s seconds running eulerian framework ---" % (time.time() - start_time))
    # start_time = time.time()
    # if args.lagrangian:
    #     lagrangian()
    # print("--- %s seconds for running lagrangian framework ---" % (time.time() - start_time))            