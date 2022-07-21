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


from metpy.units import units
from metpy.constants import Rd
from metpy.constants import Cp_d
from metpy.constants import g
from metpy.constants import Re
from metpy.calc import potential_temperature

import pandas as pd
import xarray as xr
import os
import numpy as np
import argparse
import traceback
import logging

from calc import CalcAreaAverage, CalcZonalAverage

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

class DataObject:
    def __init__(self,NetCDF_data: xr.Dataset,
                 dfVars: pd.DataFrame,
                 dfbox: pd.DataFrame=None):
        self.LonIndexer = dfVars.loc['Longitude']['Variable']
        self.LatIndexer = dfVars.loc['Latitude']['Variable']
        self.TimeIndexer = dfVars.loc['Time']['Variable']
        self.LevelIndexer = dfVars.loc['Vertical Level']['Variable']
        if dfbox is None:
            self.NetCDF_data = NetCDF_data
        else:
            self.WesternLimit = float((NetCDF_data[self.LonIndexer]
                                 [(np.abs(NetCDF_data[self.LonIndexer] - 
                                  float(dfbox.loc['min_lon']))).argmin()]))
            self.EasternLimit = float((NetCDF_data[self.LonIndexer]
                                  [(np.abs(NetCDF_data[self.LonIndexer] - 
                                   float(dfbox.loc['max_lon']))).argmin()]).values)
            self.SouthernLimit = float((NetCDF_data[self.LatIndexer]
                                   [(np.abs(NetCDF_data[self.LatIndexer] - 
                                   float(dfbox.loc['min_lat']))).argmin()]).values)
            self.NorthernLimit = float((NetCDF_data[self.LatIndexer]
                                   [(np.abs(NetCDF_data[self.LatIndexer] - 
                                   float(dfbox.loc['max_lat']))).argmin()]).values)
            self.NetCDF_data = NetCDF_data.sel(
                **{self.LatIndexer:slice(self.NorthernLimit,self.SouthernLimit),
                   self.LonIndexer: slice(self.WesternLimit,self.EasternLimit)})
        self.Temperature = self.NetCDF_data[dfVars.loc['Air Temperature']['Variable']] \
            * units(dfVars.loc['Air Temperature']['Units']).to('K')
        self.u = self.NetCDF_data[dfVars.loc['Eastward Wind Component']['Variable']] \
            * units(dfVars.loc['Eastward Wind Component']['Units']).to('m/s')
        self.v = self.NetCDF_data[dfVars.loc['Northward Wind Component']['Variable']] \
            * units(dfVars.loc['Northward Wind Component']['Units']).to('m/s')
        self.Omega = self.NetCDF_data[dfVars.loc['Omega Velocity']['Variable']] \
            * units(dfVars.loc['Omega Velocity']['Units']).to('Pa/s')
        self.GeopotHeight = self.NetCDF_data[dfVars.loc['Geopotential Height']['Variable']] \
            * units(dfVars.loc['Geopotential Height']['Units']).to('gpm')
        self.Pressure = self.NetCDF_data[self.LevelIndexer]
        self.Theta = theta = potential_temperature(
            self.Pressure,self.Temperature)
        self.dTdt = self.Temperature.differentiate(
                self.TimeIndexer,datetime_unit='s') / units('seconds')
        self.AdvHTemp = self.HorizontalTemperatureAdvection()
        self.sigma = (self.Temperature/theta) *theta.differentiate(
            self.LevelIndexer)/units.hPa
        self.ResT =  self.dTdt - self.AdvHTemp - (self.sigma * self.Omega)
        self.AdiabaticHeating = self.ResT*Cp_d
    def HorizontalTemperatureAdvection(self):
        lons,lats = self.Temperature[self.LonIndexer],\
            self.Temperature[self.LatIndexer]
        cos_lats = np.cos(np.deg2rad(lats))
        ## Horizontal advection of temperature ##
        # Differentiate temperature in respect to longitude and latitude
        dTdlambda = self.Temperature.copy(deep=True
                                          ).differentiate(self.LonIndexer)
        dTdphi = self.Temperature.copy(deep=True
                                       ).differentiate(self.LatIndexer)
        # Get the values for dx and dy in meters
        dx = np.deg2rad(lons.differentiate(self.LonIndexer))*cos_lats*Re
        dy = np.deg2rad(lats.differentiate(self.LatIndexer))*Re
        AdvHT = -1* ((self.u*dTdlambda/dx)+(self.v*dTdphi/dy)) 
        return AdvHT

def LagrangianAnalysis(LagrangianObj):
    
    # Track file
    trackfile = './track'
    track = pd.read_csv(trackfile,parse_dates=[0],delimiter=';',index_col='time')

    # Directory where results will be stored
    ResultsMainDirectory = '../CycloneThermodynamics_Results'
    # Append data limits to outfile name
    outfile_name = ''.join(infile.split('/')[-1].split('.nc'))+'_lagranigan'
    # Each dataset of results have its own directory, allowing to store results
    # from more than one experiment at each time
    ResultsSubDirectory = ResultsMainDirectory+'/'+outfile_name+'/'
    # Check if the LEC_Figures directory exists. If not, creates it
    check_create_folder(ResultsMainDirectory)
    # Check if a directory for current data exists. If not, creates it
    check_create_folder(ResultsSubDirectory)
    
    timesteps = LagrangianObj.NetCDF_data[LagrangianObj.TimeIndexer]
    
    PresLevels = LagrangianObj.Temperature[LagrangianObj.LevelIndexer].\
        metpy.convert_units('hPa').values
    
    stored_terms = ['T_AA','Omega_AA','Q_AA', # average terms
                    'T_ZE_AA','Omega_ZE_AA','Q_ZE_AA', # eddy terms (averaged)
                    'AdvHTemp','SpOmega','dTdt','ResT'] # thermodyn. eq. terms
    
    dfDict = {}
    for term in stored_terms:
        dfDict[term] = pd.DataFrame(columns=[str(t.values) for t in timesteps],
        index=[float(i) for i in PresLevels])

    
    for t in timesteps:
        # Get current time and box limits
        itime = str(t.values)
        datestr1 = pd.to_datetime(itime).strftime('%Y-%m-%d %HZ')
        datestr2 = pd.to_datetime(itime).strftime('%Y%m%d%H%M')
        min_lon, max_lon = track.loc[itime]['Lon']-7.5,track.loc[itime]['Lon']+7.5
        min_lat, max_lat = track.loc[itime]['Lat']-7.5,track.loc[itime]['Lat']+7.5
        print('\nComputing terms for '+datestr1+'...')
        print('Box limits (lon/lat): '+str(max_lon)+'/'+str(max_lat),
              ' '+str(min_lon)+'/'+str(min_lat))
        
        
        # Get closest grid point to actual track
        WesternLimit = float((NetCDF_data[LagrangianObj.LonIndexer]
                             [(np.abs(NetCDF_data[LagrangianObj.LonIndexer] - 
                              min_lon)).argmin()]))
        EasternLimit = float((NetCDF_data[LagrangianObj.LonIndexer]
                              [(np.abs(NetCDF_data[LagrangianObj.LonIndexer] - 
                              max_lon)).argmin()]).values)
        SouthernLimit = float((NetCDF_data[LagrangianObj.LatIndexer]
                               [(np.abs(NetCDF_data[LagrangianObj.LatIndexer] - 
                               min_lat)).argmin()]).values)
        NorthernLimit = float((NetCDF_data[LagrangianObj.LatIndexer]
                               [(np.abs(NetCDF_data[LagrangianObj.LatIndexer] - 
                               max_lat)).argmin()]).values)
        
        # Store area average of eddy component of temperature
        T = LagrangianObj.Temperature.sel({LagrangianObj.TimeIndexer:t}).sel(
            **{LagrangianObj.LatIndexer:slice(NorthernLimit,SouthernLimit),
               LagrangianObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        T_ZA = CalcZonalAverage(T,LagrangianObj.LonIndexer)
        T_AA = CalcAreaAverage(T_ZA,LagrangianObj.LatIndexer)
        T_ZE = T - T_ZA
        T_ZE_AA = CalcAreaAverage(T_ZE,LagrangianObj.LatIndexer,
                                  LonIndexer=LagrangianObj.LonIndexer)
        dfDict['T_AA'][itime] = T_AA  
        dfDict['T_ZE_AA'][itime] = T_ZE_AA   
        
        # Store area average of eddy component of omega
        Omega = LagrangianObj.Omega.sel({LagrangianObj.TimeIndexer:t}).sel(
            **{LagrangianObj.LatIndexer:slice(NorthernLimit,SouthernLimit),
               LagrangianObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        Omega_ZA = CalcZonalAverage(Omega,LagrangianObj.LonIndexer)
        Omega_AA = CalcAreaAverage(Omega_ZA,LagrangianObj.LatIndexer)
        Omega_ZE = Omega - Omega_ZA
        Omega_ZE_AA = CalcAreaAverage(Omega_ZE,LagrangianObj.LatIndexer,
                                  LonIndexer=LagrangianObj.LonIndexer)
        dfDict['Omega_AA'][itime] = Omega_AA
        dfDict['Omega_ZE_AA'][itime] = Omega_ZE_AA        
        
        # Store area average of eddy component of adiabatic heating
        Q = LagrangianObj.AdiabaticHeating.sel({LagrangianObj.TimeIndexer:t}).sel(
            **{LagrangianObj.LatIndexer:slice(NorthernLimit,SouthernLimit),
               LagrangianObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        Q_ZA = CalcZonalAverage(Q,LagrangianObj.LonIndexer)
        Q_AA = CalcAreaAverage(Q_ZA,LagrangianObj.LatIndexer)
        Q_ZE = Q - Q_ZA
        Q_ZE_AA = CalcAreaAverage(Q_ZE,LagrangianObj.LatIndexer,
                                  LonIndexer=LagrangianObj.LonIndexer)
        dfDict['Q_AA'][itime] = Q_AA
        dfDict['Q_ZE_AA'][itime] = Q_ZE_AA
        
        # Store area average for each term of thermodynamic equation
        AdvHTemp = LagrangianObj.AdvHTemp.sel({LagrangianObj.TimeIndexer:t}).sel(
            **{LagrangianObj.LatIndexer:slice(NorthernLimit,SouthernLimit),
               LagrangianObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        sigma = LagrangianObj.sigma.sel({LagrangianObj.TimeIndexer:t}).sel(
            **{LagrangianObj.LatIndexer:slice(NorthernLimit,SouthernLimit),
               LagrangianObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        dTdt = LagrangianObj.dTdt.sel({LagrangianObj.TimeIndexer:t}).sel(
            **{LagrangianObj.LatIndexer:slice(NorthernLimit,SouthernLimit),
               LagrangianObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        ResT = LagrangianObj.ResT.sel({LagrangianObj.TimeIndexer:t}).sel(
            **{LagrangianObj.LatIndexer:slice(NorthernLimit,SouthernLimit),
               LagrangianObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        
        dfDict['AdvHTemp'][itime] = CalcAreaAverage(AdvHTemp,
                                     LagrangianObj.LatIndexer,
                                     LonIndexer=LagrangianObj.LonIndexer).\
                                        metpy.convert_units('K/ day')
        dfDict['SpOmega'][itime] = CalcAreaAverage(sigma*Omega,
                                     LagrangianObj.LatIndexer,
                                     LonIndexer=LagrangianObj.LonIndexer).\
                                        metpy.convert_units('K/ day')
        dfDict['dTdt'][itime] = CalcAreaAverage(dTdt,
                                     LagrangianObj.LatIndexer,
                                     LonIndexer=LagrangianObj.LonIndexer).\
                                        metpy.convert_units('K/ day')
        dfDict['ResT'][itime] = CalcAreaAverage(ResT,
                                     LagrangianObj.LatIndexer,
                                     LonIndexer=LagrangianObj.LonIndexer).\
                                        metpy.convert_units('K/ day')
        
    # Save CSV files with all the terms stored above
    for term in stored_terms:
        dfDict[term].to_csv(ResultsSubDirectory+term+'.csv')
    # Make figures
    os.system("python plot_timeseries.py")
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "\
Program for analysing the thermodynamics of a cyclone. \n \
The program can use two distinct frameworks:\
    1) Lagragian framework. A box is definid in the box_lims' file and then the \
        energetics are computed for a fixed domain.\
    2) Eulerian framework. The domain is not fixed and follows the system using \
        the 'track' file.\
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
    
    varlist = './fvars'
    dfVars = pd.read_csv(varlist,sep= ';',index_col=0,header=0)
    
    infile  = args.infile
    NetCDF_data = convert_lon(xr.open_dataset(infile),
                              dfVars.loc['Longitude']['Variable'])
    
    # Run the program
    if args.eulerian:
        start_time = time.time()
        dfbox = pd.read_csv('./box_limits',header=None,delimiter=';',index_col=0)
        EulerianObj = DataObject(NetCDF_data,dfVars,dfbox)
        print('NOT IMPLEMENTED YET')
        print("--- %s seconds running eulerian framework ---" % (time.time() - start_time))
    if args.lagrangian:
        start_time = time.time()
        LagrangianObj = DataObject(NetCDF_data,dfVars)
        LagrangianAnalysis(LagrangianObj)
        print("--- %s seconds for running lagrangian framework ---" % (time.time() - start_time))            