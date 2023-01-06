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
from metpy.constants import Cp_d
from metpy.constants import Re
from metpy.calc import potential_temperature
from metpy.calc import geopotential_to_height
from metpy.calc import height_to_pressure_std
from metpy.constants import g

import pandas as pd
import xarray as xr
import os
import numpy as np
import argparse

from calc import CalcAreaAverage, CalcZonalAverage
from plot_maps import LagrangianMaps as plot_map
from select_area import draw_box_map

import time

def check_create_folder(DirName):
    """

    Check if directory exists and if not, creates it.
    
    Parameters
    ----------
    DirName : str
        directory name.

    Returns
    -------
    None. 

    """
    if not os.path.exists(DirName):
                os.makedirs(DirName)
                print(DirName+' created')
    else:
        print(DirName+' directory exists')
 
def convert_lon(xr,LonIndexer):
    """
    
    Convert longitudes from 0:360 range to -180:180

    Parameters
    ----------
    xr : xarray.DataArray 
        gridded data.
    LonIndexer : str
        corrdinate indexer used for longitude.

    Returns
    -------
    xr : xarray.DataArray 
        gridded data with longitude converted to desired format.

    """
    xr.coords[LonIndexer] = (xr.coords[LonIndexer] + 180) % 360 - 180
    xr = xr.sortby(xr[LonIndexer])
    return xr

class DataObject:
    """
    Object for storing variables from a NetCDF file on intialization.
    It also computes each term of the Quasi-Geostrophic Equation, except the
    Adiabatic Heating Term (Q) which is estimated as a residual. 
    Note that: Q = J * Cp_d
    """
    def __init__(self,NetCDF_data: xr.Dataset,
                 dfVars: pd.DataFrame,
                 dfbox: pd.DataFrame=None):
        self.LonIndexer = dfVars.loc['Longitude']['Variable']
        self.LatIndexer = dfVars.loc['Latitude']['Variable']
        self.TimeIndexer = dfVars.loc['Time']['Variable']
        self.LevelIndexer = dfVars.loc['Vertical Level']['Variable']
        # When constructing object for Stationary analysis, the data can
        # be sliced beforehand using the dfBox limits, but for the unstationary
        # analysis (dfBox not specified), we need full data and then it is
        # sliced for each timestep.
        if dfbox is None:
            self.NetCDF_data = NetCDF_data
        else:
            self.WesternLimit = float(NetCDF_data[self.LonIndexer].sel(
                {self.LonIndexer:float(dfbox.loc['min_lon'].values)},
                method='nearest'))
            self.EasternLimit =float(NetCDF_data[self.LonIndexer].sel(
                {self.LonIndexer:float(dfbox.loc['max_lon'].values)},
                method='nearest'))
            self.SouthernLimit =float(NetCDF_data[self.LatIndexer].sel(
                {self.LatIndexer:float(dfbox.loc['min_lat'].values)},
                method='nearest'))
            self.NorthernLimit =float(NetCDF_data[self.LatIndexer].sel(
                {self.LatIndexer:float(dfbox.loc['max_lat'].values)},
                method='nearest'))
            self.NetCDF_data = NetCDF_data.sel(
                **{self.LatIndexer:slice(self.SouthernLimit,self.NorthernLimit),
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
        lons,lats  = self.NetCDF_data[self.LonIndexer],\
            self.NetCDF_data[self.LatIndexer]
        cos_lats = self.NetCDF_data["coslats"]
        # Differentiate temperature in respect to longitude and latitude
        dTdlambda = self.Temperature.differentiate(self.LonIndexer)
        dTdphi = self.Temperature.differentiate(self.LatIndexer)
        # Get the values for dx and dy in meters
        dx = np.deg2rad(lons.differentiate(self.LonIndexer))*cos_lats*Re
        dy = np.deg2rad(lats.differentiate(self.LatIndexer))*Re
        AdvHT = -1* ((self.u*dTdlambda/dx)+(self.v*dTdphi/dy)) 
        return AdvHT

def UnstationaryAnalysis(UnstationaryObj):
    """
    Parameters
    ----------
    UnstationaryObj : DataObj
        Object containing meteorological data from a NetCDF file

    Returns
    -------
    For each timestep, creates a box with 15º width/height and store resuts 
    in CSV files. Also creates Figures for posterior analysis. 
    The variables stored are the Temperature, Omega and Q eddy means over the
    domain (the box created), their eddy components (value at each grid 
    point minus its zonal average) averaged over the domain and each term of
    the Quasi-Geostrophic Thermodynamic Equation.
    """
    # Track file
    if args.track:
        trackfile = './inputs/track'
        track = pd.read_csv(trackfile,parse_dates=[0],delimiter=';',index_col='time')
    
    timesteps = UnstationaryObj.NetCDF_data[UnstationaryObj.TimeIndexer]
    
    PresLevels = UnstationaryObj.Temperature[UnstationaryObj.LevelIndexer].\
        metpy.convert_units('hPa').values
    
    stored_terms = ['T_AA','Omega_AA','Q_AA', # average terms
                    'T_ZE_AA','Omega_ZE_AA','Q_ZE_AA', # eddy terms (averaged)
                    'AdvHTemp','SpOmega','dTdt','ResT'] # thermodyn. eq. terms
    
    dfDict = {}
    for term in stored_terms:
        dfDict[term] = pd.DataFrame(columns=[str(t.values) for t in timesteps],
        index=[float(i) for i in PresLevels])

    
    for t in timesteps:
        # Get current time and time strings
        itime = str(t.values)
        datestr = pd.to_datetime(itime).strftime('%Y-%m-%d %HZ')
        
        # Get current time and box limits
        if args.track:
            min_lon, max_lon = track.loc[itime]['Lon']-7.5,track.loc[itime]['Lon']+7.5
            min_lat, max_lat = track.loc[itime]['Lat']-7.5,track.loc[itime]['Lat']+7.5
        elif args.choose:
            iu_1000 = UnstationaryObj.u.sel({UnstationaryObj.TimeIndexer:t,
                                           UnstationaryObj.LevelIndexer:1000})
            iv_1000 = UnstationaryObj.u.sel({UnstationaryObj.TimeIndexer:t,
                                           UnstationaryObj.LevelIndexer:1000})
            igeop_1000 = UnstationaryObj.GeopotHeight.sel(
                {UnstationaryObj.TimeIndexer:t,
                 UnstationaryObj.LevelIndexer:1000})*g
            iheight_1000 = geopotential_to_height(igeop_1000)
            imslp = height_to_pressure_std(iheight_1000)
            draw_box_map(iu_1000, iv_1000, imslp)
            
        print('\nComputing terms for '+datestr+'...')
        print('Box limits (lon/lat): '+str(max_lon)+'/'+str(max_lat),
              ' '+str(min_lon)+'/'+str(min_lat))
        
        
        # Get closest grid point to actual track
        WesternLimit = float(NetCDF_data[UnstationaryObj.LonIndexer].sel(
            {UnstationaryObj.LonIndexer:min_lon}, method='nearest'))
        EasternLimit =float( NetCDF_data[UnstationaryObj.LonIndexer].sel(
            {UnstationaryObj.LonIndexer:max_lon}, method='nearest'))
        SouthernLimit = float(NetCDF_data[UnstationaryObj.LatIndexer].sel(
            {UnstationaryObj.LatIndexer:min_lat}, method='nearest'))
        NorthernLimit = float(NetCDF_data[UnstationaryObj.LatIndexer].sel(
            {UnstationaryObj.LatIndexer:max_lat}, method='nearest'))
    
        ## In the next code lines we will slice data for current timestep 
        # and for the specifed box -> Unstationary analysis.
        # ZA: zonal average (average through longitude circle)
        # AA: area average for the selected box
        # ZE: zonal eddy (value at each gridpoint minus its zonal average)
        # ZE_AA: zonal eddy avreaged in the selected box
        
        # Store area average and the averaged eddy component of temperature
        T = UnstationaryObj.Temperature.sel({UnstationaryObj.TimeIndexer:t}).sel(
            **{UnstationaryObj.LatIndexer:slice(SouthernLimit,NorthernLimit),
               UnstationaryObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        T_ZA = CalcZonalAverage(T)
        T_AA = CalcAreaAverage(T_ZA)
        T_ZE = T - T_ZA
        T_ZE_AA = CalcAreaAverage(T_ZE,ZonalAverage=True)
        dfDict['T_AA'][itime] = T_AA  
        dfDict['T_ZE_AA'][itime] = T_ZE_AA   
        
        
        # Store area average and the averaged eddy component of omega
        Omega = UnstationaryObj.Omega.sel({UnstationaryObj.TimeIndexer:t}).sel(
            **{UnstationaryObj.LatIndexer:slice(SouthernLimit,NorthernLimit),
               UnstationaryObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        Omega_ZA = CalcZonalAverage(Omega)
        Omega_AA = CalcAreaAverage(Omega_ZA)
        Omega_ZE = Omega - Omega_ZA
        Omega_ZE_AA = CalcAreaAverage(Omega_ZE,ZonalAverage=True)
        dfDict['Omega_AA'][itime] = Omega_AA
        dfDict['Omega_ZE_AA'][itime] = Omega_ZE_AA  
        
        # Store area average and the averaged eddy component of the adiabatic
        # heating
        Q = UnstationaryObj.AdiabaticHeating.sel({UnstationaryObj.TimeIndexer:t}).sel(
            **{UnstationaryObj.LatIndexer:slice(SouthernLimit,NorthernLimit),
               UnstationaryObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        Q_ZA = CalcZonalAverage(Q)
        Q_AA = CalcAreaAverage(Q_ZA)
        Q_ZE = Q - Q_ZA
        Q_ZE_AA = CalcAreaAverage(Q_ZE,ZonalAverage=True)
        dfDict['Q_AA'][itime] = Q_AA
        dfDict['Q_ZE_AA'][itime] = Q_ZE_AA
        
        # Plot maps of adiabatic heating term
        Q.name = "Q"
        Q['units'] = "J kg-1 s-1"
        plot_map(Q,MapsDirectory,"Q")
        Q_ZE.name = "Q'"
        Q_ZE['units'] = "J kg-1 s-1"
        
        # Plot temperature x omega (eddies)
        TO_ZE = T_ZE*Omega_ZE
        TO_ZE.name = "T'ω'"
        TO_ZE['units'] = "K Pa-1 s-1"
        
        # Plot temperature x Q (eddies)
        TQ_ZE = T_ZE*Q_ZE
        TQ_ZE.name = "T'Q'"
        TQ_ZE['units'] = "J K kg-1 s-1"
        
        if args.plot:
            plot_map(Omega,MapsDirectory,"Omega")
            plot_map(Omega_ZE,MapsDirectory,"Omega_ZE")
            plot_map(T,MapsDirectory,"T")
            plot_map(T_ZE,MapsDirectory,"T_ZE")
            plot_map(Q_ZE,MapsDirectory,"Q_ZE")
            plot_map(TO_ZE,MapsDirectory,"TOmega_ZE")
            plot_map(TO_ZE,MapsDirectory,"TOmega_ZE")
            plot_map(TQ_ZE,MapsDirectory,"TQ_ZE")

        
        # Store area average and of each thermodynamic equation
        AdvHTemp = UnstationaryObj.AdvHTemp.sel({UnstationaryObj.TimeIndexer:t}).sel(
            **{UnstationaryObj.LatIndexer:slice(SouthernLimit,NorthernLimit),
               UnstationaryObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        sigma = UnstationaryObj.sigma.sel({UnstationaryObj.TimeIndexer:t}).sel(
            **{UnstationaryObj.LatIndexer:slice(SouthernLimit,NorthernLimit),
               UnstationaryObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        dTdt = UnstationaryObj.dTdt.sel({UnstationaryObj.TimeIndexer:t}).sel(
            **{UnstationaryObj.LatIndexer:slice(SouthernLimit,NorthernLimit),
               UnstationaryObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        ResT = UnstationaryObj.ResT.sel({UnstationaryObj.TimeIndexer:t}).sel(
            **{UnstationaryObj.LatIndexer:slice(SouthernLimit,NorthernLimit),
               UnstationaryObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        
        # Store area averages for each term of the thermodynamic equation
        dfDict['AdvHTemp'][itime] = CalcAreaAverage(AdvHTemp,
                            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['SpOmega'][itime] = CalcAreaAverage(sigma*Omega,
                            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['dTdt'][itime] = CalcAreaAverage(dTdt,
                            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['ResT'][itime] = CalcAreaAverage(ResT,
                            ZonalAverage=True).metpy.convert_units('K/ day')
        
    # Save CSV files with all the terms stored above
    for term in stored_terms:
        dfDict[term].to_csv(ResultsSubDirectory+term+'.csv')
    # Make timeseries
    os.system("python plot_timeseries.py "+ResultsSubDirectory)
    
def StationaryAnalysis(StationaryObj):
    
    timesteps = StationaryObj.NetCDF_data[StationaryObj.TimeIndexer]
    
    PresLevels = StationaryObj.Temperature[StationaryObj.LevelIndexer].\
        metpy.convert_units('hPa').values
        
    stored_terms = ['T_AA','Omega_AA','Q_AA', # average terms
                    'T_ZE_AA','Omega_ZE_AA','Q_ZE_AA', # eddy terms (averaged)
                    'AdvHTemp','SpOmega','dTdt','ResT'] # thermodyn. eq. terms
    
    dfDict = {}
    for term in stored_terms:
        dfDict[term] = pd.DataFrame(columns=[str(t.values) for t in timesteps],
        index=[float(i) for i in PresLevels])
        
    # Store area average and the averaged eddy component of temperature
    T = StationaryObj.Temperature
    T_ZA = CalcZonalAverage(T)
    T_AA = CalcAreaAverage(T_ZA)
    T_ZE = T - T_ZA
    T_ZE_AA = CalcAreaAverage(T_ZE,ZonalAverage=True)
    # plot_map(T,MapsDirectory,"T")
    # plot_map(T_ZE,MapsDirectory,"T_ZE")
    
    # Store area average and the averaged eddy component of omega
    Omega = StationaryObj.Omega
    Omega_ZA = CalcZonalAverage(Omega)
    Omega_AA = CalcAreaAverage(Omega_ZA)
    Omega_ZE = Omega - Omega_ZA
    Omega_ZE_AA = CalcAreaAverage(Omega_ZE,ZonalAverage=True)
    
    # plot_map(Omega,MapsDirectory,"Omega")
    # plot_map(Omega_ZE,MapsDirectory,"Omega_ZE")
    
    # Store area average and the averaged eddy component of the adiabatic
    # heating
    Q = StationaryObj.AdiabaticHeating
    Q_ZA = CalcZonalAverage(Q)
    Q_AA = CalcAreaAverage(Q_ZA)
    Q_ZE = Q - Q_ZA
    Q_ZE_AA = CalcAreaAverage(Q_ZE,ZonalAverage=True)
            
    for t in timesteps:
        # Get current time and time strings
        itime = str(t.values)
        datestr = pd.to_datetime(itime).strftime('%Y-%m-%d %HZ')
        print("Storing results for: "+datestr)
        
        dfDict['T_AA'][itime] = T_AA.sel({StationaryObj.TimeIndexer:t})
        dfDict['T_ZE_AA'][itime] = T_ZE_AA.sel({StationaryObj.TimeIndexer:t})
        dfDict['Omega_AA'][itime] = Omega_AA.sel({StationaryObj.TimeIndexer:t})
        dfDict['Omega_ZE_AA'][itime] = Omega_ZE_AA.sel({StationaryObj.TimeIndexer:t})
        dfDict['Q_AA'][itime] = Q_AA.sel({StationaryObj.TimeIndexer:t})
        dfDict['Q_ZE_AA'][itime] = Q_ZE_AA.sel({StationaryObj.TimeIndexer:t})
          
        # Plot maps of adiabatic heating term
        Q.name = "Q"
        Q['units'] = "J kg-1 s-1"
        # plot_map(Q.sel({StationaryObj.TimeIndexer:t}),MapsDirectory,"Q")
        Q_ZE.name = "Q'"
        Q_ZE['units'] = "J kg-1 s-1"
        # plot_map(Q_ZE.sel({StationaryObj.TimeIndexer:t}),MapsDirectory,"Q_ZE")
        
        # Plot temperature x omega (eddies)
        TO_ZE = T_ZE*Omega_ZE
        TO_ZE.name = "T'ω'"
        TO_ZE['units'] = "K Pa-1 s-1"
        # plot_map(TO_ZE.sel({StationaryObj.TimeIndexer:t}),MapsDirectory,"TOmega_ZE")
        # Plot temperature x Q (eddies)
        TQ_ZE = T_ZE*Q_ZE
        TQ_ZE.name = "T'Q'"
        TQ_ZE['units'] = "J K kg-1 s-1"
        # plot_map(TQ_ZE.sel({StationaryObj.TimeIndexer:t}),MapsDirectory,"TQ_ZE")
        
        # Store area averages for each term of the thermodynamic equation
        dfDict['AdvHTemp'][itime] = CalcAreaAverage(StationaryObj.AdvHTemp.sel(
            {StationaryObj.TimeIndexer:t}),
            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['SpOmega'][itime] = CalcAreaAverage(StationaryObj.sigma.sel(
            {StationaryObj.TimeIndexer:t})*
            Omega.sel({StationaryObj.TimeIndexer:t}),
            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['dTdt'][itime] = CalcAreaAverage(StationaryObj.dTdt.sel(
            {StationaryObj.TimeIndexer:t}),
            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['ResT'][itime] = CalcAreaAverage(StationaryObj.ResT.sel(
            {StationaryObj.TimeIndexer:t}),
            ZonalAverage=True).metpy.convert_units('K/ day')
                                                 
    # Save CSV files with all the terms stored above
    for term in stored_terms:
        dfDict[term].to_csv(ResultsSubDirectory+term+'.csv')
    # Make timeseries
    os.system("python plot_timeseries.py "+ResultsSubDirectory)                                   
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "\
Program for analysing the thermodynamics of a cyclone. \
The program can use two distinct frameworks:\n\
    1) Stationary framework. A box is definid in the box_lims' file and then the\
 energetics are computed for a Stationary domain. (NOT IMPLEMENTED YET!)\
    2) Unstationary framework. The domain is not Stationary and follows the system \
using the 'track' file.\
Both frameworks can be applied at the same time, given the required files are\
provided. An auxilliary 'fvars' file is also needed for both frameworks. \
It contains the specified names used for each variable.  The results \
are stored as csv files in the 'CycloneThermodynamics_Results' directory on \
../ and it also creates figures for visualising the results.")
    parser.add_argument("infile", help = "input .nc file with temperature,\
   meridional, zonal and vertical components of the wind, in pressure levels")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--Stationary", default = False,
    action='store_true', help = "compute the energetics for a Stationary domain\
 specified by the 'inputs/box_lims' file.")
    group.add_argument("-t", "--track", default = False,
    action='store_true', help = "define the box using a track file specified \
by the 'inputs/track' file. The track indicate the central point of the system\
and a arbitraty box of 15°x15° is constructed.")
    group.add_argument("-c", "--choose", default = False,
    action='store_true', help = "For each time step, the user can choose the\
domain by clicking on the screen.")
    parser.add_argument("-m", "--multiple", default = False,
    action='store_true', help = "open multiple files at once")
    parser.add_argument("-p", "--plots", default = False,
    action='store_true', help = "create default plots")
    
    args = parser.parse_args()
    
    start_time = time.time()
    
    # Open namelist
    varlist = './inputs/fvars'
    dfVars = pd.read_csv(varlist,sep= ';',index_col=0,header=0)
    
    # Data indexers
    LonIndexer = dfVars.loc['Longitude']['Variable']
    LatIndexer = dfVars.loc['Latitude']['Variable']
    TimeIndexer = dfVars.loc['Time']['Variable']
    LevelIndexer = dfVars.loc['Vertical Level']['Variable']
    
    # Open file
    infile  = args.infile
    print('Opening file: '+infile)
    if args.multiple:
        NetCDF_data = convert_lon(xr.open_mfdataset(infile, ),
                              dfVars.loc['Longitude']['Variable'])
    else:
        NetCDF_data = convert_lon(xr.open_dataset(infile),
                              dfVars.loc['Longitude']['Variable'])
    print('Done.\nNow, running pre-processing stages...')
    # load data into memory (code optmization)
    NetCDF_data = NetCDF_data.load()
    # Assign lat and lon as radians, for calculations
    NetCDF_data = NetCDF_data.assign_coords(
        {"rlats": np.deg2rad(NetCDF_data[LatIndexer])})
    NetCDF_data = NetCDF_data.assign_coords(
        {"coslats": np.cos(np.deg2rad(NetCDF_data[LatIndexer]))})
    NetCDF_data = NetCDF_data.assign_coords(
        {"rlons": np.deg2rad(NetCDF_data[LonIndexer])})
    # Sort data coordinates as data from distinc sources might have different
    # arrangements and this might affect the results from the integrations
    NetCDF_data = NetCDF_data.sortby(LonIndexer).sortby(LevelIndexer).sortby(LatIndexer)
    print('Done.')
    
    # Directory where results will be stored
    ResultsMainDirectory = './Results/'
    # Append data method to outfile name
    outfile_name = ''.join(infile.split('/')[-1].split('.nc'))+'_lagranigan'
    # Each dataset of results have its own directory, allowing to store results
    # from more than one experiment at each time
    ResultsSubDirectory = ResultsMainDirectory+'/'+outfile_name+'/'
    # Subdirectory for saving figures and maps
    FigsDirectory = ResultsSubDirectory+'Figures/'
    MapsDirectory = FigsDirectory+'maps'
    # Check if the required directories exists. If not, creates them
    check_create_folder(ResultsMainDirectory)
    check_create_folder(ResultsSubDirectory)
    check_create_folder(FigsDirectory)
    check_create_folder(MapsDirectory)
    
    # Run the program
    if args.Stationary:
        print('Running Stationary framework.')
        dfbox = pd.read_csv('./inputs/box_limits',header=None,delimiter=';',index_col=0)
        StationaryObj = DataObject(NetCDF_data,dfVars,dfbox)
        StationaryAnalysis(StationaryObj)
        print("--- %s seconds running Stationary framework ---" % (time.time() - start_time))
    if args.track or args.choose:
        print('Running unstationary framework.')
        UnstationaryObj = DataObject(NetCDF_data,dfVars)
        UnstationaryAnalysis(UnstationaryObj)
        print("--- %s seconds for running Unstationary framework ---" % (time.time() - start_time))            