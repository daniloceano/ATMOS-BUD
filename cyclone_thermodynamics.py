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
from metpy.calc import altimeter_to_sea_level_pressure
from metpy.calc import wind_speed
from metpy.constants import g

import pandas as pd
import xarray as xr
import os
import numpy as np
import argparse

from calc import CalcAreaAverage, CalcZonalAverage
from select_area import draw_box_map
from select_area import initial_domain

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

def choose_domain_analysis(MovingObj, position, t, domain_limits):
    
    itime = str(t.values)
    
    iu_1000 = MovingObj.u.sel({MovingObj.TimeIndexer:t,
                                   MovingObj.LevelIndexer:1000})
    iv_1000 = MovingObj.v.sel({MovingObj.TimeIndexer:t,
                                   MovingObj.LevelIndexer:1000})
    # Geopotential from geopotential height
    igeop_1000 = MovingObj.GeopotHeight.sel(
        {MovingObj.TimeIndexer:t,
         MovingObj.LevelIndexer:1000})*g
    iheight_1000 = geopotential_to_height(igeop_1000)
    
    lat = iu_1000[MovingObj.LatIndexer]
    lon = iu_1000[MovingObj.LonIndexer]
    
    # # Reduce to sea level pressure
    # it_1000 = MovingObj.Temperature.sel({MovingObj.TimeIndexer:t,
    #                                MovingObj.LevelIndexer:1000})
    # ip_1000 = it_1000[MovingObj.LevelIndexer].expand_dims(
    #     {MovingObj.LatIndexer:it_1000[MovingObj.LatIndexer],
    #      MovingObj.LonIndexer:it_1000[MovingObj.LonIndexer]}
    #     )*units(it_1000[MovingObj.LevelIndexer].units)
    # imslp = altimeter_to_sea_level_pressure(ip_1000, iheight_1000, it_1000)
    
    # imslp = height_to_pressure_std(iheight_1000)
    
    imslp = iheight_1000.copy()
    
    # Select initial domain so its easier to see the systems on the map
    if not domain_limits:
        domain_limits = initial_domain(imslp, lat, lon)
    else:
        domain_limits = domain_limits
        
    iu_1000 = iu_1000.sel({MovingObj.LatIndexer:slice(
        domain_limits['min_lat'],domain_limits['max_lat'])}).sel({
        MovingObj.LonIndexer:slice(
            domain_limits['min_lon'],domain_limits['max_lon'])})
    iv_1000 = iv_1000.sel({MovingObj.LatIndexer:slice(
        domain_limits['min_lat'],domain_limits['max_lat'])}).sel({
        MovingObj.LonIndexer:slice(
            domain_limits['min_lon'],domain_limits['max_lon'])})
    imslp = imslp.sel({MovingObj.LatIndexer:slice(
        domain_limits['min_lat'],domain_limits['max_lat'])}).sel({
        MovingObj.LonIndexer:slice(
            domain_limits['min_lon'],domain_limits['max_lon'])})
    lat = iu_1000[MovingObj.LatIndexer]
    lon = iu_1000[MovingObj.LonIndexer]
    
    # Draw maps and ask user to specify corners for specifying the box
    limits = draw_box_map(iu_1000, iv_1000, imslp, lat, lon,
                          itime, domain_limits)
    
    # Store results
    time = pd.to_datetime(itime).strftime('%Y-%m-%d %H%M')
    central_lat = (limits['max_lat'] + limits['min_lat'])/2
    central_lon = (limits['max_lon'] + limits['min_lon'])/2
    dy = (limits['max_lat'] - limits['min_lat'])/2
    dx = (limits['max_lat'] - limits['min_lon'])/2
    min_slp = float(imslp.min())
    max_wind = float(wind_speed(iu_1000, iv_1000).max())
    
    keys = ['time', 'central_lat', 'central_lon', 'dy', 'dx',
            'min_slp','max_wind']
    values = [time, central_lat, central_lon, dy, dx,
              min_slp, max_wind]

    if len(position) == 0:
        for key,val in zip(keys,values):
            position[key] = [val]
    else:
        for key,val in zip(keys,values):
            position[key].append(val)
    
    return position, limits, domain_limits

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
        # When constructing object for Fixed analysis, the data can
        # be sliced beforehand using the dfBox limits, but for the Moving
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
        if 'Geopotential Height' in dfVars.index:
            self.GeopotHeight = self.NetCDF_data[dfVars.loc[
                'Geopotential Height']['Variable']] * units(dfVars.loc[
                    'Geopotential Height']['Units']).to('gpm')
        else:
            self.GeopotHeight = (self.NetCDF_data[dfVars.loc[
                'Geopotential']['Variable']] * units(dfVars.loc[
                    'Geopotential']['Units']).to('m**2/s**2'))/g
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

def MovingAnalysis(MovingObj,args):
    """
    Parameters
    ----------
    MovingObj : DataObj
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
    else:
        position = {}
        
    # Dictionary to save DataArray results and transform into nc later
    results_nc = {}
    
    timesteps = MovingObj.NetCDF_data[MovingObj.TimeIndexer]
    pres_levels = MovingObj.Temperature[MovingObj.LevelIndexer].\
        metpy.convert_units('hPa').values
    
    # Create a dictionary for saving area averages of each term
    stored_terms = ['AdvHTemp','SpOmega','dTdt','ResT'] 
    dfDict = {}
    for term in stored_terms:
        dfDict[term] = pd.DataFrame(columns=[str(t.values) for t in timesteps],
        index=[float(i) for i in pres_levels])
        results_nc[term] = []
    
    for t in timesteps:
        # Get current time and time strings
        itime = str(t.values)
        datestr = pd.to_datetime(itime).strftime('%Y-%m-%d %HZ')
        
        if args.track:
            # Get current time and box limits
            if 'dx'in track.columns:
                dx, dy = track.loc[itime]['dx'],track.loc[itime]['dy']
            else:
                dx, dy = 7.5, 7.5
            min_lon = track.loc[itime]['Lon']-dx
            max_lon = track.loc[itime]['Lon']+dx
            min_lat = track.loc[itime]['Lat']-dy
            max_lat = track.loc[itime]['Lat']+dy
            
        elif args.choose:
            if t == timesteps[0]:
                position, limits, domain_limits = choose_domain_analysis(
                MovingObj, position, t, None)
            else:
                position, limits, domain_limits = choose_domain_analysis(
                MovingObj, position, t, domain_limits)
            dx, dy = position['dx'][-1],position['dy'][-1]
            min_lon, max_lon = limits['min_lon'],  limits['max_lon']
            min_lat, max_lat = limits['min_lat'],  limits['max_lat']
            
        print('\nComputing terms for '+datestr+'...')
        print('Box min_lon, max_lon: '+str(min_lon)+'/'+str(max_lon))
        print('Box min_lat, max_lat: '+str(min_lat)+'/'+str(max_lat))
        print('Box size (longitude): '+str(dx*2))
        print('Box size (latitude): '+str(dy*2))
        
        NetCDF_data = MovingObj.NetCDF_data
        # Get closest grid point to actual track
        WesternLimit = float(NetCDF_data[MovingObj.LonIndexer].sel(
            {MovingObj.LonIndexer:min_lon}, method='nearest'))
        EasternLimit =float( NetCDF_data[MovingObj.LonIndexer].sel(
            {MovingObj.LonIndexer:max_lon}, method='nearest'))
        SouthernLimit = float(NetCDF_data[MovingObj.LatIndexer].sel(
            {MovingObj.LatIndexer:min_lat}, method='nearest'))
        NorthernLimit = float(NetCDF_data[MovingObj.LatIndexer].sel(
            {MovingObj.LatIndexer:max_lat}, method='nearest'))
    
        sigma = MovingObj.sigma.sel({MovingObj.TimeIndexer:t}).sel(
            **{MovingObj.LatIndexer:slice(SouthernLimit,NorthernLimit),
               MovingObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        omega = MovingObj.Omega.sel({MovingObj.TimeIndexer:t}).sel(
            **{MovingObj.LatIndexer:slice(SouthernLimit,NorthernLimit),
               MovingObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        
        # Compute each term of the thermodynamic equation
        AdvHTemp = MovingObj.AdvHTemp.sel({MovingObj.TimeIndexer:t}).sel(
            **{MovingObj.LatIndexer:slice(SouthernLimit,NorthernLimit),
               MovingObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        dTdt = MovingObj.dTdt.sel({MovingObj.TimeIndexer:t}).sel(
            **{MovingObj.LatIndexer:slice(SouthernLimit,NorthernLimit),
               MovingObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        ResT = MovingObj.ResT.sel({MovingObj.TimeIndexer:t}).sel(
            **{MovingObj.LatIndexer:slice(SouthernLimit,NorthernLimit),
               MovingObj.LonIndexer: slice(WesternLimit,EasternLimit)})
        SpOmega = omega * sigma
        
        term_results = [AdvHTemp,SpOmega,dTdt,ResT]
        for term,r in zip(stored_terms,term_results):
            results_nc[term].append(r.metpy.dequantify().rename(term))
        
        # Store area averages for each term of the thermodynamic equation
        dfDict['AdvHTemp'][itime] = CalcAreaAverage(AdvHTemp,
                            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['SpOmega'][itime] = CalcAreaAverage(SpOmega,
                            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['dTdt'][itime] = CalcAreaAverage(dTdt,
                            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['ResT'][itime] = CalcAreaAverage(ResT,
                            ZonalAverage=True).metpy.convert_units('K/ day')
      
    print('saving results to nc file...')
    dict_nc = {}
    for term in stored_terms:
        dict_nc[term] = xr.concat(
            results_nc[term],dim=MovingObj.TimeIndexer)
    out_nc = xr.merge([dict_nc])
    fname = ResultsSubDirectory+ outfile_name+'.nc'
    out_nc.to_netcdf(fname)
    print(fname+' created')
    
    # Save CSV files with all the terms stored above
    for term in stored_terms:
        dfDict[term].to_csv(ResultsSubDirectory+term+'.csv')
        
    if args.choose:
        # Save system position as a csv file for replicability
        track = pd.DataFrame.from_dict(position)
        track = track.rename(columns={'central_lat':'Lat','central_lon':'Lon'})
        track.to_csv(ResultsSubDirectory+outfile_name+'_track',
                     index=False, sep=";")

    # Make timeseries
    #os.system("python plot_timeseries.py "+ResultsSubDirectory)
    
def FixedAnalysis(FixedObj):
    
    timesteps = FixedObj.NetCDF_data[FixedObj.TimeIndexer]
    
    pres_levels = FixedObj.Temperature[FixedObj.LevelIndexer].\
        metpy.convert_units('hPa').values
        
    stored_terms = ['T_AA','Omega_AA','Q_AA', # average terms
                    'T_ZE_AA','Omega_ZE_AA','Q_ZE_AA', # eddy terms (averaged)
                    'AdvHTemp','SpOmega','dTdt','ResT'] # thermodyn. eq. terms
    
    dfDict = {}
    for term in stored_terms:
        dfDict[term] = pd.DataFrame(columns=[str(t.values) for t in timesteps],
        index=[float(i) for i in pres_levels])
        
    # Store area average and the averaged eddy component of temperature
    T = FixedObj.Temperature
    T_ZA = CalcZonalAverage(T)
    T_AA = CalcAreaAverage(T_ZA)
    T_ZE = T - T_ZA
    T_ZE_AA = CalcAreaAverage(T_ZE,ZonalAverage=True)
    # plot_map(T,MapsDirectory,"T")
    # plot_map(T_ZE,MapsDirectory,"T_ZE")
    
    # Store area average and the averaged eddy component of omega
    Omega = FixedObj.Omega
    Omega_ZA = CalcZonalAverage(Omega)
    Omega_AA = CalcAreaAverage(Omega_ZA)
    Omega_ZE = Omega - Omega_ZA
    Omega_ZE_AA = CalcAreaAverage(Omega_ZE,ZonalAverage=True)
    
    # plot_map(Omega,MapsDirectory,"Omega")
    # plot_map(Omega_ZE,MapsDirectory,"Omega_ZE")
    
    # Store area average and the averaged eddy component of the adiabatic
    # heating
    Q = FixedObj.AdiabaticHeating
    Q_ZA = CalcZonalAverage(Q)
    Q_AA = CalcAreaAverage(Q_ZA)
    Q_ZE = Q - Q_ZA
    Q_ZE_AA = CalcAreaAverage(Q_ZE,ZonalAverage=True)
            
    for t in timesteps:
        # Get current time and time strings
        itime = str(t.values)
        datestr = pd.to_datetime(itime).strftime('%Y-%m-%d %HZ')
        print("Storing results for: "+datestr)
        
        dfDict['T_AA'][itime] = T_AA.sel({FixedObj.TimeIndexer:t})
        dfDict['T_ZE_AA'][itime] = T_ZE_AA.sel({FixedObj.TimeIndexer:t})
        dfDict['Omega_AA'][itime] = Omega_AA.sel({FixedObj.TimeIndexer:t})
        dfDict['Omega_ZE_AA'][itime] = Omega_ZE_AA.sel({FixedObj.TimeIndexer:t})
        dfDict['Q_AA'][itime] = Q_AA.sel({FixedObj.TimeIndexer:t})
        dfDict['Q_ZE_AA'][itime] = Q_ZE_AA.sel({FixedObj.TimeIndexer:t})
          
        # Plot maps of adiabatic heating term
        Q.name = "Q"
        Q['units'] = "J kg-1 s-1"
        # plot_map(Q.sel({FixedObj.TimeIndexer:t}),MapsDirectory,"Q")
        Q_ZE.name = "Q'"
        Q_ZE['units'] = "J kg-1 s-1"
        # plot_map(Q_ZE.sel({FixedObj.TimeIndexer:t}),MapsDirectory,"Q_ZE")
        
        # Plot temperature x omega (eddies)
        TO_ZE = T_ZE*Omega_ZE
        TO_ZE.name = "T'ω'"
        TO_ZE['units'] = "K Pa-1 s-1"
        # plot_map(TO_ZE.sel({FixedObj.TimeIndexer:t}),MapsDirectory,"TOmega_ZE")
        # Plot temperature x Q (eddies)
        TQ_ZE = T_ZE*Q_ZE
        TQ_ZE.name = "T'Q'"
        TQ_ZE['units'] = "J K kg-1 s-1"
        # plot_map(TQ_ZE.sel({FixedObj.TimeIndexer:t}),MapsDirectory,"TQ_ZE")
        
        # Store area averages for each term of the thermodynamic equation
        dfDict['AdvHTemp'][itime] = CalcAreaAverage(FixedObj.AdvHTemp.sel(
            {FixedObj.TimeIndexer:t}),
            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['SpOmega'][itime] = CalcAreaAverage(FixedObj.sigma.sel(
            {FixedObj.TimeIndexer:t})*
            Omega.sel({FixedObj.TimeIndexer:t}),
            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['dTdt'][itime] = CalcAreaAverage(FixedObj.dTdt.sel(
            {FixedObj.TimeIndexer:t}),
            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['ResT'][itime] = CalcAreaAverage(FixedObj.ResT.sel(
            {FixedObj.TimeIndexer:t}),
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
    1) Fixed framework. A box is definid in the box_lims' file and then the\
 energetics are computed for a Fixed domain. (NOT IMPLEMENTED YET!)\
    2) Moving framework. The domain is not Fixed and follows the system \
using the 'track' file.\
Both frameworks can be applied at the same time, given the required files are\
provided. An auxilliary 'fvars' file is also needed for both frameworks. \
It contains the specified names used for each variable.  The results \
are stored as csv files in the 'CycloneThermodynamics_Results' directory on \
../ and it also creates figures for visualising the results.")
    parser.add_argument("infile", help = "input .nc file with temperature,\
   meridional, zonal and vertical components of the wind, in pressure levels")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fixed", default = False,
    action='store_true', help = "compute the energetics for a Fixed domain\
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
    if args.fixed:
        outfile_name = ''.join(infile.split('/')[-1].split('.nc'))+'_fixed'
    elif args.track:
        outfile_name = ''.join(infile.split('/')[-1].split('.nc'))+'_track'
    elif args.choose:
        outfile_name = ''.join(infile.split('/')[-1].split('.nc'))+'_choose'
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
    if args.fixed:
        print('Running Fixed framework.')
        dfbox = pd.read_csv('./inputs/box_limits',header=None,delimiter=';',index_col=0)
        FixedObj = DataObject(NetCDF_data,dfVars,dfbox)
        FixedAnalysis(FixedObj)
        print("--- %s seconds running Fixed framework ---" % (time.time() - start_time))
    if args.track or args.choose:
        print('Running Moving framework.')
        MovingObj = DataObject(NetCDF_data,dfVars)
        MovingAnalysis(MovingObj,args)
        print("--- %s seconds for running Moving framework ---" % (time.time() - start_time))            