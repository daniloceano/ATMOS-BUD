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
from metpy.calc import vorticity
from metpy.calc import wind_speed
from metpy.constants import g

from scipy.signal import savgol_filter    

import pandas as pd
import xarray as xr
import os
import numpy as np
import argparse

from calc import CalcAreaAverage
from select_area import draw_box_map
from select_area import initial_domain

from plot_domains import plot_fixed_domain
from plot_domains import plot_track
from plot_domains import plot_min_zeta_hgt


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

def slice_domain(NetCDF_data, args, dfVars):
    
    # Data indexers
    LonIndexer = dfVars.loc['Longitude']['Variable']
    LatIndexer = dfVars.loc['Latitude']['Variable']
    TimeIndexer = dfVars.loc['Time']['Variable']
    LevelIndexer = dfVars.loc['Vertical Level']['Variable']
    
    if args.fixed:
        method = 'fixed'
        dfbox = pd.read_csv('../inputs/box_limits',header=None,
                            delimiter=';',index_col=0)
        WesternLimit = float(NetCDF_data[LonIndexer].sel(
            {LonIndexer:float(dfbox.loc['min_lon'].values)},
            method='nearest'))
        EasternLimit =float(NetCDF_data[LonIndexer].sel(
            {LonIndexer:float(dfbox.loc['max_lon'].values)},
            method='nearest'))
        SouthernLimit =float(NetCDF_data[LatIndexer].sel(
            {LatIndexer:float(dfbox.loc['min_lat'].values)},
            method='nearest'))
        NorthernLimit =float(NetCDF_data[LatIndexer].sel(
            {LatIndexer:float(dfbox.loc['max_lat'].values)},
            method='nearest'))
        
    elif args.track:
        trackfile = '../inputs/track'
        track = pd.read_csv(trackfile,parse_dates=[0],
                            delimiter=';',index_col='time')
        if 'width' in track.columns:
            method = 'track'
            min_width, max_width = track['width'].min(), track['width'].max()
            min_length, max_length = track['length'].min(), track['length'].max()
        else:
            method = 'track-15x15'
            min_width, max_width = 15, 15
            min_length, max_length = 15, 15
            
        WesternLimit = track['Lon'].min()-(min_width/2)
        EasternLimit = track['Lon'].max()+(max_width/2)
        SouthernLimit = track['Lat'].min()-(min_length/2)
        NorthernLimit = track['Lat'].max()+(max_length/2)
        
    elif args.choose:
        method = 'choose'
        iu_850 = NetCDF_data.isel({TimeIndexer:0}).sel({LevelIndexer:850}
                        )[dfVars.loc['Eastward Wind Component']['Variable']]
        iv_850 = NetCDF_data.isel({TimeIndexer:0}).sel({LevelIndexer:850}
                        )[dfVars.loc['Northward Wind Component']['Variable']]
        zeta = vorticity(iu_850, iv_850).metpy.dequantify()
        
        # Apply filter when using high resolution gridded data
        dx = float(iv_850[LonIndexer][1]-iv_850[LonIndexer][0])
        if dx < 1:
            zeta = zeta.to_dataset(name='vorticity'
                ).apply(savgol_filter,window_length=31, polyorder=2).vorticity
        
        lat, lon = iu_850[LatIndexer], iu_850[LonIndexer]
        domain_limits = initial_domain(zeta, lat, lon)
        WesternLimit = domain_limits['min_lon']
        EasternLimit = domain_limits['max_lon']
        SouthernLimit = domain_limits['min_lat']
        NorthernLimit = domain_limits['max_lat']
        
    NetCDF_data = NetCDF_data.sel(
        **{LatIndexer:slice(SouthernLimit,NorthernLimit),
           LonIndexer: slice(WesternLimit,EasternLimit)})
    
    return NetCDF_data, method

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
        self.NetCDF_data = NetCDF_data
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
        # Get the values for width and length in meters
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
        trackfile = '../inputs/track'
        track = pd.read_csv(trackfile,parse_dates=[0],delimiter=';',index_col='time')

    # Dictionary for saving system position and attributes
    position = {}
    results_keys = ['time', 'central_lat', 'central_lon', 'length', 'width',
            'min_zeta_850','min_hgt_850','max_wind_850']
    for key in results_keys:
        position[key] =  []
        
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
        
        itime = str(t.values)
        datestr = pd.to_datetime(itime).strftime('%Y-%m-%d %HZ')
        
        iu_850 = MovingObj.u.sel({TimeIndexer:t}).sel({LevelIndexer:850})
        iv_850 = MovingObj.v.sel({TimeIndexer:t}).sel({LevelIndexer:850})
        ight_850 = MovingObj.GeopotHeight.sel({TimeIndexer:t}
                                               ).sel({LevelIndexer:850})
        zeta = vorticity(iu_850, iv_850).metpy.dequantify()
        # Apply filter when using high resolution gridded data
        dx = float(iv_850[LonIndexer][1]-iv_850[LonIndexer][0])
        if dx < 1:
            zeta = zeta.to_dataset(name='vorticity'
                ).apply(savgol_filter,window_length=31, polyorder=2).vorticity
            
        lat, lon = iu_850[MovingObj.LatIndexer], iu_850[MovingObj.LonIndexer]
        
        if args.track:
            # Get current time and box limits
            if 'width'in track.columns:
                width, length = track.loc[itime]['width'],track.loc[itime]['length']
            else:
                width, length = 15, 15
            min_lon = track.loc[itime]['Lon']-(width/2)
            max_lon = track.loc[itime]['Lon']+(width/2)
            min_lat = track.loc[itime]['Lat']-(length/2)
            max_lat = track.loc[itime]['Lat']+(length/2)
            limits = {'min_lon':min_lon,'max_lon':max_lon,
                      'min_lat':min_lat,'max_lat':max_lat}
        
        elif args.choose:
            # Draw maps and ask user to specify corners for specifying the box
            limits = draw_box_map(iu_850, iv_850, zeta, ight_850,
                                  lat, lon, itime)
            min_lon, max_lon = limits['min_lon'],  limits['max_lon']
            min_lat, max_lat = limits['min_lat'],  limits['max_lat']
        
        # Store system position and attributes
        central_lat = (limits['max_lat'] + limits['min_lat'])/2
        central_lon = (limits['max_lon'] + limits['min_lon'])/2
        length = limits['max_lat'] - limits['min_lat']
        width = limits['max_lon'] - limits['min_lon']
        min_zeta = float(zeta.min())
        min_hgt = float(ight_850.min())
        max_wind = float(wind_speed(iu_850, iv_850).max())
        
        values = [datestr, central_lat, central_lon, length, width,
                  min_zeta, min_hgt, max_wind]
        for key,val in zip(results_keys,values):
            position[key].append(val)
            
        print('\nTime: ',datestr)
        print('Box min_lon, max_lon: '+str(min_lon)+'/'+str(max_lon))
        print('Box min_lat, max_lat: '+str(min_lat)+'/'+str(max_lat))
        print('Box size (longitude): '+str(width))
        print('Box size (latitude): '+str(length))
        print('Minimum vorticity at 850 hPa:',min_zeta)
        print('Minimum geopotential height at 850 hPa:',min_hgt)
        print('Maximum wind speed at 850 hPa:',max_wind)
        
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
    term_results = [MovingObj.AdvHTemp, MovingObj.sigma*MovingObj.Omega,
                    MovingObj.dTdt,MovingObj.ResT]
    results_nc = []
    for term, da in zip(stored_terms,term_results):
        da =  da.metpy.dequantify()
        da.name = term
        results_nc.append(da)
    out_nc = xr.merge(results_nc)
    fname = ResultsSubDirectory+ outfile_name+'.nc'
    out_nc.to_netcdf(fname)
    print(fname+' created')
    
    # Save CSV files with all the terms stored above
    for term in stored_terms:
        dfDict[term].to_csv(ResultsSubDirectory+term+'.csv')
        
    # Save system position as a csv file for replicability
    track = pd.DataFrame.from_dict(position)
    track = track.rename(columns={'central_lat':'Lat','central_lon':'Lon'})
    track.to_csv(ResultsSubDirectory+outfile_name+'_track',
                 index=False, sep=";")

    plot_track(track, FigsDirectory)
    plot_min_zeta_hgt(track, FigsDirectory)
    # Make timeseries
    #os.system("python plot_timeseries.py "+ResultsSubDirectory)
    
def FixedAnalysis(FixedObj):
    
    timesteps = FixedObj.NetCDF_data[FixedObj.TimeIndexer]
    
    pres_levels = FixedObj.Temperature[FixedObj.LevelIndexer].\
        metpy.convert_units('hPa').values
        
    stored_terms = ['AdvHTemp','SpOmega','dTdt','ResT'] # thermodyn. eq. terms
    
    # Dictionary to save DataArray results and transform into nc later
    results_nc = {}
    
    dfDict = {}
    for term in stored_terms:
        dfDict[term] = pd.DataFrame(columns=[str(t.values) for t in timesteps],
        index=[float(i) for i in pres_levels])
        results_nc[term] = []
        
    for t in timesteps:
        # Get current time and time strings
        itime = str(t.values)
        datestr = pd.to_datetime(itime).strftime('%Y-%m-%d %HZ')
        print("Storing results for: "+datestr)
                
        # Store area averages for each term of the thermodynamic equation
        dfDict['AdvHTemp'][itime] = CalcAreaAverage(FixedObj.AdvHTemp.sel(
            {FixedObj.TimeIndexer:t}),
            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['SpOmega'][itime] = CalcAreaAverage(FixedObj.sigma.sel(
            {FixedObj.TimeIndexer:t})*
            FixedObj.Omega.sel({FixedObj.TimeIndexer:t}),
            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['dTdt'][itime] = CalcAreaAverage(FixedObj.dTdt.sel(
            {FixedObj.TimeIndexer:t}),
            ZonalAverage=True).metpy.convert_units('K/ day')
        dfDict['ResT'][itime] = CalcAreaAverage(FixedObj.ResT.sel(
            {FixedObj.TimeIndexer:t}),
            ZonalAverage=True).metpy.convert_units('K/ day')
    
    term_results = [FixedObj.AdvHTemp, FixedObj.Omega*FixedObj.sigma,
                    FixedObj.dTdt, FixedObj.ResT]
    for term,r in zip(stored_terms,term_results):
        results_nc[term].append(r.metpy.dequantify().rename(term))
    
    print('saving results to nc file...')
    term_results = [FixedObj.AdvHTemp, FixedObj.sigma*FixedObj.Omega,
                    FixedObj.dTdt,FixedObj.ResT]
    results_nc = []
    for term, da in zip(stored_terms,term_results):
        da =  da.metpy.dequantify()
        da.name = term
        results_nc.append(da)
    out_nc = xr.merge(results_nc)
    fname = ResultsSubDirectory+ outfile_name+'.nc'
    out_nc.to_netcdf(fname)
    print(fname+' created')
                                                 
    # Save CSV files with all the terms stored above
    for term in stored_terms:
        dfDict[term].to_csv(ResultsSubDirectory+term+'.csv')
        
    min_lon = FixedObj.NetCDF_data[FixedObj.LonIndexer].min()
    max_lon = FixedObj.NetCDF_data[FixedObj.LonIndexer].max()
    min_lat = FixedObj.NetCDF_data[FixedObj.LatIndexer].min()
    max_lat = FixedObj.NetCDF_data[FixedObj.LatIndexer].max()
    plot_fixed_domain(min_lon, max_lon, min_lat, max_lat, ResultsSubDirectory)                                
    
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
    varlist = '../inputs/fvars'
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
    # Sort data coordinates as data from distinc sources might have different
    # arrangements and this might affect the results from the integrations
    NetCDF_data = NetCDF_data.sortby(LonIndexer).sortby(LevelIndexer
                                                        ).sortby(LatIndexer)
    # Slice data so the code runs faster
    NetCDF_data, method = slice_domain(NetCDF_data, args, dfVars)
    
    # Assign lat and lon as radians, for calculations
    NetCDF_data = NetCDF_data.assign_coords(
        {"rlats": np.deg2rad(NetCDF_data[LatIndexer])})
    NetCDF_data = NetCDF_data.assign_coords(
        {"coslats": np.cos(np.deg2rad(NetCDF_data[LatIndexer]))})
    NetCDF_data = NetCDF_data.assign_coords(
        {"rlons": np.deg2rad(NetCDF_data[LonIndexer])})
    print('Done.')
    
    # Directory where results will be stored
    ResultsMainDirectory = '../Results/'
    # Append data method to outfile name
    outfile_name = ''.join(infile.split('/')[-1].split('.nc'))+'_'+method
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
    # backup of box_limits and track file for reproductability
    if args.fixed:
        os.system('cp ../inputs/box_limits '+ResultsSubDirectory)
    elif args.track:
        os.system('cp ../inputs/track '+ResultsSubDirectory)
    
    # Run the program
    if args.fixed:
        print('Running Fixed framework.')
        FixedObj = DataObject(NetCDF_data,dfVars)
        FixedAnalysis(FixedObj)
        print("--- %s seconds running Fixed framework ---" % (
            time.time() - start_time))
    if args.track or args.choose:
        print('Running Moving framework.')
        MovingObj = DataObject(NetCDF_data,dfVars)
        MovingAnalysis(MovingObj,args)
        print("--- %s seconds for running Moving framework ---" % (
            time.time() - start_time))            