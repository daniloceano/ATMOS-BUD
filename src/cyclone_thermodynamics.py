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
from metpy.calc import divergence
from metpy.calc import coriolis_parameter
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

import dask

import time

def check_create_folder(DirName, verbose=True):
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
        if verbose:
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
                 dTdt: xr.DataArray, dZdt: xr.DataArray,
                 dfVars: pd.DataFrame, args: argparse.ArgumentParser):
        self.LonIndexer = dfVars.loc['Longitude']['Variable']
        self.LatIndexer = dfVars.loc['Latitude']['Variable']
        self.TimeIndexer = dfVars.loc['Time']['Variable']
        self.LevelIndexer = dfVars.loc['Vertical Level']['Variable']
        self.NetCDF_data = NetCDF_data
        
        # When using mfdataset the variables are not DataArrays! 
        # Need to individually transform variables into DataArrays WHEN
        # performing computations!!! Who was the mf that programmed this?
        if args.gfs:
            print("loading dask array into memory, if it's too slow, complain\
with the guys from metpy")
            self.NetCDF_data = self.NetCDF_data.load()
                
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
        self.Pressure = (self.NetCDF_data[self.LevelIndexer]
                ) * units(self.NetCDF_data[self.LevelIndexer].units).to('Pa')
        
        self.lons  = self.NetCDF_data[self.LonIndexer]
        self.lats = self.NetCDF_data[self.LatIndexer]
        self.cos_lats = self.NetCDF_data["coslats"]
        # Get the values for width and length in meters
        self.dx = np.deg2rad(self.lons.differentiate(self.LonIndexer)
                        )*self.cos_lats*Re
        self.dy = np.deg2rad(self.lats.differentiate(self.LatIndexer))*Re
        self.f = coriolis_parameter(self.lats)
            
        ## Thermodynamic terms
        self.Theta = theta = potential_temperature(
            self.Pressure,self.Temperature)
        # self.dTdt = self.Temperature.differentiate(
        #         self.TimeIndexer,datetime_unit='s') / units('s')
        self.dTdt = dTdt
        self.AdvHTemp = self.HorizontalTemperatureAdvection()
        self.AdvVTemp = -1 * (self.Temperature.differentiate(self.LevelIndexer
                    ) * self.Omega) / (1 * self.Pressure.metpy.units).to('Pa')
        self.Sigma = (-1 * (self.Temperature/theta) * theta.differentiate(
            self.LevelIndexer) / units(str(self.Pressure.metpy.units))
            ).metpy.convert_units('K / Pa')
        self.ResT =  self.dTdt - self.AdvHTemp - (self.Sigma * self.Omega)
        self.AdiabaticHeating = self.ResT*Cp_d
        
        ## Vorticity terms
        self.Zeta = vorticity(self.u, self.v)
        # self.dZdt = self.Zeta.differentiate(
        #         self.TimeIndexer,datetime_unit='s') / units('s')
        self.dZdt = dZdt
        self.AdvHZeta = self.HorizontalVorticityAdvection()
        self.AdvVZeta = -1 * (self.Zeta.differentiate(self.LevelIndexer
                    ) * self.Omega) / (1 * self.Pressure.metpy.units).to('Pa')
        self.Beta = (self.f).differentiate(self.LatIndexer)/self.dy
        self.vxBeta = -1 * (self.v * self.Beta)
        self.DivH = divergence(self.u, self.v)
        self.fDivH = -1 * self.f * self.DivH
        self.ZetaDivH = -1 * (self.Zeta * self.DivH)
        self.Tilting = self.TiltingTerm()
        self.SumVorticity = self.AdvHZeta + self.AdvVZeta + self.vxBeta + \
            self.ZetaDivH + self.fDivH + self.Tilting
        self.ResZ = self.dZdt - self.SumVorticity
        
    def HorizontalTemperatureAdvection(self):
        dTdlambda = self.Temperature.differentiate(self.LonIndexer)
        dTdphi = self.Temperature.differentiate(self.LatIndexer)
        AdvHT = -1 * ((self.u*dTdlambda/self.dx)+(self.v*dTdphi/self.dy)) 
        return AdvHT
    
    def HorizontalVorticityAdvection(self):
        dZdlambda = self.Zeta.differentiate(self.LonIndexer)
        dZdphi = self.Zeta.differentiate(self.LatIndexer)  
        AdvHZeta = -1 * ((self.u*dZdlambda/self.dx) + (self.v*dZdphi/self.dy)) 
        return AdvHZeta
    
    def TiltingTerm(self):
        dOmegalambda = self.Omega.differentiate(self.LonIndexer)
        dOmegadphi = self.Omega.differentiate(self.LatIndexer)
        dOmegady = dOmegalambda/self.dy
        dOmegadx = dOmegadphi/self.dx
        dudp = self.u.differentiate(self.LevelIndexer
                                ) / (1 * self.Pressure.metpy.units).to('Pa')
        dvdp = self.v.differentiate(self.LevelIndexer
                                ) /(1 * self.Pressure.metpy.units).to('Pa')
        return (dOmegady*dudp) - (dOmegadx*dvdp)
    

def cyclone_thermodynamics(NetCDF_data, dfVars, dTdt, dZdt, args):
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
    
    # Indexers
    LonIndexer = dfVars.loc['Longitude']['Variable']
    LatIndexer = dfVars.loc['Latitude']['Variable']
    TimeIndexer = dfVars.loc['Time']['Variable']
    LevelIndexer = dfVars.loc['Vertical Level']['Variable']
    
    # timesteps = MovingObj.NetCDF_data[TimeIndexer]
    # pres_levels = MovingObj.Temperature[MovingObj.LevelIndexer].\
    #     metpy.convert_units('hPa').values
    timesteps = NetCDF_data[TimeIndexer]
    pres_levels = NetCDF_data[LevelIndexer].values
    
    # Create a dictionary for saving area averages of each term
    stored_terms = ['AdvHTemp','AdvVTemp', 'Sigma','Omega','dTdt','ResT', 
                    'Zeta', 'dZdt','AdvHZeta','AdvVZeta', 'vxBeta',
                   'ZetaDivH','fDivH', 'Tilting', 'ResZ'] 
    dfDict = {}
    for term in stored_terms:
        dfDict[term] = pd.DataFrame(columns=[str(t.values) for t in timesteps],
        index=[float(i) for i in pres_levels])
        results_nc[term] = []

    if args.fixed:
        dfbox = pd.read_csv('../inputs/box_limits',header=None,
                            delimiter=';',index_col=0)
        WesternLimit = min_lon = float(NetCDF_data[LonIndexer].sel(
            {LonIndexer:float(dfbox.loc['min_lon'].values)},
            method='nearest'))
        EasternLimit = max_lon =float(NetCDF_data[LonIndexer].sel(
            {LonIndexer:float(dfbox.loc['max_lon'].values)},
            method='nearest'))
        SouthernLimit = min_lat =float(NetCDF_data[LatIndexer].sel(
            {LatIndexer:float(dfbox.loc['min_lat'].values)},
            method='nearest'))
        NorthernLimit = max_lat =float(NetCDF_data[LatIndexer].sel(
            {LatIndexer:float(dfbox.loc['max_lat'].values)},
            method='nearest'))
    
    for t in timesteps:
        
        itime = str(t.values)
        datestr = pd.to_datetime(itime).strftime('%Y-%m-%d %HZ')
        datestr2 = pd.to_datetime(itime).strftime('%Y%m%d%H00')
        
        MovingObj = DataObject(NetCDF_data.sel({TimeIndexer:t}), dTdt=dTdt.sel({TimeIndexer:t}),
                                 dZdt=dZdt.sel({TimeIndexer:t}), dfVars=dfVars, args=args)

        # iu_850 = MovingObj.u.sel({TimeIndexer:t}).sel({LevelIndexer:850})
        # iv_850 = MovingObj.v.sel({TimeIndexer:t}).sel({LevelIndexer:850})
        # ight_850 = MovingObj.GeopotHeight.sel({TimeIndexer:t}
        #                                        ).sel({LevelIndexer:850})
        iu_850 = MovingObj.u.sel({LevelIndexer:850})
        iv_850 = MovingObj.v.sel({LevelIndexer:850})
        ight_850 = MovingObj.GeopotHeight.sel({LevelIndexer:850})
        
        # Filter doesn't work well with metpy units
        # zeta = vorticity(iu_850, iv_850).metpy.dequantify()
        zeta = MovingObj.Zeta.sel({LevelIndexer:850}).metpy.dequantify()
        # Apply filter when using high resolution gridded data
        # dx = float(iv_850[LonIndexer][1]-iv_850[LonIndexer][0])
        dx = float(NetCDF_data[LonIndexer][1]-NetCDF_data[LonIndexer][0])
        if dx < 1:
            zeta = zeta.to_dataset(name='vorticity'
                ).apply(savgol_filter,window_length=31, polyorder=2).vorticity
            
        # lat, lon = iu_850[MovingObj.LatIndexer], iu_850[MovingObj.LonIndexer]
        lat, lon = NetCDF_data[MovingObj.LatIndexer], NetCDF_data[MovingObj.LonIndexer]

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
        
        # Save figure with box used for computations
        plot_fixed_domain(min_lon, max_lon, min_lat, max_lat, ResultsSubDirectory,
                       time=datestr2, zeta=zeta, lat=zeta[LatIndexer], lon=zeta[LonIndexer], hgt=ight_850)
        
        if args.fixed:
            print("Storing results for: "+datestr)

        else:
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
            
            # NetCDF_data = MovingObj.NetCDF_data
            # Get closest grid point to actual track
            # WesternLimit = float(NetCDF_data[MovingObj.LonIndexer].sel(
            #     {MovingObj.LonIndexer:min_lon}, method='nearest'))
            # EasternLimit =float( NetCDF_data[MovingObj.LonIndexer].sel(
            #     {MovingObj.LonIndexer:max_lon}, method='nearest'))
            # SouthernLimit = float(NetCDF_data[MovingObj.LatIndexer].sel(
            #     {MovingObj.LatIndexer:min_lat}, method='nearest'))
            # NorthernLimit = float(NetCDF_data[MovingObj.LatIndexer].sel(
            #     {MovingObj.LatIndexer:max_lat}, method='nearest'))
            WesternLimit = float(NetCDF_data[LonIndexer].sel(
                {LonIndexer:min_lon}, method='nearest'))
            EasternLimit =float(NetCDF_data[LonIndexer].sel(
                {LonIndexer:max_lon}, method='nearest'))
            SouthernLimit = float(NetCDF_data[LatIndexer].sel(
                {MovingObj.LatIndexer:min_lat}, method='nearest'))
            NorthernLimit = float(NetCDF_data[LatIndexer].sel(
                {LatIndexer:max_lat}, method='nearest'))
        
        for term in stored_terms:
            # term_sliced = getattr(MovingObj,term).sel({MovingObj.TimeIndexer:t}
            #                     ).sel(**{MovingObj.LatIndexer:slice(
            #                     SouthernLimit,NorthernLimit),
            #                     MovingObj.LonIndexer: slice(
            #                             WesternLimit,EasternLimit)})
            term_sliced = getattr(MovingObj,term).sel(**{MovingObj.LatIndexer:slice(
                                SouthernLimit,NorthernLimit),
                                MovingObj.LonIndexer: slice(
                                        WesternLimit,EasternLimit)})
            dfDict[term][itime] = CalcAreaAverage(term_sliced,
                                ZonalAverage=True)
    
    print('saving results to nc file...')
    term_results = []
    for term in stored_terms:
        term_results.append(getattr(MovingObj,term))
    
    results_nc = []
    for term, da in zip(stored_terms,term_results):
        da =  da.metpy.dequantify()
        da.name = term
        results_nc.append(da)
    out_nc = xr.merge(results_nc)
    fname = ResultsSubDirectory+ outfile_name+'.nc'
    os.system('rm '+fname)
    out_nc.to_netcdf(fname, mode='w',engine="netcdf4")
    print(fname+' created')
    
    # Save CSV files with all the terms stored above
    for term in stored_terms:
        dfDict[term].to_csv(ResultsSubDirectory+term+'.csv')
        
    # Save system position as a csv file for replicability
    if not args.fixed:

        track = pd.DataFrame.from_dict(position)
        track = track.rename(columns={'central_lat':'Lat','central_lon':'Lon'})
        track.to_csv(ResultsSubDirectory+outfile_name+'_track',
                    index=False, sep=";")

        plot_track(track, FigsDirectory)
        plot_min_zeta_hgt(track, FigsDirectory)                          

    
def FixedAnalysis(FixedObj):
    
    timesteps = FixedObj.NetCDF_data[FixedObj.TimeIndexer]
    
    pres_levels = FixedObj.Temperature[FixedObj.LevelIndexer].\
        metpy.convert_units('hPa').values
        
    stored_terms = ['AdvHTemp','Sigma', 'Omega','dTdt','ResT', 
                   'Zeta', 'dZdt','AdvHZeta','AdvVZeta', 'vxBeta',
                  'ZetaDivH','fDivH', 'Tilting', 'ResZ'] 
    
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
                
        for term in stored_terms:
            term_sliced = getattr(FixedObj,term).sel({FixedObj.TimeIndexer:t})
            dfDict[term][itime] = CalcAreaAverage(term_sliced,
                                ZonalAverage=True)

    print('saving results to nc file...')
    term_results = []
    for term in stored_terms:
        term_results.append(getattr(FixedObj,term))
        
    results_nc = []
    for term, da in zip(stored_terms,term_results):
        da =  da.metpy.dequantify()
        da.name = term
        results_nc.append(da)
    out_nc = xr.merge(results_nc)
    fname = ResultsSubDirectory+ outfile_name+'.nc'
    out_nc.to_netcdf(fname, mode='w',engine="netcdf4")
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
    parser.add_argument("-g", "--gfs", default = False,
    action='store_true', help = "open multiple  GFS files at once using cfgrib\
 engine")
    parser.add_argument("-o", "--outname", default = False, type=str,
    help = "choose a name for saving results (default is\
 the same as infile)")
    
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
    infile = args.infile
    print('Opening file: ' + infile)
    if args.gfs:
        with dask.config.set(array={'slicing': {'split_large_chunks': True}}):
            NetCDF_data = convert_lon(
                xr.open_mfdataset(
                    infile,
                    engine='cfgrib',
                    parallel=True,
                    filter_by_keys={'typeOfLevel': 'isobaricInhPa'},
                    combine='nested',
                    concat_dim=TimeIndexer
                ),
                dfVars.loc['Longitude']['Variable']
            )
    else:
        with dask.config.set(array={'slicing': {'split_large_chunks': True}}):
            NetCDF_data = convert_lon(
                xr.open_dataset(infile),
                dfVars.loc['Longitude']['Variable']
            )
        
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

    # Computeterms with temporal dependency
    print('Computing zeta and temperature tendencies...')
    dTdt =  NetCDF_data[dfVars.loc['Air Temperature']['Variable']].differentiate(
                TimeIndexer,datetime_unit='s') * units('K/s')
    Zeta = vorticity(NetCDF_data[dfVars.loc['Eastward Wind Component']['Variable']],
                                  NetCDF_data[dfVars.loc['Northward Wind Component']['Variable']])          
    dZdt = Zeta.differentiate(TimeIndexer, datetime_unit='s') / units('s')
    print('Done.')
    
    # Directory where results will be stored
    ResultsMainDirectory = '../Results/'
    if args.outname:
        outfile_name = args.outname
    # Append data method to outfile name
    else:
        outfile_name = ''.join(infile.split('/')[-1].split('.nc'))+'_'+method
    # Each dataset of results have its own directory, allowing to store results
    # from more than one experiment at each time
    ResultsSubDirectory = ResultsMainDirectory+'/'+outfile_name+'/'
    # Subdirectory for saving figures
    FigsDirectory = ResultsSubDirectory+'Figures/'
    # Check if the required directories exists. If not, creates them
    check_create_folder(ResultsMainDirectory)
    check_create_folder(ResultsSubDirectory)
    check_create_folder(FigsDirectory)
    # backup of box_limits and track file for reproductability
    if args.fixed:
        os.system('cp ../inputs/box_limits '+ResultsSubDirectory)
    elif args.track:
        os.system('cp ../inputs/track '+ResultsSubDirectory)

     # Run the program
    cyclone_thermodynamics(NetCDF_data,dfVars, dTdt, dZdt, args)
    print("--- %s seconds for running the program ---" % (
            time.time() - start_time))
