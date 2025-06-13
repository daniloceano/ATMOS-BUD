# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    calculations.py                                    :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2024/02/16 16:42:55 by daniloceano       #+#    #+#              #
#    Updated: 2025/06/13 10:30:11 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import numpy as np
import pandas as pd
import xarray as xr

from metpy.calc import wind_speed

from src.data_object import DataObject    
from src.select_domain import get_domain_limits
from src.visualization import plot_fixed_domain
from src.utils import handle_track_file
from src.utils import find_extremum_coordinates
from src.utils import get_domain_extreme_values
from src.output_management import save_output_track
from src.output_management import save_results_csv
from src.output_management import save_results_netcdf
from src.visualization import hovmoller_mean_zeta

def CalcZonalAverage(VariableData):
    """
    Computates variable zonal average of some variable, for all z levels and time steps.
    
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
    xlength = VariableData['rlons'][-1] - VariableData['rlons'][0]
    return VariableData.integrate("rlons") / xlength

def CalcAreaAverage(VariableData, ZonalAverage=False):
    """
    Computates the Area Average of a function.
    
    The default is to computate the zonal average and then a meridional average.
    If the input data is already some sort of zonal quantity (average or not),
    simply set longitude_indexer to None
    
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
    return (ZA * ZA["coslats"]).integrate("rlats") / ylength

def perform_calculations(input_data, namelist_df, dTdt, dZdt, dQdt, args, app_logger, *outputs):
    """
    Performs meteorological calculations on input data and generates results and figures.

    Parameters:
    - input_data (xr.Dataset): Dataset containing meteorological data from a NetCDF file.
    - namelist_df (pd.DataFrame): DataFrame mapping variable names.
    - dTdt, dZdt, dQdt: Data arrays representing temperature, vortcity and moisture tendencies, respectively.
    - args: Command-line arguments or parameters specifying calculation options.
    - app_logger (Logger): Logger for outputting information and error messages.
    - outputs (tuple): A tuple containing paths for results and figures directories, and the output file name.

    The function processes each time step of the input data to calculate various meteorological terms,
    stores results in CSV files, and generates figures for analysis. It handles domain selection based on tracks,
    calculates area averages for specified terms, and plots diagnostic figures.
    """

    # Unpack outputs
    results_subdirectory, figures_subdirectory, outfile_name = outputs
    app_logger.info(f'Directory where results will be stored: {results_subdirectory}')
    app_logger.info(f'Directory where figures will be stored: {figures_subdirectory}')
    app_logger.info(f'Name of the output file with results: {outfile_name}')
        
    # Dictionary to save DataArray results and transform into nc later
    results_nc = {}
    
    # Indexers
    longitude_indexer = namelist_df.loc['Longitude']['Variable']
    latitude_indexer = namelist_df.loc['Latitude']['Variable']
    time_indexer = namelist_df.loc['Time']['Variable']
    vertical_level_indexer = namelist_df.loc['Vertical Level']['Variable']

    # Data time steps and pressure levels
    timesteps = pd.to_datetime(input_data[time_indexer].values)
    pres_levels = input_data[vertical_level_indexer].values

    # Load track file, if specified, and check for errors
    if args.track:
        track = handle_track_file(input_data, timesteps, longitude_indexer, latitude_indexer, app_logger)
    
    # Create a dictionary for saving area averages of each term
    stored_terms = ['AdvHTemp','AdvVTemp', 'Sigma','Omega','dTdt','ResT', 
                    'Zeta', 'dZdt','AdvHZeta','AdvVZeta', 'vxBeta',
                    'ZetaDivH','fDivH', 'Tilting', 'ResZ',
                    'dQdt', 'dQdt_integrated', 'divQ', 'divQ_integrated', 'WaterBudgetResidual', 'WaterBudgetResidual_integrated'] 
    
    # Create a dataframe for each term
    results_df_dictionary = {}
    for term in stored_terms:
        if term == 'ResQ':
            results_df_dictionary[term] = {}
        else:
            results_df_dictionary[term] = pd.DataFrame(columns=[str(t) for t in timesteps],
                                            index=[float(i) for i in pres_levels])
        results_nc[term] = []

    # Get choosen vertical level from args and convert to Pa
    plevel_Pa = int(args.level) * 100
        
    # Dictionary for saving track attributes for each timestep
    output_track_attributes = {}
    results_keys = ['time', 'central_lat', 'central_lon', 'length', 'width',
            f'{args.track_vorticity}_zeta_{args.level}', f'{args.track_geopotential}_hgt_{args.level}', f'max_wind_{args.level}']
    for key in results_keys:
        output_track_attributes[key] =  []

    for time_step in timesteps:
        
        itime = str(time_step)
        datestr = pd.to_datetime(itime).strftime('%Y-%m-%d %HZ')
        datestr2 = pd.to_datetime(itime).strftime('%Y%m%d%H00')
        app_logger.info(f'⏳ Processing time step: {datestr}')

        # Convert the timezone-aware timestamp 't' to a timezone-naive timestamp (still in UTC)
        t_naive = time_step.tz_localize(None)
        
        # Create a DataObject for the current time step, where all the calculations are performed
        MovingObj = DataObject(
            input_data.sel({time_indexer:t_naive}),
            dTdt=dTdt.sel({time_indexer:t_naive}),
            dZdt=dZdt.sel({time_indexer:t_naive}),
            dQdt=dQdt.sel({time_indexer:t_naive}),
            namelist_df=namelist_df,
            app_logger=app_logger,
            args=args)
    

        # Get variables at choosen pressure level for the current time step
        iu_plevel = MovingObj.u.sel({vertical_level_indexer:plevel_Pa})
        iv_plevel = MovingObj.v.sel({vertical_level_indexer:plevel_Pa})
        ight_plevel = MovingObj.GeopotHeight.sel({vertical_level_indexer:plevel_Pa})
        iwspd_plevel = wind_speed(iu_plevel, iv_plevel)
        zeta = MovingObj.Zeta.sel({vertical_level_indexer:plevel_Pa}).metpy.dequantify()
        lat, lon = input_data[MovingObj.latitude_indexer], input_data[MovingObj.longitude_indexer]

        # Get domain limits for the current time step
        variables_at_plevel  = [iu_plevel, iv_plevel, zeta, ight_plevel, lat, lon, itime]
        if 'track' in locals() or 'track' in globals():
            current_domain_limits = get_domain_limits(args, *variables_at_plevel, track=track)
        else:
            current_domain_limits = get_domain_limits(args, *variables_at_plevel)
        min_lat, max_lat = current_domain_limits['min_lat'], current_domain_limits['max_lat']
        min_lon, max_lon = current_domain_limits['min_lon'], current_domain_limits['max_lon']
        central_lat, central_lon = current_domain_limits['central_lat'], current_domain_limits['central_lon']
            
        # Slice variables for the domain limits
        izeta_plevel_slice = zeta.sel({latitude_indexer:slice(min_lat, max_lat), longitude_indexer:slice(min_lon, max_lon)})
        ight_plevel_slice = ight_plevel.sel({latitude_indexer:slice(min_lat, max_lat), longitude_indexer:slice(min_lon, max_lon)})
        iwspd_plevel_slice = iwspd_plevel.sel({latitude_indexer:slice(min_lat, max_lat), longitude_indexer:slice(min_lon, max_lon)})

        # Get the extreme values at choosen pressure level for the current time step
        slices_plevel = [izeta_plevel_slice, ight_plevel_slice, iwspd_plevel_slice]
        if 'track' in locals() or 'track' in globals():
            min_max_zeta, min_max_hgt, max_wind = get_domain_extreme_values(itime, args, slices_plevel, track)
        else:
            min_max_zeta, min_max_hgt, max_wind = get_domain_extreme_values(itime, args, slices_plevel)

        # Find position of the extremes
        lat_slice, lon_slice = izeta_plevel_slice[latitude_indexer], izeta_plevel_slice[longitude_indexer]
        min_max_zeta_lat, min_max_zeta_lon = find_extremum_coordinates(izeta_plevel_slice, lat_slice, lon_slice, f'{args.track_vorticity}_zeta_{args.level}', args)
        min_max_hgt_lat, min_max_hgt_lon = find_extremum_coordinates(ight_plevel_slice, lat_slice, lon_slice, f'{args.track_geopotential}_hgt_{args.level}', args)
        max_wind_lat, max_wind_lon = find_extremum_coordinates(iwspd_plevel_slice, lat_slice, lon_slice, 'max_wind', args)

        # Print results on screen
        if args.fixed:
            app_logger.info(f"Storing results for: {datestr}")
            WesternLimit = min_lon
            EasternLimit = max_lon
            NorthernLimit = max_lat
            SouthernLimit = min_lat

        else:
            # Store system position and attributes
            length, width = max_lat - min_lat, max_lon - min_lon
            domain_attributes = [datestr, central_lat, central_lon, length, width, min_max_zeta, min_max_hgt, max_wind]
            for key,val in zip(results_keys, domain_attributes):
                output_track_attributes[key].append(val)
            
            WesternLimit = float(input_data[longitude_indexer].sel({longitude_indexer: min_lon}, method='nearest'))
            EasternLimit =float(input_data[longitude_indexer].sel({longitude_indexer: max_lon}, method='nearest'))
            SouthernLimit = float(input_data[latitude_indexer].sel({latitude_indexer: min_lat}, method='nearest'))
            NorthernLimit = float(input_data[latitude_indexer].sel({latitude_indexer: max_lat}, method='nearest'))

            app_logger.info(
            f"📍 central lat: {central_lat}, central lon: {central_lon}, "
            f"size: {length} x {width}, "
            f"lon range: {min_lon} to {max_lon}, "
            f"lat range: {min_lat} to {max_lat}"
            )

            app_logger.info(
            f"📍 {int(args.level)} hPa diagnostics --> "
            f"{args.track_vorticity} ζ: {min_max_zeta:.2e}, "
            f"{args.track_geopotential} geopotential height: {min_max_hgt:.0f}, "
            f"max wind speed: {max_wind:.2f}"
            )
        
            app_logger.info(
            f"📍 {int(args.level)} hPa positions (lat/lon) --> "
            f"{args.track_vorticity} ζ: {min_max_zeta_lat:.2f}, {min_max_zeta_lon:.2f}, "
            f"{args.track_geopotential} geopotential height: {min_max_hgt_lat:.2f}, {min_max_hgt_lon:.2f}, "
            f"max wind speed: {max_wind_lat:.2f}, {max_wind_lon:.2f}")
        
        app_logger.info(f"📊 Storing results for: {datestr}")
        for term in stored_terms:
            term_sliced = getattr(MovingObj,term).sel(
                **{
                MovingObj.latitude_indexer:slice(SouthernLimit,NorthernLimit),
                MovingObj.longitude_indexer: slice(WesternLimit,EasternLimit)
                }
            )
            if term in ['divQ_integrated', 'dQdt_integrated', 'WaterBudgetResidual_integrated']:
                results_df_dictionary[term][itime] = float(CalcAreaAverage(term_sliced, ZonalAverage=True))
            else:
                results_df_dictionary[term][itime] = CalcAreaAverage(term_sliced, ZonalAverage=True)

        # Save figure with box used for computations
        dict_for_plot = {
            f'{args.track_vorticity}_zeta_{args.level}': {
                'latitude': min_max_zeta_lat,
                'longitude': min_max_zeta_lon,
                'data': zeta
            },
            f'{args.track_geopotential}_hgt_{args.level}': {
                'latitude': min_max_hgt_lat,
                'longitude': min_max_hgt_lon,
                'data': ight_plevel
            },
            'max_wind': {
                'latitude': max_wind_lat,
                'longitude': max_wind_lon,
                'data': iwspd_plevel
            },
            'lat': lat,
            'lon': lon,
        }

        app_logger.info(f"📊 Saving domain plot for {datestr2}")
        plot_fixed_domain(current_domain_limits, dict_for_plot, args, results_subdirectory, datestr2, app_logger)

        if args.save_nc_file:
            app_logger.debug(f"💾 Storing results for the NetCDF file for {datestr}")
            # Initializing out_nc on the first time step
            if time_step == timesteps[0]:
                # Create the first Dataset with the terms' variables
                results_ds = [getattr(MovingObj, term).metpy.dequantify().assign_attrs(units='').rename(term) for term in stored_terms]
                out_nc = xr.merge(results_ds)  # Merge the initial variables
                app_logger.debug(f"Created Xarray dataset for {datestr}")

            else:
                # For subsequent time steps, create results_ds and concatenate along the 'time' dimension
                results_ds = [getattr(MovingObj, term).metpy.dequantify().assign_attrs(units='').rename(term) for term in stored_terms]
                # Concatenate results_ds along the 'time' dimension
                out_nc = xr.concat([out_nc, xr.merge(results_ds)], dim='initial_time0_hours')  # Concatenate along time
                app_logger.debug(f"Appended Xarray dataset for {datestr}")


    # Save system position as a csv file for replicability
    if not args.fixed:
        save_output_track(output_track_attributes, args, results_subdirectory, figures_subdirectory, outfile_name, app_logger)
        hovmoller_mean_zeta(results_df_dictionary['Zeta'], figures_subdirectory, app_logger)
    
    app_logger.info(f"💾 Saving CSV results for {outfile_name}")
    save_results_csv(results_df_dictionary, results_subdirectory, app_logger)
    if args.save_nc_file:
        save_results_netcdf(out_nc, results_subdirectory, outfile_name, app_logger)
        

        