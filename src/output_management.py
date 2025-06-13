# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    output_management.py                               :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2024/02/16 16:13:36 by daniloceano       #+#    #+#              #
#    Updated: 2025/06/13 10:34:17 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import xarray as xr
import pandas as pd
from src.visualization import plot_track
from src.visualization import plot_min_max_zeta_hgt

def manage_output(args, infile, method): 
    """
    Manages the output directories and file names for results and figures.

    Parameters:
    - args: Parsed command-line arguments.
    - infile: The input file name.
    - method: The method used for calculations.

    Returns:
    - Tuple containing paths for results subdirectory, figures subdirectory, and the output file name.
    """   

    results_main_directory = './Results/'
    outfile_name = args.outname if args.outname else ''.join(infile.split('/')[-1].split('.nc')) + '_' + method
    results_subdirectory = os.path.join(results_main_directory, outfile_name)
    figures_subdirectory = os.path.join(results_subdirectory, 'Figures')

    # Create subdirectories for heat, moisture, and vorticity terms
    heat_subdirectory = os.path.join(results_subdirectory, 'heat_terms')
    moisture_subdirectory = os.path.join(results_subdirectory, 'moisture_terms')
    vorticity_subdirectory = os.path.join(results_subdirectory, 'vorticity_terms')

    # Create necessary directories
    for directory in [results_main_directory, results_subdirectory, figures_subdirectory,
                      heat_subdirectory, moisture_subdirectory, vorticity_subdirectory]:
        try:
            os.makedirs(directory, exist_ok=True)
        except Exception as e:
            print(f"Error creating directory {directory}: {e}")
            raise

   # Backup box_limits and track file for reproducibility
    if args.fixed:
        backup_file('inputs/box_limits', results_subdirectory)
    elif args.track:
        backup_file('inputs/track', results_subdirectory)

    return results_subdirectory, figures_subdirectory, outfile_name

def backup_file(source_path, destination_directory):
    """
    Copies a file to a destination directory for backup purposes.

    Parameters:
    - source_path: Path to the source file.
    - destination_directory: Directory where the file will be copied.
    """
    try:
        destination_path = os.path.join(destination_directory, os.path.basename(source_path))
        if not os.path.exists(destination_path):
            os.system(f'cp {source_path} {destination_path}')
    except Exception as e:
        print(f"Error backing up file from {source_path} to {destination_directory}: {e}")

def save_results_netcdf(out_nc, results_subdirectory, outfile_name, app_logger):
    """
    Saves the calculation results to a NetCDF file.

    Parameters:
    - MovingObj: The object containing calculation results.
    - stored_terms: List of terms to be saved.
    - results_subdirectory: Directory where results will be saved.
    - outfile_name: Base name for the output file.
    - app_logger: Logger for outputting information and error messages.
    """
    try:
        fname = os.path.join(results_subdirectory, f'{outfile_name}.nc')
        out_nc.to_netcdf(fname, mode='w')
        app_logger.info(f'💾 {fname} created successfully')
    except Exception as e:
        app_logger.error(f'❌ Error saving NetCDF file: {e}')

def save_results_csv(results_df_dictionary, results_subdirectory, app_logger):
    """
    Saves the calculation results to CSV files for each term.

    Parameters:
    - results_df_dictionary: Dictionary of pandas DataFrames for each term.
    - results_subdirectory: Directory where results will be saved.
    - app_logger: Logger for outputting information and error messages.
    """
    try:
        # Save CSV files for each term
        for term, df in results_df_dictionary.items():

            if term == 'ResQ':
                df = pd.DataFrame.from_dict(df, orient='index')

            if term in ['AdvHTemp','AdvVTemp', 'Sigma','Omega','dTdt','ResT']:
                budget_results_subdirectory = os.path.join(results_subdirectory, 'heat_terms')
            elif term in ['Zeta', 'dZdt','AdvHZeta','AdvVZeta', 'vxBeta',
                    'ZetaDivH','fDivH', 'Tilting', 'ResZ']:
                budget_results_subdirectory = os.path.join(results_subdirectory, 'vorticity_terms')
            elif term in ['dQdt', 'dQdt_integrated', 'divQ', 'divQ_integrated', 'WaterBudgetResidual', 'WaterBudgetResidual_integrated']:
                budget_results_subdirectory = os.path.join(results_subdirectory, 'moisture_terms')
            else:
                app_logger.error(f"Unknown term: {term}. Saving to main directory.")
                budget_results_subdirectory = results_subdirectory

            csv_file_name = os.path.join(budget_results_subdirectory, f'{term}.csv')
            df.to_csv(csv_file_name)
            app_logger.info(f'💾 {csv_file_name} created successfully')
    except Exception as e:
        app_logger.error(f'❌ Error saving CSV files: {e}')

def save_output_track(output_track_attributes, args, results_subdirectory, figures_subdirectory, outfile_name, app_logger):
    """
    Saves track data to a CSV file and generates track and min/max zeta height plots.

    Parameters:
    - output_track_attributes: Dictionary containing track attributes.
    - args: Command-line arguments or parameters specifying calculation options.
    - results_subdirectory: Directory where the track CSV will be saved.
    - figures_subdirectory: Directory where figures will be saved.
    - outfile_name: Base name for the output track file.
    - app_logger: Logger for outputting information and error messages.
    """
    try:
        track = pd.DataFrame.from_dict(output_track_attributes)
        track = track.rename(columns={'central_lat': 'Lat', 'central_lon': 'Lon'})
        track_file_path = os.path.join(results_subdirectory, f'{outfile_name}_track.csv')
        track.to_csv(track_file_path, index=False, sep=";")
        app_logger.info(f'💾 {track_file_path} created successfully')

        plot_track(track, args, figures_subdirectory, app_logger)
        plot_min_max_zeta_hgt(track.set_index('time'), args, figures_subdirectory, app_logger)
    except Exception as e:
        app_logger.error(f'❌ Error saving output track: {e}')