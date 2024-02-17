# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    output_management.py                               :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2024/02/16 16:13:36 by daniloceano       #+#    #+#              #
#    Updated: 2024/02/16 23:34:50 by daniloceano      ###   ########.fr        #
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

    results_main_directory = '../ATMOS-BUD_Results/'
    outfile_name = args.outname if args.outname else ''.join(infile.split('/')[-1].split('.nc')) + '_' + method
    results_subdirectory = os.path.join(results_main_directory, outfile_name)
    figures_subdirectory = os.path.join(results_subdirectory, 'Figures')

    # Create necessary directories
    for directory in [results_main_directory, results_subdirectory, figures_subdirectory]:
        try:
            os.makedirs(directory, exist_ok=True)
        except Exception as e:
            print(f"Error creating directory {directory}: {e}")
            raise

   # Backup box_limits and track file for reproducibility
    if args.fixed:
        backup_file('../inputs/box_limits', results_subdirectory)
    elif args.track:
        backup_file('../inputs/track', results_subdirectory)

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

def save_results(MovingObj, results_df_dictionary, stored_terms, results_subdirectory, outfile_name, app_logger):
    """
    Saves the calculation results to a NetCDF file and CSV files for each term.

    Parameters:
    - MovingObj: The object containing calculation results.
    - results_df_dictionary: Dictionary of pandas DataFrames for each term.
    - stored_terms: List of terms to be saved.
    - results_subdirectory: Directory where results will be saved.
    - outfile_name: Base name for output files.
    - app_logger: Logger for outputting information and error messages.
    """
    try:
        app_logger.info('Saving results to NetCDF file...')
        term_results = [getattr(MovingObj, term).metpy.dequantify().assign_attrs(units='').rename(term) for term in stored_terms]
        out_nc = xr.merge(term_results)
        fname = os.path.join(results_subdirectory, f'{outfile_name}.nc')
        out_nc.to_netcdf(fname, mode='w')
        app_logger.info(f'{fname} created')

        # Save CSV files for each term
        for term, df in results_df_dictionary.items():
            df.to_csv(os.path.join(results_subdirectory, f'{term}.csv'))
    except Exception as e:
        app_logger.error(f"Error saving results: {e}")

def save_output_track(output_track_attributes, results_subdirectory, figures_subdirectory, outfile_name, app_logger):
    """
    Saves track data to a CSV file and generates track and min/max zeta height plots.

    Parameters:
    - output_track_attributes: Dictionary containing track attributes.
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

        plot_track(track, figures_subdirectory)
        plot_min_max_zeta_hgt(track.set_index('time'), figures_subdirectory, app_logger)
    except Exception as e:
        app_logger.error(f"Error saving output track: {e}")