# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    atmos_bud.py                                       :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/06/18 10:09:03 by daniloceano       #+#    #+#              #
#    Updated: 2024/04/12 08:06:16 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import time
import sys
import pandas as pd
from metpy.units import units
from metpy.calc import vorticity

from src.utils import initialize_logging
from src.calculations import perform_calculations
from src.cli_interface import parse_arguments
from src.data_handling import load_data
from src.data_handling import preprocess_data
from src.output_management import manage_output
                          
    
def main():
    # Start timer
    start_time = time.time()

    if len(sys.argv) > 1:
        # Parse command line arguments
        args = parse_arguments()
    
    else:
        # For debuging:
        print('Debug mode')
        # debug_args = ['samples/sample1_ERA5.nc', '-c', '-v']
        debug_args = ['~/Documents/Programs_and_scripts/data_etc/netCDF_files/akara_teste.nc', '-f', '-v']
        args = parse_arguments(debug_args)

    # Set method
    if args.fixed:
        method = 'fixed'
    elif args.track:
        method = 'track'
    elif args.choose:
        method = 'choose'

    # Manage output directories and output files
    infile = args.infile
    outputs = manage_output(args, infile, method)

    # Initialize logger
    app_logger = initialize_logging(outputs[0], args)
    
    # Open namelist
    namelist = 'inputs/namelist'
    namelist_df = pd.read_csv(namelist, sep=';', index_col=0, header=0)
    
    # Data indexers
    longitude_indexer = namelist_df.loc['Longitude']['Variable']
    time_indexer = namelist_df.loc['Time']['Variable']
    
    # Open file
    input_data = load_data(infile, longitude_indexer, args, app_logger)
    input_data = preprocess_data(input_data, namelist_df, args, app_logger)

    # Computeterms with temporal dependency
    app_logger.info('Computing zeta and temperature tendencies...')
    dTdt =  input_data[namelist_df.loc['Air Temperature']['Variable']].differentiate(time_indexer, datetime_unit='s') * units('K/s')
    Zeta = vorticity(input_data[namelist_df.loc['Eastward Wind Component']['Variable']],
                     input_data[namelist_df.loc['Northward Wind Component']['Variable']],
                     )          
    dZdt = Zeta.differentiate(time_indexer, datetime_unit='s') / units('s')
    dQdt = input_data[namelist_df.loc['Specific Humidity']['Variable']].differentiate(time_indexer, datetime_unit='s') * units('kg/kg/s')
    app_logger.info('Done.')

    # Run the program
    perform_calculations(input_data, namelist_df, dTdt, dZdt, dQdt, args, app_logger, *outputs)
    app_logger.info("--- %s seconds for running the program ---" % (time.time() - start_time))

if __name__ == "__main__":
    main()
