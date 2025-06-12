# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    cli_interface.py                                   :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2024/02/16 15:09:06 by daniloceano       #+#    #+#              #
#    Updated: 2025/06/12 15:35:23 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import argparse

def parse_arguments(custom_args=None):
    """
    Parses command-line arguments.

    Args:
        custom_args (list, optional): A list of argument strings for debugging.

    Returns:
        args (Namespace): Namespace object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(description="""
        Program for analyzing the heat, vorticity and humidity budgets of a cyclone.
        The program offers three distinct frameworks for setting the analysis domain:
        1) Fixed framework: The energetics are computed for a domain specified by the 'inputs/box_limits' file.
        2) Moving framework: The domain follows the system, defined using the 'inputs/track' file.
        3) Interactive framework: For each time step, the user can interactively choose the domain by clicking on the map displayed.
        All frameworks require an auxiliary namelist file (inputs/namelist), specifying the variable names used.
        Results are stored as CSV files in the 'ATMOS-BUD_Results' directory (../ATMOS-BUD_Results), and figures for visualization are also created.
    """, formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("infile", 
                        help="""Input .nc file with temperature, specific humidity and meridional, zonal, and vertical components of the wind, in pressure levels.""")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fixed", action='store_true',
                       help="Compute the energetics for a Fixed domain specified by the 'inputs/box_limits' file.")
    group.add_argument("-t", "--track", action='store_true',
                       help="""Define the box using a track file specified by the 'inputs/track' file. The track indicates the central point of the system and an arbitrary box of 15°x15° is constructed.""")
    group.add_argument("-c", "--choose", action='store_true',
                       help="For each time step, the user can choose the domain by clicking on the screen.")
    parser.add_argument("-l", "--level", type=int, default=850,
                   help="Pressure level (in hPa) to show when choosing the domain (default: 850)")
    parser.add_argument("--track_vorticity", choices=["min", "max"], default="min",
                    help="Whether to track the minimum (default) or maximum vorticity.")
    parser.add_argument("--track_geopotential", choices=["min", "max"], default="min",
                    help="Whether to track the minimum (default) or maximum geopotential height.")
    parser.add_argument("-g", "--gfs", action='store_true',
                        help="Open multiple GFS files at once using cfgrib engine.")
    parser.add_argument("-o", "--outname", default=False, type=str,
                        help="""Choose a name for saving results (default is the same as infile).""")
    parser.add_argument("--save_nc_file",  default=True, type=str,
                        help="Wether to save the NetCDF file containing the results for each term, for the entire domain.")
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="Show debug messages while running.")
    
    if custom_args is not None:
        args = parser.parse_args(custom_args)
        
    else:
        args = parser.parse_args()
    
    return args