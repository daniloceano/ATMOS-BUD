# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    __init__.py                                      :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2025/08/15 by daniloceano              #+#    #+#              #
#    Updated: 2025/08/15 by daniloceano             ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
ATMOS-BUD: Atmospheric Budget Analysis Tool

A comprehensive Python package for analyzing heat, vorticity, and humidity budgets
of limited regions in the atmosphere using reanalysis data.

Version: 0.1.2
"""

__version__ = "0.1.2"
__author__ = "Danilo Couto de Souza"
__email__ = "danilo.oceano@gmail.com"
__license__ = "GPL-3.0"

# Main modules
__all__ = [
    "calculations",
    "cli_interface", 
    "data_handling",
    "data_object",
    "get_era5_data",
    "output_management",
    "select_domain",
    "utils",
    "visualization"
]
