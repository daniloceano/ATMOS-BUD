# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    data_object.py                                     :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2024/02/16 21:32:34 by daniloceano       #+#    #+#              #
#    Updated: 2025/04/22 09:08:32 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import xarray as xr
import numpy as np
import pandas as pd
import argparse
import logging

from metpy.units import units
from metpy.constants import Cp_d
from metpy.constants import Re
from metpy.calc import potential_temperature
from metpy.calc import vorticity
from metpy.calc import divergence
from metpy.calc import coriolis_parameter
from metpy.constants import g

class DataObject:
    """
    A class for processing meteorological data and computing terms of the Quasi-Geostrophic Equation.
    It calculates thermodynamic and vorticity terms, excluding Adiabatic Heating Term (Q), which is estimated as a residual.
    Note that: Q = J * Cp_d

    Attributes:
        input_data (xr.Dataset): The input dataset containing meteorological variables.
        dTdt (xr.DataArray): Data array representing the temperature tendency.
        dZdt (xr.DataArray): Data array representing the geopotential height tendency.
        dQdt (xr.DataArray): Data array representing the moisture tendency.
        namelist_df (pd.DataFrame): DataFrame mapping variable names and units.
        args (argparse.Namespace): Parsed command-line arguments.
        app_logger (logging.Logger): Logger for outputting information and error messages.
    """
    def __init__(self, input_data: xr.Dataset, dTdt: xr.DataArray, dZdt: xr.DataArray, dQdt: xr.DataArray,
                 namelist_df: pd.DataFrame, args: argparse.Namespace, app_logger: logging.Logger):
        try:
            self.extract_variables(input_data, namelist_df, args)
            self.calculate_thermodynamic_terms(dTdt)
            self.calculate_vorticity_terms(dZdt)
            self.calculate_water_budget_terms(dQdt)
        except KeyError as e:
            app_logger.error(f'Missing key in namelist DataFrame or input data: {e}')
            raise
        except Exception as e:
            app_logger.error(f'Unexpected error during DataObject initialization: {e}')
            raise

    def extract_variables(self, input_data, namelist_df, args):
        """
        Extracts variables from the input dataset using the mapping provided in the namelist DataFrame.
        Converts units, defines coordinates, and prepares variables for later computations.

        Parameters:
        - input_data (xr.Dataset): Input dataset containing raw meteorological variables.
        - namelist_df (pd.DataFrame): DataFrame mapping standardized variable names to dataset variables and units.
        - args (argparse.Namespace): Command-line arguments with configuration options.
        """
        self.input_data = input_data.load() if args.gfs else input_data

        # Extract indexers using namelist DataFrame
        self.longitude_indexer = namelist_df.loc['Longitude']['Variable']
        self.latitude_indexer = namelist_df.loc['Latitude']['Variable']
        self.time_indexer = namelist_df.loc['Time']['Variable']
        self.vertical_level_indexer = namelist_df.loc['Vertical Level']['Variable']

        # Convert units and store variables
        self.Temperature = self.convert_units('Air Temperature', namelist_df)
        self.u = self.convert_units('Eastward Wind Component', namelist_df)
        self.v = self.convert_units('Northward Wind Component', namelist_df)
        self.Omega = self.convert_units('Omega Velocity', namelist_df)
        self.GeopotHeight = self.calculate_geopotential_height(namelist_df)
        self.SpecificHumidity = self.convert_units('Specific Humidity', namelist_df)
        self.Pressure = self.input_data[self.vertical_level_indexer] * units('Pa')

        # Additional calculations for dx, dy, f, etc.
        self.calculate_additional_properties()

    def convert_units(self, var_name, namelist_df):
        """
        Converts the units of a variable from the dataset based on definitions in the namelist DataFrame.

        Parameters:
        - var_name (str): Standardized name of the variable to extract and convert.
        - namelist_df (pd.DataFrame): DataFrame with variable mappings and units.

        Returns:
        - xarray.DataArray: The variable with appropriate physical units.
        """
        return self.input_data[namelist_df.loc[var_name]['Variable']] * units(namelist_df.loc[var_name]['Units']).to(units.parse_units(namelist_df.loc[var_name]['Units']))

    def calculate_geopotential_height(self, namelist_df):
        """
        Retrieves geopotential height if available, otherwise calculates it from geopotential.

        Parameters:
        - namelist_df (pd.DataFrame): DataFrame containing variable mapping and unit definitions.

        Returns:
        - xarray.DataArray: Geopotential height (in meters).
        """
        if 'Geopotential Height' in namelist_df.index:
            return self.convert_units('Geopotential Height', namelist_df)
        else:
            geopotential = self.convert_units('Geopotential', namelist_df)
            return geopotential / g

    def calculate_additional_properties(self):
        """
        Computes grid spacing (dx, dy) in meters and Coriolis parameter (f) for each latitude.
        Also stores coordinate variables for future use.

        Sets:
        - self.dx, self.dy (xarray.DataArray): Grid spacing in x and y directions.
        - self.f (xarray.DataArray): Coriolis parameter (1/s).
        """
        self.lons = self.input_data[self.longitude_indexer]
        self.lats = self.input_data[self.latitude_indexer]
        self.cos_lats = np.cos(np.deg2rad(self.lats))
        self.dx = np.deg2rad(self.lons.differentiate(self.longitude_indexer)) * self.cos_lats * Re
        self.dy = np.deg2rad(self.lats.differentiate(self.latitude_indexer)) * Re
        self.f = coriolis_parameter(self.lats)

    def calculate_thermodynamic_terms(self, dTdt):
        """
        Calculates the atmospheric thermodynamic budget terms:
        - Horizontal and vertical advection of temperature
        - Sigma term related to adiabatic processes
        - Residual of the thermodynamic equation
        - Estimated diabatic heating (Q) from the residual

        The thermodynamic budget equation is:
            dT/dt = - (u · ∇T) - ω * ∂T/∂p + (θ/T) * ∂θ/∂p * ω + Q/Cp_d

        Where the residual is used to estimate:
            Q = Cp_d * ResT

        Parameters:
        - dTdt (xarray.DataArray): Time derivative of temperature (K/s)

        Sets:
        - self.dTdt: input field
        - self.Theta: potential temperature (K)
        - self.AdvHTemp: horizontal advection of temperature (K/s)
        - self.AdvVTemp: vertical advection of temperature (K/s)
        - self.Sigma: static stability term (K/Pa)
        - self.ResT: residual of the thermodynamic equation (K/s)
        - self.AdiabaticHeating: estimated diabatic heating (W/kg)
        """
        self.Theta = potential_temperature(self.Pressure, self.Temperature)
        self.dTdt = dTdt
        self.AdvHTemp = self.calculate_horizontal_advection(self.Temperature)
        self.AdvVTemp = -1 * (self.Temperature.differentiate(self.vertical_level_indexer) * self.Omega) / (1 * units('Pa'))
        self.Sigma = (-1 * (self.Temperature / self.Theta) * (self.Theta.differentiate(self.vertical_level_indexer) / units('Pa'))).metpy.convert_units('K / Pa')
        self.ResT =  self.dTdt - self.AdvHTemp - (self.Sigma * self.Omega)
        self.AdiabaticHeating = self.ResT * Cp_d

    def calculate_vorticity_terms(self, dZdt):
        """
        Calculates the atmospheric vorticity budget terms:
        - Horizontal and vertical advection of relative vorticity
        - Meridional advection of planetary vorticity (β effect)
        - Divergence-related terms (ζ·divV and f·divV)
        - Tilting term
        - Residual of the vorticity equation

        The vorticity budget equation is:
            dζ/dt = - (u · ∇ζ) - ω * ∂ζ/∂p - v * β - ζ * div(V) - f * div(V) + Tilting

        The residual is calculated as:
            ResZ = dZdt - [AdvHZeta + AdvVZeta + vxBeta + ZetaDivH + fDivH + Tilting]

        Parameters:
        - dZdt (xarray.DataArray): Time derivative of relative vorticity (1/s²)

        Sets:
        - self.Zeta: relative vorticity (1/s)
        - self.dZdt: input field
        - self.AdvHZeta: horizontal advection of vorticity (1/s²)
        - self.AdvVZeta: vertical advection of vorticity (1/s²)
        - self.Beta: meridional gradient of Coriolis parameter (1/m/s)
        - self.vxBeta: meridional advection of planetary vorticity (1/s²)
        - self.DivH: horizontal divergence of the wind (1/s)
        - self.ZetaDivH: term ζ·div(V) (1/s²)
        - self.fDivH: term f·div(V) (1/s²)
        - self.Tilting: tilting term (1/s²)
        - self.ResZ: residual of the vorticity budget (1/s²)
        """
        self.Zeta = vorticity(self.u, self.v)
        self.dZdt = dZdt
        self.AdvHZeta = self.calculate_horizontal_advection(self.Zeta)
        self.AdvVZeta = - 1 * (self.Zeta.differentiate(self.vertical_level_indexer) * self.Omega) / (1 * units('Pa'))
        self.Beta = (self.f).differentiate(self.latitude_indexer) / self.dy
        self.vxBeta = - 1 * (self.v * self.Beta)
        self.DivH = divergence(self.u, self.v)
        self.fDivH = - 1 * self.f * self.DivH
        self.ZetaDivH = - 1 * (self.Zeta * self.DivH)
        self.Tilting = self.tilting_term()
        self.SumVorticity = self.AdvHZeta + self.AdvVZeta + self.vxBeta + self.ZetaDivH + self.fDivH + self.Tilting
        self.ResZ = self.dZdt - self.SumVorticity

    def calculate_water_budget_terms(self, dQdt):
        """
        Calculates the atmospheric water budget terms:
        - Time derivative of specific humidity (dQdt)
        - Horizontal divergence of moisture flux (divQ)
        - Residual of the water vapor budget (WaterBudgetResidual)

        The water vapor budget equation is:
            dW/dt + div(Q) = P - E

        Parameters:
        - dQdt (xarray.DataArray): Time derivative of specific humidity (kg/kg/s)

        Sets:
        - self.dQdt: input field
        - self.dQdt_integrated: vertically integrated dQdt (kg/m^2/s)
        - self.divQ: horizontal divergence of moisture flux (1/s)
        - self.divQ_integrated: vertically integrated div(Q) (kg/m^2/s)
        - self.WaterBudgetResidual: dQdt_integrated + divQ_integrated (kg/m^2/s)
        """
        self.dQdt = dQdt
        self.dQdt_integrated = (self.dQdt.integrate(self.vertical_level_indexer) * units('Pa')) / g

        self.divQ = self.calculate_divQ()
        self.divQ_integrated = (self.divQ.integrate(self.vertical_level_indexer) * units('Pa')) / g

        self.WaterBudgetResidual = self.dQdt + self.divQ
        self.WaterBudgetResidual_integrated = self.dQdt_integrated + self.divQ_integrated

    def calculate_horizontal_advection(self, field):
        """
        Calculates horizontal advection for a given field.

        Parameters:
        - field (xarray.DataArray): The field for which to calculate the advection (e.g., temperature, vorticity, moisture).

        Returns:
        - xarray.DataArray: The horizontal advection of the given field.
        """
        dField_dlambda = field.differentiate(self.longitude_indexer)
        dField_dphi = field.differentiate(self.latitude_indexer)
        advection = -1 * ((self.u * dField_dlambda / self.dx) + (self.v * dField_dphi / self.dy))
        return advection

    def tilting_term(self):
        """
        Calculates the tilting term in the vorticity budget equation,
        related to the vertical shear of horizontal wind and horizontal gradients of Omega.

        Returns:
        - xarray.DataArray: Tilting contribution (1/s²)
        """
        dOmegalambda = self.Omega.differentiate(self.longitude_indexer)
        dOmegadphi = self.Omega.differentiate(self.latitude_indexer)
        dOmegady = dOmegalambda / self.dy
        dOmegadx = dOmegadphi / self.dx
        dudp = self.u.differentiate(self.vertical_level_indexer) / (1 * units('Pa'))
        dvdp = self.v.differentiate(self.vertical_level_indexer) / (1 * units('Pa'))
        return (dOmegady * dudp) - (dOmegadx * dvdp)
    
    def calculate_divQ(self):
        """
        Calculates the horizontal divergence of the moisture flux (q * V),
        which corresponds to the divergence of the vertically integrated
        water vapor transport in the water vapor budget equation.

        Mathematically:
            div(Q) = d(q*u)/dx + d(q*v)/dy

        Returns:
        - xarray.DataArray: Horizontal divergence of moisture flux (1/s)
        """
        q_u = self.SpecificHumidity * self.u  # zonal moisture flux
        q_v = self.SpecificHumidity * self.v  # meridional moisture flux

        dq_u_dx = q_u.differentiate(self.longitude_indexer) / self.dx
        dq_v_dy = q_v.differentiate(self.latitude_indexer) / self.dy

        divQ = dq_u_dx + dq_v_dy
        return divQ