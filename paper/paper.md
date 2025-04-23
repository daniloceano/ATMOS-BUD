---
title: 'ATMOS-BUD: A Comprehensive Python Tool for Analyzing Heat, Vorticity and Water Budgets on the Atmosphere'
tags:
  - atmospheric diagnostics
  - quasi-geostrophic
  - meteorology
  - climate
  - python
authors:
  - name: Danilo Couto de Souza
    orcid: 0000-0003-4121-7583
    affiliation: 1
  - name: Pedro Leite da Silva Dias
    orcid: 0000-0002-4051-2962
    affiliation: 2
  - name: Ricardo Hallak
    orcid: 0000-0002-8795-1700
    affiliation: 2
affiliations:
  - name:  Universidade de São Paulo (USP), Instituto de Astronomia, Geofísica e Ciências Atmosféricas (IAG), São Paulo, Brazil & Climate Risk Initiative, IRB(re), Rio de Janeiro, Brazil
    index: 1
  - name: Universidade de São Paulo (USP), Instituto de Astronomia, Geofísica e Ciências Atmosféricas (IAG), São Paulo, Brazil
    index: 2
date: 2025-04-22
---

# Summary

# Statement of need

The analysis of heat, vorticity, and water budgets is a powerful tool for diagnosing the dynamical and thermodynamical forcings that drive the development of various atmospheric systems [@shu2022vorticity;sun2024comparison], particularly cyclones, whether extratropical [@liou1987heat], subtropical [@dutra2017structure], or tropical [@raymond2011vorticity]. The heat budget describes the local temperature tendency as the result of three main processes: horizontal advection, adiabatic heating or cooling due to vertical motion, and diabatic heating processes such as latent heat release. The latter is especially important, as it is associated with convection — a key process that fuels the development of systems such as tropical cyclones. In the vorticity budget, the local change in relative vorticity is governed by several mechanisms: horizontal and vertical advection, advection of planetary vorticity, stretching due to mass divergence or convergence, tilting of horizontal vorticity into the vertical, and frictional effects. Lastly, in the water vapor budget, the local change in column-integrated moisture is primarily controlled by two processes: horizontal moisture flux convergence, and net sources or sinks due to surface evaporation and precipitation.

# Description

ATMOS-BUD is built entirely in Python using open-source libraries such as `xarray`, `metpy`, and `matplotlib`. It reads NetCDF input data and a user-defined namelist to identify variables and their units. The program calculates budget terms such as:

- Heat: `dT/dt`, advection, sigma, diabatic heating (residual)
- Vorticity: `dζ/dt`, advection, divergence effects, planetary vorticity advection, tilting
- Moisture: `dq/dt`, divergence of moisture flux, residual

The results are saved in structured folders with CSV summaries, NetCDF grids, logs, and diagnostic figures such as time series, Hovmöller diagrams, and domain maps.

# Impact

ATMOS-BUD can support academic research, operational diagnostics, and case studies in synoptic meteorology. It is particularly valuable for:

Studying tropical and extratropical cyclones
Analyzing forecast model output or reanalysis
Teaching atmospheric dynamics and budgets
The modular structure enables extensions for other vertical coordinates or additional budget terms.

# External Libraries Used

# Acknowledgements