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
  - name: IRB Brasil RE & Universidade de São Paulo (USP), São Paulo, Brazil
    index: 1
  - name: Universidade de São Paulo (USP), Instituto de Astronomia, Geofísica e Ciências Atmosféricas (IAG), São Paulo, Brazil
    index: 2
date: 2025-04-22
---

# Summary

ATMOS-BUD is a flexible and open-source Python package designed to compute and visualize the main components of atmospheric heat, vorticity, and water vapor budgets in pressure coordinates. It provides a structured and modular pipeline to process meteorological datasets, particularly reanalyses such as ERA5 or model output, and compute budget residuals using a quasi-geostrophic framework. 

The tool supports fixed and moving analysis windows (semi-Lagrangian mode), allowing users to diagnose and track synoptic and mesoscale systems such as cyclones. ATMOS-BUD facilitates both operational diagnostics and academic research, with clearly documented routines and outputs in CSV, NetCDF, and figure formats.

# Statement of need

There is a lack of streamlined and open-source tools for diagnosing atmospheric budgets with temporal continuity and multi-level integration, especially for studies involving system tracking or dynamically evolving regions. While numerical models compute these tendencies internally, researchers often lack access to component terms such as horizontal/vertical advection, diabatic heating, tilting, and divergence-related effects.

ATMOS-BUD fills this gap by offering:
- A reproducible framework for isolating and visualizing key dynamical and thermodynamical processes,
- Multi-framework support (fixed box, tracking, or interactive selection),
- Compatibility with standard reanalyses and output formats.

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