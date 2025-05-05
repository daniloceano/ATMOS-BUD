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
    affiliation: "1, 2"
  - name: Pedro Leite da Silva Dias
    orcid: 0000-0002-4051-2962
    affiliation: 1
  - name: Ricardo Hallak
    orcid: 0000-0002-8795-1700
    affiliation: 1
affiliations:
  - name:  Universidade de São Paulo (USP), Instituto de Astronomia, Geofísica e Ciências Atmosféricas (IAG), São Paulo, Brazil & 
    index: 1
  - name: Climate Risk Initiative, IRB(re), Rio de Janeiro, Brazil
    index: 2
date: 22 April 2025
bibliography: paper.bib
---

# Summary

`ATMOS-BUD` is a flexible and open-source Python package designed to compute and visualize the main components of atmospheric heat, vorticity, and water budgets in pressure coordinates. It provides a structured and modular pipeline to process meteorological datasets, particularly reanalyses such as ERA5 or model outputs, and compute budget residuals using a quasi-geostrophic framework. The tool supports both fixed and moving computational domains (Eulerian and Semi-Lagrangian frameworks), allowing users to diagnose and track synoptic and mesoscale systems such as cyclones. `ATMOS-BUD` facilitates both operational diagnostics and academic research, offering clearly documented routines and outputs in CSV, NetCDF, and figure formats.


# Background

The analysis of heat, vorticity, and water budgets is a powerful tool for diagnosing the dynamical and thermodynamical forcings that drive the development of various atmospheric systems [@shu2022vorticity] [@sun2024comparison], particularly cyclones, whether extratropical [@liou1987heat], subtropical [@dutra2017structure], or tropical [@raymond2011vorticity]. The heat budget describes the local temperature tendency as the result of three main processes: horizontal advection, adiabatic heating or cooling due to vertical motion, and diabatic heating processes such as latent heat release. The latter is especially important, as it is associated with convection — a key process that fuels the development of systems such as tropical cyclones. In the vorticity budget, the local change in relative vorticity is governed by several mechanisms: horizontal and vertical advection, advection of planetary vorticity, stretching due to mass divergence or convergence, tilting of horizontal vorticity into the vertical, and frictional effects. Lastly, in the water vapor budget, the local change in column-integrated moisture is primarily controlled by two processes: horizontal moisture flux convergence, and net sources or sinks due to surface evaporation and precipitation.

# Statement of need

`ATMOS-BUD` is a Python tool designed to compute heat, vorticity, and water budgets for limited atmospheric regions. It supports both Eulerian and Semi-Lagrangian frameworks. In the former, the user specifies a fixed computational domain for budget computations; in the latter, the domain is updated at each time step to follow the center of an atmospheric system of interest (e.g., a cyclone). This methodology mirrors that adopted in [@de2024lorenzcycletoolkit] (see their Figure 2). Additionally, computations can be performed using an interactive framework, allowing the user to manually define the domain for each time step. The program outputs budget terms averaged over the selected domain and, for moisture variables, provides both layerwise and vertically integrated diagnostics.

`ATMOS-BUD` supports academic research, operational diagnostics, and case studies in synoptic meteorology. It has also been used for several years in graduate-level courses such as *Synoptic Meteorology III* and *Tropical Meteorology* at the Department of Atmospheric Sciences at the University of São Paulo. The tool was initially developed using the Grid Analysis and Display System (GrADS) [@doty1995geophysical], but over the past three years, it has been rewritten in Python to meet pedagogical needs and facilitate integration into scientific workflows and publications.

# External Libraries Used

The input file must be in NetCDF format [@rew1990netcdf], which is processed using `xarray` [@hoyer2017xarray] and `dask` [@daniel2019data] to support efficient memory management and scalable processing of large datasets. Results are managed with `pandas` [@reback2020pandas] for tabular manipulation and exported in CSV format. Numerical computations are performed using `numpy` [@harris2020array], while `metpy` [@may2022metpy] is used to handle meteorological constants, units, and atmospheric diagnostics.


# Acknowledgements

This study was partly financed by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior – Brasil (CAPES) under Finance Code 001. The authors would like to acknowledge the Laboratory of Applied Meteorology for Regional Meteorological Systems (MASTER) for the data processing infrastructure. Special thanks are extended to Jean Peres and Djalma Vieira for their invaluable technical assistance and support in the laboratory.