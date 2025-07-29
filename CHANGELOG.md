# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2025-07-29

### Added
- Updated CDSAPI syntax in `get_era5_data.py` to support latest API version
- Added troubleshooting documentation for "request too large" errors
- Added CDO merge command example for combining multiple NetCDF files

### Changed
- Updated data download script to use modern CDSAPI syntax
- Improved documentation with additional error handling guidance

### Fixed
- Fixed deprecated CDSAPI syntax that was causing compatibility issues

## [0.0.1] - 2025-06-19

### Added
- Initial release of ATMOS-BUD
- Heat, vorticity and humidity budget analysis for atmospheric regions
- Support for ERA5 reanalysis data
- Command-line interface for budget calculations
- Documentation with getting started guide
- Example scripts for visualization
- Automated data download functionality
- Support for fixed, tracking, and choose domain frameworks

### Features
- **Data Handling**: Load and preprocess atmospheric data from NetCDF files
- **Budget Calculations**: Compute heat, vorticity, and humidity budgets
- **Domain Selection**: Support for fixed regions, storm tracking, and user-defined domains
- **Visualization**: Generate maps, vertical profiles, and Hovmöller diagrams
- **Output Management**: Save results in CSV and NetCDF formats
- **ERA5 Integration**: Direct download from Copernicus Climate Data Store

### Documentation
- Complete user guide with tutorials
- API documentation for all modules
- Examples for different analysis frameworks
- Installation and setup instructions

### Requirements
- Python >= 3.10
- NetCDF4, xarray, numpy, matplotlib
- cartopy for map projections
- cdsapi for ERA5 data download
