# Benchmarks for ATMOS-BUD

This directory contains the benchmark results for various functions used in the ATMOS-BUD project. These benchmarks are calculated using data from the NCEP-R2 dataset for the specified date, vertical level, and latitude/longitude ranges.

## How to Obtain the Data

The data used for generating the benchmarks is from the NCEP-R2 (National Centers for Environmental Prediction - Reanalysis 2) dataset. The specific data for the following parameters is required:

- **Initial time**: 2005-08-12T12
- **Vertical level**: 1000 hPa
- **Latitude range**: -22.5 to 0.0 degrees
- **Longitude range**: -100.0 to -77.5 degrees

You can access the NCEP-R2 dataset from the [NOAA ESRL website](https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.derived.html). Download the dataset for the required date and variables.

## How to Generate the Benchmarks

To generate the benchmark results, run the `generate_benchmarks.py` script. This script will calculate the benchmarks for various functions such as `CalcZonalAverage` and `CalcAreaAverage` using the NCEP-R2 data.

### Steps:
1. **Download the NCEP-R2 dataset** as mentioned above.
2. Place the dataset in the `samples/` directory.
3. Run the following command to generate the benchmarks:

```bash
python benchmarks/generate_benchmarks.py
```

## How to Use the Benchmarks

These benchmark files are used in the tests to validate that the implemented functions (like `CalcZonalAverage` and `CalcAreaAverage`) are producing correct results. During testing, the computed results are compared to these benchmarks to ensure the accuracy of the functions.

Run the tests with:

```bash
pytest
```
This will automatically compare the test results with the stored benchmarks.
