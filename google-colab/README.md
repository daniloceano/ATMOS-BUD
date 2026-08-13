# ATMOS-BUD on Google Colab

This directory contains Jupyter notebooks designed to run and support the **ATMOS-BUD** workflow using the **Google Colab** platform.

The notebooks provide a step-by-step workflow for downloading ERA5 reanalysis data, creating cyclone track files, running ATMOS-BUD using a moving domain that follows the cyclone, and optionally generating synoptic charts from ERA5 data.

---

## Workflow

The recommended order for using the notebooks is:

```text
1. Acces_to_Reanalysis_ERA5.ipynb
                  |
                  v
        Download ERA5 data
                  |
                  v
2. track_file.ipynb
                  |
                  v
   Create the cyclone track file
                  |
                  v
3. ATMOSBUD_Google_Colab.ipynb
                  |
                  v
        Run ATMOS-BUD with -t
                  |
                  v
         ATMOS-BUD results
```

The notebook

```text
Synoptic_chart_Reanalysis_ERA5.ipynb
```

is optional and can be used to generate synoptic charts from ERA5 reanalysis data.

---

## 1. Download ERA5 reanalysis data

Before running ATMOS-BUD, the meteorological fields required by the software must first be downloaded from the **ERA5 reanalysis**.

Use:

```text
Acces_to_Reanalysis_ERA5.ipynb
```

This notebook is used to retrieve the ERA5 meteorological fields and create the NetCDF file that will later be used as input for ATMOS-BUD.

The resulting file should contain the meteorological variables and atmospheric levels required by ATMOS-BUD for the period of interest.

---

## 2. Create the cyclone track file

After defining the cyclone or hurricane to be analyzed, use:

```text
track_file.ipynb
```

to create a track file with a name such as:

```text
track_<cyclone_name>.txt
```

For example:

```text
track_cyclone_catarina.txt
```

The track file contains the position and dimensions of the analysis domain for each time step used to follow the cyclone.

Its general structure is:

```text
time;Lat;Lon;length;width
YYYY-MM-DD-HHMM;latitude;longitude;length;width
```

For example:

```text
time;Lat;Lon;length;width
2004-03-20-0000;-26.5559;-40.3947;8.6842;8.3684
2004-03-20-0600;-24.4243;-43.4737;9.4737;7.8947
2004-03-20-1200;-25.5296;-45.2895;9.9474;7.8947
```

The information contained in this file is used by ATMOS-BUD to define a moving domain that follows the cyclone throughout its life cycle.

---

## Important: ERA5 and track temporal consistency

The ERA5 NetCDF file and the cyclone track file must be temporally consistent.

The dates and time steps contained in:

```text
track_<cyclone_name>.txt
```

must correspond to the dates and time steps available in the ERA5 NetCDF file used as input for ATMOS-BUD.

For example, if the track contains:

```text
2004-03-20 00 UTC
2004-03-20 06 UTC
2004-03-20 12 UTC
2004-03-20 18 UTC
```

the ERA5 NetCDF file must contain the corresponding meteorological fields for:

```text
2004-03-20 00 UTC
2004-03-20 06 UTC
2004-03-20 12 UTC
2004-03-20 18 UTC
```

This correspondence must be maintained throughout the entire analysis period.

Therefore, before running ATMOS-BUD, always verify that:

* the ERA5 period covers the complete cyclone track;
* every track time step has a corresponding ERA5 time step;
* the temporal resolution used in the track is consistent with the ERA5 data used for the analysis;
* duplicated timestamps are not present in the track file.

---

## 3. Run ATMOS-BUD on Google Colab

Once both required input files are available:

```text
ERA5 NetCDF file
```

and

```text
track_<cyclone_name>.txt
```

use:

```text
ATMOSBUD_Google_Colab.ipynb
```

to run ATMOS-BUD on Google Colab.

For cyclone-following applications, ATMOS-BUD should be executed using the track option:

```text
-t
```

The `-t` option tells ATMOS-BUD to use the track file to define a moving analysis domain that follows the cyclone or hurricane at each time step.

Conceptually, the execution is:

```text
ERA5 meteorological fields
            +
track_<cyclone_name>.txt
            |
            v
      ATMOS-BUD -t
            |
            v
Moving-domain atmospheric budget analysis
```

The Google Colab notebook provides the necessary steps to:

1. clone the ATMOS-BUD repository;
2. install the required Python packages;
3. connect Google Colab to Google Drive;
4. configure the ERA5 NetCDF input file;
5. configure the cyclone track file;
6. execute ATMOS-BUD using the `-t` option;
7. save the generated results.

---

## 4. Synoptic charts from ERA5

The notebook

```text
Synoptic_chart_Reanalysis_ERA5.ipynb
```

can be used to generate synoptic charts from ERA5 reanalysis data.

The plots are designed following the general visualization style used in the synoptic meteorology products developed by **Alicia M. Bentley**:

https://www.atmos.albany.edu/student/abentley/realtime.html

These charts can be useful for examining the large-scale atmospheric environment associated with the cyclone before, during, and after the ATMOS-BUD analysis.

This notebook is complementary to the ATMOS-BUD workflow and is not required to execute the atmospheric budget calculations.

---

## Notebook summary

| Notebook                                  | Purpose                                                               | Required     |
| ----------------------------------------- | --------------------------------------------------------------------- | ------------ |
| `Acces_to_Reanalysis_ERA5.ipynb` | Download ERA5 meteorological fields and prepare the NetCDF input file | Yes          |
| `track_file.ipynb`                        | Create `track_<cyclone_name>.txt` for the cyclone-following analysis  | Yes for `-t` |
| `ATMOSBUD_Google_Colab.ipynb`             | Run ATMOS-BUD on Google Colab                                         | Yes          |
| `Synoptic_chart_Reanalysis_ERA5.ipynb`    | Generate ERA5 synoptic charts                                         | Optional     |

---

## Recommended directory structure

A possible organization in Google Drive is:

```text
ATMOSBUD/
|
├── data/
|   ├── ERA5/
|   |   └── cyclone_era5.nc
|   |
|   └── track/
|       └── track_cyclone_name.txt
|
└── results/
```

The exact paths can be modified by the user according to their own Google Drive organization.

---

## Complete workflow

```text
Define cyclone and analysis period
              |
              v
Download ERA5 meteorological fields
Acces_to_Reanalysis_ERA5.ipynb
              |
              v
        ERA5 NetCDF file
              |
              |
              +-------------------------+
              |                         |
              v                         v
Create cyclone track          Generate synoptic charts
track_file.ipynb              Synoptic_chart_Reanalysis_ERA5.ipynb
              |                         |
              v                         |
track_<cyclone_name>.txt                |
              |                         |
              +------------+------------+
                           |
                           v
              Check temporal consistency
              ERA5 times == track times
                           |
                           v
              ATMOSBUD_Google_Colab.ipynb
                           |
                           v
                     ATMOS-BUD -t
                           |
                           v
                   Atmospheric budgets
```

---

## Notes

The Google Colab notebooks provide an additional workflow for running ATMOS-BUD in a cloud-based environment and do not replace the standard local installation of the software.

Users should consult the main ATMOS-BUD documentation and repository for additional information about the software, input variables, diagnostics, and available command-line options.

Main ATMOS-BUD repository:

https://github.com/daniloceano/ATMOS-BUD

