#!/usr/bin/env python3
"""
Interactive cyclone tracking tool.

This script reads a NetCDF file, displays a selected meteorological field
for each time step, and allows the user to draw a rectangular box around
a cyclone using the mouse.

For each selected box, the script saves:
    - time
    - centroid latitude
    - centroid longitude
    - box length in degrees longitude
    - box width in degrees latitude

The output file follows the format:

    time;Lat;Lon;length;width

Example:
    2026-04-01-0000;-35.2500;190.5000;8.0000;6.5000

The script preserves the longitude convention of the input NetCDF file:
    - 0 to 360
    - -180 to 180

Optionally, a continent shapefile can be overlaid using pyshp.
"""

import os
import argparse

import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import shapefile

# ==========================================================
# Configurações
# ==========================================================
parser = argparse.ArgumentParser(
    description=(
        "Interactively extract cyclone centroid coordinates and box "
        "dimensions from a NetCDF file."
    )
)

parser.add_argument(
    "file_nc",
    type=str,
    help="Path to the input NetCDF file."
)

parser.add_argument(
    "--var",
    type=str,
    default="msl",
    help="Variable name in the NetCDF file. Default: msl"
)

parser.add_argument(
    "--level",
    type=float,
    default=None,
    help=(
        "Optional pressure level to select if the dataset has a vertical "
        "pressure-level coordinate. Example: --level 850."
    )
)

parser.add_argument(
    "--out",
    type=str,
    default="track_ciclone.txt",
    help="Name of the output file (default: track_ciclone.txt)"
)

parser.add_argument("--shp", type=str, default=None,
                    help="Optional continent shapefile path.")

args = parser.parse_args()

file_nc = args.file_nc
outfile = args.out
varname = args.var
level_value = args.level
shp_path = args.shp

# ==========================================================
# Input checks
# ==========================================================
if not os.path.exists(file_nc):
    raise FileNotFoundError(f"File not found: {file_nc}")

if shp_path and not os.path.exists(shp_path):
    raise FileNotFoundError(f"Shapefile not found: {shp_path}")

# ==========================================================
# Open dataset
# ==========================================================
ds = xr.open_dataset(file_nc)

# ==========================================================
# Select pressure level, if requested
# ==========================================================
if level_value is not None:
    possible_level_names = [
        "level",
        "pressure_level",
        "isobaricInhPa",
        "plev",
        "lev"
    ]

    level_name = None

    for name in possible_level_names:
        if name in ds.coords or name in ds.dims:
            level_name = name
            break

    if level_name is None:
        raise ValueError(
            "The --level option was used, but no pressure-level coordinate "
            "was found in the NetCDF file."
        )

    ds = ds.sel({level_name: level_value}, method="nearest")

    print(f"Selected level: {level_name} = {float(ds[level_name].values)}")
    
# ==========================================================
# Detect coordinate names automatically
# ==========================================================
lat_name = None
lon_name = None
time_name = None

for name in ds.coords:
    lname = name.lower()
    if lname in ["lat", "latitude"]:
        lat_name = name
    elif lname in ["lon", "longitude"]:
        lon_name = name
    elif lname in ["time", "valid_time", "date"]:
        time_name = name

if lat_name is None or lon_name is None or time_name is None:
    raise ValueError(
        "Could not automatically detect latitude, longitude, or time "
        "coordinates in the NetCDF file."
    )

# ==========================================================
# Domain and longitude convention
# ==========================================================
lat_min = float(ds[lat_name].min())
lat_max = float(ds[lat_name].max())
lon_min = float(ds[lon_name].min())
lon_max = float(ds[lon_name].max())

use_lon_360 = lon_max > 180

print(
    "Detected longitude convention: "
    f"{'0-360' if use_lon_360 else '-180 to 180'}"
)

# ==========================================================
# Shapefile plotting function
# ==========================================================
def plot_shapefile(ax, 
                   shp_path, 
                   use_lon_360,
                   edgecolor="black",
                   linewidth=0.6, 
                   label=None):
    """
    Plot shapefile contours using Matplotlib.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis where the shapefile contours will be drawn.
    shp_path : str
        Path to the shapefile.
    use_lon_360 : bool
        If True, shapefile longitudes are converted from -180/180 to 0/360.
        If False, shapefile longitudes are kept or converted to -180/180.
    edgecolor : str
        Contour color.
    linewidth : float
        Contour line width.
    label : str or None
        Optional label for the legend. The label is added only once.

    Notes
    -----
    The function splits line segments when large longitude jumps are detected.
    This avoids artificial lines near the dateline.
    """

    sf = shapefile.Reader(shp_path)
    first = True

    for shp in sf.shapes():
        pts = np.array(shp.points)
        if pts.size == 0:
            continue

        x = pts[:, 0]
        y = pts[:, 1]

        if use_lon_360:
            x = np.where(x < 0, x + 360, x)
        else:
            x = np.where(x > 180, x - 360, x)

        # Split multipart geometries.
        parts = list(shp.parts) + [len(pts)]

        for i in range(len(parts) - 1):
            xs = x[parts[i]:parts[i+1]]
            ys = y[parts[i]:parts[i+1]]

            # Split segments that cross the dateline after longitude conversion.
            jumps = np.abs(np.diff(xs)) > 180

            start = 0

            for j, jump in enumerate(jumps):
                if jump:
                    if first and label is not None:
                        ax.plot(xs[start:j+1], ys[start:j+1], 
                                color=edgecolor, linewidth=linewidth, label=label)
                        first = False
                    else:
                        ax.plot(xs[start:j+1], ys[start:j+1], 
                                color=edgecolor, linewidth=linewidth)
                    start = j + 1
            
            if first and label is not None:
                ax.plot(xs[start:], ys[start:], color=edgecolor, linewidth=linewidth, label=label)
                first = False
            else:
                ax.plot(xs[start:], ys[start:], color=edgecolor, linewidth=linewidth)

# ==========================================================
# Prepare output file
# time;Lat;Lon,length;width
# ==========================================================
nt = len(ds[time_name])
da_all = ds[varname]

with open(outfile, "w") as f:
    f.write("time;Lat;Lon;length;width\n")

# ==========================================================
# Time loop
# ==========================================================
for n in range(nt):

    da = da_all.isel({time_name: n})

    # Convert Pa to hPa when the field appears to be pressure in Pa.
    if float(da.mean()) > 2000:
        da_plot = da / 100.0
    else:
        da_plot = da

    time_value = pd.to_datetime(ds[time_name].isel({time_name: n}).values)
    time_txt = time_value.strftime("%Y-%m-%d-%H%M")
    time_title = time_value.strftime("%HZ%d%b%Y").upper()

    # ======================================================
    # Create plot
    # ======================================================
    fig, ax = plt.subplots(figsize=(7, 9))
        
    vmin = np.floor(float(da_plot.min()))
    vmax = np.ceil(float(da_plot.max()))
    levels = np.arange(vmin, vmax + 4, 4)
    
    cs = ax.contour(
        ds[lon_name].values,
        ds[lat_name].values,
        da_plot.values,
        levels=levels,
        colors="black",
        linewidths=1
    )

    ax.clabel(cs, inline=True, fontsize=8, fmt="%.0f")

    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(f"Mean sea level pressure {time_title}", fontsize=14)

    ax.grid(True, linestyle="--", alpha=0.4)

    # Plot optional continent shapefile
    if shp_path is not None:
        plot_shapefile(
            ax, 
            shp_path, 
            use_lon_360,  
            edgecolor="red", 
            linewidth=1.5,
            label="Continents"
            )
    ax.legend(loc="upper right")

    print("\n======================================")
    print(f"Time step {n+1}/{nt}: {time_txt}")
    print("Left-click, drag the rectangle, and release to select the cyclone box.")
    print("Close the window to skip this time step.")
    print("======================================")

    selection = {}

    def onselect(eclick, erelease):
        """
        Callback executed after drawing the rectangle.

        The selected box limits are stored in the 'selection' dictionary.
        """
        lon_a, lat_a = eclick.xdata, eclick.ydata
        lon_b, lat_b = erelease.xdata, erelease.ydata

        if lon_a is None or lon_b is None or lat_a is None or lat_b is None:
            return

        lon_left = min(lon_a, lon_b)
        lon_right = max(lon_a, lon_b)
        lat_bottom = min(lat_a, lat_b)
        lat_top = max(lat_a, lat_b)

        selection["lon_left"] = lon_left
        selection["lon_right"] = lon_right
        selection["lat_bottom"] = lat_bottom
        selection["lat_top"] = lat_top

        plt.close(fig)

    selector = RectangleSelector(
        ax,
        onselect,
        useblit=True,
        button=[1],
        minspanx=0.1,
        minspany=0.1,
        spancoords="data",
        interactive=True
    )

    plt.show()

    if not selection:
        print("No box selected. Skipping this time step.")
        plt.close(fig)
        continue

    lon_left = selection["lon_left"]
    lon_right = selection["lon_right"]
    lat_bottom = selection["lat_bottom"]
    lat_top = selection["lat_top"]

    # Compute box centroid and dimensions in degrees.
    lon_center = (lon_left + lon_right) / 2.0
    lat_center = (lat_bottom + lat_top) / 2.0

    length = lon_right - lon_left
    width = lat_top - lat_bottom

    # Save selected box information.
    with open(outfile, "a") as f:
        f.write(
            f"{time_txt};"
            f"{lat_center:.4f};"
            f"{lon_center:.4f};"
            f"{length:.4f};"
            f"{width:.4f}\n"
        )

    print(
        f"Saved: {time_txt};"
        f"{lat_center:.4f};"
        f"{lon_center:.4f};"
        f"{length:.4f};"
        f"{width:.4f}"
    )


print("\nProcessing completed.")
print(f"Output file saved as: {outfile}")