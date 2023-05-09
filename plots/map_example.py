# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    map_example.py                                     :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>    +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/05/09 17:05:06 by daniloceano       #+#    #+#              #
#    Updated: 2023/05/09 17:05:06 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import cmocean.cm as cmo
import matplotlib.colors as colors

nc_file = '../Results/Reg1-Representative_NCEP-R2_fixed/Reg1-Representative_NCEP-R2_fixed.nc'
ds = xr.open_dataset(nc_file)

# Open longitude and latitude data
lons, lats = ds.lon_2, ds.lat_2

# Get data for a specific time and level
dTdt = ds.dTdt.sel(initial_time0_hours='2005-08-12T12').sel(lv_ISBL3=1000)


# Define the projection and extent of the map
projection = ccrs.PlateCarree()
extent = [-60, -20, -40, -15]  # [lon_min, lon_max, lat_min, lat_max]

# Create a new figure
fig = plt.figure(figsize=(8, 8))

# Add the map projection and extent to the figure
ax = plt.axes(projection=projection)
ax.set_extent(extent)

# Adjusts limits of colormap to max, min and zero
norm = colors.TwoSlopeNorm(vmin=dTdt.min(), vcenter=0, vmax=dTdt.max())

# Plot the data on the map
im = ax.contourf(lons, lats, dTdt, transform=projection, cmap=cmo.balance, norm=norm)

# Add a colorbar to the plot
cbar = plt.colorbar(im, shrink=0.5)

# Add some features to the map (coastlines, states, countries, etc.)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES)

# Set the title and axis labels
ax.set_title(r'$\frac{\partial T}{\partial t}$ 2005-08-12T12 @ 1000 hPa')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# add gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False  # turn off top labels
gl.right_labels = False  # turn off right labels
gl.xlabel_style = {'size': 12, 'color': 'gray'}  # set x-axis label style
gl.ylabel_style = {'size': 12, 'color': 'gray'}  # set y-axis label style


# Show the plot
plt.show()