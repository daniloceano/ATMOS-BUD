# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    hovmoller_example.py                               :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>    +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/05/09 15:10:06 by daniloceano       #+#    #+#              #
#    Updated: 2023/05/09 15:10:06 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.dates import DateFormatter
import cmocean.cm as cmo
import matplotlib.colors as colors
import matplotlib.dates as mdates

# Dummy file
file = '../Results/Reg1-Representative_NCEP-R2_fixed/AdvHTemp.csv'

# Open data in a Dataframe
df = pd.read_csv(file, index_col=[0])
# Convert column names to datetime objects
df.columns = pd.to_datetime(df.columns)
# Format column names to 'YYYY-MM-DDTHH'
df.columns = df.columns.strftime('%Y-%m-%dT%H')

# Adjusts limits of colormap to max, min and zero
norm = colors.TwoSlopeNorm(vmin=df.to_numpy().min(), vcenter=0, vmax=df.to_numpy().max())

# plot hovmoller diagram
fig, ax = plt.subplots(figsize=(10,10))
im = ax.contourf(df.columns, df.index, df.values, cmap=cmo.balance, levels=7)
cbar = fig.colorbar(im)

# invert y-axis
ax.invert_yaxis()

# format x-axis as date and rotate tick labels
ax.xaxis.set_major_locator(mdates.AutoDateLocator())
ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
ax.tick_params(axis='x', labelrotation=45)

# Set title and axis labels
ax.set_title('AdvHTemp', fontsize=16)
ax.set_ylabel('Pressure (hPa)', fontsize=14)

# edit colorbar legend
cbar.set_label('[K / s]', fontsize=12)
cbar.ax.tick_params(labelsize=10)

plt.show()