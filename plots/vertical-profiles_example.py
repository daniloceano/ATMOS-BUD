#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 19:16:52 2023

@author: daniloceano
"""

import matplotlib.pyplot as plt
import pandas as pd
import os 

###### Inputs ######
start_date = '2005-08-09T12'
end_date = '2005-08-10T18'
budget = 'vorticity' # Can be 'heat', 'vorticity' or 'moisture'
# Note that for the NCEP-R2 fixed domain example, there are missing levels for moisture budget terms

# Set colors
colors = [
    "#81d8f7", "#1a94d0", "#730000", "#db0000", "#f47920", "#f9b233",
    "#00ff00", "#3fae2a", "#00c1a1", "#138d75", "#ef7d90", "#cfcfcf", 
    "#f6a6e1", "#b23aff", "#8377d1", "#8a9aff", "#b4e3ce",
]

# Set path
path = f'./Results/Reg1-Representative_NCEP-R2_fixed/{budget}_terms/'

# Open all csv files in the directory
files = [f for f in os.listdir(path) if f.endswith('.csv')]

fig, ax = plt.subplots(figsize=(10, 10))

i = 0
for file in files:
        variable = file.split('.')[0]

        if (variable == 'vxBeta') or (variable == 'Zeta'):
                continue

        if (variable == 'Omega') or (variable == 'Sigma'):
                continue

        if 'integrated' in variable:
                continue

        # Open data in a Dataframe
        df = pd.read_csv(os.path.join(path, file), index_col=[0])

        # Convert column names to datetime objects
        df.columns = pd.to_datetime(df.columns)

        # Format column names to 'YYYY-MM-DDTHH'
        df.columns = df.columns.strftime('%Y-%m-%dT%H')

        # Convert index from Pa to hPa
        df.index = df.index / 100

        # Slice for a specific range and compute the mean
        selected_data = df.loc[:,start_date:end_date]
        selected_data_mean = selected_data.mean(axis=1)

        # Plot variable on x axis and pressure levels on y axis
        ax.plot(selected_data_mean, selected_data_mean.index, color=colors[i],
                 linewidth=2, linestyle='solid', label=variable)
        ax.axvline(0, c = 'k', linewidth=0.75)
        i += 1

# Reverse y-axis
plt.gca().invert_yaxis()

# Set y limits
plt.ylim(1000, 10)

# Add horizontal grid lines
ax.grid(axis='y', linestyle='--', color='gray', alpha=0.7, linewidth=0.5)

ax.legend(fontsize=14)

# Set title and axis labels
ax.set_title(f'{variable} from {start_date} to {end_date}', fontsize=16)
ax.set_ylabel('Pressure (hPa)', fontsize=14)

# Show plot
plt.tight_layout()
plt.savefig(f'./vertical-profiles_example_{budget}_terms.png')
plt.show()