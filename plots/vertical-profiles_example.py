#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 19:16:52 2023

@author: daniloceano
"""

import matplotlib.pyplot as plt
import pandas as pd

# Dummy file
file = '../Results/Reg1-Representative_NCEP-R2_fixed/AdvHTemp.csv'

# Open data in a Dataframe
df = pd.read_csv(file, index_col=[0])
# Convert column names to datetime objects
df.columns = pd.to_datetime(df.columns)
# Format column names to 'YYYY-MM-DDTHH'
df.columns = df.columns.strftime('%Y-%m-%dT%H')

# Slice for a specific range and compute the mean
date1 = '2005-08-09T12'
date2 = '2005-08-10T18'
selected_data = df.loc[:,date1:date2]
selected_data_mean = selected_data.mean(axis=1)

# Plot variable on x axis and pressure levels on y axis
fig, ax = plt.subplots(figsize=(10, 10))
ax.plot(selected_data_mean, df.index, label = 'AdvHTemp',
        color='#453F78', linewidth=2, linestyle='solid')
ax.axvline(0, c = 'k', linewidth=0.75)

# Reverse y-axis
plt.gca().invert_yaxis()

# Set y limits
plt.ylim(1000, 10)

# Add horizontal grid lines
ax.grid(axis='y', linestyle='--', color='gray', alpha=0.7, linewidth=0.5)

ax.legend(fontsize=14)

# Set title and axis labels
ax.set_title('2005-08-09T12 to 2005-08-10T18', fontsize=16)
ax.set_xlabel('[K / s]', fontsize=14)
ax.set_ylabel('Pressure (hPa)', fontsize=14)

# Show plot
plt.show()