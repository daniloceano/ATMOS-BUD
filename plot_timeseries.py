#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 09:43:42 2022

@author: danilocoutodsouza
"""

import matplotlib.pyplot as plt
import matplotlib.dates
import numpy as np
from main import check_create_folder
import glob

def time_series(fname,data2D_1,unit_label1,
                data2D_2=None,unit_label2=None):    
    times = full_data[TimeName].values
    dates = []
    for t in times:
        dates.append(str(t)[:19])
    dates = matplotlib.dates.date2num(times)
    plt.close('all')
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(14,10))
    pres = [1000,850,700,500,300,200]
    max1,min1 = np.amax(data2D_1), np.amin(data2D_1)
    if data2D_2 is not None:
        max2,min2 = np.amax(data2D_2), np.amin(data2D_2)
    i = 0
    for row in range(2):
        for col in range(3):
            ax = axs[row,col]
            p = pres[i]
            ax = axs[row,col]
            ax.plot(times,data2D_1.sel({VerticalCoordIndexer:p}),
                    color='#3B95BF',linewidth=2)
            ax.axvline(dates[itime],zorder=0,c='#383838',alpha=0.2,linewidth=2,
                        label='')
            ax.axhline(0,zorder=0,c='#383838',alpha=.8,linewidth=.5,
                        label='')
            ax.tick_params(axis='x',rotation=45)
            ax.tick_params(axis='both',which='major',labelsize=10,
                           labelcolor='#383838')
            ax.set_ylim(min1*1.2,max1*1.2)
            ax.set_title(str(p)+' hPa ',c='#383838',fontsize=12)
            ax.grid(c='#383838',linewidth=0.25,linestyle='-',alpha=.5)
            if col == 0:
                ax.set_ylabel(unit_label1,fontsize=12, c='#3B95BF')
            else:
                ax.set_yticklabels([])
            if data2D_2 is not None:
                ax2 = ax.twinx()
                ax2.plot(times,data2D_2.sel({VerticalCoordIndexer:p}),
                    '--', color='#BF3D3B',linewidth=2)
                ax2.set_ylim(min2*1.2,max2*1.2)
                ax2.tick_params(axis='both',which='major',labelsize=10,
                           labelcolor='#383838')
                if col == 2:
                    ax2.set_ylabel(unit_label2,fontsize=12, c='#BF3D3B')
                else:
                    ax2.set_yticklabels([])
            i += 1
            
    plt.gcf().autofmt_xdate()
    # ax.legend(fontsize=12)
    plt.savefig(outdir+'timeseries_'+fname,bbox_inches='tight') 
    print(outdir+'timeseries_'+fname+' created!')
    
if __name__ == "__main__":
    
    ResultsSubDirectory = '../CycloneThermodynamics_Results/Reg1-Representative_NCEP-R2_lagranigan/'
    # Directory for saving Figures
    FigsSubDirectory = ResultsSubDirectory+'/Figures/'
    # Check if the LEC_Figures directory exists. If not, creates it
    check_create_folder(FigsSubDirectory)
    # Check if a directory for current data exists. If not, creates it
    check_create_folder(FigsSubDirectory)
    
    # csvFiles = glob.glob(ResultsSubDirectory+'/*.csv')
    # for file in csvFiles:
        
    