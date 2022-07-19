#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 09:43:42 2022

@author: danilocoutodsouza
"""

import matplotlib.pyplot as plt
import numpy as np
from main import check_create_folder
import pandas as pd
import matplotlib.dates as mdates

def time_series(fname,df1,label1,
                df2=None,label2=None):    
    times = df1.columns
    pdtime = pd.to_datetime(times) 
    plt.close('all')
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(14,10))
    pres = [1000,850,700,500,300,200]
    max1,min1 = np.amax(np.amax(df1)), np.amin(np.amin(df1))
    if df2 is not None:
        max2,min2 = np.amax(np.amax(df2)), np.amin(np.amin(df2))
    i = 0
    for row in range(2):
        for col in range(3):
            ax = axs[row,col]
            p = pres[i]
            ax = axs[row,col]
            ax.plot(pdtime,df1.loc[p],
                    color='#3B95BF',linewidth=2)
            ax.axhline(0,zorder=0,c='#383838',alpha=.8,linewidth=.5,
                        label='')
            ax.tick_params(axis='x',rotation=45)
            ax.tick_params(axis='both',which='major',labelsize=10,
                           labelcolor='#383838')
            ax.set_ylim(min1*1.2,max1*1.2)
            ax.set_title(str(p)+' hPa ',c='#383838',fontsize=12)
            ax.grid(c='#383838',linewidth=0.25,linestyle='-',alpha=.5)
            if col == 0:
                ax.set_ylabel(label1,fontsize=12, c='#3B95BF')
            else:
                ax.set_yticklabels([])
            if df2 is not None:
                ax2 = ax.twinx()
                ax2.plot(pdtime,df2.loc[p],
                    '--', color='#BF3D3B',linewidth=2)
                ax2.set_ylim(min2*1.2,max2*1.2)
                ax2.tick_params(axis='both',which='major',labelsize=10,
                           labelcolor='#383838')
                if col == 2:
                    ax2.set_ylabel(label2,fontsize=12, c='#BF3D3B')
                else:
                    ax2.set_yticklabels([])
            i += 1
            
    plt.gcf().autofmt_xdate()
    # ax.legend(fontsize=12)
    plt.savefig(FigsSubDirectory+'timeseries_'+fname,bbox_inches='tight') 
    print(FigsSubDirectory+'timeseries_'+fname+' created!')
    
if __name__ == "__main__":
    
    ResultsSubDirectory = '../CycloneThermodynamics_Results/Reg1-Representative_NCEP-R2_lagranigan/'
    # Directory for saving Figures
    FigsSubDirectory = ResultsSubDirectory+'/Figures/'
    # Check if the LEC_Figures directory exists. If not, creates it
    check_create_folder(FigsSubDirectory)
    # Check if a directory for current data exists. If not, creates it
    check_create_folder(FigsSubDirectory)
    
    dfQ = pd.read_csv(ResultsSubDirectory+'Q_ZE_AA.csv',index_col=0)
    dfT = pd.read_csv(ResultsSubDirectory+'T_ZE_AA.csv',index_col=0)
    dfO = pd.read_csv(ResultsSubDirectory+'Omega_ZE_AA.csv',index_col=0)
    
    time_series('T_ZE_AA',dfT,"[T'] (K)",df2=None,label2=None)
    time_series('TQ_ZE_AA',dfT,"[T'] (K)",df2=dfQ,label2="[Q'] (J Kg-1 s-1)")
    time_series('TOmega_ZE_AA',dfT,"[T'] (K)",df2=dfO,label2="[ω'] (Pa s-1)")
    time_series('OmegaQ_ZE_AA',dfO,"[ω'] (Pa s-1)",df2=dfQ,label2="[Q'] (J Kg-1 s-1)")
        
    