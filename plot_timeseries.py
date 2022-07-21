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
import cmocean as cmo
import matplotlib.colors as colors
from matplotlib import cm

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
                    color='#BF3D3B',linewidth=2)
            if label1 != "T_AA":
                ax.axhline(0,zorder=0,c='#383838',alpha=.8,linewidth=.5,
                            label='')
            ax.tick_params(axis='x',rotation=45)
            ax.tick_params(axis='both',which='major',labelsize=10,
                           labelcolor='#383838')
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
            if label1 != "T_AA":
                ax.set_ylim(min1*1.2,max1*1.2)
            ax.set_title(str(p)+' hPa ',c='#383838',fontsize=12)
            ax.grid(c='#383838',linewidth=0.25,linestyle='-',alpha=.5)
            if col == 0:
                ax.set_ylabel(label1,fontsize=12, c='#BF3D3B')
            elif label1 == "T_AA":
                pass
            else:
                ax.set_yticklabels([])
            if df2 is not None:
                ax2 = ax.twinx()
                ax2.plot(pdtime,df2.loc[p],
                    '--', color='#3B95BF',linewidth=2)
                if label2 == "T_AA":
                    ax.set_ylim(200,360)
                else:
                    ax2.set_ylim(min2*1.05,max2*1.05)
                ax2.tick_params(axis='both',which='major',labelsize=10,
                           labelcolor='#383838')
                if col == 2:
                    ax2.set_ylabel(label2,fontsize=12, c='#3B95BF')
                elif label2 == "T_AA":
                    pass
                else:
                    ax2.set_yticklabels([])
            i += 1
    plt.savefig(FigsSubDirectory+'timeseries_'+fname,bbox_inches='tight') 
    print(FigsSubDirectory+'timeseries_'+fname+' created!')
            
def time_series_thermodyn(ThermDict):    
    times = ThermDict['AdvHTemp'].columns
    pdtime = pd.to_datetime(times) 
    plt.close('all')
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(14,10))
    pres = [1000,850,700,500,300,200]
    cols = ['#3B95BF','#87BF4B','#BFAB37','#BF3D3B']
    # Data range
    maxv = minv = []
    for term in ThermDict.keys():
        maxv.append(np.amax(ThermDict[term]))
        minv.append(np.amin(ThermDict[term]))
    max1, min1 = np.amax(maxv), np.amin(minv)
    # iterator
    i = 0
    for row in range(2):
        for col in range(3):
            ax = axs[row,col]
            p = pres[i]
            ax = axs[row,col]
            for term,c in zip(ThermDict.keys(),cols):
                ax.plot(pdtime,ThermDict[term].loc[p],
                        color=c,linewidth=2,label=term)
            ax.axhline(0,zorder=0,c='#383838',alpha=.8,linewidth=.5,
                        label='')
            ax.tick_params(axis='x',rotation=45)
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
            ax.tick_params(axis='both',which='major',labelsize=10,
                           labelcolor='#383838')
            ax.set_ylim(min1,max1)
            ax.set_title(str(p)+' hPa ',c='#383838',fontsize=12)
            ax.grid(c='#383838',linewidth=0.25,linestyle='-',alpha=.5)
            if col == 0:
                ax.set_ylabel('K day-1',fontsize=12, c='#383838')
            else:
                ax.set_yticklabels([])
            if col == 2 and row == 0:
                    ax.legend(fontsize=10,ncol=1,
                              bbox_to_anchor=(0.95, -0.2, 0.5, 1))
            i += 1
            
    plt.gcf().autofmt_xdate()
    # ax.legend(fontsize=12)
    plt.savefig(FigsSubDirectory+'timeseries_thermodynamics',bbox_inches='tight') 
    print(FigsSubDirectory+'timeseries_thermodynamics created!')
    
def plot_Hovmoller(df,units,fname):
    times = df.columns
    pdtime = pd.to_datetime(times) 
    plt.close('all')
    fig, ax = plt.subplots(figsize=(14,10))
    max1,min1 = np.amax(np.amax(df)), np.amin(np.amin(df))
    if min1 < 0:
        norm = colors.TwoSlopeNorm(vmin=min1, vcenter=0, vmax=max1)
        cmap='cmo.curl'
    else:
        norm = cm.colors.Normalize(vmax=max1, vmin=min1)
        cmap='cmo.tarn'
    levs = sorted(df.index.values, reverse=True)
    cf = ax.contourf(pdtime,levs,df,cmap=cmap, extend='both',norm=norm)
    ax.contour(pdtime,levs,df,colors='k', extend='both')
    ax.tick_params(axis='x',labelrotation=20)
    ax.tick_params(size=12)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    ax.set_title(fname,fontdict={'fontsize':14})
    ax.set_ylim(levs[0],levs[-1])
    # plt.yscale('log')
    # ax.set_yticks(levs)
    # colorbar
    cb_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(cf, cax=cb_ax)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel(units, rotation=270,fontsize=12)
    for t in cbar.ax.get_yticklabels():
         t.set_fontsize(10)
    plt.savefig(FigsSubDirectory+'hovmoller_'+fname,bbox_inches='tight') 
    print(FigsSubDirectory+'hovmoller_'+fname+' created!')
    
if __name__ == "__main__":
    
    ResultsSubDirectory = '../CycloneThermodynamics_Results/Reg1-Representative_NCEP-R2_lagranigan/'
    # Directory for saving Figures
    FigsSubDirectory = ResultsSubDirectory+'/Figures/'
    # Check if the LEC_Figures directory exists. If not, creates it
    check_create_folder(FigsSubDirectory)
    # Check if a directory for current data exists. If not, creates it
    check_create_folder(FigsSubDirectory)
    
    # Plot averaged terms
    terms = ['T_AA','Omega_AA','Q_AA',
                    'T_ZE_AA','Omega_ZE_AA','Q_ZE_AA']
    labels = ["T (K)","ω (Pa s-1)","Q (J Kg-1 s-1)",
              "[T'] (K)","[ω'] (Pa s-1)","[Q'] (J Kg-1 s-1)"]
    for term,label in zip(terms,labels):
        df = pd.read_csv(ResultsSubDirectory+term+'.csv',index_col=0)
        time_series(term,df,label,df2=None,label2=None)
        plot_Hovmoller(df,label,term)
    # compare eddy area averages 
    for term1,term2,label1,label2 in zip(terms[:3],terms[:3],
                                         labels[:3],labels[:3]):
        if term1 == term2:
            pass
        else:
            df1 = pd.read_csv(ResultsSubDirectory+term1+'.csv',index_col=0)
            df2 = pd.read_csv(ResultsSubDirectory+term2+'.csv',index_col=0)
            time_series(term,df1,label1,df2=df2,label2=term2)
    # compare area averages       
    for term1,term2,label1,label2 in zip(terms[3:],terms[3:],
                                         labels[3:],labels[3:]):
        if term1 == term2:
            pass
        else:
            df1 = pd.read_csv(ResultsSubDirectory+term1+'.csv',index_col=0)
            df2 = pd.read_csv(ResultsSubDirectory+term2+'.csv',index_col=0)
            time_series(term,df1,label1,df2=df2,label2=term2)
            
    # # Plot terms of the thermodynamic equation
    ThermDict = {}
    terms = ['AdvHTemp','SpOmega','dTdt','ResT']
    for term in terms:
        ThermDict[term] = pd.read_csv(ResultsSubDirectory+term+'.csv',index_col=0)
        time_series_thermodyn(ThermDict)
        plot_Hovmoller(ThermDict[term],'K day-1',term)