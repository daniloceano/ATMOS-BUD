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
import cmocean
import matplotlib.colors as colors
from matplotlib import cm
import sys
import os
import matplotlib.colors

def white_cmap(cmap):
    n=35
    x = 0.5
    lower = cmap(np.linspace(0, x, n))
    white = np.ones((80-2*n,4))
    upper = cmap(np.linspace(1-x, 1, n))
    colors = np.vstack((lower, white, upper))
    tmap = matplotlib.colors.LinearSegmentedColormap.from_list('map_white', colors)
    return tmap


def time_series(fname,df1,label1,
                df2=None,label2=None):  
    """

    Parameters
    ----------
    fname : str
        Name to append to outfile.
    df1 : pandas.DataFrame
        Dataframe containing data to be ploted. Indexes are vertical levels and
        Columns are timesteps
    label1 : str
        label for y-axis.
    df2 : andas.DataFrame, optional
        same as in df1. The default is None.
    label2 : str, optional
        same as in label1. The default is None.

    Returns
    -------
    Plot time series for the input data vertical levels closest to 1000, 850, 
    700, 500, 300 and 200 hPa and save figures. 

    """
    times = df1.columns
    pdtime = pd.to_datetime(times) 
    plt.close('all')
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(14,10))
    # Find vertical levels that closely matches the desired levels for plotting
    MatchingLevels = []
    for pres in [1000,850,700,500,300,200]:
        MatchingLevels.append(min(df1.index, key=lambda x:abs(x-pres))) 
    max1,min1 = np.amax(np.amax(df1)), np.amin(np.amin(df1))
    if df2 is not None:
        max2,min2 = np.amax(np.amax(df2)), np.amin(np.amin(df2))
    i = 0
    for row in range(2):
        for col in range(3):
            ax = axs[row,col]
            p = MatchingLevels[i]
            ax = axs[row,col]
            ax.plot(pdtime,df1.loc[p],
                    color='#BF3D3B',linewidth=2)
            # If plotting temperature, do not plot horizontal line for 0.
            if label1 != "T (K)" and df2 is None:
                ax.axhline(0,zorder=0,c='#383838',alpha=.8,linewidth=.5,
                            label='')
            ax.tick_params(axis='x',rotation=45)
            ax.tick_params(axis='both',which='major',labelsize=10,
                           labelcolor='#383838')
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
            # If plotting temperature, avoid setting a constant limit for the
            # y-axis as for each vertical level, the values greatly differ.
            if label1 != "T (K)":
                ax.set_ylim(min1*1.2,max1*1.2)
            ax.set_title(str(p)+' hPa ',c='#383838',fontsize=12)
            ax.grid(c='#383838',linewidth=0.25,linestyle='-',alpha=.5)
            # Hide axis except for the first column
            if col == 0:
                ax.set_ylabel(label1,fontsize=12, c='#BF3D3B')
            # ...but as when plotting temperature there is not a fixed limit,
            # it is needed to show the axis values, in this case
            elif label1 == "T (K)":
                pass
            #  we can hide the axsis values for the other cases
            else:
                ax.set_yticklabels([])
            # Plot the second data when it's provided
            if df2 is not None:
                ax2 = ax.twinx()
                ax2.plot(pdtime,df2.loc[p],
                    '--', color='#3B95BF',linewidth=2)
                if label2 != "T (K)":
                    ax2.set_ylim(min2*1.05,max2*1.05)
                ax2.tick_params(axis='both',which='major',labelsize=10,
                           labelcolor='#383838')
                ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
                if col == 2:
                    ax2.set_ylabel(label2,fontsize=12, c='#3B95BF')
                elif label2 == "T (K)":
                    pass
                else:
                    ax2.set_yticklabels([])
            i += 1
    outdir = FigsSubDirectory+'/timeseries/'; os.makedirs(
        outdir, exist_ok=True)
    if df2 is not None:
        outfilename = outdir+'timeseries_compare_'+fname
    else:
        outfilename = outdir+'timeseries_'+fname
    plt.savefig(outfilename,bbox_inches='tight') 
    print(outfilename+' created!')
            
def time_series_thermodyn(ThermDict):    
    """
    
    Parameters
    ----------
    ThermDict : dict
        Dictonary containing DataFrames with each term of the quasi-geostrphic
        thermodynamic equation. For each Dataframe, Indexes are vertical levels 
        and Columns are timesteps

    Returns
    -------
    Plot time series for each term vertical levels closest to 1000, 850, 700,
    500, 300 and 200 hPa and save figures.

    """
    times = ThermDict['AdvHTemp'].columns
    pdtime = pd.to_datetime(times) 
    plt.close('all')
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(14,10))
    # Find vertical levels that closely matches the desired levels for plotting
    MatchingLevels = []
    for pres in [1000,850,700,500,300,200]:
        MatchingLevels.append(min(ThermDict['AdvHTemp'].index, key=lambda x:abs(x-pres)))   
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
            p = MatchingLevels[i]
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
    outdir = FigsSubDirectory+'/thermodynamics/'; os.makedirs(
        outdir, exist_ok=True)
    plt.savefig(outdir+'timeseries_thermodynamics',bbox_inches='tight') 
    print(FigsSubDirectory+'timeseries_thermodynamics created!')
    
def plot_Hovmoller(df,units,fname):
    """
    
    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing data to be plotted. Indexes are vertical levels
        and Columns are timesteps.
    units : str
        label for y-axis.
    fname : str
        Name to append to outfile.

    Returns
    -------
    Create Hovmoller Plots (pressure x time) and store figures.

    """
    times = df.columns
    pdtime = pd.to_datetime(times) 
    plt.close('all')
    fig, ax = plt.subplots(figsize=(14,10))
    levs = sorted(df.index.values, reverse=True)
    max1,min1 = np.amax(np.amax(df)), np.amin(np.amin(df))
    interval = 16
    if min1 < 0:
        norm = colors.TwoSlopeNorm(vmin=min1, vcenter=0, vmax=max1)
        cmap=white_cmap(cmocean.cm.balance)
        if np.abs(max1) > np.abs(min1):
            clev = np.linspace(-max1,max1,interval)
        else:
            clev = np.linspace(min1,-min1,interval)
    else:
        norm = cm.colors.Normalize(vmax=max1, vmin=min1)
        cmap=cmocean.cm.tarn
        clev = np.linspace(min1,max1,interval)
    cf = ax.contourf(pdtime,levs,df,
                     cmap=cmap, extend='both',norm=norm,
                     levels=clev)
    ct = ax.contour(pdtime,levs,df,
                    colors='#383838', extend='both',norm=norm,
                    levels=clev, linewidths=0.5)
    ax.tick_params(axis='x',labelrotation=20)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    ax.set_title(units,fontdict={'fontsize':14})
    ax.set_ylim(levs[0],levs[-1])
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16) 
    # colorbar
    cb_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(cf, cax=cb_ax)
    cbar.add_lines(ct)
    cbar.ax.get_yaxis().labelpad = 15
    for t in cbar.ax.get_yticklabels():
         t.set_fontsize(14)
    outdir = FigsSubDirectory+'/hovmoller/'; os.makedirs(
        outdir, exist_ok=True)
    plt.savefig(outdir+'hovmoller_'+fname,bbox_inches='tight') 
    print(FigsSubDirectory+'hovmoller_'+fname+' created!')
    
def plot_periods_vertical(ThermDict):
    
    linecolors = ['#3B95BF','#87BF4B','#BFAB37','#BF3D3B']
    markers = ['s','o','^','v','<','>']     
    markercolor =  ['#59c0f0','#b0fa61','#f0d643','#f75452','#f07243','#bc6ff7']   
    
    terms = ThermDict.keys()
    periods = pd.read_csv('./periods',sep= ';',header=0)
    
    for i in range(len(periods)):
        fig = plt.figure(figsize=(10, 12))
        ax = fig.add_subplot()
        start,end = periods.iloc[i]['start'],periods.iloc[i]['end']
        period =  periods.iloc[i]['Period']
        
        for term,j in zip(terms,range(len(terms))):
            t = ThermDict[term]
            levs = t.index.values
            # Need to adjust columns to dateformat
            dates = pd.to_datetime(t.columns.values)
            t.columns = dates
            # Get data for selected periods
            selected_dates = t.columns[(t.columns >= start) & (t.columns <= end)]
            data_period = t[selected_dates].transpose().mean()
            # plot
            ax.plot(data_period,levs,label=term,
                    c = linecolors[j], marker=markers[i],linewidth=2,
                    markerfacecolor=markercolor[j])
        # plot cosmedics
        plt.grid(visible=True,c='gray',linewidth=0.25,linestyle='dashdot')
        ax.axvline(0,c='#383838',linewidth=0.5,zorder=1)
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        ax.set_xlabel('(K day-1)',fontsize=18)
        ax.set_ylabel('Pressure (hPa)',fontsize=18)
        plt.legend(prop={'size': 18})
        plt.title(period,fontsize=20)
        plt.ylim(levs[-1],levs[0])
        # save
        outdir = FigsSubDirectory+'/thermodynamics/'; os.makedirs(
            outdir, exist_ok=True)
        plt.savefig(outdir+'vertical_periods'+period,
                    bbox_inches='tight')
    
if __name__ == "__main__":
    
    ResultsSubDirectory = sys.argv[1]
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
        # time_series(term,df,label,df2=None,label2=None)
        plot_Hovmoller(df,label,term)
    # compare eddy area averages 
    for term1,label1 in zip(terms[:3],labels[:3]):
        for term2,label2 in zip(terms[:3],labels[:3]):
            if term1 == term2:
                # do not plot if the same term appears in both axes
                pass
            else:
                df1 = pd.read_csv(ResultsSubDirectory+term1+'.csv',index_col=0)
                df2 = pd.read_csv(ResultsSubDirectory+term2+'.csv',index_col=0)
                time_series(term1[0]+term2,df1,label1,df2=df2,label2=label2)
    # compare area averages       
    for term1,label1 in zip(terms[3:],labels[3:]):
        for term2,label2 in zip(terms[3:],labels[3:]):
            if term1 == term2:
                # do not plot if the same term appears in both axes
                pass
            else:
                df1 = pd.read_csv(ResultsSubDirectory+term1+'.csv',index_col=0)
                df2 = pd.read_csv(ResultsSubDirectory+term2+'.csv',index_col=0)
                # time_series(term1[0]+term2,df1,label1,df2=df2,label2=label2)
            
    # Plot each term of the thermodynamic equation
    ThermDict = {}
    terms = ['AdvHTemp','SpOmega','dTdt','ResT']
    for term in terms:
        ThermDict[term] = pd.read_csv(ResultsSubDirectory+term+'.csv',index_col=0)
        time_series_thermodyn(ThermDict)
        plot_Hovmoller(ThermDict[term],term+' [K day-1]',term)
    
    # Plot vertical profiles for the terms, for each period of the system
    # life cycle
    plot_periods_vertical(ThermDict)