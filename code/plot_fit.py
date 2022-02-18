import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
import numpy as np
sn.set()

from compute_5PL import r_squared, PL5, fit_standard
from collect_luminex import get_standard

def plot_standard(ax, s, conc, fit, color="red"):
    _ = ax.scatter(conc, fit, color=color, s=.5, alpha=0.5)
    _ = ax.scatter(s['conc'], s['FI'], color=color, s=45)


def fit_curve(params, zero_value=0.1):
    '''
    returns the fit curve for plotting
    '''

    conc = np.power(10, np.linspace(np.log10(zero_value), 5, 10000))
    fit = PL5(conc, params)

    return conc, fit



def fit_curve(params, zero_value=0.1, ymax=100000):
    '''
    returns the fit curve for plotting
    '''

    conc = np.power(10, np.linspace(np.log10(zero_value), np.log10(ymax*1.2)+1, 10000))
    fit = PL5(conc, params)
    return conc, fit
    

def plot_standard(ax, st_df, params, color="blue", ymax=10000, zero_value=0.1, **kwargs):
    '''
    subplot the standard curve
    '''
    
    # fit the curve
    conc, fit = fit_curve(params, zero_value=zero_value, ymax=ymax)
    _ = ax.scatter(conc, fit, color=color, s=.1, alpha=0.5)
    _ = ax.scatter(st_df['conc'], st_df['FI'], **kwargs)

    
def plot_values(ax, data_df, **kwargs):
    '''
    subplot the standard curve
    '''
    
    # fit the curve
    _ = ax.scatter(data_df['conc'], data_df['FI'], **kwargs)

    
def plot_fitting(st_df, sample_df=pd.DataFrame(), control_df=pd.DataFrame(), params=[], R="", zero_value=0.1, show_params=True):
    '''
    calculate samples for a given set of 
    '''
    
    # make sure, params is a list
    params = list(params)
    
    
    # init the plot
    fig, ax = plt.subplots(figsize=(12,12))
    
    # calculate R if R is not given
    if not R:
        R = r_squared(params, st_df) 
    ymax = st_df['FI'].max()
    
    if not sample_df.empty:
        plot_values(ax, sample_df, color="green", s=60, alpha=0.8, marker="x")
        ymax = max(ymax, sample_df['FI'].max())

    if not control_df.empty:
        plot_values(ax, control_df, color="yellow", s=150)
        ymax = max(ymax, control_df['FI'].max())

    plot_standard(ax, st_df, params, color="blue", s=50)
    # plots the params info

    x_pos = st_df['conc'].min()
    if show_params:
        for i,p in enumerate("ABCDE"):
            _ = ax.text(x_pos,ymax - i*ymax/10, f"{p}: {round(params[i],3)}", size=30, color="darkgray")
    
    # other settings for plot
    _ = ax.set_ylim(-200,ymax*1.2)
    _ = ax.set_xscale('log')
    _ = ax.set_xlabel('conc', fontsize=20)
    _ = ax.set_ylabel('FI', fontsize=20)
    _ = plt.xticks(fontsize=20)
    _ = plt.yticks(fontsize=20)

    gene = st_df['Gene'].iloc[0]
    plt.title(f"{gene} | R={round(R,5)}", fontsize=30)
    
    return fig, ax


def plot_multi_s(data_df, col_df, gene="M-CSF", runs=["20211021", "20211222"]):
    
    fig, ax = plt.subplots(figsize=(12,12))
    
    c = {"20211021":"red", "20211222":"green"}
    
    ymax = 0
    for r in runs:
        s = get_standard(data_df, col_df, run=r, gene=gene)
        # fit the curve
        params, _R = fit_standard(s)
        conc, fit = fit_curve(params)
        plot_standard(ax, s, conc, fit, color=c[r])
        ymax = max(ymax, s['FI'].max())
    _ = ax.set_ylim(0,ymax*1.1)
    _ = ax.set_xscale('log')
    _ = ax.set_xlabel('conc', fontsize=20)
    _ = ax.set_ylabel('FI', fontsize=20)
    _ = plt.xticks(fontsize=20)
    _ = plt.yticks(fontsize=20)
    _ = ax.set_xticks([10**i for i in range(-3,6)])
    plt.title(f"{gene} | {'vs '.join(runs)}", fontsize=30)
    
    return