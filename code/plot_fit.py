import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
import numpy as np
sn.set()

from compute_5PL import r_squared, PL5


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
    _ = ax.scatter(st_df['conc'], st_df['FI'], color=color, **kwargs)

    
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


def plot_multi(data_dict_list=[dict(
        Run="",
        Gene="",
        st=pd.DataFrame(),
        ctrl=pd.DataFrame(),
        data=pd.DataFrame(),
        params=[],
        R=0,
        color="black"
    )], zero_value=0.1,
        show_info=True
              ):
    '''
    calculate samples for a given set of data_dictionaries
    '''

    
    # init the plot
    fig, ax = plt.subplots(figsize=(12,12))
    # init ymax and xpos and i for text info
    ymax,x_pos,i = (0,np.Inf,0)

    for data in data_dict_list:
        if not data['data'].empty:
            plot_values(ax, data['data'], color=data['color'], s=80, alpha=0.8, marker="x")
            ymax = max(ymax, data['data']['FI'].max())
        if not data['ctrl'].empty:
            plot_values(ax, data['ctrl'], color="yellow", s=150)
            ymax = max(ymax, data['ctrl']['FI'].max())

        plot_standard(ax, data['st'], data['params'], color=data['color'], s=50)
        # plots the params info
        x_pos = min(x_pos, data['st']['conc'].min())
        if show_info:
            i += 1
            print(x_pos,ymax - i*ymax/10)
            _ = ax.text(x_pos,ymax - i*ymax/10, f"Run {data['Run']}: R={round(data['R'],5)}", size=30, color=data['color'])
    
    # other settings for plot
    _ = ax.set_ylim(-200,ymax*1.2)
    _ = ax.set_xscale('log')
    _ = ax.set_xlabel('conc', fontsize=20)
    _ = ax.set_ylabel('FI', fontsize=20)
    _ = plt.xticks(fontsize=20)
    _ = plt.yticks(fontsize=20)

    gene = data_dict_list[0]['Gene']
    plt.title(f"{gene}", fontsize=30)
    
    return fig, ax