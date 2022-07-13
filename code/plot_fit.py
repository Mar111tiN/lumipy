import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
import numpy as np
sn.set()

from compute_5PL import r_squared, PL5
from script_utils import show_output

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


def plot_controls(ax, data_df, **kwargs):
    '''
    subplot the standard curve
    '''
    
    # fit the curve
    _ = ax.scatter(data_df['Cmin'], data_df['FI'], **kwargs, s=150) #  - np.abs(data_df['Coff']) * 74)
    _ = ax.scatter(data_df['Cmax'], data_df['FI'], **kwargs, s=150) #  - np.abs(data_df['Coff']) * 74)


def plot_values(ax, data_df, **kwargs):
    '''
    subplot the standard curve
    '''
    
    # fit the curve
    _ = ax.scatter(data_df['conc'], data_df['FI'], **kwargs)

    
def plot_fitting(standard_row, zero_value=0.1, plot_folder="", 
                 show_params=True, 
                 plot_type='pdf',
                 plot_font_size=20, 
                 control_color="yellow",
                 sample_color="green",
                 standard_color="blue",
                 verbose=True,
                 **kwargs
                ):
    '''
    calculate samples for a given standard_row (containing all pertaining data sets as ['ss', 'sc', 'data'])
    '''
    
    # extract the data fields
    params, R, st_df, control_df, sample_df = list(standard_row.loc[['params', 'R^2', 'ss', 'sc', 'data']])
    params = [float(p) for p in params.split(" | ")]
    # adjust the zero_value
    st_df.loc[st_df['conc'] == 0, 'conc'] = zero_value
    fig, ax = plt.subplots(figsize=(12,12))
    
    # calculate R if R is not given
    if not R:
        R = r_squared(params, st_df)

    ymax = st_df['FI'].max()
    
    if not sample_df.empty:
        plot_values(ax, sample_df, color=sample_color, s=60, alpha=0.8, marker="x")
        ymax = max(ymax, sample_df['FI'].max())

    if not control_df.empty:
        plot_controls(ax, control_df, color=control_color)
        ymax = max(ymax, control_df['FI'].max())

    plot_standard(ax, st_df, params, color=standard_color, s=50)
    # plots the params info

    x_pos = st_df['conc'].min()
    if show_params:
        for i,p in enumerate("ABCDE"):
            _ = ax.text(x_pos,ymax - i*ymax/10, f"{p}: {round(params[i],3)}", size=30, color="darkgray")
    
    # other settings for plot
    _ = ax.set_ylim(-200,ymax*1.2)
    _ = ax.set_xscale('log')
    _ = ax.set_xlabel('conc', fontsize=plot_font_size)
    _ = ax.set_ylabel('FI', fontsize=plot_font_size)
    _ = plt.xticks(fontsize=plot_font_size)
    _ = plt.yticks(fontsize=plot_font_size)

    plt.title(f"{standard_row['Protein']} | {standard_row['Run']} | R={round(R,5)} | StQ={standard_row['StMax']}", fontsize=plot_font_size*1.2)

    if plot_folder:
        # set (and create if neccessary) the fig_plot_path ( = plot_folder/Run) 
        if not os.path.isdir(fig_plot_path := os.path.join(plot_folder, standard_row['Run'])):
            os.makedirs(fig_plot_path)
        fig_file_path = os.path.join(fig_plot_path, f"{standard_row['Protein']}.{plot_type}")
        fig.savefig(fig_file_path)
        plt.close()
        if verbose:
            show_output(f"Saving fitted data to {'/'.join(fig_file_path.split('/')[-3:])}")
    return fig, ax



def plot_multi(data_dict_list=[dict(
        Run="",
        Gene="",
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

    for data_dict in data_dict_list:
        s = data_dict['data']
        # extract the control, standard and samples from the data_df
        st_df = s.loc[s['Type'].str.match(r"^S[1-8]$"),:]
        sample_df = s.loc[~s['Type'].str.match("^[SC][1-8]$"), :]
        control_df = s.loc[s['Type'].str.match("^[C][12]$"), :]
        if not sample_df.empty:
            plot_values(ax, sample_df, color=data_dict['color'], s=80, alpha=0.8, marker="x")
            ymax = max(ymax, sample_df['FI'].max())
        if not control_df.empty:
            plot_values(ax, control_df, color="yellow", s=150)
            ymax = max(ymax, control_df['FI'].max())
        if not st_df.empty:
            plot_standard(ax, st_df, data_dict['params'], color=data_dict['color'], s=50)
            ymax = max(ymax, st_df['FI'].max())
            # plots the params info
        x_pos = min(x_pos, st_df['conc'].min())

    # add this list only after ymax has been determined
    if show_info:
        for data_dict in data_dict_list:
            R_string = data_dict['R'] if isinstance(data_dict['R'], str) else f"R={round(data_dict['R'],5)}"
            _ = ax.text(x_pos,ymax*1.09 - i*ymax/12, f"Run {data_dict['Run']}: {R_string}", size=18, color=data_dict['color'])
            i += 1
    
    # other settings for plot
    _ = ax.set_ylim(-200,ymax*1.2)
    _ = ax.set_xscale('log')
    _ = ax.set_xlabel('conc', fontsize=20)
    _ = ax.set_ylabel('FI', fontsize=20)
    _ = plt.xticks(fontsize=20)
    _ = plt.yticks(fontsize=20)

    protein = data_dict_list[0]['Protein']
    plt.title(f"{protein}", fontsize=30)
    
    return fig, ax