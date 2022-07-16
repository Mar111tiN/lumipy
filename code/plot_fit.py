import os
import random as rnd
import matplotlib.pyplot as plt
import seaborn
import numpy as np
from matplotlib.patches import Rectangle as Rect
import matplotlib as mpl
seaborn.set()

from compute_5PL import r_squared, PL5
from script_utils import show_output

def r(normalize=True):
    return round(rnd.randint(0,255) / 255, 3) if normalize else rnd.randint(0,255)


def random_color(alpha=None, normalize=True):
    '''
    create a random color with a given 
    '''
    
    color = (r(normalize), r(normalize), r(normalize))
    return color + (alpha,) if alpha else color


def fit_curve(params, plot_zero=0.1, ymax=100000):
    '''
    returns the fit curve for plotting
    '''

    conc = np.power(10, np.linspace(np.log10(plot_zero), np.log10(ymax*1.2)+1, 10000))
    fit = PL5(conc, params)
    return conc, fit
    

def plot_standard(
    ax, standard_row,
    s=50,
    ymax=10000, 
    plot_zero=0.1, 
    **kwargs
    ):
    '''
    subplot the standard curve
    colors are passed through
    '''


    # extract data from standard_row
    params, ss = list(standard_row.loc[['params', 'ss']])
    params = [float(p) for p in params.split(" | ")]
    # adjust the zero_value
    ss.loc[ss['conc'] < plot_zero, 'conc'] = plot_zero

    # fit the curve
    conc, fit = fit_curve(params, plot_zero=plot_zero, ymax=ymax)
    _ = ax.scatter(conc, fit, s=.1, alpha=0.5, **kwargs)
    _ = ax.scatter(ss['conc'], ss['FI'],  s=s, **kwargs)
    # return maximum y value
    return ss['FI'].max(), ss['conc'].max()


def plot_controls(ax, standard_row, 
    s=150,
    lw=5,
    **kwargs
    ):
    '''
    subplot the standard curve
    color in kwargs
    '''
    sc = standard_row['sc']
    if len(sc):
        # fit the curve
        for i in range(len(sc)):
            _ = ax.plot(sc.loc[:,['Cmin', 'Cmax']].iloc[i], sc.loc[:,['FI', 'FI']].iloc[i],
                lw=lw,
                **kwargs
            )
        for c in ['Cmin', 'Cmax']:
            _ = ax.scatter(sc[c], sc['FI'], 
            s=s, 
            **kwargs
            ) #  - np.abs(data_df['Coff']) * 74)
        return sc['FI'].max(), sc['conc'].max()
    return 0, 0


def plot_values(ax, standard_row,
        sample_point_size=80,
        sample_marker="o",
        plot_zero=0.1,
        sample_color="green",
        sample_alpha=1,
        sample_off_marker="x",
        sample_off_size=50,
        sample_off_color="gray",
        sample_off_alpha=0.8,
        **kwargs
        ):
    '''
    subplot the standard curve
    '''

    data_df, Fmin, Fmax = list(standard_row.loc[['data', 'Fmin', 'Fmax']])
    if len(data_df):
        in_range = (data_df['FI']>= Fmin) & (data_df['FI']<= Fmax)
        out_off_range = ~in_range & (data_df['conc'] > plot_zero)
        # plot all computed values 
        _ = ax.scatter(data_df.loc[in_range, 'conc'], data_df.loc[in_range, 'FI'], 
            marker=sample_marker,
            fc='none',
            color=sample_color,
            s=sample_point_size,
            alpha=sample_alpha,
            )
        _ = ax.scatter(data_df.loc[out_off_range, 'conc'], data_df.loc[out_off_range, 'FI'], 
            marker=sample_off_marker, 
            color=sample_off_color,
            s=sample_off_size,
            alpha=sample_off_alpha,
            )
        # return maximum y value
        return data_df['FI'].max(), data_df['conc'].max()
    return 0, 0


def plot_confidence_rect(ax, standard_row, conf_rect_color=[0.701, 0.764, 0.858], conf_rect_alpha=0.5, **kwargs):
    '''
    plots the confidence interval for a given standard_row
    '''

    # extract the needed params
    Fmin, Fmax,Cmin, Cmax = list(standard_row.loc[['Fmin', 'Fmax', 'ConcMin', 'ConcMax']])
    fc = tuple(conf_rect_color) + (conf_rect_alpha,)
    confidence_rect = Rect((Cmin, Fmin), Cmax-Cmin, Fmax-Fmin, fc=fc)
    ax.add_patch(confidence_rect)
    

def plot_info(
    ax, standard_row, 
    show_fit_params=False,
    plot_zero=0.1,
    ypos=5000,
    info_font_size=30,
    info_font_color="darkgray",
    **kwargs
    ):
    '''
    plot the significant values for the 
    '''
    
    xpos = plot_zero * 2
    # extract the required data fields
    params, R, StQ, C1fit, C2fit = list(standard_row.loc[['params', 'R^2', 'StMax', 'C1fit', 'C2fit']])
    params = [float(p) for p in params.split(" | ")]

    # either show params or info field
    if show_fit_params:
        for i,p in enumerate("ABCDE"):
            _ = ax.text(xpos,ypos - i*ypos/10, f"{p}: {round(params[i],3)}", size=info_font_size, color=info_font_color)
    else:
        #### show other marker
        for i, p in enumerate({'R':R, 'StQ':StQ, 'C1fit': C1fit, 'C2fit': C2fit}.items()):
            _ = ax.text(xpos,ypos - i*ypos/10, f"{p[0]}: {round(p[1],5)}", size=info_font_size, color=info_font_color)


def plot_fitting(
    standard_row,
    plot_confidence=True,
    control_color="yellow",
    control_lw=5,
    sc_point_size=150,
    standard_color="blue",
    ss_point_size=50,
    plot_zero=0.1, 
    plot_folder="",
    plot_type='pdf',
    plot_font_size=20,
    verbose=True,
    **kwargs
    ):
    '''
    calculate samples for a given standard_row (containing all pertaining data sets as ['ss', 'sc', 'data'])
    '''
    
    fig, ax = plt.subplots(figsize=(12,12))
    
    if plot_confidence:
        plot_confidence_rect(ax, standard_row, **kwargs)

    ymax, xmax = 0, 0
    ymax_control, xmax_control = plot_controls(
        ax, standard_row, 
        s=sc_point_size, 
        lw=control_lw, 
        color=control_color
        )
    ymax = max(ymax, ymax_control)
    xmax = max(xmax, xmax_control)
    
    ymax_standards, xmax_standards = plot_standard(
        ax, standard_row, 
        color=standard_color, 
        s=ss_point_size
        )
    ymax = max(ymax, ymax_standards)
    xmax = max(xmax, xmax_standards)
    
    ymax_samples, xmax_samples = plot_values(
        ax, standard_row, 
        **kwargs
        )
    ymax = max(ymax, ymax_samples)
    xmax = max(xmax, xmax_samples)
    
    plot_info(
        ax, standard_row,
        ypos=ymax-500,
        plot_zero=plot_zero,
        **kwargs
        )

    # other settings for plot
    _ = ax.set_ylim(-200,ymax*1.2)
    _ = ax.set_xlim(plot_zero - 0.01, xmax * 3)
    _ = ax.set_xscale('log')
    _ = ax.set_xlabel('conc [pg/ml]', fontsize=plot_font_size)
    _ = ax.set_ylabel('FI', fontsize=plot_font_size)
    _ = ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    _ = plt.xticks(fontsize=plot_font_size)
    _ = plt.yticks(fontsize=plot_font_size)

    plt.title(f"{standard_row['Protein']} | {standard_row['Run']}", fontsize=plot_font_size*1.5)

    if plot_folder:
        # set (and create if neccessary) the fig_plot_path ( = plot_folder/Run) 
        if not os.path.isdir(fig_plot_path := os.path.join(plot_folder, standard_row['Run'])):
            os.makedirs(fig_plot_path)
        fig_file_path = os.path.join(fig_plot_path, f"{standard_row['Protein']}.{plot_type}")
        fig.savefig(fig_file_path)
        # plt.close()
        if verbose:
            show_output(f"Saving fitted data to {'/'.join(fig_file_path.split('/')[-3:])}")
    return fig, ax