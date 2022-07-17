import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn
seaborn.set()

from plot_utils import (
    plot_standard, 
    plot_controls, 
    plot_info, 
    plot_confidence_rect, 
    plot_values,
    plot_info_multi, 
    plot_external, 
    log_tick_x
)
from script_utils import show_output


def plot_fitting(
    standard_row,
    plot_confidence=True,
    control_color="yellow",
    control_lw=5,
    figsize=(10,10),
    sc_point_size=150,
    standard_color="blue",
    ss_point_size=50,
    plot_folder="",
    plot_type='pdf',
    plot_font_size=20,
    hide=True,
    verbose=True,
    **kwargs
    ):
    '''
    calculate samples for a given standard_row (containing all pertaining data sets as ['ss', 'sc', 'data'])
    '''
    
    ####### INIT ############################
    fig, ax = plt.subplots(figsize=figsize)
    ymax, xmax, xmin = 0, 0, np.inf


    ####### PLOT CONFIDENCE RECT ############
    if plot_confidence:
        plot_confidence_rect(ax, standard_row, ec='none')

    ####### PLOT CONTROLS ###################
    ymax_control, xmax_control = plot_controls(
        ax, standard_row, 
        s=sc_point_size, 
        lw=control_lw, 
        color=control_color
        )
    # update the coords
    ymax = max(ymax, ymax_control)
    xmax = max(xmax, xmax_control)
    
    ####### PLOT STANDARD S#################
    ymax_standards, xmax_standards, conc_zero = plot_standard(
        ax, standard_row, 
        color=standard_color, 
        s=ss_point_size
        )
    # update the coords
    xmin = min(xmin, conc_zero)
    ymax = max(ymax, ymax_standards)
    xmax = max(xmax, xmax_standards)

    ####### PLOT SAMPLE VALUES #############
    ymax_samples, xmax_samples = plot_values(
        ax, standard_row, 
        **kwargs
        )
    # update the coords
    ymax = max(ymax, ymax_samples)
    xmax = max(xmax, xmax_samples)
    
    ####### PLOT INFO ######################
    plot_info(
        ax, standard_row,
        ypos=ymax*1.1,
        xpos=xmin,
        **kwargs
        )

    ####### SCALING AND LABELS #############
    # other settings for plot
    _ = ax.set_ylim(-200,ymax*1.2)
    _ = ax.set_xlim(xmin / 2, xmax * 3)
    log_tick_x(ax, xmin, xmax * 2, **kwargs)

    _ = ax.set_xlabel('conc [pg/ml]', fontsize=plot_font_size)
    _ = ax.set_ylabel('FI', fontsize=plot_font_size)
    _ = plt.xticks(fontsize=plot_font_size)
    _ = plt.yticks(fontsize=plot_font_size)

    plt.title(f"{standard_row['Protein']} | {standard_row['Run']}", fontsize=plot_font_size*1.5)

    ####### OUTPUT ########################
    if plot_folder:
        # set (and create if neccessary) the fig_plot_path ( = plot_folder/Run) 
        if not os.path.isdir(fig_plot_path := os.path.join(plot_folder, standard_row['Run'])):
            os.makedirs(fig_plot_path)
        fig_file_path = os.path.join(fig_plot_path, f"{standard_row['Protein']}.{plot_type}")
        fig.savefig(fig_file_path)
        if hide:
            plt.close()
        if verbose:
            show_output(f"Saving fitted data to {'/'.join(fig_file_path.split('/')[-3:])}")
    return fig, ax


def plot_multi(standard_df,data_df,
    protein="ABC",
    plot_folder="",
    plot_type='pdf',
    figsize=(10,10),
    plot_font_size=20,
    verbose=True,
    plot_confidence=True,
    control_lw=5,
    sc_point_size=150,
    ss_point_size=50,
    run_colors={},
    hide=True,
    **kwargs
    ):
    '''
    calculate samples for a given standard_row (containing all pertaining data sets as ['ss', 'sc', 'data'])
    '''
    
    ####### INIT ############################
    # get the protein
    prot_standard = standard_df.query("Protein == @protein")
    fig, ax = plt.subplots(figsize=figsize)
    ymax, xmax, xmin = 0, 0, np.inf
    
    ##### cycle through all the plot routines
    # has to be staggered to maintained layers (overlap of confidence rect etc)

    ####### PLOT CONFIDENCE RECT ############
    for _, standard_row in prot_standard.iterrows():
        if plot_confidence:
            plot_confidence_rect(ax, standard_row, ec=run_colors[standard_row['Run']])

    ####### PLOT CONTROLS ###################
    for _, standard_row in prot_standard.iterrows():        
        ymax_control, xmax_control = plot_controls(
            ax, standard_row, 
            s=sc_point_size, 
            lw=control_lw, 
            color=run_colors[standard_row['Run']]
            )
        # update the coords
        ymax = max(ymax, ymax_control)
        xmax = max(xmax, xmax_control)

    ####### PLOT STANDARD S#################
    for _, standard_row in prot_standard.iterrows():
        ymax_standards, xmax_standards, conc_zero = plot_standard(
            ax, standard_row, 
            color=run_colors[standard_row['Run']], 
            s=ss_point_size
            )
        
        # update the coords
        xmin = min(xmin, conc_zero)
        ymax = max(ymax, ymax_standards)
        xmax = max(xmax, xmax_standards)

    ####### PLOT SAMPLE VALUES #############
    for _, standard_row in prot_standard.iterrows():
        ymax_samples, xmax_samples = plot_values(
            ax, standard_row, 
            multi=True,
            **kwargs
            )
        # update the coords
        ymax = max(ymax, ymax_samples)
        xmax = max(xmax, xmax_samples)

    ####### PLOT SAMPLE VALUES #############
    plot_info_multi(
        ax, prot_standard, 
        xpos=xmin,
        ypos=ymax * 1.1,
        run_colors=run_colors,
        **kwargs
        )
    
    ymax_ext, xmax_ext = plot_external(
        ax, prot_standard, data_df,
        protein=protein,
        **kwargs
    )
    ymax = max(ymax, ymax_ext)
    xmax = max(xmax, xmax_ext)
    
    
    ####### SCALING AND LABELS #############
    _ = ax.set_ylim(-200, ymax*1.2)
    _ = ax.set_xlim(xmin / 2, xmax * 3)
    log_tick_x(ax, xmin, xmax * 2, **kwargs)
    _ = ax.set_xlabel('conc [pg/ml]', fontsize=plot_font_size)
    _ = ax.set_ylabel('FI', fontsize=plot_font_size)
    
    _ = plt.xticks(fontsize=plot_font_size)
    _ = plt.yticks(fontsize=plot_font_size)

    plt.title(f"Luminex - {standard_row['Protein']}", fontsize=plot_font_size*1.5)

    ####### OUTPUT ########################
    if plot_folder:
        # set (and create if neccessary) the fig_plot_path ( = plot_folder/Run) 
        if not os.path.isdir(fig_plot_path := os.path.join(plot_folder, "ProtPlots")):
            os.makedirs(fig_plot_path)
        fig_file_path = os.path.join(fig_plot_path, f"{standard_row['Protein']}.{plot_type}")
        fig.savefig(fig_file_path)
        if hide:
            plt.close()
        if verbose:
            show_output(f"Saving combined plot for {protein} to {'/'.join(fig_file_path.split('/')[-3:])}")
    return fig, ax