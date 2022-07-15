import os
import matplotlib.pyplot as plt
from plot_fit import plot_confidence_rect, plot_controls, plot_standard, plot_values
from script_utils import show_output
import matplotlib.ticker as mtick

def plot_external(
        ax, prot_standard, data_full_df,
        protein="ABC",
        external_point_size=5,
        external_marker=".",
        plot_zero=0.1,
        external_color="purple",
        external_alpha=1,
        external_connect_lw=1,
        external_connect_alpha=0.3,
        sample_off_marker="x",
        sample_off_size=50,
        sample_off_color="gray",
        sample_off_alpha=0.8,
        **kwargs
    ):
    '''
    print all other values that have been extrapolated
    '''

    # get the data from standard_df and data_full_df
    # get the common Fmin and Fmax
    Fmin = prot_standard['Fmin'].max()
    Fmax = prot_standard['Fmax'].min()
    prot_df = data_full_df.query('Protein == @protein and conc != conc')

    if len(prot_df):
        in_range = (prot_df['FI']>= Fmin) & (prot_df['FI']<= Fmax)
        out_off_range = ~in_range & (prot_df['concMean'] > plot_zero)
        
        _ = ax.scatter(prot_df.loc[in_range, 'concMean'], prot_df.loc[in_range, 'FI'], 
            marker=external_marker,
            fc=external_color,
            color=external_color,
            s=external_point_size,
            alpha=external_alpha,
            )
        _ = ax.scatter(prot_df.loc[out_off_range, 'concMean'], prot_df.loc[out_off_range, 'FI'], 
            marker=sample_off_marker, 
            color=sample_off_color,
            s=sample_off_size,
            alpha=sample_off_alpha,
            )
            # return maximum y value

        # plot the connecting lines
        conc_cols = [f"conc{run}" for run in prot_standard['Run'].unique()]
        # reduce the inrange data to applicable concentrations and get min and max values
        # FI  min  max
        line_df = prot_df.loc[in_range, ['FI'] + conc_cols].set_index('FI').agg(
            ["min", "max"], axis=1
            ).reset_index()
        for i, row in line_df.iterrows():
            _ = ax.plot([row['min'], row['max']], [row['FI'], row['FI']],
            lw=external_connect_lw,
            alpha=external_connect_alpha,
            color=external_color
            )

        return prot_df['FI'].max(), prot_df['concMean'].max()
    return 0, 0

def plot_info_multi(
    ax, standard_df, 
    show_fit_params=False,
    plot_zero=0.1,
    ypos=5000,
    info_font_size=30,
    run_colors={},
    **kwargs
    ):
    '''
    plot the significant values for the 
    '''
    xpos = plot_zero * 2
    i = 0
    for _, standard_row in standard_df.reset_index(drop=True).iterrows():
        Run, R, StQ, C1fit, C2fit = list(standard_row.loc[['Run', 'R^2', 'StMax', 'C1fit', 'C2fit']])
        _ = ax.text(xpos, ypos * (1-i/20), Run, size=info_font_size, color=run_colors[standard_row['Run']])
        info_text = f"R={round(R, 3)} | StQ={round(StQ, 2)}"
        _ = ax.text(xpos, ypos * (1 - (i+1)/20), info_text, size=info_font_size * 0.84, color=run_colors[standard_row['Run']])
        info_text = f"C1|2-fit=[{round(C1fit,3)}|{round(C2fit, 3)}]"
        _ = ax.text(xpos, ypos * (1 - (i+2)/20), info_text, size=info_font_size * 0.84, color=run_colors[standard_row['Run']])
        i += 4
        

def plot_multi(standard_df,data_df,
    protein="ABC",
    plot_folder="",
    plot_type='pdf',
    plot_font_size=20,
    verbose=True,
    plot_confidence=True,
    control_lw=5,
    sc_point_size=150,
    ss_point_size=50,
    plot_zero=0.1, 
    run_colors={},
    **kwargs
    ):
    '''
    calculate samples for a given standard_row (containing all pertaining data sets as ['ss', 'sc', 'data'])
    '''
    
    ######## INIT ##################
    # get the protein
    prot_standard = standard_df.query("Protein == @protein")
    fig, ax = plt.subplots(figsize=(12,12))
    ymax, xmax = 0, 0
    
    for _, standard_row in prot_standard.iterrows():
        if plot_confidence:
            plot_confidence_rect(ax, standard_row, **kwargs)

    for _, standard_row in prot_standard.iterrows():        
        ymax_control, xmax_control = plot_controls(
            ax, standard_row, 
            s=sc_point_size, 
            lw=control_lw, 
            color=run_colors[standard_row['Run']]
            )
        ymax = max(ymax, ymax_control)
        xmax = max(xmax, xmax_control)

    for _, standard_row in prot_standard.iterrows():
        ymax_standards, xmax_standards = plot_standard(
            ax, standard_row, 
            color=run_colors[standard_row['Run']], 
            s=ss_point_size
            )
        ymax = max(ymax, ymax_standards)
        xmax = max(xmax, xmax_standards)
        
    for _, standard_row in prot_standard.iterrows():
        ymax_samples, xmax_samples = plot_values(
            ax, standard_row, 
            **kwargs
            )
        ymax = max(ymax, ymax_samples)
        xmax = max(xmax, xmax_samples)

    plot_info_multi(ax, prot_standard, 
                    plot_zero=plot_zero,
                    ypos=ymax*1.2-500,
                    run_colors=run_colors,
                    **kwargs
                   )
    
    ymax_ext, xmax_ext = plot_external(
        ax, prot_standard, data_df,
        protein=protein,
        plot_zero=plot_zero,
        **kwargs
    )
    ymax = max(ymax, ymax_ext)
    xmax = max(xmax, xmax_ext)
    
    
    # other settings for plot
    _ = ax.set_ylim(-200, ymax*1.2)
    _ = ax.set_xlim(plot_zero - 0.01, xmax * 4)
    _ = ax.set_xscale('log')
    _ = ax.set_xlabel('conc [pg/ml]', fontsize=plot_font_size)
    _ = ax.set_ylabel('FI', fontsize=plot_font_size)
    _ = ax.get_xaxis().set_major_formatter(mtick.ScalarFormatter())
    _ = plt.xticks(fontsize=plot_font_size)
    _ = plt.yticks(fontsize=plot_font_size)

    plt.title(f"{standard_row['Protein']}", fontsize=plot_font_size*1.5)

    if plot_folder:
        # set (and create if neccessary) the fig_plot_path ( = plot_folder/Run) 
        if not os.path.isdir(fig_plot_path := os.path.join(plot_folder, "ProtPlots")):
            os.makedirs(fig_plot_path)
        fig_file_path = os.path.join(fig_plot_path, f"{standard_row['Protein']}.{plot_type}")
        fig.savefig(fig_file_path)
        plt.close()
        if verbose:
            show_output(f"Saving combined plot for {protein} to {'/'.join(fig_file_path.split('/')[-3:])}")
    return fig, ax