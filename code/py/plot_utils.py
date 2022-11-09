import numpy as np
import random as rnd
import math
import matplotlib.ticker as mtick
from matplotlib.patches import Rectangle as Rect

from compute_5PL import r_squared, PL5


def r(normalize=True):
    return round(rnd.randint(0,255) / 255, 3) if normalize else rnd.randint(0,255)


def random_color(alpha=None, normalize=True):
    '''
    create a random color with a given 
    '''
    
    color = (r(normalize), r(normalize), r(normalize))
    return color + (alpha,) if alpha else color


def fit_curve(params, xmin=0.1, ymax=100000):
    '''
    returns the fit curve for plotting
    '''

    conc = np.power(10, np.linspace(np.log10(xmin), np.log10(ymax), 10000))
    fit = PL5(conc, params)
    return conc, fit
    

def plot_standard(
    ax, standard_row,
    s=50,
    zero_log_dist=.5,
    alpha=1,
    **kwargs
    ):
    '''
    subplot the standard curve
    colors are passed through
    '''

    # extract data from standard_row
    params, ss = list(standard_row.loc[['params', 'ss']])
    # copy to keep standard_df['ss'] immutable
    ss = ss.copy()
    params = [float(p) for p in params.split(" | ")]
    ### adjust the zero_value
    # get the minimum conc above the blank
    min_conc = ss.drop_duplicates('Type').reset_index()['conc'].iloc[-2]
    # take the next lower 10x level (with a dist)
    zero_conc = math.pow(10,math.floor(np.log10(min_conc)-zero_log_dist))
    ss.loc[ss['conc'] == 0, 'conc'] = zero_conc
    max_conc = ss['conc'].max()
    max_FI = ss['FI'].max()
    ### fit the curve
    # plot one decade more than needed
    conc, fit = fit_curve(params, xmin=zero_conc, ymax=max_conc*10)
    _ = ax.scatter(conc, fit, s=.1, alpha=0.5, **kwargs)
    _ = ax.scatter(ss['conc'], ss['FI'],  s=s, alpha=alpha, **kwargs)
    # plot the zero-plot a bit nicer
    _ = ax.scatter(
        ss.loc[ss['conc'] == zero_conc, 'conc'],
        ss.loc[ss['conc'] == zero_conc, 'FI'],
        s=s*15, alpha=alpha*0.1,**kwargs
    )
    # return maximum y value
    return max_FI, max_conc, zero_conc


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
        multi=False,
        **kwargs
        ):
    '''
    subplot the standard curve
    '''

    if multi:
        sample_off_color=sample_color

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


def plot_confidence_rect(
    ax, standard_row, 
    conf_rect_color=[0.701, 0.764, 0.858],
    conf_rect_alpha=0.5,
    **kwargs
    ):
    '''
    plots the confidence interval for a given standard_row
    '''

    # extract the needed params
    Fmin, Fmax,Cmin, Cmax = list(standard_row.loc[['Fmin', 'Fmax', 'ConcMin', 'ConcMax']])
    fc = tuple(conf_rect_color) + (conf_rect_alpha,)
    confidence_rect = Rect((Cmin, Fmin), Cmax-Cmin, Fmax-Fmin, fc=fc, **kwargs)
    ax.add_patch(confidence_rect)
    

def plot_info(
    ax, standard_row, 
    show_fit_params=False,
    xpos=0.1,
    ypos=5000,
    info_font_size=30,
    info_font_color="darkgray",
    **kwargs
    ):
    '''
    plot the significant values for the 
    '''
    
    # extract the required data fields
    params, R, StQ, C1fit, C2fit = list(standard_row.loc[['params', 'R^2', 'StMax', 'C1fit', 'C2fit']])
    params = [float(p) for p in params.split(" | ")]

    # either show params or info field
    if show_fit_params:
        for i,p in enumerate("ABCDE"):
            _ = ax.text(xpos,ypos - i*ypos/15, f"{p}: {round(params[i],3)}", size=info_font_size, color=info_font_color)
    else:
        #### show other marker
        for i, p in enumerate({'R':R, 'StQ':StQ, 'C1fit': C1fit, 'C2fit': C2fit}.items()):
            _ = ax.text(xpos,ypos - i*ypos/15, f"{p[0]}: {round(p[1],5)}", size=info_font_size, color=info_font_color)


def plot_external(
        ax, prot_standard, data_full_df,
        protein="ABC",
        external_point_size=5,
        external_marker=".",
        xmin=0.1,
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
        out_off_range = ~in_range & (prot_df['concMean'] > xmin)

        _ = ax.scatter(prot_df.loc[out_off_range, 'concMean'], prot_df.loc[out_off_range, 'FI'], 
            marker=sample_off_marker, 
            color=external_color,
            s=sample_off_size,
            alpha=sample_off_alpha,
            )
            # return maximum y value
        # only plot, if there is something in range left
        if sum(in_range):
            _ = ax.scatter(prot_df.loc[in_range, 'concMean'], prot_df.loc[in_range, 'FI'], 
                marker=external_marker,
                fc=external_color,
                color=external_color,
                s=external_point_size,
                alpha=external_alpha,
                )
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
    xpos=0.1,
    ypos=5000,
    info_font_size=30,
    run_colors={},
    **kwargs
    ):
    '''
    plot the significant values for each standard_row
    '''

    i = 0
    for _, standard_row in standard_df.reset_index(drop=True).iterrows():
        Run, R, StQ, C1fit, C2fit = list(standard_row.loc[['Run', 'R^2', 'StMax', 'C1fit', 'C2fit']])
        _ = ax.text(xpos, ypos * (1-i/20), Run, size=info_font_size, color=run_colors[standard_row['Run']])
        info_text = f"R={round(R, 3)} | StQ={round(StQ, 2)}"
        _ = ax.text(xpos, ypos * (1 - (i+1)/20), info_text, size=info_font_size * 0.84, color=run_colors[standard_row['Run']])
        info_text = f"C1|2-fit=[{round(C1fit,3)}|{round(C2fit, 3)}]"
        _ = ax.text(xpos, ypos * (1 - (i+2)/20), info_text, size=info_font_size * 0.84, color=run_colors[standard_row['Run']])
        i += 4


def fmt(num):
    if num < 1:
        return str(num)
    else:
        return str(int(num))
    
def log_tick_steps(xmin=1, xmax=1000, step=1):
    '''
    return the required positions for log-ticks
    '''
    
    logmin=math.floor(np.log10(xmin))
    logmax=math.ceil(np.log10(xmax))

    major_ticks = [np.power(10.,l) for l in range(logmin, logmax)]
    minor_ticks = [x*mt for mt in major_ticks for x in range(step,10,step)]
    # now add the last major_tick
    major_ticks.append(np.power(10.,logmax))
    major_tick_labels = [fmt(label) for label in major_ticks]

    return major_ticks, minor_ticks, major_tick_labels
    

def log_tick_x(
    ax, xmin=1, xmax=1000, log_tick_step=1,
    minor_grid_lw= 0.5,
    minor_grid_alpha= 0.5,
    **kwargs
):
    '''
    sets the log_ticks for x-coord from min and max values
    '''
    ax.set_xscale('log')
    ax.xaxis.grid(True, which="both")
    major_ticks, minor_ticks, major_tick_labels = log_tick_steps(xmin,xmax,log_tick_step)
    ax.get_xaxis().set_ticks(major_ticks, which="major")
    ax.get_xaxis().set_ticks(minor_ticks, minor=True)
    ax.tick_params(which="major", width=1, length=8, bottom=True, left=True)
    ax.tick_params(which="minor", width=1, length=4, bottom=True, 
                   grid_linewidth=minor_grid_lw, grid_alpha=minor_grid_alpha, grid_linestyle="-")
    # ax.get_xaxis().set_tick_params(which="major", bottom=True)
    # ax.get_xaxis().set_tick_params(which="minor", top=True)
    if xmax < 500000:
        ax.get_xaxis().set_major_formatter(mtick.FixedFormatter(major_tick_labels))
    return ax