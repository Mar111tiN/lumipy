import os
import numpy as np
import pandas as pd
from script_utils import show_output

def get_data_types(data_df, col_df, run="20211021", protein='M-CSF', dilution=4, zero_value=0.1):
    '''
    returns the standard for that run and the respective gene
    '''
    
    # retrieve only the controls from the raw data
    run = int(run)
    s = data_df.query('Run == @run and Protein == @protein')
    if s.empty:
        show_output(f"No data for Run {run} and Protein {protein}!", color="warning")
        return
    s.loc[:, 'Type'] = s['Type'].fillna("")

    # retrieve standards and control from data_df
    ss = s.loc[s['Type'].str.match(r"^S[1-8]$"),:]
    sc = s.loc[s['Type'].str.match("^C[12]$"),:]
    sx = s.loc[~s['Type'].str.match("^[SC][1-8]$"), :]
    # get the starting concentration for that dilution
    conc = col_df.query('Run == @run and Protein == @protein')['S1'].iloc[0]
    # fill the dilution series with the last being 0
    ss.loc[:, 'conc'] = conc / np.power(dilution, ss['Type'].str.extract("S([1-8])", expand=False).astype(int)-1)
    ss.loc[ss['Type'] == "S8", 'conc'] = zero_value

    return ss, sc, sx



def get_standard(data_df, col_df, run="20211021", protein='M-CSF', dilution=4, zero_value=0.1):
    '''
    returns the standard for that run and the respective protein
    '''
    
    # retrieve only the controls from the raw data
    control_df = data_df.loc[data_df['Type'].fillna("").str.match(r"^[SC][1-9]?"),:]
    
    run = int(run)
    # retrieve standards and control from control_df
    s = control_df.query('Run == @run and Protein == @protein')
    if s.empty:
        show_output(f"No data for Run {run} and Protein {protein}!", color="warning")
        return
    ss = s.loc[s['Type'].str.match(r"^S[1-8]$"),:]
    sc = s.loc[s['Type'].str.match("^C[12]$"),:]
    
    # get the starting concentration for that dilution
    conc = col_df.query('Run == @run and Protein == @protein')['S1'].iloc[0]
    # fill the dilution series with the last being 0
    ss.loc[:, 'conc'] = conc / np.power(dilution, ss.loc[:, 'Type'].str.extract("S([1-8])", expand=False).astype(int)-1)
    ss.loc[ss['Type'] == "S8", 'conc'] = zero_value
    return ss


def load_controls(sc, col_df):
    '''
    loads the range of known control concentrations
    and adds as cols Cmean and Cbound (upper lower spread)
    '''
    
    range_df = col_df.merge(sc.loc[:, ['Run', 'Protein']]).drop_duplicates().loc[:, ['C1', 'C2']].T[0].str.extract(r"(?P<Cmin>[0-9]+(?:\.[0-9])?) ?[-–] ?(?P<Cmax>[0-9]+(?:\.[0-9])?)").astype(float)
    # after the first round, Cmin and Cmax will have been added and need to be remove lest (Cmin_x, Cmin_y) be created
    sc = sc.drop(['Cmin', 'Cmax'], axis=1, errors="ignore").merge(range_df.loc[:, ['Cmin', 'Cmax']], left_on="Type", right_index=True)
    return sc



def fit4protein(data_df, col_df, protein="M-CSF", run="20211021", apply_standard_from_run=0, fig_path="", zero_value=0.1, dilution=4, confidence=0.9, plot=True):
    '''
    master function for computing and applying fits
    if fig_path is set, a figure for the standard curve is plotted
    '''

    # retrieve all the relevant data from the DB
    ss, sc, sample_df = get_data_types(data_df, col_df, run=run, protein=protein, zero_value=zero_value, dilution=dilution)
    # load other standard if apply_standard_from_run is set to a value
    if apply_standard_from_run:
        ss = get_standard(data_df, col_df, run=apply_standard_from_run, protein=protein, zero_value=zero_value, dilution=dilution)

    # fit the params
    params, R = fit_standard(ss)
    
    # apply the data point to the col_df if not used from other run
    if not apply_standard_from_run:
        col_df.loc[(col_df['Run'] == run) & (col_df['Protein'] == protein), ['params', 'R^2']] = [" | ".join([str(round(p,3)) for p in params]), R]

    # load the values for the controls
    sc = load_controls(sc, col_df)
    # get the Cmin and Cmax from the params for sample_df
    sample_df.loc[:, ["Cmin", "Cmax"]] = get_confidence(params, fraction=confidence)
    ss.loc[:, ["ConcMin", "ConcMax"]] = get_confidence(params)
    sample_df = pd.concat([sample_df, sc, ss])
    sample_df = compute_samples(sample_df, params)
    # apply the computed values to the data_df
    data_df.loc[(data_df['Run'] == run) & (data_df['Protein'] == protein), ['conc', 'Cmin', 'Cmax', 'Coff']]  = np.round(sample_df.loc[:, ['conc', 'Cmin', 'Cmax',  'Coff']],2)
    # apply the maximum Coff from the standards to the col_df as SCoffMax
    # this is a measure of the reach of the maximal standard concentrations
    # SCoffMax < 0.6 mean the sigmoidal curve is largely extrapolated 
    col_df.loc[(col_df['Run'] == run) & (col_df['Protein'] == protein), "SCoffMax"] = sample_df.loc[sample_df['Type'].str.match(r"[eE]?[S][1-8]"), 'Coff'].max().round(2)
    col_df.loc[(col_df['Run']  == run) & (col_df['Protein'] == protein), ["C1CoffMean", "C2CoffMean"]] = list(sample_df.loc[sample_df['Type'].str.startswith("C"), :].groupby('Type')['Coff'].mean().round(2))
    # create the data_dict for export to multiplotting
    data_dict = dict(
        Run=run,
        Protein=protein,
        data=sample_df,
        params=params,
        R=R
    )
    
    if plot:
        fig, ax = plot_fitting(data_dict)
        # save output
        fig.savefig(fig_path)
        plt.close()
        show_output(f"Saving fitted data to {'/'.join(fig_path.split('/')[-3:])}")
    return data_df, col_df, data_dict



def get_params_from_col_df(col_df, run="", protein=""):
    '''
    retrieves the params (5PL and R) from the standard_df
    '''
    run = int(run)
    params, R = col_df.query("Run == @run and Protein == @protein").iloc[0].fillna("").loc[['params','R^2']]
    params = [float(p) for p in params.split(" | ")] if params else []
    if not R:
        R = "No standard used!"
    return params,R



def fit4all(data_df, col_df, plate_df, output_path=".", zero_value=0.1, dilution=4, confidence=0.9, plot_fit=True, plot_combined_graphs=True, run_colors={}, apply_standard_from_run=0, fig_type="svg"):
    '''
    takes care of visualization for all the proteins for all the runs
    compute the fit and save the standard curves
    '''
    
    # create figure folder
    if plot_fit or plot_combined_graphs:
        fig_base_folder = os.path.join(output_path, "curves")
        if not os.path.isdir(fig_base_folder):
            os.makedirs(fig_base_folder)    
        
    # init ddld ("data_dict_list_dict")
    # ddld is the dictionary of genes and the respective lists of data_dicts for multiplotting
    # will be used afterwards for the plotting function as an iterable
    ddld= {protein:[] for protein in data_df['Protein'].unique()}
    for run in (runs:= data_df['Run'].unique()):
        if plot_fit:
            fig_folder = os.path.join(fig_base_folder, str(run))
            if not os.path.isdir(fig_folder):
                os.makedirs(fig_folder)

        # check if run has standard
        if plate_df.query('Run == @run')['hasStandard'].any():
            fallback_run = 0
            add_string = ""
        else:
            # if not use the last run
            fallback_run = apply_standard_from_run
            add_string = f" using standards from run {fallback_run}"
        for protein in data_df.query('Run == @run')['Protein'].unique():
            show_output(f"Computing params for {protein} in run {run}{add_string}")
            fig_path = os.path.join(fig_folder, f"{run}_{protein}.{fig_type}") if plot_fit else ""
            data_df, col_df, data_dict = fit4protein(data_df, col_df, protein=protein, run=run, apply_standard_from_run=fallback_run, fig_path=fig_path, zero_value=zero_value, dilution=dilution, confidence=confidence, plot=plot_fit)
            # add the gene-run data to the ddld
            data_dict['color'] = run_color[run]
            if fallback_run: 
                data_dict['R'] = f"Standard from {fallback_run} used!"
            ddld[protein].append(data_dict)
            
    # plot the combined plots
    if plot_combined_graphs:
        for protein in data_df['Protein'].unique():
            # check for the runs that actually contain data for this data point
            fig, ax = plot_multi(ddld[protein], show_info=True)
            fig_path = os.path.join(fig_base_folder, f"{protein}.{fig_type}")
            fig.savefig(fig_path)
            plt.close()
            
    return data_df, col_df