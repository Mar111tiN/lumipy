import pandas as pd
import os
import matplotlib.pyplot as plt

from script_utils import show_output
from lumipy_utils import *
from compute_5PL import fit_standard, retro_5PL # only needed for placeholder function load_controls
from plot_fit import plot_fitting

data_cols = ['Run', 'Plex', 'Well', 'Type', 'Gene']


def read_raw_plate(file):
    '''
    reads a Luminex raw data file
    autodetects format
    returns plate_data as series and well_raw data a df
    '''
    
    # get info from filename
    run, plex = get_run_plex(file)
    
    ### header
    # read file depending on extension 
    is_excel = file.split(".")[-1].startswith("xls")
    
    # read header and add 
    header = read_header(file, is_excel)
    header['Run'] = run
    header['Plex'] = plex
    header = header.loc[:, ['Run', 'Plex'] + list(header.columns[:-2])]
    # read the raw data body
    data = pd.read_excel(file, skiprows=7) if is_excel else pd.read_csv(file, skiprows=7, sep=";", encoding = "ISO-8859-1")
    data = data.rename({'Sampling Errors':'SamplingErrors'}, axis=1)

    # adjust the Gene headers
    # get the genes and headers into the col_df
    col_df = pd.DataFrame(data.columns[2:-1])[0].str.extract(r"([^(/]+)/?([^(/]+)? \(([0-9]+)\)").rename({0:"Gene", 1:"altGene",2:"col"}, axis=1)
    cols = col_df.columns
    col_df['Run'] = run
    col_df['Plex'] = plex
    col_df = col_df.loc[:, ['Run', 'Plex'] + list(cols)]
    
    
    # apply the cleaned gene names to the column names
    data.columns = list(data.columns[:2]) + list(col_df['Gene']) + list(data.columns[-1:])
    # stack the data
    data = data.melt(id_vars=['Well', 'Type', 'SamplingErrors'], var_name="Gene", value_name="FI")
    data.loc[:, 'FI'] = data['FI'].str.replace(",", ".").str.replace(r"***", "0", regex=False).astype(float)
    # add run as id
    data['Run'] = run
    data['Plex'] = plex

    raw_data_cols = data_cols + ['FI', 'SamplingErrors']
    data = data.loc[:, raw_data_cols]
    
    # detect if standard has been used
    has_standard = len(data['Type'].str.extract(r"^(S[1-8])$").dropna()[0].unique()) == 8
    header['hasStandard'] = has_standard
    data = data.sort_values(['Run','Gene', 'Plex', 'Well'])
    return header, data, col_df



def read_conc_plate(file):
    '''
    reads from a checkimmune output excel file both the expected concentrations of the respective plex --> col_df
    and the computed values --> conc_df
    '''
    
    # read the plex info
    col_df = read_standard_from_conc(file)
    
    # read the concentration
    conc_df = pd.read_excel(file, skiprows=7, sheet_name="Obs Conc").iloc[1:, :].reset_index(drop=True)
    
    # apply the cleaned gene names to the column names
    conc_df.columns = ['Type', 'Well'] + list(col_df['Gene'])
    # keep only data columns
    conc_df = conc_df.query('Well == Well')
    conc_df = conc_df.loc[~(conc_df['Type'].str.match(r"[eE]?[SC][1-8]")), :].reset_index(drop=True)
    conc_df = conc_df.melt(id_vars=['Well', 'Type'], var_name="Gene", value_name="conc")

    # add run as id
    run, plex = get_run_plex(file)
    
    conc_df = convert2float(conc_df.set_index(['Type', 'Well'])).reset_index()
    cols = conc_df.columns
    conc_df['Run'] = run
    conc_df['Plex'] = plex
    conc_cols = data_cols + ['conc']
    conc_df = conc_df.loc[:, conc_cols]

    # sort_values
    conc_df = conc_df.sort_values(['Run','Gene', 'Plex', 'Well'])
    return conc_df, col_df


def compute_fit(col_df, data_df, gene="M-CSF", run="20211021", apply_standard_from_run=0, fig_path="", zero_value=0.1, dilution=4):
    '''
    master function for computing and applying fits
    if fig_path is set, a figure for the standard curve is plotted
    '''

    # retrieve all the relevant data from the DB
    ss, sc, sample_df = get_data_types(data_df, col_df, run=run, gene=gene, zero_value=zero_value, dilution=dilution)
    
    # load other standard if apply_standard_from_run is set to a value
    if apply_standard_from_run:
        ss = get_standard(data_df, col_df, run=apply_standard_from_run, gene=gene, zero_value=zero_value, dilution=dilution)

    # fit the params
    params, R = fit_standard(ss)
    
    # apply the data point to the col_df if not used from other run
    if not apply_standard_from_run:
        col_df.loc[(col_df['Run'] == run) & (col_df['Gene'] == gene), ['params', 'R^2']] = [" | ".join([str(round(p,3)) for p in params]), R]

    # load the values for the controls
    sc.loc[:, 'conc'] = load_controls(sc, col_df, params)
    
    # compute the values for the samples
    sample_df.loc[:, 'conc'] = retro_5PL(sample_df['FI'], params)

    # apply the computed values to the conc_df
    data_df.loc[(data_df['Run'] == run) & (data_df['Gene'] == gene), 'conc']  = np.round(pd.concat([sample_df, sc, ss]).loc[:, 'conc'],1)
    # apply the computed values to the conc_df for the controls
    # data_df.loc[(data_df['Run'] == run) & (data_df['Gene'] == gene)), 'conc']  = sc.loc[:, 'conc']


    if fig_path:
        fig, ax = plot_fitting(ss, control_df=sc, sample_df=sample_df, params=params, R=R)
        # save output
        fig.savefig(fig_path)
        plt.close()
        show_output(f"Saving fitted data to {'/'.join(fig_path.split('/')[-3:])}")
    return data_df, col_df


def read_luminex(data_folder, raw_pattern="Rawdata", conc_pattern="ISA_conc", output_path=".", excel_out="luminexel", default_run = 20211021, fig_type="svg"):
    '''
    reads all luminex data from one folder
    '''
    
    #### file reading
    # init the file lists
    raw_file_list = []
    
    conc_file_list = []
    
    for f in [folder for folder in os.walk(data_folder)]:
        folder = f[0]
        raw_files = [os.path.join(folder, file) for file in f[2] if raw_pattern in file and not os.path.basename(file).startswith("~$")]
        raw_file_list += raw_files
        
        conc_files = [os.path.join(folder, file) for file in f[2] if conc_pattern in file and not os.path.basename(file).startswith("~$")]
        conc_file_list += conc_files

    # ######### raw files
    # load in all the raw_files and store data in dfs
    raw_col_dfs = []
    data_dfs = []
    plate_dfs = []
    # cycle through raw files
    for file in raw_file_list:
        show_output(f"Loading raw data file {file}")
        plate_df, data_df, raw_col_df = read_raw_plate(file)
        raw_col_dfs.append(raw_col_df)
        plate_dfs.append(plate_df)
        data_dfs.append(data_df)
    
    plate_df = pd.concat(plate_dfs).sort_values(['Run', 'Plex'])
    raw_col_df = pd.concat(raw_col_dfs).sort_values(['Run', 'col'])
    data_df = pd.concat(data_dfs).sort_values(['Run','Gene', 'Plex', 'Type']).reset_index(drop=True)
    
    # ######## conc_files
    conc_col_dfs = []
    conc_dfs = []

    # cycle through conc files
    for file in conc_file_list:
        show_output(f"Loading concentration file {file}")
        conc_df, col_df = read_conc_plate(file)
        conc_col_dfs.append(col_df)
        conc_dfs.append(conc_df)
    
    conc_col_df = pd.concat(conc_col_dfs).sort_values(['Run','col']).loc[:, ['Run', 'Gene', 'col', 'S1']]
    conc_df = pd.concat(conc_dfs).sort_values(['Run','Gene', 'Plex', 'Well']).reset_index(drop=True)    
    
    # check for consistency beween the plexes in raw and conc files
    col_df = raw_col_df.merge(conc_col_df, how="outer")
    col_df.loc[:, 'Run'] = col_df['Run'].astype(int)
    data_df.loc[:, 'Run'] = data_df['Run'].astype(int)
    data_df.loc[~data_df['Type'].str.match(r"^[SC][1-9]$"), 'Type'] = "X"
    
    plate_df.loc[:, 'Run'] = plate_df['Run'].astype(int)
    conc_df.loc[:, 'Run'] = conc_df['Run'].astype(int)
    conc_df = conc_df.rename({'conc':'conc_CI'}, axis=1)
    
    # merge data_df and conc_df for single output
    # expand the wells for duplicates
    conc_df = conc_df.set_index(['Run', 'Plex', 'Gene', 'conc_CI'])['Well'].str.extractall(r"(?P<Well>[A-H][0-9]+)").reset_index().drop('match', axis=1)

    data_df = data_df.merge(conc_df, how="left")
    
    
    #### compute the fit and save the standard curves
    fig_base_folder = os.path.join(output_path, "curves")
    if not os.path.isdir(fig_base_folder):
        os.makedirs(fig_base_folder)

    for run in (runs:= data_df['Run'].unique()):
        fig_folder = os.path.join(fig_base_folder, str(run))
        if not os.path.isdir(fig_folder):
            os.makedirs(fig_folder)
    
        # check if run has standard
        if plate_df.query('Run == @run')['hasStandard'].any():
            fallback_run = 0
            add_string = ""
        else:
            # if not use the last run
            fallback_run = default_run
            add_string = f" using standards from run {fallback_run}"
        for gene in data_df.query('Run == @run')['Gene'].unique():
            show_output(f"Computing params for {gene} in run {run}{add_string}")
            fig_path = os.path.join(fig_folder, f"{run}_{gene}.{fig_type}")
            data_df, col_df = compute_fit(col_df, data_df, gene=gene, run=run, apply_standard_from_run=fallback_run, fig_path=fig_path, zero_value=0.1, dilution=4)

        # plot the combined plots

    # ##### output
    # 
    if excel_out:
        excel_file = os.path.join(output_path, f"{excel_out}.xlsx")
        show_output(f"Writing excel output to {excel_file}")
        with pd.ExcelWriter(excel_file, mode="w") as writer:
            plate_df.to_excel(writer, sheet_name="Plates", index=False)
            col_df.to_excel(writer, sheet_name="Plexes", index=False)
            data_df.to_excel(writer, sheet_name="RawData", index=False)
    
    show_output(f"Finished collecting Luminex data for folder {data_folder}", color="success")
    return plate_df, col_df, data_df