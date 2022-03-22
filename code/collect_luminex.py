import pandas as pd
import os
import matplotlib.pyplot as plt

from script_utils import show_output
from lumipy_utils import *
from compute_5PL import fit_standard, get_confidence, compute_samples
from plot_fit import plot_fitting, plot_multi

data_cols = ['Run', 'Plex', 'Well', 'Type', 'Gene']
run_color={20211021:"green", 20211102:"orange", 20211222:"brown"}

def read_raw_plate(plate, control_df):
    '''
    reads a Luminex raw data file
    autodetects format
    returns plate_data as series and well_raw data a df
    '''
    
    # get info from filename
    run, plex = plate['run'], plate['plex']
    
    ### header
    # read file depending on extension 
    is_excel = plate['raw'].split(".")[-1].startswith("xls")
    show_output(f"Run {plate['run']} {plate['plex']}: Loading raw data file.")
    # read header and add 
    header = read_header(plate['raw'], is_excel)
    header['Run'] = run
    header['Plex'] = plex
    header['FolderPath'] = plate['raw']
    header = header.loc[:, ['Run', 'Plex', 'FolderPath'] + list(header.columns[:-3])]
    # read the raw data body
    data_df = pd.read_excel(plate['raw'], skiprows=7) if is_excel else pd.read_csv(plate['raw'], skiprows=7, sep=";", encoding = "ISO-8859-1")
    data_df = data_df.rename({'Sampling Errors':'SamplingErrors'}, axis=1)

    # adjust the Gene headers
    # get the Plex setup into the col_df and merge with control_df containing the Luminex params
    col_df = pd.DataFrame(data_df.columns[2:-1]).rename({0:"PlexName"}, axis=1).merge(control_df, how="left")
    
    # check for consistency in Luminex Params
    if (len(miss_prots := col_df.loc[col_df['Protein'] != col_df['Protein'], 'PlexName'])):
        miss_list = '; '.join([f'<{prot}>' for prot in miss_prots])
        show_output(f"The following Proteins were not found in the Luminex Params [{miss_list}]", color="warning")
        return
    
    cols = col_df.columns
    col_df['Run'] = run
    col_df['Plex'] = plex
    col_df = col_df.loc[:, ['Run', 'Plex'] + list(cols)]
    
    # apply the cleaned gene names to the column names
    data_df.columns = list(data_df.columns[:2]) + list(col_df['Protein']) + list(data_df.columns[-1:])
    data_df = data_df.melt(id_vars=['Well', 'Type', 'SamplingErrors'], var_name="Protein", value_name="FI")
    data_df.loc[:, 'FI'] = data_df['FI'].str.replace(",", ".").str.replace(r"***", "0", regex=False).astype(float)
    
    # add info and adjust columns
    data_cols = list(data_df.columns)
    # add run as id
    data_df['Run'] = run
    data_df['Plex'] = plex
    raw_data_cols = ['Run','Plex'] + data_cols
    data_df = data_df.loc[:, raw_data_cols]
    
    # detect if standard has been used
    has_standard = len(data_df['Type'].str.extract(r"^(S[1-8])$").dropna()[0].unique()) == 8
    header['hasStandard'] = has_standard
    data_df = data_df.sort_values(['Run','Plex', 'Type', 'Protein', 'Well'])
    return header, data_df, col_df


def read_conc_plate(plate, col_df):
    '''
    reads from a checkimmune output excel file the computed values --> conc_df
    '''

    # output empty dataframe with respective columns for compatibility
    if plate['conc']:
        # read the concentration
        show_output(f"Run {plate['run']} {plate['plex']}: Loading precomputed concentrations from {plate['conc']}")
        conc_df = pd.read_excel(plate['conc'], skiprows=7, sheet_name="Obs Conc").iloc[1:, :].reset_index(drop=True)
    else:
        # return the empty df
        show_output(f"Run {plate['run']} {plate['plex']}: No concentration file found!", color="warning")
        return pd.DataFrame(columns=['Run','Plex', 'Protein', 'Well', 'conc'])
    
    # apply the cleaned gene names to the column names
    conc_df.columns = ['Type', 'Well'] + list(col_df['Protein'])
    # keep only data columns
    conc_df = conc_df.query('Well == Well')
    # remove external standards (eS1...)
    conc_df = conc_df.loc[~(conc_df['Type'].str.match(r"[eE][SC][1-8]")), :].reset_index(drop=True)
    conc_df = conc_df.melt(id_vars=['Well', 'Type'], var_name="Protein", value_name="conc")

    # add run as id
    run, plex = plate['run'], plate['plex']
    
    conc_df = convert2float(conc_df)
    org_cols = list(conc_df.columns)
    conc_df['Run'] = run
    conc_df['Plex'] = plex
    conc_cols = ['Run','Plex'] + org_cols
    conc_df = conc_df.loc[:, conc_cols]

    # sort_values
    conc_df = conc_df.sort_values(['Run','Plex', 'Protein', 'Well'])
    return conc_df


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
    ss.loc[:, ["Cmin", "Cmax"]] = get_confidence(params)
    sample_df = pd.concat([sample_df, sc, ss])
    sample_df = compute_samples(sample_df, params)
    # apply the computed values to the data_df
    data_df.loc[(data_df['Run'] == run) & (data_df['Protein'] == protein), ['conc', 'Cmin', 'Cmax', 'Coff']]  = np.round(sample_df.loc[:, ['conc', 'Cmin', 'Cmax',  'Coff']],2)
    # apply the maximum Coff from the standards to the col_df as SCoffMax
    # this is a measure of the reach of the maximal standard concentrations
    # SCoffMax < 0.6 mean the sigmoidal curve is largely extrapolated 
    col_df.loc[(col_df['Run'] == run) & (col_df['Protein'] == protein), "SCoffMax"] = sample_df.loc[sample_df['Type'].str.match(r"[eE]?[S][1-8]"), 'Coff'].max().round(2)
    col_df.loc[(col_df['Run'] == run) & (col_df['Protein'] == protein), ["C1CoffMean", "C2CoffMean"]] = list(sample_df.loc[sample_df['Type'].str.startswith("C"), :].groupby('Type')['Coff'].mean().round(2))
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


def read_luminex_folder(data_folder, raw_pattern="Rawdata", conc_pattern="ISA_conc", params_file="", output_path=".", excel_out="luminexel", default_run=20211021, fig_type="svg", plot_fit=True, plot_combined_graphs=True):
    '''
    reads all luminex data from one folder
    '''
    # !!!!
    # adjust run colors more dynamically
    
    run_color={20211021:"green", 20211102:"orange", 20211222:"brown"}
    
    # load in the Luminex params
    if params_file == "":
        show_output("You have to provide a Luminex params file [params_file=""]!!! Aborting", color="warning")
        return
    control_df = pd.read_excel(params_file, sheet_name="Controls").merge(pd.read_excel(params_file, sheet_name="Proteins").loc[:, ['PlexName', 'Protein']]).loc[:, ['PlexName', 'Protein', 'C1', 'C2', 'S1']]
    
    # check if excel_file already exists and 
    if excel_out:
        excel_file = os.path.join(output_path, f'{excel_out.replace(".xlsx", "").replace("xls", "")}.xlsx')
        # create has_excel flag for use in plate_list and below for saving to file
        use_file = excel_file if (has_excel := os.path.isfile(excel_file)) else ""
    #### file reading
    plate_list = get_luminex_plates(data_folder, use_file=use_file, raw_pattern=raw_pattern, conc_pattern=conc_pattern)
    if not(len(plate_list)):
        show_output(f"No new data found in {data_folder}!", color="warning")
        if excel_out:
            return read_luminexcel(excel_file)
        else:
            return pd.DataFrame(), pd.DataFrame(),pd.DataFrame()
    # ######### raw files
    # load in all the raw_files and store data in dfs
    col_dfs = []
    data_dfs = []
    plate_dfs = []

    # cycle through raw files
    for plate in plate_list:
        plate_df, data_df, col_df = read_raw_plate(plate, control_df)
        conc_df = read_conc_plate(plate, col_df)
        # expand the wells for duplicates
        conc_df = conc_df.set_index(['Run', 'Plex', 'Protein', 'conc'])['Well'].str.extractall(r"(?P<Well>[A-H][0-9]+)").reset_index().drop('match', axis=1)
        # merge data_df and conc_df for single output
        data_df = data_df.merge(conc_df.loc[:, ['Well', 'Protein', 'conc']], how="left", on=['Well', 'Protein'])
        
        col_dfs.append(col_df)
        plate_dfs.append(plate_df)
        data_dfs.append(data_df)
         
    # combine to dfs
    plate_df = pd.concat(plate_dfs).sort_values(['Run', 'Plex'])
    # remove the absolute path from the FolderPath for the raw_plates to keep it compatible to moving the data
    plate_df.loc[:, 'FolderPath'] = plate_df['FolderPath'].str.replace(f"{data_folder}/", "")

    col_df = pd.concat(col_dfs).sort_values(['Run', 'Plex', 'Protein']).drop_duplicates(['Run', 'Plex', 'Protein'])
    data_df = pd.concat(data_dfs).sort_values(['Run', 'Plex', 'Protein', 'Type', 'Well']).reset_index(drop=True)

    for df in [col_df, data_df, plate_df]:
        df.loc[:, 'Run'] = df['Run'].astype(int)
    
    data_df.loc[~data_df['Type'].str.match(r"^[SC][1-9]$"), 'Type'] = "X"
    
    # assign concentration as coming from CheckImmune as "CI"
    data_df = data_df.rename({'conc':'conc_CI'}, axis=1)
    data_df, col_df = fit4all(data_df, col_df, plate_df, output_path=output_path, apply_standard_from_run=default_run, fig_type=fig_type, plot_fit=plot_fit, plot_combined_graphs=plot_combined_graphs)
    
           
    # ##### output
    if excel_out:
        # check whether the file already exists
        if has_excel:
            # add the new stuff to the old sheets and sort again
            show_output(f"Adding new data to {excel_file}")
            plate_df_old, col_df_old, data_df_old = read_luminexcel(excel_file)
            plate_df = pd.concat([plate_df_old, plate_df]).sort_values(['Run', 'Plex'])
            col_df = pd.concat([col_df_old, col_df]).sort_values(['Run', 'Plex', 'Protein']).drop_duplicates(['Run', 'Plex', 'Protein'])
            data_df = pd.concat([data_df_old, data_df]).sort_values(['Run', 'Plex', 'Protein', 'Type', 'Well']).drop_duplicates(['Run', 'Plex', 'Protein', 'Well'])
        else:
            show_output(f"Writing excel output to {excel_file}")
        with pd.ExcelWriter(excel_file, mode="w") as writer:
            plate_df.to_excel(writer, sheet_name="Plates", index=False)
            col_df.to_excel(writer, sheet_name="Plexes", index=False)
            # Cmin and Cmax should be dropped for actual output
            data_df.to_excel(writer, sheet_name="RawData", index=False)       
    show_output(f"Finished collecting Luminex data for folder {data_folder}", color="success")
    return plate_df, col_df, data_df          

