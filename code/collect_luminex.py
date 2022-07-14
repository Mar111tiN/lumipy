import pandas as pd
import os
import matplotlib.pyplot as plt

from script_utils import show_output
from lumipy_utils import *
from compute_5PL import *
from plot_fit import plot_fitting, plot_multi


run_color={20211021:"green", 20211102:"orange", 20211222:"brown"}


def fit_standard_row(standard_row, data_df=pd.DataFrame(), **fit_config):
    '''
    for every row of the standards the fit_params and other coefficients
    pass kwargs for standards and plotting
    '''
    
    # filter for the protein
    s = data_df.loc[data_df['Protein'] == standard_row['Protein'], :]  
    
    standard_row = analyse_standard(standard_row, s, **fit_config)
    # get the confidence information about the control samples
    standard_row = analyse_control(standard_row, s)
    # get the samples
    sx = s.loc[s['Type'] == "X", :]
    standard_row['data'] = compute_conc(sx, standard_row)
    
    return standard_row


def apply_external_standards(data_df, standard_df, fit_config):
    '''
    computes concentrations for all samples for all available standards
    '''

    for _, standard_row in standard_df.iterrows():
        data_df = compute_external_fit(standard_row, df=data_df)
    
    # now get the summary statistics
    sum_cols = ['concMean', 'concStd', 'FposMean']
    data_df.loc[:, sum_cols] = data_df.apply(mean_row, standard_df=standard_df, axis=1, minFpos=fit_config['minFpos'])
    return data_df


def make_protein_summary(data_df, minFpos=0.1, minFI=2000):
    '''
    creates summary statistics for proteins
    '''

    sum_df = data_df.groupby("Protein", as_index=False).agg(
        FposMean=pd.NamedAgg("Fpos", "mean"),
        goodFposfrac=pd.NamedAgg("Fpos", lambda prot: np.sum(prot > minFpos) / len(prot)),
        FIMean=pd.NamedAgg("FI", "mean"),
        goodFIfrac=pd.NamedAgg("FI", lambda prot: np.sum(prot > minFI) / len(prot)),
        FposExtMean=pd.NamedAgg("FposMean", "mean"),
        goodFposExtfrac=pd.NamedAgg("FposMean", lambda prot: np.sum(prot > minFpos) / len(prot))
    )

    return sum_df


def read_raw_plate(plate, control_df, config={}):
    '''
    reads a Luminex raw data file
    autodetects format
    returns: 
        plate_info as series with information about the plate
        standard_df with computed fit params
        raw_df containing raw FI and computed concentrations
    '''

    # set the raw_file with the data_path
    raw_file = os.path.join(config['data_path'], plate['rawPath'])
    ### GET HEADER ##########################
    # read file depending on extension 
    is_excel = raw_file.split(".")[-1].startswith("xls")
    show_output(f"Run {plate['Run']} {plate['Plex']} Plate{plate['Plate']}: Loading raw data file {plate['rawPath']}.")
    # read header and add 
    plate = pd.concat([plate,read_header(raw_file, is_excel)])
    # read the raw data body
    if is_excel:
        data_df = pd.read_excel(raw_file, skiprows=7)
    else:
        data_df = pd.read_csv(raw_file, skiprows=7, sep=";", encoding = "ISO-8859-1")
    data_df = data_df.rename({'Sampling Errors':'SE'}, axis=1)
    # clean the Type from 
    data_df.loc[:, 'Type'] = data_df['Type'].str.strip(" ").fillna("")
     # detect if standard has been used
    # has_standard = len(raw_df['Type'].str.extract(r"^(S[1-8])$").dropna()[0].unique()) > 2
    # !!! some samples have Type S1 to S20 or so which are NOT standards --> so here is the fix:
    standard_length = len(data_df['Type'].str.extract(r"^(S[1-9]+)$").dropna()[0])
    plate['hasStandard'] = (standard_length > 2 and standard_length < 17)
    # assign every non-standard Type (not control or standard or blank) an Type="X"
    if plate['hasStandard']:
        data_df.loc[~data_df['Type'].str.match(r"^[BSC][1-9]?$"), 'Type'] = "X"
    else:
         data_df.loc[~data_df['Type'].str.match(r"^[BC][12]?$"), 'Type'] = "X"
 
    # hard fix for inconsistency in PDFG-AA / PDGFAA
    data_df.columns = [col.replace("PDGF-", "PDGF") for col in data_df.columns]

    # adjust the Gene headers
    # get the Plex setup into the standard_df and merge with control_df containing the Luminex params
    standard_df = pd.DataFrame(data_df.columns[2:-1]).rename({0:"PlexName"}, axis=1).merge(control_df, how="left")
    
    # check for consistency in Luminex Params
    if (len(miss_prots := standard_df.loc[standard_df['Protein'] != standard_df['Protein'], 'PlexName'])):
        miss_list = '; '.join([f'<{prot}>' for prot in miss_prots])
        show_output(f"The following Proteins were not found in the Luminex Params [{miss_list}]", color="warning")
        return

    # apply the cleaned gene names to the column names
    data_df.columns = list(data_df.columns[:2]) + list(standard_df['Protein']) + list(data_df.columns[-1:])



    # tidy the data into Well-Type-SamplingErrors
    raw_df = data_df.melt(id_vars=['Well', 'Type', 'SE'], var_name="Protein", value_name="FI")
    raw_df.loc[:, 'FI'] = raw_df['FI'].str.replace(",", ".").str.replace(r"***", "0", regex=False).astype(float)
    # add RunPlexPlate and reorder cols
    data_cols = list(raw_df.columns)
    base_cols = ['Run', 'Plex', 'Plate']
    raw_df.loc[:, base_cols] = list(plate.loc[base_cols])
    raw_df = raw_df.loc[:, base_cols + data_cols]

    # detect if control has been included
    plate['hasControl'] = len(raw_df.loc[raw_df['Type'].str.match(r"^C[12]$"), :].index) > 0
    # detect if there is any data
    plate['DataWells'] = (has_data := len(data_df.query('Type == "X"')))
    if not has_data:
        if config['verbose']:
            show_output(f"Plate {plate['rawPath']} has no data!", color = "warning")


    # has_standard = False
    if plate['hasStandard']:
        # add RunPlexPlate and reorder cols
        # else no output to standard_df is needed
        standard_cols = list(standard_df.columns)
        standard_df.loc[:, base_cols] = list(plate.loc[base_cols])
        standard_df = standard_df.loc[:, base_cols + standard_cols]
        # calculate the standard fit and add to standard_df
        standard_df = standard_df.apply(fit_standard_row, data_df=raw_df, axis=1, **config['fitting'])
        if config['plot_fit']:
            plot_config = config['plotting']
            standard_df.apply(plot_fitting, axis=1, **plot_config, verbose=config['verbose'])
            
        # extract the data from the standard_df into raw_df
        raw_df = pd.concat(list(standard_df['data']))
        standard_df = standard_df.drop(['ss', 'sc', 'data'], axis=1)
    else:
        # return no standard_df if there is no standard_data
        # standard_df = None
        if config['verbose']:
            show_output(f"Plate {plate['rawPath']} has no standards!", color = "warning")
        if has_data:
            # at least get the data - nothing else to do, already stored in raw_df
            pass
        standard_df = pd.DataFrame()

    raw_df = raw_df.sort_values(['Run','Plex', 'Type', 'Protein', 'Well'])
    
    return plate, standard_df, raw_df


def read_conc_plate(plate, control_df, config={}):
    '''
    reads from a checkimmune output excel file the computed values
    and adds them to the data_df in the standard_row
    '''

    # output empty dataframe with respective columns for compatibility
    conc_file = os.path.join(config['data_path'], plate['concPath'])
    if conc_file:
        # read the concentration
        show_output(f"Run {plate['Run']} {plate['Plex']} Plate{plate['Plate']}: Loading precomputed concentrations from {plate['concPath']}")
        conc_df = pd.read_excel(conc_file, skiprows=7, sheet_name="Obs Conc").iloc[1:, :].reset_index(drop=True).rename(
            {'Unnamed: 0': 'Type', 'Unnamed: 1': 'Well'}, axis=1).dropna(subset="Well")
    else:
        # return an empty df if there is no conc file
        show_output(f"Run {plate['Run']} {plate['Plex']} Plate{plate['Plate']}: No concentration file found!", color="warning")
        return pd.DataFrame()
    # remove external standards (eS1...)
    conc_df = conc_df.loc[~(conc_df['Type'].str.match(r"[eE][SC][1-8]")), :].reset_index(drop=True)
    
    # pivot the data
    conc_df = conc_df.melt(id_vars=['Well', 'Type'], var_name="PlexName", value_name="concCI")
    # adjust the column names from control_df
    conc_df = conc_df.merge(control_df.loc[:, ['PlexName', 'Protein']]).loc[:, ['Well', 'Protein', 'concCI']]
    # convert all the marked values to allow float-cast
    conc_df = convert2float(conc_df, cols=['concCI'])
    org_cols = list(conc_df.columns)
    base_cols = ['Run', 'Plex', 'Plate']
    conc_df.loc[:, base_cols] = list(plate.loc[base_cols])

    return conc_df.loc[:, base_cols + org_cols]


def read_luminex_folder(analysis_name="results", config_file={}, **kwargs):
    '''
    read all the luminex data from one folder
    
    '''
    
    base_cols = ['Run', 'Plex', 'Plate']
    data_cols = base_cols + ['Well', 'Type', 'Protein']
    
    
    # load the config_file and pass extra kwargs
    show_output(f"Loading Luminex configs from {config_file}.")
    config = load_lumi_config(config_file=config_file, analysis_name=analysis_name, **kwargs)

    show_output(f"Starting analysis of Luminex data folder {config['data_path']}.", color="success")
    # load in the Luminex params
    params_file = config['params_file']
    if not params_file:
        show_output("You have to provide a Luminex params file [params_file=""]!!! Aborting", color="warning")
        return
    if not os.path.isfile(params_file):
        show_output(f"Luminex params file {params_file} cannot be found!!! Aborting", color="warning")
        return
    
    control_df = pd.read_excel(params_file, sheet_name="Controls").merge(pd.read_excel(params_file, sheet_name="Proteins").loc[:, ['PlexName', 'Protein']]).loc[:, ['PlexName', 'Protein', 'C1', 'C2', 'S1']]
    
    # create the output file and check for existing
    excel_file = os.path.join(config['analysis_folder'], "luminexcel.xlsx")
    csv_file = os.path.join(config['analysis_folder'], "luminextern.csv.gz")
    # create this empty old_data marker
    use_old = 0
    
    if (old_file :=config['use_file']):
        if os.path.isfile(old_file):
            old_data = load_existing(old_file)
            use_old = 1
        else:
            show_output(f"{old_file} for appending new data cannot be found! Please check!", color="warning")
            return
    else:
        if os.path.isfile(excel_file) and config['use_existing']:
            old_data = load_existing(excel_file)
            use_old = 2

    # load the plates and remove the duplicates from old runs
    plate_df = get_luminex_plates(config['data_path'], raw_pattern=config['raw_pattern'])
    
    if use_old:
        # only use the new plates (drop duplicates)
        old_plates = old_data['Plates'].loc[:, base_cols + ['rawPath', 'concPath']]
        plate_df = pd.concat([plate_df, old_plates]).drop_duplicates(keep=False).sort_values(base_cols)

        ############ Pseudo-duplicates ##############
        # this stuff is quite complicated as it allows duplicate data to be created which is kinda hard to solve later
        # get the count for pseudo_duplicates
        plate_df.loc[:, ['count']] = plate_df.groupby(base_cols)['rawPath'].transform(lambda x: x.count())
        # extract the removers_df containing base_cols of pseudoduplicates
        remover_df = plate_df.loc[plate_df['count'] >1, base_cols].drop_duplicates()
        old_data['Plates'] = old_data['Plates'].merge(remover_df, how="outer", indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)
        old_data['Standards'] = old_data['Standards'].merge(remover_df, how="outer", indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)
        old_data['tidyData'] = old_data['tidyData'].merge(remover_df, how="outer", indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)

        # adjust pseudo_duplicate lines if only one sample has been added (collect data for both entries)
        plate_df = plate_df.fillna("").groupby(base_cols).agg({'rawPath': max, 'concPath':max}).reset_index()

        if not len(plate_df.index):
            show_output(f"No new data found in {config['data_path']}. Exiting!", color="success")
            # only recompute the protein stats in case something has been added in code
            sum_df = make_protein_summary(old_data['tidyData'], **config['summary'])
            return old_data['Plates'], old_data['Standards'], sum_df, old_data['tidyData']

    ########### READ IN RAW DATA ##############################################
    # load in all the raw_files and store data in dfs
    plate_rows = []
    standard_dfs = []
    data_dfs = []
    # cycle through files using df.iterrows
    for _, plate in plate_df.iterrows():
        if not plate['rawPath']:
            show_output(f"Run {plate['Run']} {plate['Plex']} Plate{plate['Plate']}: No raw data file detected. Skipping {plate['concPath']}", color="warning")
            continue
        plate, standard_df, data_df = read_raw_plate(plate, control_df, config=config)
        if plate['concPath']:
            conc_df = read_conc_plate(plate, control_df, config=config)
            data_df = data_df.merge(conc_df, how="left")
        plate_rows.append(plate)
        standard_dfs.append(standard_df)
        data_dfs.append(data_df)
    
    # combine to dfs
    if len(plate_rows):
        plate_df = pd.DataFrame(plate_rows)
    else:
        show_output(f"No new data found in {config['data_path']}. Exiting!", color="success")
        if use_old:
            print("HELLO")
            return old_data['Plates'], old_data['Standards'], old_data['ProteinStats'], old_data['tidyData']
        # really nothing there
        else:
            return [pd.DataFrame()] * 4
    # maybe no new standard has been added
    if len(standard_df.index):
        standard_df = pd.concat(standard_dfs).sort_values(base_cols + ['Protein']).drop_duplicates(base_cols + ['Protein'])
    data_df = pd.concat(data_dfs).sort_values(['Run', 'Plex', 'Protein', 'Type', 'Well']).reset_index(drop=True)


    ############ ADD EXTERNAL STANDARDS  ################################
    show_output("Computing concentrations from external standards")
    data_cols = list(data_df.columns)
    sum_cols = ['concMean', 'concStd', 'FposMean']
    data_df = apply_external_standards(data_df, standard_df, config['fitting'])
    # reduce the data_df to fewer output
    data_full = data_df.copy()
    data_df = data_df.loc[:, data_cols + sum_cols]

    ############ OUTPUT #################################################
    # ##### output
    if config['write_excel']:
        if use_old:
            # add the new stuff to the old sheets and sort again
            if use_old == 1:
                show_output(f"Combining preexisting data from {config['use_file']} and new data to {excel_file}")
            if use_old == 2:
                show_output(f"Adding new data to {excel_file}")
            plate_df = pd.concat([old_data['Plates'], plate_df]).sort_values(base_cols)
            standard_df = pd.concat([old_data['Standards'], standard_df]).sort_values(base_cols + ['Protein']).drop_duplicates()
            data_df = pd.concat([old_data['tidyData'], data_df]).sort_values(data_cols).drop_duplicates()
        else:
            show_output(f"Writing excel output to {excel_file}")

        ############ PROTEIN SUMMARY ################
        # protein summary should be performed on combined data
        sum_df = make_protein_summary(data_df, **config['summary'])

        with pd.ExcelWriter(excel_file, mode="w") as writer:
            plate_df.to_excel(writer, sheet_name="Plates", index=False)
            standard_df.to_excel(writer, sheet_name="Standards", index=False)
            sum_df.to_excel(writer, sheet_name="ProteinStats", index=False)
            data_df.to_excel(writer, sheet_name="tidyData", index=False)
            if config['output_untidy']:
                for col in ['FI', 'conc', 'concCI', 'Fpos']:
                    set_cols = ['Run', 'Plex', 'Plate', 'Well', 'Type', 'SE']
                    pivot_df = data_df.set_index(set_cols).pivot(columns="Protein", values=col).loc[:, list(control_df['Protein'])].dropna(how="all").reset_index(drop=False)
                    pivot_df.to_excel(writer, sheet_name=col, index=False)
    show_output(f"Writing complete external conc file output to {csv_file}")
    data_full.to_csv(csv_file, index=False, sep="\t", compression="gzip")
    show_output(f"Finished collecting Luminex data for folder {config['data_path']}", color="success")
    return plate_df, standard_df, sum_df, data_df        