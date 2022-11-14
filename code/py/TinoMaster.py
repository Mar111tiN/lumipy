import pandas as pd
from script_utils import show_output
from TinoMastercols import *


#### Convert TINO data
def get_patient_code(df, org_col="sample code", target_col="PatientCode", remove_org_col=True, prefix="P"):
    '''
    takes that 14.111 like code and transforms it into a valid string
    '''
    
    # flatten
    df.loc[:, target_col] = df[org_col].str.replace(".", "", regex=False).astype(int)
    df.loc[:, target_col] = prefix + ((df[target_col] - df[target_col] % 1000) / 1000).astype(int).astype(str) + "." + (df[target_col] % 1000).astype(str)
    if remove_org_col:
        df = df.drop(org_col, axis=1)
    return df

####
def get_subdata_index(data, sheet, sep="Run"):
    '''
    finds the structure in the data table for the incremental import
    retrieve the coords of the subtables in "all_martin" by detecting the "Run" (or sep) column header
    '''
    # load the df's first column
    df = pd.read_excel(data, sheet_name=sheet, header=None, usecols="A").rename({0:"col"}, axis=1)
    indices = list(df.query('col == @sep').index) + [len(df.index)]
    return indices

def read_subtable(tino_file, all_sheet="", **args):
    '''
    imports all the subfiles
    takes positional args (skiprows, nrows) from consumption of the subdata_index and returns one full data table with some data edits
    '''
    
    # read the df with the respective 
    df2 = pd.read_excel(tino_file, sheet_name=all_sheet, **args).rename({'weight [g]': 'Weight', 'Platte':'Plate'}, axis=1)
    # print(list(df2.columns)[:14])
    df2 = df2.loc[:, [col for col in df2.columns if not col.startswith("Unnamed")]]
    # stack the proteins using melt
    df2 = df2.melt(id_vars=list(df2.columns[:13]), var_name="Protein", value_name="conc").dropna(subset="conc", axis=0)
    # extract the Proteins
    df2.loc[:, ['Gene', 'altGene', 'PlexCol']] = df2['Protein'].str.extract(r"([^(/]+)(?:/([^(/]+))?(?: \(([0-9]+)\))").rename({0:"Gene", 1:"altGene",2:"PlexCol"}, axis=1)
    df2 = df2.drop("Protein", axis=1)    
    return df2
    
    
####### Aggregate TINO DATA from back then
def get_patients_clean(patients_file, project_name="NAC"):
    '''
    get the patients file and extract relevant columns
    '''
    # rename column headers and prepare sample code for transfer
    df = pd.read_excel(patients_file, sheet_name="all", dtype={"sample code": str}).rename({'patient code': "PatientCodeAlt", "Date of Death": "DOT"}, axis=1)
    df.loc[:, 'Project'] = project_name
    
    ### transcode the PatientCode
    #  column sample code is actually a patient code
    # now, it is a mix of formats
    # use the prefix P to make sure its format is adjusted
    df = get_patient_code(df, remove_org_col=True)
    
    # adjust the sex column
    df.loc[df['Female'] == 1, "Sex"] = "F"
    df.loc[df['Male'] == 1, "Sex"] = "M"
    df = df.drop(['Female', 'Male'], axis=1)
    
    ####### SAMPLE NAME #####################
    # edit sample type info into SampleName
    # remove spaces
    for col in ['sample type', 'Serum drawn', 'sentinel/metast', 'tumor locus']:
        df.loc[:, col] = df[col].str.strip(" ")

    # condense sample data into sample name
    df.loc[df['sample type'] == "serum", 'suff'] = "Se"
    df.loc[df['sample type'] == "node", 'suff'] = "P" + df['node #'].fillna(0).astype(int).astype(str)
    df.loc[df['sample type'] == "tumor", 'suff'] = df['tumor locus'].str.replace("na", "T")
    df.loc[df['sample type'] == "normal", 'suff'] = "N"
    df.loc[:, "SampleName"] = df['PatientCode'] + "_" + df['suff']

    # condense extra info (node #, Serum drawn, etc. into SampleType)
    df.loc[df['sample type'] == "serum", 'SampleType'] = "serum " + df['Serum drawn'].fillna("?")
    df.loc[df['sample type'] == "node", 'SampleType'] = "node" + df['node #'].fillna(0).astype(int).astype(str) + " " + df['sentinel/metast'].str.lower().str.replace("non", "non-sentinel")
    df.loc[df['sample type'] == "tumor", 'SampleType'] = "tumor " + df['tumor locus'].str.replace("na", "unknown")
    df.loc[df['sample type'] == "normal", 'SampleType'] = df['sample type']
    return df
    

def import_tino_data(tino_file, sheet=""):
    '''
    import data incrementally and wrangle some data
    cycles through the subtable_index and imports and concatenates the data and wrangles the data a bit
    
    '''
    
    # detect the row coords of the sub tables
    ins = get_subdata_index(tino_file, sheet=sheet)
    show_output(f"Detected header rows in data at following lines: {', '.join([str(i+1) for i in ins[:-1]])}")

    
    # go through all the subtables and import
    dfs = []
    for i in range(len(ins)-1):
        df = read_subtable(tino_file, all_sheet=sheet, skiprows=ins[i], nrows=ins[i+1] - ins[i]-1)
        dfs.append(df)
    
    df = pd.concat(dfs).reset_index(drop=True)
    # print(list(df2.columns))
    # data wrangling

    # streamline the patient data
    df.loc[:, ["Patient","Type", "Extra"]] = df['SampleName'].str.split(" ", expand=True).rename({0:"Patient", 1:"Type", 2:"Extra"}, axis=1)
    df = get_patient_code(df, org_col="Patient").sort_values(['Run', 'PatientCode', 'Type'])
    # make the node names congruent: 5 --> P5 etc. 
    df.loc[df['Type'].str.match(r"^[0-9]+"), 'Type'] = "P" + df['Type'].astype(str)

    # create the SampleName ID
    df.loc[:, "SampleName"] = df['PatientCode'] + "_" + df['Type']
    # add the Project
    df['Project'] = "NAC"
    # store that non-scalar value from Weight in SampleType
    df.loc[df['Weight'].astype(str).str.match(r"[a-zA-Z ]+"), "SampleType"] = df['Weight']
    df.loc[df['SampleType'] == "only blood", "SampleType"] = "blood"
    # make the weight 0 for float conversion
    df.loc[df['Weight'].astype(str).str.match(r"[a-zA-Z ]+"), "Weight"] = "0"
    df.loc[:, "Dilution"] = df['Dilution'].fillna(1).astype(int)
    df.loc[:, "Run"] = df['Run'] - 20000000
    return df.reset_index(drop=True)
    
def makeTinoMaster(tino_patient_excel, lumi_param_excel="", tino_data_excel="", tino_data_sheet="", excel_out=""):
    '''
    wrangle all the data from Tinos data tables to create the TinoMaster file
    '''
    
    # load the tino patients and other required data
    clean_df = get_patients_clean(tino_patient_excel)
    plex_df = pd.read_excel(lumi_param_excel, sheet_name="Proteins").loc[:, ['Protein', 'Plex']]
    
    
    # retrieve the patient columns
    patient_df = clean_df.loc[:, tino_pat_cols].drop_duplicates()
    
    # retrieve the case columns
    case_df = clean_df.loc[:, tino_case_cols].sort_values(['PatientCode', 'tumor surgery']).groupby("PatientCode").first().reset_index()
    case_df.loc[:, 'tumor surgery'] = case_df['tumor surgery'].fillna("unknown")
    case_df = case_df.drop_duplicates().sort_values('PatientCode')  # .drop_duplicates(['PatientCode', 'Date of TURB'])
    
    # retrieve the sample columns
    sample_df = clean_df.loc[:, tino_sample_cols].sort_values(['PatientCode', 'SampleType'])
    
    # import the Tino Data Tables

    df = import_tino_data(tino_data_excel, sheet=tino_data_sheet) # .rename({'SampleName': 'SampleNameAlt', 'Type': 'TypeAlt'}, axis=1)
    
    s2_df = sample_df.merge(df.loc[:, ['SampleName', 'Weight', 'SampleType']].drop_duplicates(['SampleName', 'Weight']).rename({'SampleType': 'SerumDrawn'}, axis=1), on="SampleName", how="outer", indicator=True)
    # message incongruities
    if (len(incon_df := s2_df.query("_merge == 'right_only'").drop_duplicates())):
        show_output(f"There is some incongruity between Patient Data and the data tables!\n {incon_df}", color="warning")
    
    # modify and rename data
    s2_df.loc[s2_df['SerumDrawn'].isin(["nodes", "blood"]), "SampleType"] = s2_df["SampleType"] + s2_df['SerumDrawn']
    sample_df = s2_df.loc[:, [
        'Project',
        'PatientCode',
        'SampleName',
        'SampleType',
        'tumor for analysis',
        'Messrunde',
        'SerumDrawn',
        '_merge'
    ]].rename({
        'tumor for analysis': 'included',
        '_merge':'missingInSamples'
    }, axis=1)

    # modify indicator columns 
    sample_df.loc[:, 'included'] = sample_df['included'].fillna(0).astype(int)
    sample_df.loc[:, 'missingInSamples'] = sample_df['missingInSamples'].astype(str)
    sample_df.loc[sample_df['missingInSamples'] != "right_only", "missingInSamples"] = "0"
    sample_df.loc[sample_df['missingInSamples'] == "right_only", "missingInSamples"] = "1"
    sample_df.loc[:, 'missingInSamples'] = sample_df['missingInSamples'].fillna(0).astype(int)
    
    # assign the proper plex from the plex data
    df = df.drop("Plex", axis=1).rename({'Gene':'Protein'}, axis=1).merge(plex_df)
    well_df = df.loc[:, tino_well_cols].drop_duplicates().sort_values(['Run', 'Plate', 'Plex', 'Well']).reset_index(drop=True)
    # get conc from all conc data
    lumi_df = df.loc[:, [
        'Run',
        'Plate',
        'Plex',
        'Well',
        'Protein',
        'conc'
    ]].sort_values(['Run', 'Plate', 'Plex', 'Well']).reset_index(drop=True)

    # get the qPCR data
    qPCR_df = clean_df.loc[:, qPCR_cols].dropna(subset="CD3")

    # last edits
    patient_df = patient_df.rename({"DOT":"DOD"}, axis=1)
    patient_df['Note'] = ""
    # write to file
    
    if excel_out:
        show_output(f"Writing output to {excel_out}")
        with pd.ExcelWriter(excel_out, mode="w") as writer:
            patient_df.to_excel(writer, sheet_name="Patients", index=False)
            case_df.to_excel(writer, sheet_name="Cases", index=False)
            sample_df.to_excel(writer, sheet_name="Samples", index=False)
            well_df.to_excel(writer, sheet_name="Wells", index=False)
            lumi_df.to_excel(writer, sheet_name="Luminex", index=False)
            qPCR_df.to_excel(writer, sheet_name="qPCR", index=False)
    show_output("Finished data aggregation for TinoData!", color="success")
    return patient_df, case_df, sample_df, well_df, lumi_df, qPCR_df