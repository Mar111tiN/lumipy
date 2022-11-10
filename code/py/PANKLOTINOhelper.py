import pandas as pd
from script_utils import show_output
from columns import *
#### utility functions to include all the data from the different data sources to create a 
# combined database file


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



############### 
def expand_plexes(df, plexes=[]):
    '''
    if there are no plexes given, duplicate the well data to all plexes
    '''
    dfs = []
    for plex in plexes:
        plex_df = df.copy()
        plex_df.loc[:, "Plex"] = f"{plex}-Plex"
        dfs.append(plex_df)
    df = pd.concat(dfs).sort_values(['Run', 'Plex', 'SampleName'])
    return df


def lumi21_wells(master_info_excel, plexes21=[3,11,23,38]):
    '''
    turns the master info wells into a df and converts the Run to 6 digit years
    '''
    df = pd.read_excel(master_info_excel, sheet_name="LumiWells2021").rename({'sample_name':'SampleName', 'Weights':"Weight"}, axis=1)
    # fix the date to short form
    df.loc[:, "Run"] = df['Run'] - 20000000
    # add all the plexes
    df = expand_plexes(df, plexes21)
    df.loc[:, ['Plate', 'RunDesc']] = [1, 'LuminexValidation2021']
    # reshuffle columns
    return df.loc[:, ['Project', 'Run', 'RunDesc', 'Plate', 'Plex', 'SampleName', 'DOX', 'Well', 'Weight', 'Dilution', 'Note']].reset_index(drop=True)


def read_excel_96plate(excel_file, anchor="A1", orient="wide", **kwargs):
    '''
    returns for an excel file with a plate setup the well setup in optional tidy format
    only works within the A-X range
    '''
    
    # set skiprows and usecols based on anchor
    l = "ABCDEFGHIJKLMNOPQRSTUVWXY"
    
    # separate col and row in anchor
    row = int(anchor[1:])
    col = anchor[0]
        
    width, height = 12,8
    if orient!="wide":
        width, height = height, width
    # set the "B:N"-like col_range
    col_range = f"{col}:{l[l.find(col) + width]}"
    # sheet_name is passed as kwargs
    df = pd.read_excel(excel_file, **kwargs, skiprows=row-1, nrows=height, usecols=col_range)
    # reset and extract wells
    df = df.set_index(df.columns[0]).stack().reset_index().rename({'Unnamed: 1': "L", 'level_1':"N", 0:"SampleName"}, axis=1)
    df['Well'] = df['L'] + df['N'].astype(str)
    return df.loc[:, ['Well', 'SampleName']]


def lumi22_wells(excel_file, anchor="B4"):
    '''
    returns the plate 
    '''
    dfs = []

    for plate in [1,2]:
        df = read_excel_96plate(excel_file, sheet_name=f"Platte{plate}", anchor=anchor)
        df.loc[:, ['Run', 'RunDesc', 'Plate', 'Plex', 'Note']] = [220906, "Nachmessung2022-20Plex", plate, "21-Plex", ""]
        dfs.append(df)
    df = pd.concat(dfs)
    # clean standards and controls
    df = df.loc[~df['SampleName'].str.match("^[SC][1-9]?$"), :]
    # autodetect project
    df.loc[df['SampleName'].str.startswith("CXCR3 "), 'Project'] = "PankreasCharite"
    df.loc[df['SampleName'].str.startswith("MSH "), 'Project'] = "LO"
    df.loc[df['Project'] != df['Project'], 'Project'] = "NAC"
    df = df.reset_index(drop=True)
    df.loc[df['Project'] == 'NAC', 'SampleName'] = df['SampleName'].str.replace(" [AB]$", "", regex=True).str.replace(" prep1", "").str.replace(" TX", "")
    return df.loc[:, ['Project', 'Run', 'RunDesc', 'Plate', 'Plex', 'SampleName', 'Well', 'Note']]


def load_tino_master(tino_master_excel):
    '''
    load the TinoMasterData and do some data wrangling
    '''

    tino_patients = pd.read_excel(tino_master_excel, sheet_name="Patients").loc[:, ['Project', 'PatientCode', 'PatientCodeAlt', 'DOD', 'Sex', 'Note']]
    # some patients have several entries (with and without extended data --> only use one)
    tino_patients = tino_patients.sort_values(['PatientCode', 'DOD', 'Sex', 'Note']).groupby(['PatientCode', 'PatientCodeAlt']).first().reset_index()


    tino_cases = pd.read_excel(tino_master_excel, sheet_name="Cases").rename({
        "Age at TURB": "Age",
        'tumor surgery': 'TumorSurgery',
        ' Tumor state TURB': 'TumorStateTURB',
        'Tumor grading': 'TumorGrading',
        ' Tumor state RC': 'TumorStateRC',
        "response on NAC": "NACresponse"
        }, axis=1)

    # add the disease to the cases
    tino_cases.loc[:, 'Disease'] = "bladder tumor"
    # create a trans_dictionary for column rename
    trans_dict = {col:col.replace("Date of ", "Do").lstrip(" ").replace(" ", "_") for col in tino_cases.columns if " " in col}
    tino_cases = tino_cases.rename(trans_dict, axis=1)
    # condense the metast data into one
    for met in ['skeletal', 'nodal', 'visceral']:
        tino_cases.loc[tino_cases[f"{met}_Metast"] == 1, 'Metast'] = met
        tino_cases = tino_cases.drop(f"{met}_Metast", axis=1)
    # DoTURB --> DOX
    tino_cases.loc[:, 'DOX'] = tino_cases['DoTURB']
    tino_cases.loc[tino_cases['NAC'] == "NAC", 'Therapy'] = 'NAC:' + tino_cases['NAC_regime'] + "_" + tino_cases['NAC_cycles'].fillna(0).astype(int).astype(str) + "x"
    tino_cases.loc[tino_cases['NAC'] == "NO-NAC", 'Therapy'] = 'NO_NAC'
    tino_cases = tino_cases.drop(['NAC_regime', 'NAC_cycles', 'NAC'], axis=1)

    tino_samples = pd.read_excel(tino_master_excel, sheet_name="Samples")

    tino_wells = pd.read_excel(tino_master_excel, sheet_name="Wells")
    tino_wells.loc[:, "Project"] = "NAC"
    tino_wells.loc[:, "Weight"] = tino_wells['Weight'] * 1000

    # add Note field if not existing
    if not "Note" in tino_patients.columns:
        tino_patients.loc[:, 'Note'] = ""
    if not "Note" in tino_cases.columns:
        tino_cases.loc[:, 'Note'] = ""
    if not "Note" in tino_samples.columns:
        tino_samples.loc[:, 'Note'] = ""
    return tino_patients, tino_cases, tino_samples, tino_wells


def format_TinoSampleName(df, pat_df):
    '''
    reformat all sample names in order to merge them properly
    extracts Tino-specific columns from Tino Sample Names
    '''
    # replace extra spaces 
    df['SampleName'] = df['SampleName'].str.replace(r"[ ]+", " ", regex=True) \
    .str.replace("/", ".") \
    .str.lstrip("0") \
    .str.replace("CP2", "CP") \
    .str.replace("^([0-9]+)", r"P\1", regex=True) \
    .str.replace(r" [0-9]+\.[0-9]$", "", regex=True) \
    .str.replace(" K[0-9]", "", regex=True) \
    .str.replace("Node ", "P") # change P to node
    # extract PatientCode, Type and node from SampleName
    df.loc[df['Project'] == "NAC", ['PatientCode', 'suff', 'node']] = df['SampleName'].str.extract(r"(?P<PatientCode>^[A-Z][-A-Z0-9.]+) (?P<suff>[PTIFC]+[0-9]*)(?: (?P<node>[SN]+))?")
    
    # merge in the patientname
    df = df.merge(pat_df.rename({'PatientCodeAlt':'PatientCode', 'PatientCode':'PatientCodeAlt'}, axis=1).iloc[:, :2], how='left')
    
    df.loc[df['PatientCodeAlt'] == df['PatientCodeAlt'], 'PatientCode'] = df['PatientCodeAlt']
    df.loc[df['Project'] == "NAC", "SampleName"] = df['PatientCode'] + "_" + df['suff']
    return df.drop('PatientCodeAlt', axis=1)


def convert2data(s):
    '''
    converts a data field in numeric format to date-type
    '''

    s = s.fillna(0).astype(int).astype(str)
    s.loc[s != "0"] = s.str[:2] + "-" + s.str[2:4] + "-" + s.str[4:]
    s.loc[s == "0"] = "00-00-00"
    return pd.to_datetime(s, errors="coerce", yearfirst=True)


def merge_2122(df21, df22, pat_df):
    '''
    bring together the weirdly formatted data tables from 2021 and 2022 and merge the 2021 Weight,DOX data into the 2022 data
    '''

    df21 = format_TinoSampleName(df21, pat_df)
    # format SampleName for merging data from 2021 into the 2022 data
    df22 = format_TinoSampleName(df22, pat_df).drop('Note', axis=1).merge(df21.loc[:, ['SampleName', 'DOX', 'Weight', 'Dilution', 'Note']].drop_duplicates())
    df2122 =  pd.concat([df21, df22]).reset_index(drop=True)
    df2122['DOX'] = convert2data(df2122['DOX'])
    return df2122


def flatten_notes(df):
    '''
    takes a df with Note_x and Note_y and combines them into one Note
    '''
    for col in ["Note_x", "Note_y"]:
        df[col] = df[col].fillna("")
    df.loc[(df['Note_x'] ==  "") & (df['Note_y'] ==  ""), 'Note'] = ""
    df.loc[(df['Note'] != "") & (df['Note_x'] ==  df['Note_y']), 'Note'] = df['Note_x']
    df.loc[(df['Note'] != "") & (df['Note_x'] !=  df['Note_y']), 'Note'] = df['Note_x'].fillna("").astype(str) + "|" + df['Note_y'].fillna("").astype(str) .fillna("")
    df.loc[:, 'Note'] = df['Note'].str.strip("|")
    df = df.drop(['Note_x', 'Note_y'], axis=1)
    return df

def mark_missing(df, text="|missing_in_TinoMaster"):
    '''
    detects the files that could not be found in TinoMaster file
    '''
    df.loc[df['_merge'] == "right_only", "Note"] = df['Note'] + text
    df.loc[:, 'Note'] = df['Note'].str.strip("|")
    return df.drop("_merge", axis=1)


def get_tino_pat_codes(tino_samples, tino_patients):
    '''
    creates a patient_code_df correlating the old to the new patient codes
    it finds the ambiguous patient codes and creates a translator
        from ambiguous pat_codes to used patient codes

    >>>
    + get all PatientCodes from the Samples (including the ones I added (pos_df) and keep 'included Field (which biopsy was used for analysis))
    + get all PatientCodes (with old named codes from back then) 
    + merge and check what is missing into Notes --> pat_df
    + add what is missing
    + for Patients with different PatientCodes, use the one that has the "included" tag
    + BEL-36 as an example has two PatientCodes but P17.113 has the "Analyzed" tag and will be used further on
    '''
    # get the deduplicated PatientCodes from all samples
    pos_df = tino_samples.loc[:, ['Project', 'PatientCode', 'included', 'Note']].sort_values(['PatientCode', 'included'], ascending=False).groupby("PatientCode").first().reset_index()
    pos_df.loc[:, 'included'] = pos_df['included'].fillna(0).astype(int)
    # get the PatientCodeAlt (and Notes) from tino_patients and combine into pat_df
    pat_df = tino_patients.loc[:, ['PatientCode', 'PatientCodeAlt', 'Note']].merge(pos_df, how="outer", on='PatientCode', indicator=True)
    pat_df = flatten_notes(pat_df)
    pat_df.loc[pat_df['_merge'] == "right_only", "Note"] = pat_df['Note'] + " | MissingInTinoPatientData!"
    pat_df = pat_df.loc[:, ['Project', 'PatientCode', 'PatientCodeAlt', 'Note', 'included']].sort_values(['PatientCodeAlt', 'included'], ascending=[True, False])
    # merge the PatientCodes onto the PatientCodes deduplicated on PatientCodeAlt (only taking the ones that have the 'included' flag )
    tino_pat_codes = pat_df.merge(pat_df.groupby('PatientCodeAlt').agg({
        'PatientCode':'first'
    }).reset_index().rename({
        'PatientCode':'PatientCodeUsed'
    }, axis=1), how="left")

    # fill NA values for PatientCodeUsed with actual PatientCode
    tino_pat_codes.loc[tino_pat_codes['PatientCodeAlt'] != tino_pat_codes['PatientCodeAlt'], "PatientCodeUsed"] = tino_pat_codes['PatientCode']
    tino_pat_codes = tino_pat_codes.sort_values(['PatientCode']).reset_index(drop=True)

    # fix the columns
    tino_pat_codes = tino_pat_codes.loc[:, patient_code_cols]
    return tino_pat_codes
    

def clean_PatientCode(df, tino_pats):
    '''
    translates the old PatientCode into the disambiguated PatientCode
    '''
    df = df.merge(tino_pats.loc[:, ['PatientCode', 'PatientCodeUsed']]).sort_values("PatientCodeUsed").reset_index(drop=True)
    df.loc[df['PatientCodeUsed'] == df['PatientCodeUsed'], "PatientCode"] = df['PatientCodeUsed']
    return df.drop("PatientCodeUsed", axis=1)


def get_tino_all(lumi21_excel, lumi22_excel, tino_master_excel, plexes21=[3,11,23,38]):
    '''
    get the tino data from TinoMaster and from 2021/22 and merge them
    bring together all tino data and pass the remaining df2122 to the pancreas and LO data analyser
    '''
    # get the 2021/2022 sample info
    df21 = lumi21_wells(lumi21_excel, plexes21=plexes21)
    df22 = lumi22_wells(lumi22_excel)
    tino_patients, tino_cases, tino_samples, tino_wells = load_tino_master(tino_master_excel)

    df2122 = merge_2122(df21, df22, tino_patients)

    tino2122 = df2122.query('Project == "NAC"').sort_values(['PatientCode', 'SampleName', 'Note'])
    # expand SampleType info
    tino2122.loc[tino2122['suff'].str.startswith("P"), "SampleType"] = tino2122['suff'].str.replace("P", "node")
    tino2122.loc[tino2122['node'] == "SN", "SampleType"] = tino2122['SampleType'] + " sentinel"
    tino2122.loc[tino2122['node'] == "NSN", "SampleType"] = tino2122['SampleType'] + " non-sentinel"
    tino2122.loc[tino2122['suff'] == "T", "SampleType"] = "tumor unknown"
    tino2122.loc[tino2122['suff'] == "IF", "SampleType"] = "tumor IF"
    tino2122.loc[tino2122['suff'] == "CP", "SampleType"] = "tumor CP"
    
    #### merge the patients
    # merge the patients
    tino_patients = tino_patients.merge(tino2122.loc[:, ['Project', 'PatientCode', 'Note']].drop_duplicates(['PatientCode']), on=['Project', 'PatientCode'], how="outer", indicator=True)
    
    
    # merge the cases
    tino_cases = tino_cases.merge(tino2122.loc[:, ['Project', 'PatientCode', 'Note']].drop_duplicates(['PatientCode']), on=['Project', 'PatientCode'], how="outer", indicator=True).sort_values(['PatientCode']).reset_index(drop=True)
    
    # merge the samples
    tino_samples = tino_samples.merge(tino2122.loc[:, [
        'Project', 
        'PatientCode',
        'SampleName',
        'SampleType',
        'Note'
    ]].drop_duplicates([
        'PatientCode',
        'SampleName'
    ]
    ), on=[
        'Project',
        'PatientCode',
        'SampleName'
    ], how="outer", indicator=True).sort_values(['PatientCode', 'SampleName']).reset_index(drop=True).rename({'SampleType_x':'SampleType'}, axis=1)
    # tino_samples.loc[:, ]
    
    # get node info from 2122 into tino_samples in case they are missing in tino_samples
    tino_samples.loc[tino_samples['SampleType'] != tino_samples['SampleType'], 'SampleType'] = tino_samples['SampleType_y']
    
    # unify the merged Notes and mark supposed new patients
    tino_patients, tino_cases, tino_samples = [mark_missing(flatten_notes(df)) for df in [tino_patients, tino_cases, tino_samples]]

    # set the columns
    # tino_patients = tino_patients.loc[:, patient_cols]
    tino_cases = tino_cases.loc[:, simple_case_cols + [
       'Metast', 'Therapy', 'Radiotherapy', 'TumorSurgery', 'DoTURB', 'DoRC', 'muscle-invasive', 'NACresponse',
        'TumorStateTURB', 'TumorGrading', 'TumorStateRC'
    ]]
    # create a Tissue field from SampleType
    tino_samples.loc[:, 'Tissue'] = tino_samples['SampleType'].str.split(" ").str[0].str.extract("([A-Za-z]+)", expand=False).str.capitalize().str.replace(r"Tumor|Normal", 'Gall bladder', regex=True)
  
    ### translate all the PatientCodes for disambiguation
    # also, get rid of the PatientCodeAlt
    # extract the patient_code_df from samples and patients
    tino_pat_codes = get_tino_pat_codes(tino_samples, tino_patients)

    # ! The PatientCodes have been unified but I left the SampleNames in original Patient nomenclatur for finding cryosamples
    tino_patients, tino_cases, tino_samples = [clean_PatientCode(df, tino_pat_codes) for df in [tino_patients, tino_cases, tino_samples]]

    # tino_patients have to be condensed for unique patients --> groupby makes collecting duplicate data easier
    tino_patients = tino_patients.sort_values([
        'PatientCode',
        'DOD',
        'Sex']
        ).groupby(['Project', 'PatientCode']).agg({'Sex':'first', 'DOD':'first', 'Note':'sum'}).reset_index(drop=False).loc[:, patient_cols]
    # tino_cases have to be condensed as well --> groupby
    tino_cases = tino_cases.sort_values([
        'PatientCode',
        'TumorSurgery'
    ]).groupby(['Project', 'PatientCode']).first().reset_index(drop=False)
    
    tino_samples = tino_samples.loc[:, sample_biopsy_cols]
    # keep the luminex DOX data for the samples
    swc = sample_well_cols + ["DOX"]
    # concat the wells because they are unique
    tino_wells = pd.concat([tino_wells, tino2122.loc[:, swc]]).loc[:, swc]
    # remove other tino_data from df2122
    df2122 = df2122.query('Project != "NAC"').drop(['PatientCode', 'suff', 'node'], axis=1)
    return df2122, tino_pat_codes, tino_patients, tino_cases, tino_samples, tino_wells

    

############# PANKREAS DATA ########################################

def get_pankreas(df2122, pancreas_excel):
    '''
    + extract pankreas wells from df
    + load data from LuminexMasterInfo
    + check integrity (luminex wells match?)
    + add data to sample_df
    '''
    pank_wells = df2122.query('Project == "PankreasCharite"').sort_values("SampleName")
    samples_pank = pd.read_excel(pancreas_excel, sheet_name="PankreasSamples").rename({'PatientID':'PatientCodeAlt'}, axis=1).drop(["Run", "Well"], axis=1)
    samples_pank['SampleName'] = samples_pank['SampleName'].str.replace("  ", " ")
    ### check integrity
    # compare sample names
    # merge and look for indicator
    if (len(pm := pank_wells.merge(samples_pank, how="outer", indicator=True).query('_merge != "both"').index)):
        show_output("PankWells shows incoherency:", color="warning")
        show_output(pm, color="warning")
        del pm
    else:
        show_output("Pankreas data> Integrity check: Pass!", color="success")

    # merge and extract some info about these samples
    pank_df = pank_wells.merge(samples_pank, how="outer")
    ############################################
    #### do some data extraction and addition of needed info fields

    # remove redundant info from SampleName
    pank_df.loc[:, 'SampleName'] = pank_df['SampleName'].str.replace("CXCR3 ", "").str.replace(" Pankreas", "").str.replace(" Leber", "")
    # store Luminex extraction in LumiDOX
    pank_df.loc[:, "LumiDOX"] = pank_df['DOX']
    # extract the data (of extraction?)
    pank_df.loc[:, 'DOX'] = pd.to_datetime(pank_df['SampleName'].str.extract(r"(2[0-9]{3}-[0-9]{2}-[0-9]{2})")[0])
    pank_df.loc[:, 'SampleName'] = pank_df['SampleName'].str.replace(r"(2[0-9]{3}-[0-9]{2}-[0-9]{2})", "", regex=True)
    # extract patient codes
    pank_df.loc[:, "PatientCode"] = pank_df['SampleName'].str.extract(r"-(?P<PatientCode>[0-9]+)")
    # pank_df.loc[:, 'SampleName'] = pank_df['SampleName'].str.replace(r"-[0-9]+", "", regex=True)
    pank_df.loc[:, "PatientCode"] = "PP" + pank_df['PatientCode'].str.zfill(4)
    pank_df.loc[:, 'PatientCodeUsed'] = pank_df['PatientCode']

    pank_df.loc[:, 'SampleName'] = pank_df['SampleName'].str.strip(" ").str.replace(r"[ ]+", "_", regex=True)
    pank_df = pank_df.sort_values("PatientCodeAlt")
    # extract a sample_type
    pank_df.loc[pank_df['SampleName'].str.find("NG") > -1, 'SampleType'] = 'normal ' +  pank_df['Tissue']
    pank_df.loc[pank_df['SampleName'].str.find("TG") > -1, 'SampleType'] = 'tumor ' +  pank_df['Tissue']

    # add extra data
    pank_df.loc[pank_df['Tissue'] == "Leber", "Tissue"] = "Liver"
    pank_df.loc[:, ['Sex', 'Age', 'DOD']] = 'Unknown', -1, -1
    pank_df['DOD'] = pd.to_datetime(pank_df['DOD'].astype(str), errors="coerce")
      
    # split the data into patients, cases, samples, wells
    pank_pat_codes = pank_df.loc[:, patient_code_cols].drop_duplicates(['PatientCode', 'PatientCodeAlt']).reset_index(drop=True)
    pank_patients = pank_df.copy().loc[:, patient_cols].drop_duplicates('PatientCode').reset_index(drop=True)
    pank_cases = pank_df.loc[:, simple_case_cols].query('Disease != "Normal"').drop_duplicates(['PatientCode', 'DOX']).sort_values("PatientCode").reset_index(drop=True)
    pank_samples = pank_df.loc[:, sample_biopsy_cols].drop_duplicates(['PatientCode', 'SampleName'])

    # keep the DOX data for sample/wells
    pank_wells = pank_df.loc[:, sample_well_cols + ['LumiDOX']].rename({'LumiDOX':"DOX"}, axis=1)
    return pank_pat_codes, pank_patients, pank_cases, pank_samples, pank_wells


def get_LO_data(df2122, LO_excel):
    '''
    get the LungOrganoid data and combine with the well data
    '''
    
    # get the patient/case-relevant data
    LO_patient_df = pd.read_excel(LO_excel, sheet_name="PatientData").rename({'Operation': 'TumorSurgery'}, axis=1)
    # extract and add data
    LO_patient_df['PatientCode'] = LO_patient_df['PatientCode'].str.replace("_", "-")
    LO_patient_df.loc[:, ['Project', 'DOD']] = ["LO", -1]
    LO_patient_df['DOD'] = pd.to_datetime(LO_patient_df['DOD'].astype(str), errors="coerce")
    
    
    for col in [ 'PatientCodeUsed', 'PatientCodeAlt']:
        LO_patient_df.loc[:, col] = LO_patient_df['PatientCode']
    LO_patient_df.loc[:, "Sex"] = LO_patient_df['Sex'].str.upper().str.replace("W", "F")

    # extract the disease from the Pathology
    LO_patient_df.loc[:, ['Disease', 'Pathology']] = LO_patient_df['Pathology'].str.extract(r"(?P<Disease>[A-Z][A-Za-z ]+)(?:\((?P<Pathology>.*))?")
        # squeeze out the pat_codes from patients
    LO_pat_codes = LO_patient_df.copy().loc[:, patient_code_cols]
    LO_patients = LO_patient_df.copy().loc[:, patient_cols]
    LO_cases = LO_patient_df.loc[:, simple_case_cols + ['Therapy', 'Pathology', 'TumorSurgery']]

    
    ## read the data tables for sample meta data
    LO_sample_df = pd.read_excel(LO_excel, sheet_name="dataTable").rename({
        'Tumor from patient': 'PatientCode',
        'name for Checkimmune': 'SampleName',
        'mg in sample':'amount',
        'tissue [mg]': "Weight"
        }, axis=1).query('SampleName != "no serum"')
    # remove space from sample_name
    LO_sample_df['SampleName'] = LO_sample_df['SampleName'].str.replace(r" ([0-9])", r"\1", regex=True).str.replace(" ", "_")
    LO_sample_df['SampleType'] = 'Lung ' + LO_sample_df['Type']
    LO_sample_df.loc[LO_sample_df['Type'] == "serum", 'Weight'] = 0

    LO_sample_df.loc[:, ["Tissue", 'Project']] = ["Lung", "LO"]
    # extract the samples
    LO_samples = LO_sample_df.loc[:, sample_biopsy_cols]
    # load the wells from master info but remove
    LO_wells = df2122.query('Project == "LO"').loc[:, sample_well_cols + ['DOX']]
    LO_wells.loc[:, 'SampleName'] = LO_wells['SampleName'].str.replace(r" ([0-9])", r"\1", regex=True).str.replace(" ", "_")
    # merge the weights into LO_wells
    LO_wells = LO_wells.drop("Weight", axis=1).merge(LO_sample_df.loc[:, ['SampleName', 'Weight']], how="outer").query('Project == Project')

    LO_wells = LO_wells.loc[:, sample_well_cols + ['DOX']]
    
    return LO_pat_codes, LO_patients, LO_cases, LO_samples, LO_wells


def gather_PANKLOTINO(
    lumi21_excel="", 
    lumi22_excel="",
    tino_master_excel="",
    pancreas_excel="",
    LO_excel="",
    plexes21=[3,11,23,38],
    excel_out=""
    ):
    '''
    bring all data together and return individual dfs
    '''


    # collect all the tino and 2021/22 data
    df2122, tino_pat_codes, tino_patients, tino_cases, tino_samples, tino_wells = get_tino_all(lumi21_excel, lumi22_excel, tino_master_excel, plexes21=plexes21)

    pank_pat_codes, pank_patients, pank_cases, pank_samples, pank_wells = get_pankreas(df2122, pancreas_excel)

    LO_pat_codes, LO_patients, LO_cases, LO_samples, LO_wells = get_LO_data(df2122, LO_excel)

    # concat the data
    pat_code_df = pd.concat([tino_pat_codes, pank_pat_codes, LO_pat_codes]).reset_index(drop=True)
    patient_df = pd.concat([tino_patients, pank_patients, LO_patients]).sort_values(['Project', 'PatientCode']).reset_index(drop=True)
    case_df = pd.concat([tino_cases, pank_cases, LO_cases]).sort_values(['Project', 'PatientCode']).reset_index(drop=True)
    sample_df = pd.concat([tino_samples, pank_samples, LO_samples]).sort_values(['Project', 'PatientCode', 'SampleType', 'SampleName']).reset_index(drop=True)
    sample_df.loc[:, 'Tissue'] = sample_df['Tissue'].str.replace('node', 'LymphNode')
    well_df = pd.concat([tino_wells, pank_wells, LO_wells]).sort_values(['Project', 'Run', 'SampleName', 'Well']).reset_index(drop=True)
    for col in ['Plate', 'Run']:
        well_df.loc[:, col] = well_df[col].fillna(1).astype(int)

    ####### create the biopsy
    # turn SampleName into BiopsyName and rename sample_df to biopsy_df
    # add counter value for cumsumming
    # sort by Run data and remove duplicates for SampleName and Weight (to keep the earliest)
    # place the cumsum into that df as indicator
    biopsy_df = sample_df.copy().rename({'SampleName': 'BiopsyName', 'SampleType': 'BioType'}, axis=1)

    sample_df = well_df.copy()
    sample_df.loc[:, "counter"] = 1
    sample_df = sample_df.sort_values("Run").drop_duplicates(['SampleName', 'Weight']).reset_index(drop=True)
    # get the count_info by cumsumming the counter per SampleName (de-duplicated for SampleName+Weight) --> one number per weight
    sample_df['cum'] = sample_df.groupby(['SampleName'])['counter'].cumsum()
    # create BiopsyName field
    sample_df.loc[:, "BiopsyName"] = sample_df['SampleName']
    sample_df.loc[:, "SampleName"] = sample_df['BiopsyName'] + "_#" + sample_df['cum'].astype(str)
    ### split sample_df into well_df and the actual (brand_new) sample_df
    # move the Sample_name into the well_df
    # merge the SampleNames into the well_df based on BiopsyName and Weight
    well_df = well_df.rename({'SampleName': 'BiopsyName'}, axis=1).merge(sample_df.loc[:, ['Weight', 'BiopsyName', 'SampleName']], how="inner")
    well_df = well_df.loc[:, well_cols]

    # create a Method field (can be extended with method ID) to mark the generation from Biopsy (or Sample) to Sample
    sample_df.loc[:, 'Method'] = "ProteinExtraction_Beatblaster"
    sample_df.loc[:, 'SourceUnit'] = "mg"
    # use the Run field as Date of Extration (or first use in that case) as an artificial extration time marker
    sample_df.loc[sample_df['DOX'] != sample_df['DOX'], "DOX"] = convert2data(sample_df['Run'])

    sample_df = sample_df.rename({'Weight':'SourceAmount'}, axis=1).loc[:, sample_cols]

    # make the Note columns explicit
    well_df = well_df.rename({'Note':'Note_Well'}, axis=1).sort_values(['Run', 'Plate','Project', 'SampleName', 'Plex'])
    sample_df = sample_df.rename({'Note':'Note_Sample'}, axis=1).sort_values(['Project', 'BiopsyName', 'SampleName'])
    biopsy_df = biopsy_df.rename({'Note':'Note_Biopsy'}, axis=1).sort_values(['Project', 'PatientCode', 'BiopsyName'])
    case_df = case_df.rename({'Note':'Note_Case'}, axis=1).sort_values(['Project', 'PatientCode', 'DOX'])
    patient_df = patient_df.rename({'Note':'Note_Patient'}, axis=1).sort_values(['Project', 'PatientCode'])
    pat_code_df = pat_code_df.rename({'Note':'Note_PatCode'}, axis=1).sort_values(['Project', 'PatientCode'])

    # convert datetime to date
    patient_df['DOD'] = pd.to_datetime(patient_df['DOD']).dt.date
    for col in ['DOX', 'DoTURB', 'DoRC']:
        case_df[col] = pd.to_datetime(case_df[col]).dt.date
    sample_df['DOX'] = pd.to_datetime(sample_df['DOX']).dt.date
    ### write to file
    if excel_out:
        show_output(f"Data written to {excel_out}!", color="success")
        with pd.ExcelWriter(excel_out, mode="w") as writer:
            pat_code_df.to_excel(writer, sheet_name="PatientCodes", index=False)
            patient_df.to_excel(writer, sheet_name="Patients", index=False)
            case_df.to_excel(writer, sheet_name="Cases", index=False)
            biopsy_df.to_excel(writer, sheet_name="Biopsies", index=False)
            sample_df.to_excel(writer, sheet_name="Samples", index=False)
            well_df.to_excel(writer, sheet_name="Wells", index=False)

    return pat_code_df, patient_df, case_df, biopsy_df, sample_df, well_df




####### Aggregate TINO DATA from back then
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