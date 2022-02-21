import pandas as pd
import os
from compute_5PL import retro_5PL  # only needed for placeholder function for control values
import numpy as np
from script_utils import show_output
import re


def read_csv_header(csv_file):
    '''
    reads the plate_info from a csv raw data file
    returns series with info
    '''
    
    plate_info = pd.read_csv(csv_file, nrows=6, names=['info'], sep="\t", encoding = "ISO-8859-1")
    plate_info['PlateID'] = plate_info['info'].str.split(": ").str[0]
    plate_info['data'] = plate_info['info'].str.split(": ").str[1].str.rstrip(";")
    plate_info = plate_info.drop('info', axis=1).set_index('PlateID')
    return plate_info['data']


def read_excel_header(excel_file):
    '''
    reads the plate_info from a csv raw data file
    returns series with info
    '''
    
    info = pd.read_excel(excel_file, nrows=6, header=None)
    plate_info = pd.DataFrame()
    plate_info['PlateID'] = info[0].str.split(": ").str[0]
    plate_info['data'] = info[0].str.split(": ").str[1].str.rstrip(";")
    plate_info = plate_info.set_index('PlateID')
    
    return plate_info['data']


def read_header(file, is_excel=False):
    '''
    reads all the info from a luminex header (excel or csv)
    '''

    # read basic data based on extension
    plate_info = read_excel_header(file) if is_excel else read_csv_header(file)
        
    # wrangle the data and add run and plex from file name
    plate_df = plate_info.rename({
        "Plate ID": "orgPlateID",
        "File Name": "RawdataPath",
        "Acquisition Date": "AcquisitionTime",
        "Reader Serial Number": "ReaderID", 
        "RP1 PMT (Volts)": "RP1_PMT",
        "RP1 Target": "RP1_Target"
    }).to_frame().T.reset_index(drop="True")
    
    return plate_df

def get_run_plex(file):
    '''
    retrieve run and flex from the file name
    '''
    s = os.path.basename(file).split("_")
    run = s[0]
    plex = [d for d in s if d.endswith("Plex")][0]
    return run, plex


def read_standard_from_conc(file):
    '''
    reads the expected start concentration for the standards
    '''
    # read_out the standard concentrations
    df = pd.read_excel(file, skiprows=7, sheet_name="Exp Conc")
    # create a Gene df from the columns
    col_df = pd.DataFrame(df.columns[2:])[0].str.extract(r"([^(/]+)/?([^(/]+)? \(([0-9]+)\)").rename({0:"Gene", 1:"altGene",2:"col"}, axis=1)
    col_df['S1'] = df.iloc[1:2,2:].T.reset_index().iloc[:,1].str.replace(",", ".").astype(float)
    cols = col_df.columns
    run, plex = get_run_plex(file)
    col_df['Run'] = run
    col_df['Plex'] = plex
    col_df = col_df.loc[:, ['Run', 'Plex'] + list(cols)]
    return col_df


def convert2float(df):

    for col in ['conc']:
        df.loc[:, col] = df[col].str.replace("---", "-1")
        df.loc[:, col] = df[col].str.replace(",", ".", regex=False)
        df.loc[:, col] = df[col].str.replace("OOR <", "-2", regex=False).str.replace("OOR >", "-1", regex=False)
        df.loc[:, col] = df[col].str.replace("***", "-3", regex=False).str.replace("*", "", regex=False)
        df.loc[:, col] = df[col].astype(float)
    return df


def get_standard(data_df, col_df, run="20211021", gene='M-CSF', dilution=4, zero_value=0.1):
    '''
    returns the standard for that run and the respective gene
    '''
    
    # retrieve only the controls from the raw data
    control_df = data_df.loc[data_df['Type'].fillna("").str.match(r"^[SC][1-9]?"),:]
    
    run = int(run)
    # retrieve standards and control from control_df
    s = control_df.query('Run == @run and Gene == @gene')
    if s.empty:
        show_output(f"No data for Run {run} and Gene {gene}!", color="warning")
        return
    ss = s.loc[s['Type'].str.match(r"^S[1-8]$"),:]
    sc = s.loc[s['Type'].str.match("^C[12]$"),:]
    
    # get the starting concentration for that dilution
    conc = col_df.query('Run == @run and Gene == @gene')['S1'].iloc[0]
    # fill the dilution series with the last being 0
    ss.loc[:, 'conc'] = conc / np.power(dilution, ss.loc[:, 'Type'].str.extract("S([1-8])", expand=False).astype(int)-1)
    ss.loc[ss['Type'] == "S8", 'conc'] = zero_value
    return ss.loc[:, ['Run', 'Gene', 'conc', 'FI']]


def get_data_types(data_df, col_df, run="20211021", gene='M-CSF', dilution=4, zero_value=0.1):
    '''
    returns the standard for that run and the respective gene
    '''
    
    # retrieve only the controls from the raw data
    run = int(run)
    s = data_df.query('Run == @run and Gene == @gene')
    if s.empty:
        show_output(f"No data for Run {run} and Gene {gene}!", color="warning")
        return
    s.loc[:, 'Type'] = s['Type'].fillna("")

    # retrieve standards and control from data_df
    ss = s.loc[s['Type'].str.match(r"^S[1-8]$"),:]
    sc = s.loc[s['Type'].str.match("^C[12]$"),:]
    sx = s.loc[~s['Type'].str.match("^[SC][1-8]$"), :]
    # get the starting concentration for that dilution
    conc = col_df.query('Run == @run and Gene == @gene')['S1'].iloc[0]
    # fill the dilution series with the last being 0
    ss.loc[:, 'conc'] = conc / np.power(dilution, ss['Type'].str.extract("S([1-8])", expand=False).astype(int)-1)
    ss.loc[ss['Type'] == "S8", 'conc'] = zero_value
    
    ss = ss.loc[:, ['Run', 'Gene', 'conc', 'FI']]
    return ss, sc, sx


def load_controls(sc, col_df, params):
    '''
    is supposed to load the known control concentrations
    here used as a place holder with computed values
    '''
    return retro_5PL(sc['FI'], params)


def get_params_from_string(string):
    '''
    little helper
    extracts the params from the luminex curve fit string
    "Std. Curve: FI = 27,863 + (7268,64 - 27,863) / ((1 + (Conc / 101819)^-1,51961))^0,62292"
    --> [27.863, 7268.64, 101819.0, -1.51961, 0.62292]
    auto-detects if params are 4PL and adds a 1 for the 5th param for conformity
    "Std. Curve: FI = 48,6354 + (12224,4 - 48,6354) / (1 + (Conc / 915,213)^-1,26429)"
    --> [48.6354, 12224.4, 915.213, -1.26429, 1]
    '''
    
    nums = [float(f) for f in re.findall(r"-?[0-9]+(?:\.[0-9]+)?", string.replace(",", "."))]
    if len(nums) == 6:
        A,B,_,_,C,D = nums
        E = 1
    elif len(nums) == 7:
        A,B,_,_,C,D,E = nums
    return [A,B,C,D,E]


def get_params_from_col_df(col_df, run="", gene=""):
    '''
    retrieves the params (5PL and R) from the col_df
    '''
    run = int(run)
    params, R = col_df.query("Run == @run and Gene == @gene").iloc[0].loc[['params','R^2']]
    params = [float(p) for p in params.split(" | ")]
    return params,R
    

def get_data_dict(data_df, col_df, run="20211021", gene='M-CSF', dilution=4, zero_value=0.1):
    '''
    provides all the data needed for multi_plotting
    '''
    
    # get the data
    s, c, x = get_data_types(data_df, col_df, run=run, gene=gene, dilution=dilution, zero_value=zero_value)
    
    # get the params from col_df
    params, R = get_params_from_col_df(col_df, run=run, gene=gene)

    # store in dictionary
    data_dict = dict(
        Run=run,
        Gene=gene,
        st=s,
        ctrl=c,
        data=x,
        params=list(params),
        R=R
    )
    return data_dict
