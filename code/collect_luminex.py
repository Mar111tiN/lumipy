import pandas as pd
import numpy as np
import os

from script_utils import show_output

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
    cols = data.columns
    # add run as id
    data['Run'] = run
    data['Plex'] = plex
    data = data.loc[:, ['Run', 'Plex'] + list(cols)]
    
    # detect if standard has been used
    has_standard = len(data['Type'].str.extract(r"^(S[1-8])$").dropna()[0].unique()) == 8
    header['hasStandard'] = has_standard
    
    return header, data, col_df

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
    conc_df = conc_df.loc[:, ['Run', 'Plex'] + list(cols)] 
    return conc_df, col_df