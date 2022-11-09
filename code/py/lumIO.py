import pandas as pd
from script_utils import show_output
from columns import *


def show_dups(df):
    '''
    return 
    '''
    if (l:=len(dup_df:=df.loc[df.duplicated(keep=False), :])):
        show_output(f"{l} duplicates detected in {len(df)} entries!!", color="warning")
        return dup_df
    else:
        show_output(f"No dups detected in {len(df)} entries!", color="success")


def read_lumi_data(lumi_file):
    '''
    read all the luminex conc data into lumi_dict
    '''

    ld = {}
    for data in ['Plates', 'Standards', 'ProteinStats', 'tidyData', 'tidyDataFull']:
        ld[data] = pd.read_excel(lumi_file, sheet_name=data)
    return ld


def read_MSH_DB(DB_file):
    '''
    read everything into a data dictionary
    '''

    dd = {}
    for data in ['PatientCodes', 'Patients', 'Cases', 'Biopsies', 'Samples', 'Wells']:
        dd[data] = pd.read_excel(DB_file, sheet_name=data)
    return dd


def flatten_notes(df):
    
    note_cols = [col for col in df.columns if col.startswith("Note")]
    new_cols = [col + "_1" for col in note_cols]
    df = df.rename({col:col + "_1" for col in note_cols}, axis=1)
    df['Note'] = ""
    for note_col in new_cols:
        df[note_col] = df[note_col].fillna("")
        df.loc[df['Note'] != df[note_col], 'Note'] = df['Note'] + "|" + df[note_col].fillna("")
        df['Note'] = df['Note'].str.strip("|")
    df = df.drop(new_cols, axis=1)
    return df


def flatten_DB(DB_file):
    '''
    take all levels of the MSHDB and and return the files for the different projects
    '''

    pass

    # load the MSHDB as a data_dictionary
    db_dict = read_MSH_DB(DB_file)

    # merge wells and samples and flatten notes
    dbdf = db_dict['Wells'].merge(db_dict['Samples'], on=["SampleName", "Project"], how="left")
    dbdf = flatten_notes(dbdf.loc[:, merge_col1]).rename({'DOX': 'LumiDOX'}, axis=1)
    # merge with the Biopsies and with the Cases
    dbdf = flatten_notes(dbdf.merge(db_dict['Biopsies'], how="left"))
    dbdf = flatten_notes(dbdf.merge(db_dict['Cases'], how="left"))
    dbdf = flatten_notes(dbdf.merge(db_dict['Patients'], how="left"))