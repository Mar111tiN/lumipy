import pandas as pd

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