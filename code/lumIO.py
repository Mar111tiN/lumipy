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
