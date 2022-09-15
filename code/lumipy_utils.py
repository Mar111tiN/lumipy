import os
import re
import numpy as np
import pandas as pd

from script_utils import show_output, load_config


def load_lumi_config(analysis_name="results", config_file="", create_folders = True, **kwargs):
    '''
    loads the config from a yaml file and updates with keyword arguments
    flattens the paths to absolute values
    sets the experiment folder and creates the folders for it
    '''
    
    config = load_config(config_file)
    # set the base_path
    base_path = config['paths']['base_path']
    config['paths']['base_path'] = base_path if base_path.startswith("/") else os.path.join(os.environ['HOME'], base_path)

    # flatten the paths and add base_path if not abs_path
    for key, value in config['paths'].items():
        if not key == "base_path":
            config[key] = value if value.startswith("/") else os.path.join(config['paths']['base_path'], value)
    del config['paths']
    config['analysis_folder'] = os.path.join(config['output_path'], analysis_name)
    # create plotting folder
    p = config['plotting']['plot_folder_name']
    # leave plot_folder as "" if nothing is in there
    config['plotting']['plot_folder'] = os.path.join(config['analysis_folder'], p) if p else ""
    # load in the kwargs to overwrite
    config.update(kwargs)
    # create the folders
    if create_folders:
        # check for "" and then skip the folder building
        if config['analysis_folder']:
            if not os.path.isdir(config['analysis_folder']):
                show_output(f"Creating analysis folder {config['analysis_folder']}")
                os.makedirs(config['analysis_folder'])
        if (plot_folder := config['plotting']['plot_folder']) and config['plot_fit']:
            if not os.path.isdir(plot_folder):
                os.makedirs(plot_folder)
    return config
    

def load_existing(luminexcel_file):
    '''
    preload an existing file into a data dict
    '''
    
    old_data = {}
    for sheet in ['Plates', 'Standards', 'ProteinStats', 'tidyData']:
        df = pd.read_excel(luminexcel_file, sheet_name=sheet)
        if "Run" in df.columns:
            df.loc[:, 'Run'] = df['Run'].astype(str)
        old_data[sheet] = df
        
    return old_data


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
        "File Name": "sourcePath",
        "Acquisition Date": "AcquisitionTime",
        "Reader Serial Number": "ReaderID", 
        "RP1 PMT (Volts)": "RP1_PMT",
        "RP1 Target": "RP1_Target"
    })
    
    return plate_df

def get_run_plex(file, plex_pattern="Plex", plate_pattern="plate", **kwargs):
    '''
    retrieve run and flex from the file name using regex parts
    '''
    plate = 1
    for s in os.path.basename(file).split(".")[0].split("_"):
        if re.match(f"[0-9]+-{plex_pattern}", s):
            plex = s
        elif re.match(r"[0-9]{8}", s):
                run = str(int(s) - 20000000)
        elif re.match(f"{plate_pattern}([0-9]+)", s):
            plate = int(s.replace(plate_pattern, ""))
    return (run, plex, plate)


def isMatch(conc_file, run, plex, plate, **kwargs):
    '''
    checks if a conc_file fits to a run-plex-plate combo
    '''
    
    crun, cplex, cplate = get_run_plex(conc_file, **kwargs)
    
    return (crun == run) & (cplex == plex) & (cplate == plate)


def get_luminex_plates(*, data_path, raw_pattern="rawdata", conc_pattern="conc", **kwargs):
    '''
    create a data_list containing the raw_data and conc excel files
    for a given folder (recursively)
    output:
    list of
    {
        'run': 20211021,
        'plex': '11-Plex',
        'raw': <path_to_raw_data_file>,
        'conc': <path_to_conc_excel_file> or None
    }
    '''
    
    # init the file lists
    raw_file_list = []
    conc_file_list = []
    
    # stop if no pattern is given
    if not raw_pattern and not conc_pattern:
        return "No patterns!"

    # find recursively all relevant files
    for f in [folder for folder in os.walk(data_path)]:
        folder = f[0]
        # exclude temp files of open excel files
        cand_files = [os.path.join(folder, file) for file in f[2] if not os.path.basename(file).startswith("~$") and not os.path.basename(file).startswith(".") and os.path.splitext(file)[1] in [".csv", ".xls", ".xlsx"]]
        if raw_pattern:
            raw_files = [file for file in cand_files if raw_pattern in file.lower()]
        else: # if raw_pattern == "" it will be set by absence of conc pattern
            raw_files = [file for file in cand_files if not conc_pattern in file.lower()]
        raw_file_list += raw_files
        
        # conc files are defined by not containing the raw_pattern
        if conc_pattern:
            conc_files = [file for file in cand_files if conc_pattern in file.lower()]
        else:
            conc_files = [file for file in cand_files if not raw_pattern in file.lower()]
        conc_file_list += conc_files


    # find the matching conc files
    plate_list = []                
    for raw_file in raw_file_list:
        short_file = raw_file.replace(f"{data_path}/", "")
        run, plex, plate = get_run_plex(raw_file, **kwargs)
        conc_files = [f for f in conc_file_list if isMatch(f, run, plex, plate, **kwargs)]
        conc_file = conc_files[0].replace(f"{data_path}/", "") if len(conc_files) else None
        plate_list.append(dict(Run=run, Plex=plex, Plate=plate, rawPath=short_file, concPath=conc_file))       
    # check for isolated conc-files
    for conc_file in conc_file_list:
        run, plex, plate = get_run_plex(conc_file, **kwargs)
        raw_files = [f for f in raw_file_list if isMatch(f, run, plex, plate, **kwargs)]
        if not len(raw_files):
            plate_list.append(dict(Run=run, Plex=plex, Plate=plate, rawPath=None, concPath=conc_file.replace(f"{data_path}/", "")))
    
    # convert into df
    plates_df = pd.DataFrame(plate_list).sort_values(['Run', 'Plex', 'Plate']).reset_index(drop=True)
    return plates_df


def convert2float(df, cols=['conc']):

    for col in cols:
        df.loc[:, col] = df[col].str.replace("---", "-1")
        df.loc[:, col] = df[col].str.replace(",", ".", regex=False)
        df.loc[:, col] = df[col].str.replace("OOR <", "-2", regex=False).str.replace("OOR >", "-1", regex=False)
        df.loc[:, col] = df[col].str.replace("***", "-3", regex=False).str.replace("*", "", regex=False)
        df.loc[:, col] = df[col].astype(float)
    return df


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



def read_luminexcel(excel_path):
    '''
    load all the data from the luminexcel file
    '''
    plate_df = pd.read_excel(excel_path, sheet_name="Plates")
    col_df = pd.read_excel(excel_path, sheet_name="Plexes")
    data_df = pd.read_excel(excel_path, sheet_name="RawData")
    return plate_df, col_df, data_df