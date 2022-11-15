import pandas as pd
import numpy as np
from script_utils import show_output
from DBcols import *


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
    '''
    aggregates multiple Note-fields (starting with "Note") into one "Note" field
    removes duplicate info
    '''
    note_cols = [col for col in df.columns if col.startswith("Note")]
    new_cols = [col + "_1" for col in note_cols]
    df = df.rename({col:col + "_1" for col in note_cols}, axis=1)
    df['Note'] = ""
    for note_col in new_cols:
        df[note_col] = df[note_col].fillna("")
        df.loc[df['Note'] != df[note_col], 'Note'] = df['Note'] + "|" + df[note_col].fillna("")

    # add leading and trailing "|" for deduplication    
    df['Note'] = ("|" + df['Note'] +  "|").str.replace('\|+', '|', regex=True).str.replace("^\|+$", "", regex=True)
    # do the deduplication
    df.loc[df['Note'] != "", 'Note'] = df['Note'].str.replace("(\|[^|]+)(\|.*)?\\1", "\\1\\2", regex=True)
    df = df.drop(new_cols, axis=1)
    return df


def merge_DBsamples(DB_file):
    '''
    take all levels of the MSHDB and and return the files for the different projects
    '''

    # load the MSHDB as a data_dictionary
    db_dict = read_MSH_DB(DB_file)

    # merge wells and samples and flatten notes
    dbdf = db_dict['Samples'].merge(db_dict['Wells'], on=["SampleName", "Project"], how="right")
    show_output("Merging Wells and Samples")
    show_dups(dbdf)

    dbdf = flatten_notes(dbdf.loc[:, DBmerge_cols])
    # merge with the Biopsies and with the Cases
    dbdf = flatten_notes(db_dict['Biopsies'].merge(dbdf, how="right"))
    show_output("Merging Biopsies")
    show_dups(dbdf)
    return dbdf, db_dict


def compute_resConc(df):
    '''
    make separate concentration measurements for serum and tissue
    '''

    # set the Serum unit to ""
    df.loc[df['Tissue'] == "Serum", "SourceUnit"] = ""

    # set separate units and calculations for serum and tissue
    df.loc[df['SourceUnit'] == "mg", 'ResUnit'] = "pg Protein / g Tissue"
    df.loc[df['SourceUnit'] == "", 'ResUnit'] = "pg Protein / ml Serum"
    df.loc[df['SourceUnit'] == "mg", 'ResConc'] = df['concMean'] * df['ExtractVolume'] / df['SourceAmount']
    df.loc[df['SourceUnit'] == "", 'ResConc'] = df['concMean']

    # pivot the sample using all informative columns as index
    index_cols = [col for col in df.columns if not col in lumi_cols]
    df = df.set_index(index_cols).pivot(columns="Protein", values="ResConc")
    return df


def get_measurements(df, agg_dict):
    '''
    compute the doublicate measurements per plate
    '''

    measure_df = df.drop(['Well',
       'Type', 'SE'], axis=1).groupby(sample_cols + ['Run', 'RunGroup', 'Plex', 'Plate']).agg(agg_dict).reset_index()
    
    # var_df is not really needed - replicats are predicted to be close and are 
    var_df = df.drop(['Well', 'Note',
       'Type', 'SE'], axis=1).groupby(sample_cols + ['Run', 'RunGroup', 'Plex', 'Plate']).agg(np.var).fillna(-1).reset_index()
    return measure_df


def make_sample_df(measure_df, agg_dict):
    '''
    identifies the samples with 
    '''

    # get a run_df for detecting for all samples in which runs they were measured
    run_cols = ['Project', 'BiopsyName', 'SampleName', 'LumiDOX', 'RunGroup']
    run_df = measure_df.set_index(run_cols).pivot(columns="Plex", values="Run")  # .reset_index()
    run_df = run_df.fillna(0).astype(int).reset_index()
    
    # get an indicator "repeatPlexes" for the samples with repeat measurements"
    s = run_df.set_index(run_cols).fillna(0).astype(bool).astype(int).groupby([
        'Project', 'BiopsyName', 'SampleName',
    ]).agg('sum')
    s['repeatPlexes'] = (np.sum((s > 1).astype(int), axis=1) > 0)
    
    # split by merging with indicator s according to repeat or not repeat
    measure_df = measure_df.merge(s.reset_index().loc[:, ['SampleName', 'repeatPlexes']], how="outer")
    repeat_measure_df = measure_df.query('repeatPlexes == True')
    single_measure_df = measure_df.query('repeatPlexes == False')

    # for single measurements (no repeats in any plex) everything can be grouped by sample_name
    single_sample_df = single_measure_df.drop(['Plex', 'Plate', 'repeatPlexes', 'Run'], axis=1).groupby(sample_cols).agg(agg_dict).reset_index()
    single_sample_df['RunGroup'] = "Messung"
    
    # for repeat measurements, we need the RunGroup reduced to Messung or Nachmessung for grouping
    # unify the RunGroups to Messung or Nachmessung
    repeat_measure_df['RunGroup']  = repeat_measure_df['RunGroup'].str.replace('essung.*', 'essung', regex=True)
    # perform the grouping
    repeat_sample_df = repeat_measure_df.drop(['Plex', 'Plate', 'repeatPlexes', 'Run'], axis=1).groupby(sample_cols + ['RunGroup']).agg(agg_dict).reset_index()
    # also use simplified RunGroup for run_df for correct merging
    run_df['RunGroup']  = run_df['RunGroup'].str.replace('essung.*', 'essung', regex=True)
    # recombine the 

    run_df = run_df.groupby(run_cols).agg(np.sum).reset_index()
    sample_df = pd.concat([repeat_sample_df, single_sample_df]).merge(run_df.drop('RunGroup', axis=1)).reset_index(drop=True).sort_values(['Project', 'PatientCode', 'BiopsyName', 'SampleName'])
    return sample_df


def mergeDB2Lumi(dbdf, luminexcel="", verbose=False):
    '''
    load in the luminex data and merge in the respective Luminexdata
    '''
    
    # load the computed Luminex data per well
    lumi_dict = read_lumi_data(luminexcel)
    lumi_df = lumi_dict['tidyData']
    
    lumi_df.loc[lumi_df['Plex'] == "21-Plex-StandardOnly", 'Plex'] = "21-Plex"
    
    # split off the std_df containing all the standards
    std_df = lumi_df.loc[lumi_df['Type'].str.match("^[BCS]"), :].sort_values(['Run', 'Plex', 'Protein', 'Type'])
    
    # keep on with the lumi_data without std data
    lumi_df = lumi_df.loc[~lumi_df['Type'].str.match("^[BCS]"), :]
    
    # CONSISTENCY CHECK
    # merge the database and luminex for consistency check
    merge_df = dbdf.loc[:, ['Project', 'Run', 'Plate', 'Plex', 'Well']].merge(lumi_df.loc[:, ['Run', 'Plate', 'Plex', 'Well']].drop_duplicates(), how="outer", indicator=True)
    # check for wells missing in luminex
    if len(ddbdf := merge_df.query('_merge == "left_only"')):
        show_output(f"{len(ddbdf)} database wells were not found in luminex data!", color="warning")
        print(ddbdf)
    else:
        show_output("All database wells have been found in Luminex data!", color="success")
    # check for wells missing in database using the adhoc merge
    missing_wells = False
    missing_length = 0
    for run in merge_df['Run'].unique():
        for plate in merge_df.query('Run == @run')['Plate'].unique():
            df = merge_df.query('Run == @run and Plate == @plate and _merge == "right_only"')
            if (l := len(df)):
                missing_wells = True
                missing_length += l
                if verbose:
                    show_output(f"{l} wells missing from Run {run} - Plate {plate}", color="warning")
                    show_output(list(df['Well'].unique()), color="warning")
    if missing_wells:
        if not verbose:
            show_output(f"{missing_length} luminex wells were missing in database", color="warning")
    else:
        show_output("All luminex wells have been found in database")
    
    #### MERGE
    # do the merge with the lumi data
    merge_df = dbdf.merge(lumi_df, how="left").reset_index(drop=True)
    
    # keep the duplicates in a separate df
    dup_df = pd.concat([
        # get dups from standards
        std_df.loc[std_df.duplicated(['Run', 'Plex', 'Plate', 'Type', 'Protein'], keep=False), :],
        merge_df.loc[merge_df.duplicated(['Run', 'Plate', 'Plex', 'SampleName', 'Protein'], keep=False), :]
    ]).loc[:, dup_cols] 
    
    # compute the weight-integrated resConc
    # returns a df indexed for all metadata columns 
    conc_df = compute_resConc(merge_df)
    # keep the protein_cols and index_cols for flexible grouping in the next step to collate the different plexes per RunGroup
    # all info except for the proteins is in the index
    protein_cols = list(conc_df.columns)
    # reset index for basic well_df (protein measurements per well)
    well_df = conc_df.reset_index()

    # aggregate the duplicates
    # make an agg_dict that aggregates all proteins with np.mean and note with "first"
    agg_dict = {'Note':"first"}
    agg_dict.update({prot:np.mean for prot in protein_cols})
    measure_df = get_measurements(well_df, agg_dict)
    # sample_df has to be created differently    
    sample_df = make_sample_df(measure_df, agg_dict).loc[:, sample_cols + ['RunGroup', 'Note'] + plex_cols + protein_cols]
    return merge_df, std_df, dup_df, well_df, measure_df, sample_df


def merge_DBcases(sample_df, db_dict):
    '''
    take all levels of the MSHDB and and return the files for the different projects
    '''

    dbdf = flatten_notes(db_dict['Cases'].merge(sample_df, how="left"))
    show_output("Merging Cases")
    show_dups(dbdf)
    dbdf = flatten_notes(db_dict['Patients'].merge(dbdf, how="left"))
    show_output("Merging Patients")
    show_dups(dbdf)

    # dbdf = dbdf.loc[:, db_cols]
    return dbdf


def collect_DB_data(DB_file="", luminexcel="", verbose=False, excel_out=""):
    '''
    collect all data 
    '''
    # get the sample metadata for merging with luminexcel
    dbdf, db_dict = merge_DBsamples(DB_file)
    
    # combine DB data with luminex data
    merge_df, std_df, dup_df, well_df, measure_df, sample_df = mergeDB2Lumi(dbdf, luminexcel, verbose=verbose)
    
    # merge the clinical data
    sample_df = flatten_notes(merge_DBcases(sample_df, db_dict))
    # sample_df.loc[nac_samples['Note'] != "", :]['Note'] # .str.replace("(\|[^|]+)(\|.*)?\\1", "\\1\\2", regex=True)
    measure_df = flatten_notes(merge_DBcases(measure_df, db_dict))

    # convert datetime to date
    sample_df['DOD'] = pd.to_datetime(sample_df['DOD']).dt.date
    for df in [well_df, measure_df, sample_df]:
        for col in ['DOX', 'LumiDOX', 'DoTURB', 'DoRC']:
            if col in df.columns:
                df[col] = pd.to_datetime(df[col]).dt.date


    # transform date values
    ### write to file
    if excel_out:
        show_output(f"Data written to {excel_out}!", color="success")
        with pd.ExcelWriter(excel_out, mode="w") as writer:
            # merge_df.to_excel(writer, sheet_name="Data", index=False)
            std_df.to_excel(writer, sheet_name="Controls", index=False)
            dup_df.to_excel(writer, sheet_name="Duplicates", index=False)
            well_df.to_excel(writer, sheet_name="Wells", index=False)
            measure_df.to_excel(writer, sheet_name="Measurements", index=False)
            sample_df.to_excel(writer, sheet_name="Samples", index=False)
    
    return merge_df, std_df, dup_df, well_df, measure_df, sample_df


############ ANALYSE #####################################

def getNAC(measure_df, sample_df):
    '''
    retrieve the NAC specific data
    '''
    nac_mdf = measure_df.query('Project == "NAC"').loc[:, lambda x: ~x.isna().all()]
    for col in ['Run', 'Plate']:
        nac_mdf[col] = nac_mdf[col].astype(int)
    nac_sdf = sample_df.query('Project == "NAC"').loc[:, lambda x: ~x.isna().all()]
    
    # get the statistics for the responders based on plates
    nac_mdf.loc[:, ['CR', 'NR', 'PA', 'PROG']] = pd.get_dummies(nac_mdf['NACresponse'])
    nacc = nac_mdf.loc[:, ['PatientCode', 'Sex', 'Age', 'muscle-invasive', 'SampleName', 'Run', 'Plex', 'Plate', 'CR', 'NR', 'PA', 'PROG']]
    plate_df = nacc.drop(['PatientCode', 'Sex', 'Age', 'muscle-invasive', 'SampleName'], axis=1).groupby(['Run', 'Plex', 'Plate']).agg(np.sum).loc[lambda x: x.sum(axis=1)>0, :]# 
    plate_df['ratioCR'] = plate_df['CR'] / plate_df.sum(axis=1)
    plate_df = plate_df.sort_values(['Plex', 'ratioCR']).reset_index()
    return nac_mdf, nac_sdf, plate_df
    nac_samples = sample_df.query('Project == "NAC"').loc[:, lambda x: ~x.isna().all()]