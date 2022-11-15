# COlUMNS FOR DBmerging

# all columns regarding any measurement

well_cols = [
    'Run',
    'RunGroup',
    'Plate',
    'Plex',
    'Well'
    ]


DBmerge_cols = [
    'Project',
    'BiopsyName',
    'SampleName',
    'LumiDOX',
    'Method',
    'SourceAmount',
    'SourceUnit',
    'ExtractVolume',
    'Note_Sample'] + well_cols + [
    'Note_Well'
    ]
    
lumi_cols = [
    'Protein',
    'FI',
    'conc',
    'Fpos',
    'concCI',
    'concMean',
    'concStd',
    'FposMean',
    'ResConc'
    ]

sample_cols = [
    'Project',
    'PatientCode',
    'BiopsyName',
    'BioType',
    'Tissue',
    'SampleName',
    'LumiDOX',
    'Method',
    'SourceAmount',
    'SourceUnit',
    'ExtractVolume',
    'ResUnit'
    ]



# columns for the duplicate df
dup_cols = [
    'Run',
    'Plex',
    'Plate',
    'Well',
    'Type',
    'SampleName',
    'SE',
    'Protein',
    'FI',
    'conc',
    'Fpos',
    'concCI',
    'concMean',
    'concStd', 
    'FposMean'
    ]

plex_cols = [f"{str(i)}-Plex" for i in [3,11,21,23,38]]


#################################
# 
case_cols = [
    'Disease',
    'Pathology',
    'Metast',
    'Therapy',
    'Radiotherapy',
    'TumorSurgery',
    'DoTURB',
    'DoRC',
    'muscle-invasive',
    'NACresponse',
    'TumorStateTURB',
    'TumorGrading',
    'TumorStateRC'
    ]

patient_cols = [
    'Project',
    'PatientCode',
    'Sex',
    'DOD',
    'DOX',
    'Age'
]

db_cols = patient_cols + \
    case_cols + \
        [col for col in sample_cols if not col in patient_cols] + \
            ['Note']