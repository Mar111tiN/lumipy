# COlUMNS FOR DBmerging


merge_col1 = [
    'Project',
    'Run',
    'RunDesc',
    'Plate',
    'Plex',
    'SampleName',
    'Dilution',
    'Well',
    'BiopsyName',
    'DOX',
    'Method',
    'SourceAmount',
    'SourceUnit',
    'Note_Sample',
    'Note_Well'
    ]
 
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

db_cols = [
    'Project',
    'PatientCode',
    'Sex',
    'DOD',
    'DOX',
    'Age'
] + case_cols + [
    'BiopsyName',
    'BioType',
    'Tissue',
    'SampleName',
    'Method',
    'LumiDOX',
    'SourceAmount',
    'SourceUnit',
    'Dilution',
    'Run',
    'RunDesc',
    'Plate',
    'Plex',
    'Well',
    'Note'
]


## CONSTANTS for PANKLOTINOhelper
patient_code_cols = [
    'Project', 
    'PatientCode', 
    'PatientCodeAlt', 
    'PatientCodeUsed', 
    'Note'
    ]

patient_cols = [
    'Project',
    'PatientCode', 
    'Sex',
    'DOD',
    'Note'
    ]

simple_case_cols = [
    'Project',
    'PatientCode',
    'DOX',
    'Age',
    'Disease',
    'Note'
    ]

sample_biopsy_cols = [
    'Project',
    'PatientCode',
    'SampleName',
    'SampleType',
    'Tissue',
    'Note'
    ]

biopsy_cols = [
    'Project',
    'PatientCode',
    'BiopsyName',
    'BiopsyType',
    'Tissue',
    'Note'
]

# well_cols combined with weights
sample_well_cols = [
    'Project',
    'Run',
    'RunDesc',
    'Plate',
    'Plex',
    'SampleName',
    'Note',
    'Weight',
    'Dilution',
    'Well'
]

well_cols = [
    'Project', 
    'Run', 
    'RunDesc', 
    'Plate', 
    'Plex', 
    'SampleName', 
    'Dilution', 
    'Well', 
    'Note'
]
# the new 
sample_cols = [
    'Project',
    'BiopsyName',
    'SampleName',
    'Method',
    'DOX',
    'SourceAmount',
    'SourceUnit',
    'Note'
]


### Columns for TinoMasterFile

tino_pat_cols = [
    'PatientCode',
    'Project',
    'PatientCodeAlt',
    'DOT',
    'Sex'
    ]

tino_case_cols = [
    'PatientCode',
    'Project',
    'tumor surgery',
    'Date of TURB',
    'Age at TURB',
    'Date of RC',
    'DOT',
    'muscle-invasive',
    'NAC',
    'Radiotherapy',
    'NAC cycles',
    'NAC regime',
    'response on NAC',
    ' Tumor state TURB',
    'Tumor grading',
    ' Tumor state RC',
    'nodal Metast',
    'visceral Metast',
    'skeletal Metast'
    ]

tino_sample_cols = [
    'Project',
    'PatientCode',
    'SampleName',
    'SampleType',
    'tumor for analysis',
    'Messrunde'
    ]


tino_well_cols = [
    'Run',
    'RunDesc',
    'Plate',
    'Plex',
    'SampleName',
    'Note',
    'Weight',
    'Dilution',
    'SName',
    'Well'
    ]

qPCR_cols = [
    'SampleName',
    'SampleType',
    'CXCR3alt_1',
    'CXCR3alt_2',
    'CXCR3A_1',
    'CXCR3A_2',
    'CXCR3B_1',
    'CXCR3B_2',
    'CD3_1',
    'CD3_2',
    'IPO8_1',
    'IPO8_2',
    'CDKN1B_1',
    'CDKN1B_2',
    'CXCR3alt',
    'CXCR3A',
    'CXCR3B',
    'CD3'
    ]
