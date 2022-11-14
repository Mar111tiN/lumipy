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
    'rep',
    'Well'
]

well_cols = [
    'Project', 
    'Run', 
    'RunDesc', 
    'Plate', 
    'Plex', 
    'SampleName',
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
    'Dilution',
    'SourceAmount',
    'SourceUnit',
    'Note'
]