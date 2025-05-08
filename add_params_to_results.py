from math import log10, pi
import pandas as pd

# "HazardStudy" db table - index = StudyID
haz_studies_df = pd.read_csv("data/csvs/haz_studies.csv")

df = pd.read_csv('results_complete.csv')

def get_data_from_pandas(elem):
    if isinstance(elem, pd.Series):
        elem = elem.values[0]
    return elem

params_to_find = ['StoragePressure', 'ExitPressure', 'CatastrophicRelease', 'DischargeDiameter', 'ReleaseOrientation']
param_title_to_store = ['storage_pressure', 'exit_pressure', 'catastrophic_release', 'discharge_diameter', 'release_orientation']

output = []
for _, row in df.iterrows():
    row_for_iter = row.to_dict()
    row_out = row.to_dict()
    study_id = row_for_iter['study_id']
    study_id = get_data_from_pandas(study_id)
    study_row = haz_studies_df[haz_studies_df['StudyID'] == study_id]
    for param, title in zip(params_to_find, param_title_to_store):
        val = study_row[param]
        val = get_data_from_pandas(val)
        row_out[title] = val
    row_out['pressure'] = row_out['storage_pressure']
    if pd.isna(row_out['pressure']):
        row_out['pressure'] = row_out['exit_pressure']
    del row_out['storage_pressure']
    del row_out['exit_pressure']

    output.append(row_out)

df_out = pd.DataFrame(output)
df_out.to_csv('results_complete_take_2.csv')

apple = 1