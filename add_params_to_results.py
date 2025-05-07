from math import log10
import pandas as pd

# "HazardStudy" db table - index = StudyID
haz_studies_df = pd.read_csv("data/csvs/haz_studies.csv")

df = pd.read_csv('results_complete.csv')

def get_data_from_pandas(elem):
    if isinstance(elem, pd.Series):
        elem = elem.values[0]
    return elem

output = []
for _, row in df.iterrows():
    row_for_iter = row.to_dict()
    row_out = row.to_dict()
    study_id = row_for_iter['study_id']
    study_id = get_data_from_pandas(study_id)
    study_row = haz_studies_df[haz_studies_df['StudyID'] == study_id]
    press_psig = get_data_from_pandas(study_row['StoragePressure'])
    row_out['press_psig'] = press_psig
    output.append(row_out)

df_out = pd.DataFrame(output)
df_out.to_csv('results_complete_take_2.csv')

apple = 1