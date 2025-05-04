import json
import pandas as pd
from py_lopa.model_interface import Model_Interface
from py_lopa.calcs import helpers

# "HazardStudy" db table - index = StudyID
haz_studies_df = pd.read_csv("data/csvs/haz_studies.csv")

# "StudyChemicalComposition" table - index = HazardStudy
chem_comps_df = pd.read_csv("data/csvs/chem_comps.csv")

last_completed_study_id = -1
try:
    with open('curr_idx.txt') as idx_file:
        last_completed_study_id = idx_file.read()
except:
    pass

haz_studies_df = haz_studies_df[haz_studies_df['StudyID'] > last_completed_study_id]

for idx, row in haz_studies_df.iterrows():
    study_id = helpers.get_data_from_pandas_series_element(row['StudyID'])
    if study_id <= last_completed_study_id:
        continue
    chems_rows = chem_comps_df[chem_comps_df['HazardStudy'] == study_id]
    if isinstance(chems_rows, pd.Series):
        chems_info = [chems_rows.to_dict()]
    elif isinstance(chems_rows, pd.DataFrame):
        chems_info = chems_rows.to_dict(orient='records')
    else:
        print(f'issue with study id {study_id}')
        apple = 1

    data = {
        'PrimaryInputs': row.to_dict(),
        'ChemicalComponents': chems_info,
    }
    
    
    json_data = json.dumps(data)

    m_io = Model_Interface()
    m_io.set_inputs_from_json(json_data=json_data)
    m_io.run()
    apple = 1

#simulating for now until back on emn connection
# m_io = main()

apple = 1
# store shi at process and storage conditions (based on key comp for toxic and overall for flammable)

# store footprints for lfl / 50% lfl / 25% lfl / erpg-1 / erpg-2 / erpg-3 / 10xerpg-3
