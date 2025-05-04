import json
import pandas as pd
from py_lopa.model_interface import Model_Interface
from py_lopa.calcs import helpers
from py_lopa.calcs.consts import Consts

# "HazardStudy" db table - index = StudyID
haz_studies_df = pd.read_csv("data/csvs/haz_studies.csv")

# "StudyChemicalComposition" table - index = HazardStudy
chem_comps_df = pd.read_csv("data/csvs/chem_comps.csv")

last_completed_study_id = -1
try:
    with open('curr_idx.txt', 'r') as idx_file:
        last_completed_study_id = idx_file.read()
except:
    print('could not get last completed study id.')
    pass

haz_studies_df = haz_studies_df[haz_studies_df['StudyID'] > last_completed_study_id]
haz_studies_df.sort_values(['StudyID'])

results = []

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
    shi_data = m_io.mc.chems.shi_analysis_data
    data_for_output = {
        'ave_nbp_deg_c': shi_data['ave_nbp_deg_c'],
    }
    for condition in ['storage', 'discharge']:
        key_component_idx = shi_data[condition]['shi_idx']
        key_component = m_io.mc.mi.CHEM_MIX[key_component_idx]
        shi_value = shi_data[condition]['shi_tox']
        lfl = shi_data[condition]['lfl']
        for data_item in ['temp_deg_c', 'total_vapor_moles', 'component_vapor_moles', 'total_moles', 'component_total_moles']:
            elem = shi_data[condition][data_item]
            if isinstance(elem, list):
                elem = elem[key_component_idx]
            data_for_output[f'{condition}_{data_item}'] = elem
            if data_item != 'temp_deg_c':
                data_for_output[f'{condition}_shi_times_{data_item}'] = elem * shi_value
    p_disp = m_io.phast_disp
    conc_idx = 1
    for rel_dur_sec, p_disp_at_dur in p_disp.items():
        for wx, p_disp_at_dur_at_wx in p_disp_at_dur.items():
            for haz, p_disp_at_dur_at_wx_at_haz in p_disp_at_dur_at_wx.items():
                areas_m2 = p_disp_at_dur_at_wx_at_haz.areas_m2
                conc_idx = 1
                for area_data in areas_m2:
                    data_for_output[f'conc_{conc_idx}'] = area_data['conc_ppm']
                    data_for_output[f'area_{conc_idx}'] = area_data['area_m2']
                    conc_idx += 1


    apple = 1
    last_completed_study_id = study_id
    try:
        with open('curr_idx.txt', 'w') as idx_file_update:
            idx_file_update.write(last_completed_study_id)
    except Exception as e:
        print('could not update last completed study id in text file')


#simulating for now until back on emn connection
# m_io = main()

apple = 1
# store shi at process and storage conditions (based on key comp for toxic and overall for flammable)

# store footprints for lfl / 50% lfl / 25% lfl / erpg-1 / erpg-2 / erpg-3 / 10xerpg-3
