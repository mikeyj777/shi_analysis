import csv
import copy
import json
import pandas as pd
from datetime import datetime as dt
from py_lopa.model_interface import Model_Interface
from py_lopa.calcs import helpers
from py_lopa.calcs.consts import Consts

# "HazardStudy" db table - index = StudyID
haz_studies_df = pd.read_csv("data/csvs/haz_studies.csv")

# "StudyChemicalComposition" table - index = HazardStudy
chem_comps_df = pd.read_csv("data/csvs/chem_comps.csv")

last_completed_study_id = -1
header_written = False
try:
    with open('curr_idx.txt', 'r') as idx_file:
        last_completed_study_id = idx_file.read()
        header_written = True
except:
    print('could not get last completed study id.')
    pass

haz_studies_df.sort_values(['StudyID'])
haz_studies_df = haz_studies_df[haz_studies_df['StudyID'] > last_completed_study_id]


results = []
t0_0 = dt.now()
models_completed = 0
for idx, row in haz_studies_df.iterrows():
    t0 = dt.now()
    study_id = helpers.get_data_from_pandas_series_element(row['StudyID'])
    if study_id <= last_completed_study_id:
        continue
    chems_rows = chem_comps_df[chem_comps_df['HazardStudy'] == study_id]
    if isinstance(chems_rows, pd.Series):
        chems_info = [chems_rows.to_dict()]
    elif isinstance(chems_rows, pd.DataFrame):
        chems_info = chems_rows.to_dict(orient='records')
    else:
        print(f'issue getting chems from study id {study_id}')
        last_completed_study_id = study_id
        continue
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
        data_for_output[f'{condition}_shi'] = shi_value
        lfl = shi_data[condition]['lfl']
        data_for_output[f'{condition}_lfl'] = lfl
        for data_item in ['temp_deg_c', 'total_vapor_moles', 'component_vapor_moles', 'total_moles', 'component_total_moles']:
            elem = shi_data[condition][data_item]
            if isinstance(elem, list):
                elem = elem[key_component_idx]
            data_for_output[f'{condition}_{data_item}'] = elem
            if data_item != 'temp_deg_c':
                data_for_output[f'{condition}_shi_times_{data_item}'] = elem * shi_value
    p_disp = m_io.phast_disp
    if len(p_disp) == 0:
        print(f'issue running study id {study_id}')
        last_completed_study_id = study_id
        continue

    for rel_dur_sec, p_disp_at_dur in p_disp.items():
        data_for_output['release_duration_sec'] = rel_dur_sec
        for wx, p_disp_at_dur_at_wx in p_disp_at_dur.items():
            data_for_output['weather'] = wx
            for haz, p_disp_at_dur_at_wx_at_haz in p_disp_at_dur_at_wx.items():
                data_for_output['hazard_type'] = haz
                areas_m2 = p_disp_at_dur_at_wx_at_haz.areas_m2
                for area_data in areas_m2:
                    data_for_output[f'conc_ppm'] = area_data['conc_ppm']
                    data_for_output[f'area_m2'] = area_data['area_m2']
                    results.append(copy.deepcopy(data_for_output))

    try:
        with open('resuls.csv', 'a') as results_out:
            fieldnames = list(data_for_output.keys())
            writer = csv.DictWriter(results_out, fieldnames=fieldnames)
            if not header_written:
                writer.writeheader()
                header_written = True
            writer.writerows(results)
            results = []
    except Exception as e:
        print(f'\n\n\n***********\n\n\ncould not store current results to file.  will try again after next model completes.\nmessage provided: {e}')
    
    t1 = dt.now()
    t_delta_model = t1 - t0
    t_delta_model_secs = t_delta_model.total_seconds()
    t_delta_total = t1 - t0_0
    t_delta_total_secs = t_delta_total.total_seconds()
    models_completed += 1
    models_remaining = len(haz_studies_df) - idx
    rate = 0
    if t_delta_total_secs > 0:
        rate = models_completed / t_delta_total_secs
    duration_remaining_sec = models_remaining * rate
    duration_remaining_min = duration_remaining_sec / 60

    print(f'\n\n\n\n****************\n\n\ncurrent study id just completed: {study_id}.  this is model {idx} out of {len(haz_studies_df)}.  time_to_run_model: {t_delta_model_secs} sec | time remaining: {duration_remaining_min} min')

    last_completed_study_id = study_id
    try:
        with open('curr_idx.txt', 'w') as idx_file_update:
            idx_file_update.write(last_completed_study_id)
    except Exception as e:
        print('could not update last completed study id in text file')

    apple = 1

apple = 1
