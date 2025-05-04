import csv
import copy
import json
import numpy as np
import pandas as pd
from datetime import datetime as dt
from datetime import timedelta
from py_lopa.model_interface import Model_Interface
from py_lopa.calcs import helpers
from py_lopa.calcs.consts import Consts

# "HazardStudy" db table - index = StudyID
haz_studies_df = pd.read_csv("data/csvs/haz_studies.csv")

# "StudyChemicalComposition" table - index = HazardStudy
chem_comps_df = pd.read_csv("data/csvs/chem_comps.csv")

last_completed_study_id = np.inf
header_written = False
try:
    with open('curr_idx.txt', 'r') as idx_file:
        last_completed_study_id = idx_file.read()
        if last_completed_study_id is not None and last_completed_study_id != '':
            last_completed_study_id = int(last_completed_study_id)
            header_written = True
        else:
            last_completed_study_id = np.inf
except:
    print('could not get last completed study id.')
    pass

haz_studies_df = haz_studies_df.sort_values(['StudyID'], ascending=[False])
haz_studies_df = haz_studies_df[haz_studies_df['StudyID'] < last_completed_study_id]
haz_studies_df = haz_studies_df[haz_studies_df['ReleaseElevation'] <= 6]


results = []
t0_0 = dt.now()
models_completed = 0
for idx, row in haz_studies_df.iterrows():
    t0 = dt.now()
    study_id = helpers.get_data_from_pandas_series_element(row['StudyID'])
    if study_id >= last_completed_study_id:
        continue
    print(f'Initiating Study ID {study_id}')
    chems_rows = chem_comps_df[chem_comps_df['HazardStudy'] == study_id]
    if isinstance(chems_rows, pd.Series):
        chems_info = [chems_rows.to_dict()]
    elif isinstance(chems_rows, pd.DataFrame):
        chems_info = chems_rows.to_dict(orient='records')
    else:
        print(f'issue getting chems from study id {study_id}')
        last_completed_study_id = study_id
        continue
    elev_m = helpers.get_data_from_pandas_series_element(row['ReleaseElevation'])
    if elev_m > 6:
        continue
    data = {
        'PrimaryInputs': row.to_dict(),
        'ChemicalComponents': chems_info,
    }
    json_data = json.dumps(data)
    m_io = Model_Interface()
    m_io.set_inputs_from_json(json_data=json_data)
    try:
        m_io.run()
    except:
        print(f'issue with model run for study id {study_id}')
        continue
    shi_data = m_io.mc.chems.shi_analysis_data
    data_for_output = {
        'study_id': study_id,
        'elev_m': elev_m,
        'ave_nbp_deg_c': shi_data['ave_nbp_deg_c'],
    }
    for condition in ['storage', 'discharge']:
        data_for_output[f'{condition}_key_component'] = None
        if 'shi_idx' in shi_data[condition]:
            if shi_data[condition]['shi_idx'] >= 0:
                key_component_idx = shi_data[condition]['shi_idx']
                key_component = m_io.mc.mi.CHEM_MIX[key_component_idx]
                data_for_output[f'{condition}_key_component'] = key_component
        data_for_output[f'{condition}_shi'] = None
        if 'shi_tox' in shi_data[condition]:
            shi_value = shi_data[condition]['shi_tox']
            data_for_output[f'{condition}_shi'] = shi_value
        data_for_output[f'{condition}_lfl'] = None
        if 'lfl' in shi_data[condition]:
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
                if len(areas_m2) == 0:
                    continue
                conc_types = ['low', 'high']
                for area_data, conc_type in zip(areas_m2, conc_types):
                    data_for_output[f'conc_{conc_type}_ppm'] = area_data['conc_ppm']
                    data_for_output[f'area_conc_{conc_type}_m2'] = area_data['area_m2']
                results.append(copy.deepcopy(data_for_output))
    
    if len(results) > 0:
        try:
            
            with open('resuls.csv', 'a', newline='') as results_out:
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
    models_remaining = len(haz_studies_df) - models_completed
    rate = 0
    if t_delta_total_secs > 0:
        rate = models_completed / t_delta_total_secs
    duration_remaining_sec = 0
    if rate > 0:
        duration_remaining_sec = models_remaining / rate
    duration_remaining_min = duration_remaining_sec / 60
    final_time = t1 + timedelta(minutes=duration_remaining_min)
    

    print(f'\n\n\n\n****************\n\n\nstudy id completed: {study_id}.  model {models_completed} / {idx}.  model runtime: {t_delta_model_secs:.2f} sec | est time remaining: {duration_remaining_min:.2f} min | est completion time: {final_time.strftime("%Y-%m-%d %H:%M:%S")}')

    last_completed_study_id = study_id
    try:
        with open('/mapping_shi_to_conseq_target/curr_idx.txt', 'w') as idx_file_update:
            idx_file_update.write(str(last_completed_study_id))
    except Exception as e:
        print('could not update last completed study id in text file')

    apple = 1

apple = 1
