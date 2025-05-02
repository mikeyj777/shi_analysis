import os
import json
import datetime
import pandas as pd
from datetime import datetime as dt

import tkinter as tk
from tkinter import filedialog as fd

from pypws.entities import ResultCode

from py_lopa.calcs import helpers
from py_lopa.calcs.get_phys_props import Get_Phys_Props
from py_lopa.model_interface import Model_Interface
from py_lopa.kml_processing.kml_generator import Kml_Generator
import summarize_model_runs_from_json

def run_gpp(m_io):

    chem_mix = m_io.inputs['chemical_mix']
    chem_composition = m_io.inputs['composition']
    molar_basis = m_io.inputs['comp_is_moles']
    temp_C = m_io.inputs['temp_deg_c']
    press_psig = m_io.inputs['pressure_psig']

    gpp = Get_Phys_Props(chem_mix=chem_mix, chem_composition=chem_composition, molar_basis=molar_basis, temp_C=temp_C, press_psig=press_psig)
    gpp.run()
    print(gpp.data) 
    print('\n\n\n\n\n---------------------\n\n\n\n\n')

def run_m_io(m_io:Model_Interface, path_to_file, include_google_earth_output = False, summarize_files_when_complete=False):
    status = m_io.run()

    if status != ResultCode.SUCCESS:
        print('no bueno')
        return
    
    model_description = path_to_file.split('/')[-1]
    model_description = helpers.replace_chars_in_string_prior_to_naming_windows_text_file(model_description)
    data = m_io.get_data()
    print(data)
    if summarize_files_when_complete:
        output_folder = os.path.join(os.getcwd(), 'model_runs_from_json_output')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        helpers.output_variable_to_textfile(data, file_nm = f'{output_folder}/{model_description}')
    
    # launch kml generator

    if summarize_files_when_complete:
        output_folder = os.path.join(os.getcwd(), 'google_earth_output')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    run_info = m_io.time_stamp_utc
    if include_google_earth_output:
        kml_handler = m_io.mi.KML_HANDLER
        if not summarize_files_when_complete:
            kml_handler("NOTE: kml will be generated, but per specficiation, it will not be output to a file.")
        kml_gen = Kml_Generator(dispersion_dicts=m_io.dispersion_dicts, mi=m_io.mi, save_all_disps_locally=False, save_merged_output_locally=summarize_files_when_complete, output_folder=output_folder, run_info=run_info, send_to_kml_hander=False, kml_handler=kml_handler)
        kml_gen.run()

def update_m_io_inputs_with_mods_dict(m_io, mods):
    for k, v in mods.items():
        m_io.inputs[k] = v
    
    return m_io

def main(down_with_gpp = False, generate_google_earth_output = False, summarize_files_when_complete=False, modifications_dict = {}):

    root = tk.Tk()
    root.withdraw()

    files_to_run = fd.askopenfilenames(parent=root, title='Choose files to run')

    for path_to_file in files_to_run:

        if path_to_file[-5:] != '.json':
            continue

        m_io = Model_Interface()
        m_io.set_inputs_from_json(path_to_json_file=path_to_file)

        print(f'{path_to_file}\n\n')
        
        if down_with_gpp:
            run_gpp(m_io)

        m_io = update_m_io_inputs_with_mods_dict(m_io=m_io, mods=modifications_dict)

        run_m_io(m_io, path_to_file=path_to_file, include_google_earth_output=generate_google_earth_output, summarize_files_when_complete=summarize_files_when_complete)

        t1 = dt.now(datetime.UTC)
        delta_t = t1 - m_io.t0
        print(f'Start Time {m_io.t0}  End time {t1}  time to run {delta_t}')

    if summarize_files_when_complete:
        summarize_model_runs_from_json.run()
    
    return m_io

if __name__ == '__main__':
    down_with_gpp = True
    generate_google_earth_output = False
    summarize_files_when_complete = True

    m_io = main(down_with_gpp=down_with_gpp, generate_google_earth_output=generate_google_earth_output, summarize_files_when_complete=summarize_files_when_complete)