import sys
import copy
import json
import pandas as pd
import datetime
from datetime import datetime as dt
from pypws.enums import ResultCode, Phase, FluidSpec, VesselConditions, MixtureModelling
from pypws.entities import FlashResult, State

from py_lopa.data.exception_enum import Exception_Enum
from py_lopa.classes.plotting_and_debug_data import Plotting_and_Debug_Data
from py_lopa.classes.conseq_assess import Conseq_Assess
from py_lopa.classes.all_class_results import All_Class_Results
from py_lopa.calcs import helpers, thermo_pio
from py_lopa.calcs.consts import Consts, Wx_Enum
from py_lopa.phast_io import phast_prep
from py_lopa.phast_io.phast_discharge import Phast_Discharge
# import for phast dispersion is at the bottom.  workaround for recursive calls
from py_lopa.model_work.model_inputs import Model_Inputs
from py_lopa.phast_io import phast_prep
from py_lopa.classes.interim_results_to_log_handler import Interim_Results_to_Log_Handler
from py_lopa.calcs.consts import Consts
from py_lopa.classes.vce import VCE
from py_lopa.phast_io.dnv_test_results import DNV_Test_Results
from py_lopa.calcs.pv_burst_blast_calculation import Pv_Burst_Blast_Calc
from py_lopa.calcs.flattening import Flattening

cd = Consts().CONSEQUENCE_DATA

class Model_Controller:

    def __init__(self, inputs, release_durations_to_model_sec = []):

        self.dose_response_exclusion_zones_df = pd.DataFrame()
        self.inputs = inputs
        self.resultCode = ResultCode.FAIL_EXECUTION
        self.data = ''
        self.discharges_arr = []
        self.dispersion_dicts = []
        self.release_durations_to_model_sec = release_durations_to_model_sec
        if len(self.release_durations_to_model_sec) == 0:
            self.release_durations_to_model_sec = copy.deepcopy(Consts.RELEASE_DURATIONS_TO_MODEL_SEC)

        # including the material, discharge and dispersions as properties.  
        # these can be used in a header analysis to quickly reference already completed dispersions/discharges for 
        # redundant models.
        self.phast_discharge = {}
        self.phast_disp = {}
        self.material = None
        self.flashresult = None
        self.vce_data = {}
        self.dnv_benchmark_results = None
        self.output_list = []
        self.debug_data = None
        self.timings = []
        self.dispersion_plots = []
        self.impact_areas_m2_dict = {}
        

    def run(self):

        #get all data from input sources (release parameters, chem info data, etc.).
        #parse inputs

        mi:Model_Inputs = Model_Inputs(self.inputs)
        self.mi = mi
        self.chems = mi.chems

         # if no hazard identified by the model, it will return a json string
        # with null for all consequence values.  this code develops the empty
        # table that will be used.

        default_data_list = helpers.prep_default_list(mi=mi)

        if not mi.INHALATION and not mi.FLASH_FIRE:
            mi.LOG_HANDLER('identified release does not contain hazardous material')
            self.data = helpers.list_to_json(default_data_list)
            return ResultCode.FAIL_VALIDATION

        #info on the chemicals being released and chem property data tables.
        chems = mi.chems

        t0 = dt.now(datetime.UTC)
        mi.LOG_HANDLER(f'\n********\n\nInitiating model sequence.  Current time in UTC {t0}\n\n')

        material = self.material
        if material is None:
            material = phast_prep.prep_material(chems)

        #get initial state using input temp and pressure
        lf = None
        if mi.EXIT_VAPOR_MASS_FRACTION is not None:
            lf = 1 - mi.EXIT_VAPOR_MASS_FRACTION
        state = phast_prep.prep_state(
            press_pa=mi.PRESS_PA, 
            temp_K=mi.TEMP_K, 
            lf=lf, 
            use_multicomponent_modeling=mi.USE_MULTICOMPONENT_METHOD
        )

        flashresult = self.flashresult
        
        if flashresult is None:
            flashresult = phast_prep.flash_calc(state, material)

        if state.flash_flag != FluidSpec.TP:
            mi.PRESS_PA = flashresult.pressure
            mi.TEMP_K = flashresult.temperature
            mi.EXIT_VAPOR_MASS_FRACTION = 1 - flashresult.liquid_mass_fraction

        # if the flash result returns a liquid, but the vapor relief flag is true,
        # then the vessel will need to be padded with nitrogen until the thermo allows for a
        # vapor release analysis.

        # determine the phase that will be released, pad as needed.  
        # store final composition in "chems.releasing_stream_chem_data['mole_fraction']"

        # iterate on final temp such that system temp is just above sat temp with n2 pad

        resp = chems.final_system_checks(flashresult = flashresult, state = state, material=material)

        state = resp['state']
        flashresult:FlashResult = resp['flashresult']
        mi = resp['model_inputs']
        material = resp['material']

        self.material = material
        self.flashresult = flashresult

        # initial iteration of vessel size.
        # based on 60 min release duration. 
        # discharge rate will be calculated using bernoulli, then improved upon.
        
        vessel_details = phast_prep.prep_vessel(
            mi = mi, 
            state = state, 
            material = material, 
            chems = chems, 
            flashresult = flashresult, 
            release_duration_sec = self.release_durations_to_model_sec[-1]
        )
        
        vessel = vessel_details['vessel']

        vessel_size = vessel_details['vessel_size']
        
        # if state.flashFlag != FluidSpec.TP:
        #     vessel.vesselConditions = VesselConditions.PRESSURIZED_LIQUID_VESSEL

        self.rho_kgm3 = vessel_size.rho_kgm3

        leak_height_fract_of_vessel = vessel_size.leak_height
        stop_after_this_release = vessel_size.stop_after_this_release
        mass_released_kg = vessel_size.mass_released_kg
        release_duration_sec = vessel_size.release_duration_sec
        
        #  self.height_m = ht_m 
        #     self.diam_m = diam_m
        #     self.liquid_fill = liquid_fill_fract
        #     self.vessel_condition = vc
        #     self.leak_height = leak_height_fract_of_vessel
        #     self.stop_after_this_release = False
        #     self.mass_released_kg = mass_released_kg

        leak = phast_prep.prep_leak(mi, leak_height_fract_of_vessel)
        discharge = None

        if release_duration_sec in self.phast_discharge:
            discharge = self.phast_discharge[release_duration_sec]

        flash_at_orifice_allowed = True

        phase:Phase = flashresult.fluid_phase
        if isinstance(phase, int):
            phase = Phase(phase).name
        print(f'mc.  after system checks.  fluid phase is {phase}.')
        print(f'vessel definition:  {vessel.__dict__}')

        if discharge is None:
            
            # homogenous 2-phase releases state needing a p_lf flash.
            mixture_method = MixtureModelling.PC
            if mi.USE_MULTICOMPONENT_METHOD:
                mixture_method = MixtureModelling.MC_SINGLE_AEROSOL
            if (vessel.vessel_conditions == VesselConditions.HOMOGENEOUS_VESSEL_TYPE) or (vessel.vessel_conditions == VesselConditions.STRATIFIED_TWO_PHASE_VESSEL):
                state = State(pressure=mi.PRESS_PA, temperature=None, liquid_fraction=flashresult.liquid_mole_fraction, flash_flag=FluidSpec.PLF, mixture_modelling=mixture_method)
                flash_at_orifice_allowed = False
            
            discharge = Phast_Discharge(
                leak = leak, 
                state = state,
                material = material,
                vessel = vessel,
                chems = chems,
                mi = mi,
                flashresult=flashresult,
                flash_at_orifice_allowed=flash_at_orifice_allowed)
            
            discharge.run()

            self.phast_discharge[release_duration_sec] = discharge

        disch_diag = helpers.discharge_diagnosis(discharge)

        if disch_diag != ResultCode.SUCCESS:
            
            self.data = helpers.list_to_json(default_data_list)
            return disch_diag
        
        # if mi.EXIT_VAPOR_MASS_FRACTION is not None:
        #     lmf = 1 - mi.EXIT_VAPOR_MASS_FRACTION
        #     t_p_diam_dict = thermo_pio.check_for_diff_bet_target_liquid_mass_fraction_at_rate_and_pws_liquid_mass_fraction_return_adjusted_temperature_deg_k_press_pa_and_diameter_m(pd=discharge, target_liquid_mass_fraction=lmf)
        #     mi.TEMP_K = t_p_diam_dict['temp_k']
        #     mi.PRESS_PA = t_p_diam_dict['press_pa']
        #     # mi.MAX_HOLE_SZ_M = t_p_diam_dict['diam_m']

        if mi.RELIEF_CALC_AVAILABLE:
            if len(self.phast_discharge) > 0:
                min_stored_discharge_duration = min(self.phast_discharge.keys())
                discharge = self.phast_discharge[min_stored_discharge_duration]
            discharge_mass_flow = discharge.vesselLeakCalculation.discharge_records[0].mass_flow
            targ_mass_flow = mi.RELIEF_CALC_RATE_KG_S
            if abs(discharge_mass_flow - targ_mass_flow) / targ_mass_flow > Consts.TOLERABLE_RELIEF_RATE_DEVIATION_FROM_VESSEL_CALCS:
                discharge = thermo_pio.match_relief_rate(discharge = discharge, mi = mi, relief_rate=mi.RELIEF_CALC_RATE_KG_S)
                disch_diag = helpers.discharge_diagnosis(discharge)
                if disch_diag != ResultCode.SUCCESS:
                    self.data = helpers.list_to_json(default_data_list)
                    return disch_diag
                leak = discharge.leak
                mi.MAX_HOLE_SZ_M = leak.hole_diameter

        # if mi.CATASTROPHIC_VESSEL_FAILURE:
        #     rel_sec = Consts.CATASTROPHIC_RELEASE_DURATIONS_TO_MODEL_SEC[0]
        #     if rel_sec <= 0:
        #         rel_sec = 10
        #     mass_flow_kg_s = mi.STORAGE_MASS_KG / rel_sec
        #     discharge = thermo_pio.match_relief_rate(discharge = discharge, mi = mi, relief_rate=mass_flow_kg_s, rel_duration_sec=rel_sec)
        #     disch_diag = helpers.discharge_diagnosis(discharge)
        #     if disch_diag != ResultCode.SUCCESS:
        #         self.data = helpers.list_to_json(default_data_list)
        #         return disch_diag
        #     mi.MAX_HOLE_SZ_M = discharge.leak.holeDiameter
        
        diam_str = '{:.4}'.format(mi.MAX_HOLE_SZ_M * 3.28084 * 12)
        mi.id_string += f'hole size: {diam_str} in\n'

        mass_flow_kg_s = discharge.vesselLeakCalculation.discharge_records[0].mass_flow

        output_list = []
        debug_output_list = []
        dispersion_plot_list = []
        
        if mi.CATASTROPHIC_VESSEL_FAILURE:
            self.release_durations_to_model_sec = copy.deepcopy(Consts.CATASTROPHIC_RELEASE_DURATIONS_TO_MODEL_SEC)

        stop_after_this_release = False

        vlc = discharge.vesselLeakCalculation
        data_from_discharge = thermo_pio.get_vapor_phase_comp_and_flash_calc_from_discharge_result(vessel_leak_calc=vlc, chems=chems, mi=mi)

        mi.ESTIMATE_MASS_FLOW_KG_S = mass_flow_kg_s
        mi.ESTIMATE_MASS_FLOW_LB_HR = mass_flow_kg_s * 3600 * Consts.LB_PER_KG
        mi.DOWNSTREAM_ESTIMATE_DISCHARGE_IDEAL_FLASH_CALC_ON_EXIT_MATERIAL = data_from_discharge
        mi.UPSTREAM_FLUID_PHASE_FROM_PWS_FLASH_CALC  = flashresult.fluid_phase
        mi.UPSTREAM_FLASH_RESULT = flashresult.__dict__
        mi.DOWNSTREAM_PWS_DISCHARGE_RESULT_INITIAL = discharge.vesselLeakCalculation.discharge_records[0].__dict__
        mi.DOWNSTREAM_PWS_DISCHARGE_RESULT_INITIAL_STATE = discharge.vesselLeakCalculation.discharge_records[0].final_state.__dict__

        fin_state = discharge.vesselLeakCalculation.discharge_records[0].final_state
        lf = fin_state.liquid_fraction

        downstream_phase = 'LIQUID'

        if lf == 0:
            downstream_phase = 'VAPOR'

        if lf < 1 and lf > 0:
            downstream_phase = 'TWO-PHASE'

        mi.DATA_FOR_FLOW_BASED_EVALUATION = {
            'discharge_temp_deg_c' : data_from_discharge['discharge_temp_c'],
            'discharge_press_psig' : data_from_discharge['discharge_press_psig'],
            'estimate_mass_flow_lb_hr': mi.ESTIMATE_MASS_FLOW_LB_HR,
            'discharge_phase': downstream_phase
        }

        # sending info around the process prior to initiating main loop dispersion modeling and analysis
        if mi.DISPLAY_DIAGNOSTIC_INFORMATION:
            helpers.send_model_inputs_to_log_handler(mi)

        self.indoor_modeling_objects = []

        dose_response_exclusion_zones = []

        if mi.QUICK_TEST:
            # for models identified as quick tests, only model one release duration.  this method caps the release duration at 14 minutes.
            if not mi.CATASTROPHIC_VESSEL_FAILURE:
                mi.STORAGE_MASS_KG = min(mi.ESTIMATE_MASS_FLOW_KG_S * 14 * 60, mi.STORAGE_MASS_KG)
                mi.STORAGE_VOLUME_M3 = mi.STORAGE_MASS_KG / self.rho_kgm3
        
        if mi.GET_PHAST_DISCHARGE_ONLY:
            release_duration_sec = self.release_durations_to_model_sec[-1]
            self.phast_discharge[release_duration_sec] = discharge
            return ResultCode.SUCCESS

        # only model the one-hour release size (largest)
        if mi.VALID_HAZARDS[cd.HAZARD_TYPE_VCE] or mi.VALID_HAZARDS[cd.HAZARD_TYPE_PV_BURST]:
            self.release_durations_to_model_sec = [self.release_durations_to_model_sec[-1]]
        
        if mi.MAX_RELEASE_DURATION_MIN is not None:
            self.release_durations_to_model_sec = [mi.MAX_RELEASE_DURATION_MIN]
        
        if mi.VALID_HAZARDS[cd.HAZARD_TYPE_PV_BURST]:
            pv_burst = Pv_Burst_Blast_Calc(phast_discharge = discharge)
            pv_burst.run()
            
            bldgs = self.mi.bldgs

            self.mi.LOG_HANDLER("Pressure Vessel Burst Effects on Buildings")
            for bldg in bldgs:
                if not hasattr(bldg, 'pv_burst_overpressure_psi'):
                    continue
                pv_burst_overpressure_psi = bldg.pv_burst_overpressure_psi
                dist_m = bldg.dist_m
                occupancy = bldg.occupancy
                self.mi.LOG_HANDLER(f'distance: {dist_m} m | occ: {occupancy} | overpressure: {pv_burst_overpressure_psi} psi')
            
            if not mi.VALID_HAZARDS[cd.HAZARD_TYPE_VCE]:
                return ResultCode.SUCCESS

        for release_duration_sec in self.release_durations_to_model_sec:

            # debug xxx apple
            # for this analysis, we'll only consider the amount which can be released in 1 hr.

            release_duration_sec = self.release_durations_to_model_sec[-1]

            # for indoor releases, one indoor release model will be conducted per 
            # release duration.  the exfiltration will be studied in the "indoor_modeling"
            # class.  while that class launches its own model_controller, that model_controller
            # is studying the exfiltrated material.  it sets the indoor release flag to false.

            # the "flag_to_jump_to_next_time_step" is used to jump from the completed phast dispersion
            # used for the exfiltration modeling up to the next time step.  it has to break out of two loops, 
            # so using a flag to do that. 
            
            flag_to_jump_to_next_time_step = False
            # create a reference to "actual_release duration".  
            # This will be used in case the true release duration is shorter than the current studied value, calculated below.

            actual_release_duration_sec = release_duration_sec

            if not mi.CATASTROPHIC_VESSEL_FAILURE:
                # vessel will be sized to fit material such that it will be released over the target duration.
                vessel_details = phast_prep.prep_vessel(
                    mi = mi, 
                    state = state, 
                    material = material, 
                    chems = chems, 
                    flashresult = flashresult, 
                    release_duration_sec = actual_release_duration_sec,
                    mass_flow_kg_s = mass_flow_kg_s
                )
                
                vessel = vessel_details['vessel']

                vessel_size = vessel_details['vessel_size']
                
                self.rho_kgm3 = vessel_size.rho_kgm3

                leak_height_fract_of_vessel = vessel_size.leak_height
                stop_after_this_release = vessel_size.stop_after_this_release
                mass_released_kg = vessel_size.mass_released_kg
                actual_release_duration_sec = vessel_size.release_duration_sec

                leak.hole_height_fraction = leak_height_fract_of_vessel # = phast_prep.prep_leak(mi = None, leak_height_fract_of_vessel = leak_height_fract_of_vessel, leak = leak)

                discharge = None
                if release_duration_sec in self.phast_discharge:
                    discharge = self.phast_discharge[release_duration_sec]
                
                if discharge is None:
                    discharge = Phast_Discharge(
                        leak = leak, 
                        state = state, 
                        material = material, 
                        vessel = vessel, 
                        chems = chems,
                        mi = mi,
                        flashresult = flashresult
                    )

                    discharge.run()

                disch_diag = helpers.discharge_diagnosis(discharge)
                if disch_diag != ResultCode.SUCCESS:
                    self.data = helpers.list_to_json(default_data_list)
                    out_list = copy.deepcopy(default_data_list)
                    def_out_dict = out_list[0]
                    def_out_dict[cd.OUTPUT_DICT_RELEASE_DURATION_SEC] = max(1,int(actual_release_duration_sec))
                    def_out_dict[cd.OUTPUT_DICT_MASS_RELEASED_KG] = max(1,int(mass_released_kg))
                    output_list.extend(out_list)
                    self.data = helpers.list_to_json(output_list)
                    return disch_diag

                if mi.RELIEF_CALC_AVAILABLE:

                    # the diameter sizing to match relief rate will be based on the delta P from the liquid height to the leak point.  
                    # This should not be a significant change between 15/30/etc. minute runs.  
                    # However, if it does get to be over a 5% deviation in the mass rate calc and the target relief rate,
                    # run the "match relief rate" method again to get a better estimate on diameter.

                    mass_flow_kg_s = discharge.vesselLeakCalculation.discharge_records[0].mass_flow
                    if abs(mass_flow_kg_s - mi.RELIEF_CALC_RATE_KG_S) / mi.RELIEF_CALC_RATE_KG_S > Consts.TOLERABLE_RELIEF_RATE_DEVIATION_FROM_VESSEL_CALCS:
                        discharge = thermo_pio.match_relief_rate(discharge = discharge, mi = mi, relief_rate=mi.RELIEF_CALC_RATE_KG_S, rel_duration_sec = actual_release_duration_sec)

                        if disch_diag != ResultCode.SUCCESS:
                            self.data = helpers.list_to_json(default_data_list)
                            out_list = copy.deepcopy(default_data_list)
                            def_out_dict = out_list[0]
                            def_out_dict[cd.OUTPUT_DICT_RELEASE_DURATION_SEC] = max(1,int(actual_release_duration_sec))
                            def_out_dict[cd.OUTPUT_DICT_MASS_RELEASED_KG] = max(1,int(mass_released_kg))
                            output_list.extend(out_list)
                            self.data = helpers.list_to_json(output_list)
                            return disch_diag
                    
                    mass_flow_kg_s = discharge.vesselLeakCalculation.discharge_records[0].mass_flow

                self.phast_discharge[release_duration_sec] = discharge


            self.discharges_arr.append({
                cd.OUTPUT_DICT_RELEASE_DURATION_SEC: actual_release_duration_sec,
                'discharge_model': copy.deepcopy(discharge)
            })

            # determine toxic or flammable limits to use in model.  
            # use exit material from discharge to set locs
            chems.set_limits_of_concern(
                mi = mi,
                flash_results=flashresult, 
                vessel_leak_calc=discharge.vesselLeakCalculation
            )

            if not True in mi.VALID_HAZARDS.values():
                mi.LOG_HANDLER('insufficient hazardous concentration in vapor phase')
                self.data = helpers.list_to_json(default_data_list)
                return ResultCode.FAIL_VALIDATION

            output_dict = {}

            ca = Conseq_Assess(chems=chems, mi=mi)
            acr = All_Class_Results(mi)
            
            if mi.USE_ONE_MET_CONDITION_WORST_CASE:
                wx_to_eval = [cd.WX_WORST_CASE]
            else:
                wx_to_eval = cd.WX_ALL_TYPES

            # while nighttime conditions will spread a cloud further, 
            # daytime conditions have the chance to affect more people.
            # the difference in flammable zone between day and night should be negligible, and 
            # styding potential impact to a greater number of people should be the greater concern.
                
            if mi.VALID_HAZARDS[cd.HAZARD_TYPE_VCE] or mi.VALID_HAZARDS[cd.HAZARD_TYPE_PV_BURST]:
                wx_to_eval = [cd.WX_DAY]

            # debug xxx apple
            # for this analysis, we'll maintain all "ca" objects. they contain the footprint results for all models.
            # we will also only analyze the nighttime condition to reduce runtime
            
            self.impact_areas_m2_dict[release_duration_sec] = {}
            wx_to_eval = [cd.WX_WORST_CASE]
            
            for wx in wx_to_eval:

                if flag_to_jump_to_next_time_step:
                    break
                
                self.impact_areas_m2_dict[release_duration_sec][wx] = {}
                for haz in cd.HAZARD_ALL_TYPES:
                    if flag_to_jump_to_next_time_step:
                        break
                    if not mi.VALID_HAZARDS[haz]:
                        continue
                    # vce models are dispersed using the flash fire hazard type and studied below.
                    if haz == cd.HAZARD_TYPE_VCE:
                        continue

                    phast_dispersion = None

                    if release_duration_sec in self.phast_disp:
                        p_disp_for_all_wx_condits = self.phast_disp[release_duration_sec]
                        if wx in p_disp_for_all_wx_condits:
                            phast_dispersion_for_all_hazards = p_disp_for_all_wx_condits[wx]
                            if haz in phast_dispersion_for_all_hazards:
                                phast_dispersion = phast_dispersion_for_all_hazards[haz]
                        else:
                            self.phast_disp[release_duration_sec][wx] = {}
                    else:
                        self.phast_disp[release_duration_sec] = {}
                        self.phast_disp[release_duration_sec][wx] = {}


                    if phast_dispersion is None:
                        phast_dispersion = Phast_Dispersion(
                            phast_discharge = discharge,
                            chems = chems,
                            mi = mi,
                            flashresult= flashresult,
                            release_duration_sec=actual_release_duration_sec,
                            hazard_type=haz,
                            wx_enum = wx,
                            mc = self
                        )

                        phast_dispersion.run()

                        self.phast_disp[release_duration_sec][wx][haz] = phast_dispersion

                    if not phast_dispersion.has_vapor:

                        out_list = copy.deepcopy(default_data_list)
                        def_out_dict = out_list[0]
                        def_out_dict[cd.OUTPUT_DICT_RELEASE_DURATION_SEC] = max(1,int(actual_release_duration_sec))
                        def_out_dict[cd.OUTPUT_DICT_MASS_RELEASED_KG] = max(1,int(mass_released_kg))
                        output_list.extend(out_list)
                        self.data = helpers.list_to_json(output_list)

                        return ResultCode.FAIL_VALIDATION
                    
                    self.impact_areas_m2_dict[release_duration_sec][wx][haz] = copy.deepcopy(phast_dispersion.areas_m2)

                    continue

                    if mi.VALID_HAZARDS[cd.HAZARD_TYPE_VCE]:
                        p_disp = self.phast_disp[release_duration_sec][cd.WX_DAY][cd.HAZARD_TYPE_FLASH_FIRE]
                        vce = VCE(phast_dispersion=p_disp)
                        self.vce_data = vce.get_overall_flammable_envelope_and_maximum_downwind_extent()
                        self.vce = vce
                        return ResultCode.SUCCESS


                    if mi.RELEASE_INDOORS and not mi.VALID_HAZARDS[cd.HAZARD_TYPE_VCE]:
                        # indoor release modeling is handled within the phast_dispersion object.
                        # once the phast dispersion object has returned control here, this can capture the data
                        # and jump to the next time step.

                        # indoor model is not run in conjunction with VCE.
                        # indoor modeling is for exfiltration of vapor from the building, not analyzing consequence within the building.
                        # for indoor models where VCE is a concern, the model will skip the exfiltration analysis and will allow you to perform VCE modeling at a later point.
                        
                        flag_to_jump_to_next_time_step = True
                        # store data from phast dispersion around model results from indoor model
                        self.indoor_modeling_objects.append(phast_dispersion.indoor_model)
                        break

                    self.timings.extend(phast_dispersion.timings)

                    self.dispersion_dicts.append({
                        cd.WEATHER_CONDITION: wx,
                        cd.KEYS_TARG_AND_TYPE_FLAM_OR_INHAL: haz,
                        cd.OUTPUT_DICT_RELEASE_DURATION_SEC:  actual_release_duration_sec,
                        cd.OUTPUT_DISPERSION_OBJECT: phast_dispersion,
                    })

                    ca.assessment(
                        dispersion=phast_dispersion,
                        mass_released_kg=max(1,int(mass_released_kg)),
                        wx_enum = wx,
                        hazard_type=haz,
                        duration_sec=release_duration_sec
                    )

                    if mi.DNV_BENCHMARK:
                            # for testing against benchmarks from upstream-fed models.  
                            # the benchmarks are configured as quick tests, so there will only be one dispersion model
                            
                        dnv_benchmark_specs = self.mi.DNV_BENCHMARK_SPECS
                        if len(dnv_benchmark_specs) > 0:
                            if 'get_results' in dnv_benchmark_specs:
                                get_results = dnv_benchmark_specs['get_results']
                                if get_results:
                                    self.dnv_benchmark_results = DNV_Test_Results(phast_disp = phast_dispersion, conseq_assess = ca)
                                    self.dnv_benchmark_results.run()

                        # debug apple xxx - shorter test. only need the concentration-based output
                        return ResultCode.SUCCESS

                    if ca.dose_probit_data_for_kml_zxy is not None:
                        dose_response_exclusion_zones.append({
                            cd.OUTPUT_DICT_RELEASE_DURATION_SEC: actual_release_duration_sec,
                            cd.WEATHER_CONDITION: wx,
                            cd.KEYS_TARG_AND_TYPE_FLAM_OR_INHAL: haz,
                            cd.DOSE_PROBIT_DATA_FOR_KML_ZXY: ca.dose_probit_data_for_kml_zxy
                        })

                    acr.assessment(ca=ca, haz_type = haz, wx_enum = wx)

                    interim_results = Interim_Results_to_Log_Handler(
                        release_duration_sec=actual_release_duration_sec, 
                        mass_released_kg=mass_released_kg, 
                        wx_enum = wx,
                        results=acr.results,
                        haz_type=haz,
                        id_string=mi.id_string,
                        log_handler=mi.LOG_HANDLER)
                    interim_results.output_to_handler()

                    pnd = phast_dispersion.pnd

                    debug_output_dict = {
                        cd.OUTPUT_DICT_RELEASE_DURATION_SEC: max(1,int(actual_release_duration_sec)),
                        cd.WEATHER_CONDITION: wx,
                        'conc_recs': pnd.debug_conc_records,
                        'analysis_df': pnd.debug_analysis_df
                    }

                    dispersion_plot_dict = {
                        cd.OUTPUT_DICT_RELEASE_DURATION_SEC: max(1,int(actual_release_duration_sec)),
                        cd.WEATHER_CONDITION: wx,
                        'disp_plot': pnd.disp_plot
                    }

                    debug_output_list.append(
                        copy.deepcopy(debug_output_dict)
                    )

                    dispersion_plot_list.append(
                        copy.deepcopy(dispersion_plot_dict)
                    )

            output_dict = {
                cd.OUTPUT_DICT_RELEASE_DURATION_SEC: max(1,int(actual_release_duration_sec)),
                cd.OUTPUT_DICT_MASS_RELEASED_KG: max(1,int(mass_released_kg)),
                cd.OUTPUT_DICT_HAZARD_RECS: acr.results,
            }

            output_list.append(copy.deepcopy(output_dict))
            
            
            
            if stop_after_this_release:
                break
        
        if len(dose_response_exclusion_zones) > 0:
            self.dose_response_exclusion_zones_df = pd.DataFrame(dose_response_exclusion_zones)

        if mi.RELEASE_INDOORS:
            self.parse_indoor_results()
            if len(self.data_not_json) > 0:
                return ResultCode.SUCCESS
            out_list = copy.deepcopy(default_data_list)
            def_out_dict = out_list[0]
            def_out_dict[cd.OUTPUT_DICT_RELEASE_DURATION_SEC] = Consts.RELEASE_DURATIONS_TO_MODEL_SEC[0]
            def_out_dict[cd.OUTPUT_DICT_MASS_RELEASED_KG] = max(1,mi.STORAGE_MASS_KG)
            self.data = helpers.list_to_json([def_out_dict])
            return ResultCode.FAIL_VALIDATION

        self.timings = pd.DataFrame(self.timings)
        self.data = helpers.list_to_json(output_list)
        self.debug_data = helpers.list_to_json(debug_output_list)
        self.dispersion_plots = helpers.list_to_json(dispersion_plot_list)
        self.output_list = output_list

        return ResultCode.SUCCESS

    def parse_indoor_results(self):
        timings = []
        data = []
        debug_data = []
        dispersion_plots = []
        dispersion_dicts = []
        for im in self.indoor_modeling_objects:
            if im.resultCode != ResultCode.SUCCESS:
                continue
            timings.append(im.mc.timings)
            curr_data = json.loads(im.mc.data)
            data.append(curr_data)
            curr_debug_data = json.loads(im.mc.debug_data)
            debug_data.append(curr_debug_data)
            curr_dispersion_plots = json.loads(im.mc.dispersion_plots)
            dispersion_plots.extend(curr_dispersion_plots)
            dispersion_dicts.extend(im.mc.dispersion_dicts)
        timings = copy.deepcopy(timings)
        data = copy.deepcopy(data)
        debug_data = copy.deepcopy(debug_data)
        dispersion_plots = copy.deepcopy(dispersion_plots)
        self.timings = pd.concat(timings)
        self.data = helpers.list_to_json(data)
        self.debug_data = helpers.list_to_json(debug_data)
        self.dispersion_dicts = copy.deepcopy(dispersion_dicts)
        
        self.data_not_json = data

        self.dispersion_plots = []
        try:
            self.dispersion_plots = helpers.list_to_json(dispersion_plots)
        except:
            self.mi.LOG_HANDLER('there was an issue storing the dispersion plot output from this run')
        
from py_lopa.phast_io.phast_dispersion import Phast_Dispersion







