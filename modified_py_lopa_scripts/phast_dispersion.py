import sys
import copy
import jsons

import datetime
from datetime import datetime as dt
import pandas as pd
import numpy as np

from pypws.entities import DispersionOutputConfig, FlashResult, DispersionParameters, MixtureModelling, ContourType
from pypws.calculations import DispersionCalculation, _DispersionCalculationRequest, _DispersionCalculationRequestSchema, VesselLeakCalculation
from pypws.calculations import MaxConcDistanceCalculation, MaxConcFootprintCalculation, DistancesAndFootprintsToConcentrationLevelsCalculation
from pypws.calculations import _MaxConcDistanceCalculationRequest, _MaxConcFootprintCalculationRequest, _MaxConcDistanceCalculationRequestSchema, _MaxConcFootprintCalculationRequestSchema
from pypws.calculations import _DistancesAndFootprintsToConcentrationLevelsCalculationRequest, _DistancesAndFootprintsToConcentrationLevelsCalculationRequestSchema
from pypws.enums import ResultCode, SpecialConcentration, Resolution

from py_lopa.calcs import helpers
from py_lopa.data.exception_enum import Exception_Enum
from py_lopa.calcs.consts import Consts, Wx_Enum
from py_lopa.classes.disp_analysis_prep import Disp_Analysis_Prep
from py_lopa.classes.bldg_heights_to_eval import Bldg_Heights_To_Evaluate
from py_lopa.phast_io import phast_prep
from py_lopa.classes.plotting_and_debug_data import Plotting_and_Debug_Data

from py_lopa.calcs.flattening import Flattening

cd = Consts().CONSEQUENCE_DATA

class Phast_Dispersion:

    # evaluate impact at different elevations to find worst case profile.
    
    def __init__(self, phast_discharge, chems = None, mi = None, flashresult:FlashResult = None, release_duration_sec = 3600, hazard_type = cd.HAZARD_TYPE_FLASH_FIRE, wx_enum = cd.WX_WORST_CASE, mc = None):
        self.testing = []
        self.timings = []
        self.resolution = Resolution.MEDIUM
        self.phast_discharge = phast_discharge
        self.chems = chems
        self.mi = mi
        self.vce = mi.VAPOR_CLOUD_EXPLOSION
        self.wx_enum = wx_enum
        self.weather = phast_prep.prep_weather(wx_enum=wx_enum)
        self.substrate = phast_prep.prep_substrate(mi=self.mi)
        self.flashresult = flashresult
        self.min_z = 1e6
        self.release_duration_sec = release_duration_sec
        self.hazard_type = hazard_type
        self.mc = mc
        self.zxyc_df = {}
        
        self.disp_recs = {}
        self.conc_profiles = []
        self.haz_cat_conc_footprints_df = {}
        self.distances_and_footprints_results = []
        self.analysis_df = {}
        self.has_vapor = True
        bh = Bldg_Heights_To_Evaluate(mi = self.mi)
        self.study_bldg_hts = bh.calc_elevs()
        self.elevs_m = []
        self.dispersion_calc_request = (
            _DispersionCalculationRequest,
            _DispersionCalculationRequestSchema,
        )

        self.max_conc_distance_calc_request = (
            _MaxConcDistanceCalculationRequest,
            _MaxConcDistanceCalculationRequestSchema,
        )

        self.max_conc_footprint_calc_request = (
            _MaxConcFootprintCalculationRequest,
            _MaxConcFootprintCalculationRequestSchema,
        )

        self.distances_and_footprints_to_concentration_levels_calculation_request = (
            _DistancesAndFootprintsToConcentrationLevelsCalculationRequest,
            _DistancesAndFootprintsToConcentrationLevelsCalculationRequestSchema,
        )
        self.cat_and_targ_concs_dict = None

        self.requests = []
        self.footprints_conc_elev_z_x_y_list = []
        self.footprints_conc_elev_z_x_y_df = {}
        self.areas_m2 = {}
    
    def run(self):
        
        self.prep_disp_params_set_averaging_time_and_concentration_endpoint()
        self.prep_dispersion()
        self.prep_dispersion_analysis()
        self.phast_request(calc_obj = self.dispersionCalculation, calc_name='Dispersion', calc_request=self.dispersion_calc_request)
        if self.dispersionCalculation.result_code != ResultCode.SUCCESS:
            raise Exception(Exception_Enum.DISPERSION_CALCULATION_ERROR)
        if self.dispersionCalculation.scalar_udm_outputs.record_count <= 0:
            self.has_vapor = False
            if self.mi is not None:
                self.mi.LOG_HANDLER('no appreciable vapor formation from the release.')
            return
        self.store_pool_records()
        
        if self.mi.RELEASE_INDOORS:
            self.indoor_model = Indoor_Modeling_3(phast_disp=self)
            self.indoor_model.run()
            return
        
        if self.vce:
            return

        self.get_elevs()

        self.define_categories_and_target_concentrations()

        self.batch_call_footprints()

        # self.max_conc_at_distance()

        # self.max_conc_footprint()

        apple = 1

    def define_categories_and_target_concentrations(self):
        haz_type = self.hazard_type
        self.cat_and_targ_concs_dict = self.get_targ_concs(haz_type)

    def batch_call_footprints(self):

        # get footprints at all required concs and elevations
        
        targ_concs = list(self.cat_and_targ_concs_dict.values())
        
        self.distancesAndFootprintsCalc = DistancesAndFootprintsToConcentrationLevelsCalculation(
            scalar_udm_outputs = self.dispersionCalculation.scalar_udm_outputs, 
            weather = self.weather, 
            dispersion_records = self.dispersionCalculation.dispersion_records, 
            dispersion_record_count = len(self.dispersionCalculation.dispersion_records), 
            substrate = self.substrate, 
            dispersion_output_configs = [], 
            dispersion_output_config_count = None, 
            dispersion_parameters = self.dispersionCalculation.dispersion_parameters, 
            material  = self.phast_discharge.vesselLeakCalculation.exit_material
        )

        curr_conc = -1
        targ_concs = list(set(targ_concs))
        targ_concs = sorted(targ_concs)
        targ_concs = targ_concs[-2:]
        for targ_conc in targ_concs:
            for elev in self.elevs_m:
                dispOutputCfg = DispersionOutputConfig()
                dispOutputCfg.resolution = Resolution.LOW
                # if targ_conc == min(targ_concs):
                #     # the low conc is used to identify impacts to buildings.  this needs to be high res.
                #     dispOutputCfg.resolution = Resolution.HIGH
                # if targ_conc == max(targ_concs):
                #     # serious+ needs medium resolution for accuracy
                #     dispOutputCfg.resolution = Resolution.MEDIUM
                dispOutputCfg.downwind_distance = np.inf
                dispOutputCfg.special_concentration = SpecialConcentration.NOT_DEFINED
                dispOutputCfg.concentration = targ_conc
                dispOutputCfg.elevation = elev
                dispOutputCfg.contour_type = ContourType.FOOTPRINT
                self.distancesAndFootprintsCalc.dispersion_output_configs.append(dispOutputCfg)

        self.mi.LOG_HANDLER('\n***\n\nInitiating Model:  Footprint Analysis at Elevation')
        t0 = dt.now(datetime.UTC)
        self.distancesAndFootprintsCalc.dispersion_output_config_count = len(self.distancesAndFootprintsCalc.dispersion_output_configs)
        res = self.distancesAndFootprintsCalc.run()

        log_msg = f'Model run successful.  run time: {dt.now(datetime.UTC) - t0} sec'
        if res != ResultCode.SUCCESS:
            log_msg = f'\n\nIssue with flammable envelope calc.  error messages:  {self.distancesAndFootprintsCalc.messages}'
        
        self.mi.LOG_HANDLER(log_msg)

        self.parse_batch_call_footprints()

        return res
        
    def parse_batch_call_footprints(self):
        dists_and_footprints_calc = self.distancesAndFootprintsCalc
        curr_pt = 0
        self.footprints_conc_elev_z_x_y_list = []
        curr_disp_output_config = -1
        for num_pts in dists_and_footprints_calc.n_contour_points:
            curr_disp_output_config += 1
            disp_output_config = dists_and_footprints_calc.dispersion_output_configs[curr_disp_output_config]
            if num_pts == 0:
                continue
            cps = dists_and_footprints_calc.contour_points[curr_pt:curr_pt + num_pts]
            for cp in cps:
                if abs(cp.x) > 1e10:
                    continue
                self.footprints_conc_elev_z_x_y_list.append({
                    'conc_ppm': disp_output_config.concentration * 1e6,
                    'elev_m': disp_output_config.elevation,
                    'x': cp.x,
                    'y': cp.y,
                    'z': cp.z,
                })
            curr_pt += num_pts
        
        if len(self.footprints_conc_elev_z_x_y_list) == 0:
            return
        self.footprints_conc_elev_z_x_y_df = pd.DataFrame(self.footprints_conc_elev_z_x_y_list)
        df_for_areas = self.footprints_conc_elev_z_x_y_df[self.footprints_conc_elev_z_x_y_df['z'] <= 6]
        grouped = df_for_areas.groupby('conc_ppm')
        self.areas_m2 = []
        flattening = Flattening()
        for conc_ppm, df in grouped:
            zxy = df[['z', 'x', 'y']].values
            self.areas_m2.append({
                'conc_ppm': conc_ppm,
                'area_m2': flattening.get_zxy_area(zxy=zxy)
            })
        apple = 1





    def prep_disp_params_set_averaging_time_and_concentration_endpoint(self):
        
        self.ep_conc = -1

        averaging_time_sec = Consts.AVERAGING_TIME_SEC[cd.HAZARD_TYPE_FLASH_FIRE]

        haz = self.hazard_type

        if haz == cd.HAZARD_TYPE_INHALATION:
            if self.mi.VALID_HAZARDS[cd.HAZARD_TYPE_INHALATION]:
                self.ep_conc = self.chems.inhalation_locs[1]
                averaging_time_sec = Consts.AVERAGING_TIME_SEC[cd.HAZARD_TYPE_INHALATION]

        if haz == cd.HAZARD_TYPE_FLASH_FIRE:
            if self.mi.VALID_HAZARDS[cd.HAZARD_TYPE_FLASH_FIRE]:
                self.ep_conc = self.chems.flam_loc_fracts[1]
            if self.vce:
                self.ep_conc = self.chems.flam_loc_fracts[3]

        self.dispParams = DispersionParameters()
        self.dispParams.averaging_time = averaging_time_sec

        # per Maria Fernandez (DNV - 2023 12 13) the dispersion relativeTolerance can cause 
        # solvers to not converge if set incorrectly.  for one case, it caused the target
        # tables to crash.  will test at 1e-4 per their guidance.  will adjust as needed 
        # post-testing.

        self.dispParams.relative_tolerance = 1e-4

        # there is a current limitation in discharge calculations such that any subcooled system cannot be modeled using a multicomponent method.  
        # for subcooled systems, the discharge calculation reverts to pseudocomponent modeling.
        # per Maria Fernandez (DNV) 9APR2024, once the discharge calculation is complete, the discharge record object contains a property for mixture modeling 
        # where the pseudocomponent spec can be overrided multicomponent.

        # in testing, the mc method returns many more "not a concern" results for liquid releases.  per Maria, this is a much more accurate method.
        # 19AUG2024 - DNV reported that the "not a concern" records were due to an error with multicomponent modeling.
        # reverting to PC until fix is available.
        
        if self.mi.USE_MULTICOMPONENT_METHOD:
            for rec in self.phast_discharge.vesselLeakCalculation.discharge_records:
                fin_state = rec.final_state
                setattr(fin_state, 'mixture_modelling', MixtureModelling.MC_SINGLE_AEROSOL)
        
        if self.ep_conc <= 0:
            raise Exception(Exception_Enum.HAZ_MATS_NOT_SIGNIFICANT_IN_VAPOR_PHASE)


    def prep_dispersion(self):

        phast_disch = self.phast_discharge

        self.dispersionCalculation = DispersionCalculation(material=None, substrate=None, discharge_result=None, discharge_records=None, discharge_record_count=None, weather=None, dispersion_parameters=None, end_point_concentration=None)
        self.dispersionCalculation.discharge_records = phast_disch.vesselLeakCalculation.discharge_records
        self.dispersionCalculation.discharge_record_count = len(self.dispersionCalculation.discharge_records)
        self.dispersionCalculation.discharge_result = phast_disch.vesselLeakCalculation.discharge_result
        self.dispersionCalculation.material = phast_disch.vesselLeakCalculation.exit_material
        self.dispersionCalculation.substrate = self.substrate
        self.dispersionCalculation.weather = self.weather
        self.dispersionCalculation.end_point_concentration = self.ep_conc

        self.dispersionCalculation.dispersion_parameters = self.dispParams

    def get_elevs(self):

        bldg_elevs_m = self.study_bldg_hts
        onsite_elevs_m = copy.deepcopy(Consts.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS)
        if self.mi is not None:
            onsite_elevs_m = self.mi.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS

        vert_tol_m = 2
        elevs_not_needed = []
        for b_e in bldg_elevs_m:
            for os_e in onsite_elevs_m:
                if abs(b_e - os_e) <= vert_tol_m:
                    elevs_not_needed.append(b_e)
        elevs_m = list(set(bldg_elevs_m) - set(elevs_not_needed))
        elevs_m.extend(onsite_elevs_m)
        elevs_m.sort()

        self.elevs_m = elevs_m

    def phast_request(self, calc_obj, calc_name = '', calc_request = None, calc_request_set = False):

        t0 = dt.now(datetime.UTC)
        if not isinstance(calc_name, str):
            calc_name = ''

        if calc_request is not None:
            calc_request_obj = calc_request[0]
            calc_request_schema = calc_request[1]
            req_obj = None
            if calc_request_obj is _DispersionCalculationRequest:
                
                req_obj = calc_request_obj(
                    material = calc_obj.material, 
                    substrate = calc_obj.substrate,
                    discharge_result = calc_obj.discharge_result, 
                    discharge_records = calc_obj.discharge_records, 
                    discharge_record_count = calc_obj.discharge_record_count, 
                    weather = calc_obj.weather, 
                    dispersion_parameters = calc_obj.dispersion_parameters, 
                    end_point_concentration = calc_obj.end_point_concentration
                )

            elif calc_request_obj is _DistancesAndFootprintsToConcentrationLevelsCalculationRequest:
                
                req_obj = calc_request_obj(
                    scalar_udm_outputs = calc_obj.scalar_udm_outputs, 
                    weather = calc_obj.weather, 
                    dispersion_records = calc_obj.dispersion_records, 
                    dispersion_record_count = calc_obj.dispersion_record_count, 
                    substrate = calc_obj.substrate, 
                    material = calc_obj.material, 
                    dispersion_parameters = calc_obj.dispersion_parameters,
                    dispersion_output_configs = calc_obj.dispersion_output_configs,
                    dispersion_output_config_count = calc_obj.dispersion_output_config_count
                )
            
            else:
                # falls back to max conc dist and max conc footprint requests
                
                req_obj = calc_request_obj(
                    scalar_udm_outputs = calc_obj.scalar_udm_outputs, 
                    weather = calc_obj.weather, 
                    dispersion_records = calc_obj.dispersion_records, 
                    dispersion_record_count = calc_obj.dispersion_record_count, 
                    substrate = calc_obj.substrate, 
                    material = calc_obj.material, 
                    dispersion_parameters = calc_obj.dispersion_parameters,
                    dispersion_output_config = calc_obj.dispersion_output_config,
                )

            req_schema = calc_request_schema()

            request_json = req_schema.dumps(req_obj)
            calc_obj.request_json_tt = request_json

            self.requests.append(request_json)
        
        # for any benchamrking where only requests are needed, we can return after storing the requests.  
        # however, the dispersion model will need to be run so that the post-processing model requests 
        # can also be generated.
        if calc_request_obj is not _DispersionCalculationRequest:
            get_results = helpers.get_spec_from_benchmark_specs_dict(mi=self.mi, spec = 'get_results')
            if get_results is not None:
                if not get_results:
                    return

        resultCode = ResultCode.FAIL_EXECUTION

        # if self.mi.QUICK_TEST and calc_request_obj is not _DispersionCalculationRequest:
        #     return

        if calc_name != '':
            self.mi.LOG_HANDLER(f'\n***\n\nInitiating Model:  {calc_name}')
        for i in range(5):
            try:
                resultCode = calc_obj.run()
            except:
                continue
            break
        
        # this will output to local for any calculations that have messages.  
        # This may be duplicated for some models, with some data going to log handler.
        # it will not show in the target table web output
        if hasattr(calc_obj, 'messages'):
            if calc_obj.messages is not None:
                if len(calc_obj.messages) > 0:
                    print(f'model:  {calc_name}.  messages:  {calc_obj.messages}')
        
        if resultCode != ResultCode.SUCCESS:

            return

        if self.mi is not None:
            t1 = dt.now(datetime.UTC)
            runtime = t1-t0
            runtime = runtime.total_seconds()
            rel_dur_min = max(1, int(self.release_duration_sec / 60))
            self.timings.append({
                cd.KEYS_TARG_AND_TYPE_FLAM_OR_INHAL:  self.hazard_type,
                'description': calc_name,
                'release duration':  rel_dur_min,
                'runtime': runtime
            })
            self.mi.LOG_HANDLER(f'***Weather Condition:  {self.wx_enum}.  Hazard Type: {self.hazard_type}.  Model: {calc_name}.  Release Duration: {rel_dur_min} min. Runtime: {runtime} sec. Current time: {t1} UTC. Status: Completed OK.\n')


    def phast_request_wrapper(self, calc, calc_descr, conc_targ, output_list_of_dicts, stored_data_descr, elevs_m, resolution = Resolution.MEDIUM, calc_request=None):
        phast_disch = self.phast_discharge
        #(scalar_udm_outputs=None, weather=None, dispersion_records=None, dispersion_record_count=None, substrate=None, dispersion_output_config=None, material=None, dispersion_parameters=None)
        calc.weather = self.weather
        calc.dispersion_records = self.dispersionCalculation.dispersion_records
        calc.dispersion_record_count = len(calc.dispersion_records)
        calc.substrate = self.substrate
        calc.material = phast_disch.vesselLeakCalculation.exit_material
        calc.scalar_udm_outputs = self.dispersionCalculation.scalar_udm_outputs
        calc.dispersion_parameters = self.dispersionCalculation.dispersion_parameters

        dispOutputCfg = DispersionOutputConfig()
        dispOutputCfg.resolution = resolution
        dispOutputCfg.downwind_distance = np.inf
        dispOutputCfg.special_concentration = SpecialConcentration.NOT_DEFINED
        dispOutputCfg.concentration = conc_targ

        for elev_m in elevs_m:
            dispOutputCfg.elevation = elev_m
            calc.dispersion_output_config = dispOutputCfg
            conc_ppm = conc_targ * 1e6
            calc_name = f'{calc_descr} near {int(elev_m)} m.  Conc: {conc_ppm} ppm'
            self.phast_request(calc_obj = calc, calc_name= calc_name, calc_request=calc_request)
            # if self.mi.QUICK_TEST:
            #     output_list_of_dicts.append({
            #         cd.CONC_CALC_ELEV_M: elev_m, 
            #         stored_data_descr: copy.deepcopy(calc)
            #     })
            #     break
            calc_copy = copy.deepcopy(calc)
            output_list_of_dicts.append({
                cd.CONC_CALC_ELEV_M: elev_m, 
                stored_data_descr: calc_copy,
                'successful_run': True
            })

            if calc.result_code == ResultCode.SUCCESS:
                # create a dataframe of concentrations at points downwind at each studied elevation.
                if isinstance(calc_copy, MaxConcDistanceCalculation):
                    conc_pfls = output_list_of_dicts
                    if len(conc_pfls) == 0:
                        continue
                    flat = Flattening(conc_pfls=conc_pfls)
                    zxyc_np = flat.get_zxyc_array_from_conc_pfls_bet_min_and_max_elevation(min_ht_m=0, max_ht_m=np.inf)
                    calc_copy.zxyc_df = pd.DataFrame([], columns=['z', 'x', 'y', 'c'])
                    if len(zxyc_np) == 0:
                        continue
                    calc_copy.zxyc_df = pd.DataFrame(zxyc_np, columns=['z', 'x', 'y', 'c'])

                    apple = 1
                    
            else:
                output_list_of_dicts[-1]['successful_run'] = False
                self.mi.LOG_HANDLER('\n\n----------------Warning returned during Post-Processing of Dispersion Model, but proceeded normally.\n')
                if hasattr(calc, 'messages'):
                    if calc.messages is not None:
                        if len(calc.messages) > 0:
                            self.mi.LOG_HANDLER(f'Phast Web Services returned the following warning: "{calc.messages[0]}"\n')
                            self.mi.LOG_HANDLER('----------------\n\n\n')
    
    def run_footprint_models_for_vce(self, vce_targ_concs):
        
        self.distancesAndFootprintsCalc = DistancesAndFootprintsToConcentrationLevelsCalculation(
            scalar_udm_outputs = self.dispersionCalculation.scalar_udm_outputs, 
            weather = self.weather, 
            dispersion_records = self.dispersionCalculation.dispersion_records, 
            dispersion_record_count = len(self.dispersionCalculation.dispersion_records), 
            substrate = self.substrate, 
            dispersion_output_configs = [], 
            dispersion_output_config_count = None, 
            dispersion_parameters = self.dispersionCalculation.dispersion_parameters, 
            material  = self.phast_discharge.vesselLeakCalculation.exit_material
        )

        for elev in range(51):
            for targ_conc in vce_targ_concs:
                dispOutputCfg = DispersionOutputConfig()
                dispOutputCfg.resolution = Resolution.LOW
                dispOutputCfg.downwind_distance = np.inf
                dispOutputCfg.special_concentration = SpecialConcentration.NOT_DEFINED
                dispOutputCfg.concentration = targ_conc
                dispOutputCfg.elevation = elev
                dispOutputCfg.contour_type = ContourType.FOOTPRINT
                self.distancesAndFootprintsCalc.dispersion_output_configs.append(dispOutputCfg)

        self.mi.LOG_HANDLER('\n***\n\nInitiating Model:  VCE Flammable Envelope')
        t0 = dt.now(datetime.UTC)
        self.distancesAndFootprintsCalc.dispersion_output_config_count = len(self.distancesAndFootprintsCalc.dispersion_output_configs)
        res = self.distancesAndFootprintsCalc.run()

        log_msg = f'Model run successful.  run time: {dt.now(datetime.UTC) - t0} sec'
        if res != ResultCode.SUCCESS:
            log_msg = f'\n\nIssue with flammable envelope calc.  error messages:  {self.distancesAndFootprintsCalc.messages}'
        
        self.mi.LOG_HANDLER(log_msg)

        return res


    def max_conc_at_distance(self):
        
        self.maxConcDistCalc = MaxConcDistanceCalculation(scalar_udm_outputs=None, weather=None, dispersion_records=None, dispersion_record_count=None, substrate=None, dispersion_output_config=None, material=None, dispersion_parameters=None)
        
        self.phast_request_wrapper(calc = self.maxConcDistCalc, calc_descr='Max Conc Distance at Elevation', conc_targ=self.ep_conc, output_list_of_dicts=self.conc_profiles, stored_data_descr=cd.CONC_CALC_CONC_PFL_CALC, elevs_m = self.elevs_m, resolution=Resolution.HIGH, calc_request=self.max_conc_distance_calc_request)

        # if len(self.conc_profiles) == 0:
        #     raise Exception(Exception_Enum.DISPERSION_CALCULATION_ERROR)

    def prep_dispersion_analysis(self):
        
        if len(self.analysis_df) == 0:
            self.disp_analysis_prep = Disp_Analysis_Prep(mi = self.mi, chems = self.chems)
            self.analysis_df = self.disp_analysis_prep.analysis_df
    
    def store_pool_records(self):

        self.pool_data_df = {}
        pool_recs = []
        if hasattr(self.dispersionCalculation, 'pool_records'):
            pool_recs = self.dispersionCalculation.pool_records

        if len(pool_recs) > 0:
            pool_data = [{'time_sec': x.time, 'spill_rate_kg_per_sec': x.spill_rate, 'vap_rate_kg_per_sec':x.vapourisation_rate, 'temp_k': x.temperature, 'mass_spilled_kg': x.mass_spilt} for x in pool_recs]
                
            self.pool_data_df = pd.DataFrame(pool_data)

    def max_conc_footprint_for_kml(self, elevs_m):

        # use the loc-1 / loc-2 / loc-3 and, if applicable, 10xloc-3.
        # all potential concs are in analysis_df 'conc_volf' column.

        df = self.analysis_df
        haz_df = df[df[cd.KEYS_TARG_AND_TYPE_FLAM_OR_INHAL] == self.hazard_type]
        concs = haz_df[cd.CONC_VOLF_TITLE]
        concs_unique = concs.unique()
        concs_unique = concs_unique[~pd.isna(concs_unique)]
        self.maxConcFootprintCalc = MaxConcFootprintCalculation(scalar_udm_outputs=None, weather=None, dispersion_records=None, dispersion_record_count=None, substrate=None, dispersion_output_config=None, material=None, dispersion_parameters=None)
        self.maxConcFootprintCalc.dispersion_parameters = self.dispParams
        conc_footprint_dicts = []
        haz_type = self.hazard_type
        for targ in concs_unique:
            self.curr_conc_footprints = []
            self.phast_request_wrapper(calc = self.maxConcFootprintCalc, calc_descr='Max Conc Footprint at Elevation', conc_targ=targ, output_list_of_dicts=self.curr_conc_footprints, stored_data_descr=cd.CONC_CALC_CONC_FOOTPRINT, elevs_m=elevs_m, resolution=Resolution.LOW, calc_request=self.max_conc_footprint_calc_request)
            # if len(self.curr_conc_footprints) == 0:
                #raise Exception(Exception_Enum.DISPERSION_CALCULATION_ERROR)
            conc_footprint_dicts.append({
                    cd.CONC_VOLF_TITLE: targ,
                    cd.CONC_CALC_CONC_FOOTPRINT: copy.deepcopy(self.curr_conc_footprints)
                })
        curr_targ = targ
        
        self.curr_conc_footprints = []
        self.conc_footprints_df = pd.DataFrame(conc_footprint_dicts)

    def max_conc_footprint(self, elevs_m = [], get_only_higher_level_conseq_targs = False):
        
        # get the coordinates of the max footprint of impact to target concentrations
        # the concentrations of interest are the top 3 classes
        # to determine if the impact is serious, critical or catastrophic, the total impacted
        # area to all heights below elev range (currently 6 m) is compared to the provided population

        # study the following elevations:
        # ground level
        # max elevation range (Consts.ELEVATION_RANGE_FOR_CONC_EVAL_M)
        # mean of above two values
        if len(elevs_m) == 0:
            elevs_m = copy.deepcopy(Consts.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS)
            if self.mi is not None:
                elevs_m = self.mi.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS
            if self.vce:
                elevs_m = self.elevs_m

        if self.vce:
            # will restrict the number of models run to only the lfl boundary
            get_only_higher_level_conseq_targs = True

        haz_cat_conc_footprints = []
        self.curr_conc_footprints = []
        curr_targ = -1
        haz_type = self.hazard_type
        cat_and_targ_concs_dict = self.get_targ_concs(haz_type, get_only_higher_level_conseq_targs)
        for cat, targ in cat_and_targ_concs_dict.items():
            if curr_targ == targ:
                haz_cat_conc_footprints.append(
                {
                    cd.KEYS_TARG_AND_TYPE_FLAM_OR_INHAL: haz_type,
                    cd.CAT_TITLE: cat,
                    cd.CONC_CALC_CONC_FOOTPRINT: copy.deepcopy(self.curr_conc_footprints)
                })
                continue

            self.maxConcFootprintCalc = MaxConcFootprintCalculation(scalar_udm_outputs=None, weather=None, dispersion_records=None, dispersion_record_count=None, substrate=None, dispersion_output_config=None, material=None, dispersion_parameters=None)
            self.maxConcFootprintCalc.dispersion_parameters = self.dispParams
            self.curr_conc_footprints = []
            resolution = Resolution.LOW
            # resolution = Resolution.MEDIUM
            # if cat == cd.CAT_MINOR or cat == cd.CAT_MODERATE or self.vce:
            #     resolution = Resolution.LOW
            self.phast_request_wrapper(calc = self.maxConcFootprintCalc, calc_descr='Max Conc Footprint at Elevation', conc_targ=targ, output_list_of_dicts=self.curr_conc_footprints, stored_data_descr=cd.CONC_CALC_CONC_FOOTPRINT, elevs_m=elevs_m, resolution=resolution, calc_request = self.max_conc_footprint_calc_request)
            # if len(self.curr_conc_footprints) == 0:
            #     raise Exception(Exception_Enum.DISPERSION_CALCULATION_ERROR)
            
            haz_cat_conc_footprints.append({
                    cd.KEYS_TARG_AND_TYPE_FLAM_OR_INHAL: haz_type,
                    cd.CAT_TITLE: cat,
                    cd.CONC_CALC_CONC_FOOTPRINT: copy.deepcopy(self.curr_conc_footprints)
                })
            curr_targ = targ

        
        self.curr_conc_footprints = []
        self.haz_cat_conc_footprints_df = pd.DataFrame(haz_cat_conc_footprints)

    def get_targ_concs(self, haz_type, get_only_higher_level_conseq_targs = False):
        df = self.analysis_df
        haz_df = df[df[cd.KEYS_TARG_AND_TYPE_FLAM_OR_INHAL] == haz_type]
        haz_onsite_df = haz_df[haz_df[cd.CLASS_TITLE] == cd.CLASS_ONSITE]
        if get_only_higher_level_conseq_targs:
            haz_onsite_df = haz_onsite_df[haz_onsite_df[cd.KEYS_TARG_AND_TYPE_CATEGORY].isin(cd.CAT_HIGHER_LEVEL_IMPACT)]
        # if self.mi.QUICK_TEST:
        #     # for quick tests, only get requests for moderate and serious for onsite
        #     targs = [cd.CAT_MODERATE, cd.CAT_SERIOUS]
        #     haz_onsite_df = haz_onsite_df[haz_onsite_df[cd.KEYS_TARG_AND_TYPE_CATEGORY].isin(targs)]
        cat_and_targ_concs_df = haz_onsite_df[[cd.KEYS_TARG_AND_TYPE_CATEGORY, cd.KEYS_TARG_AND_TYPE_CONC_VOLF]]
        cat_and_targ_concs_dict = {k: v for (k, v) in zip(cat_and_targ_concs_df[cd.KEYS_TARG_AND_TYPE_CATEGORY], cat_and_targ_concs_df[cd.KEYS_TARG_AND_TYPE_CONC_VOLF])}

        return cat_and_targ_concs_dict
    
    def create_zxyc_df(self):
        conc_pfls = self.conc_profiles
        if len(conc_pfls) == 0:
            return
        flat = Flattening(conc_pfls=conc_pfls)
        zxyc_np = flat.get_zxyc_array_from_conc_pfls_bet_min_and_max_elevation(min_ht_m=0, max_ht_m=np.inf)
        if len(zxyc_np) > 0:
            self.zxyc_df = pd.DataFrame(zxyc_np, columns = ['z', 'x', 'y', 'c'])
    
    def export_zxyc_df(self):
        if len(self.zxyc_df) == 0:
            return
        
        rel_dur_sec = helpers.replace_chars_in_string_prior_to_naming_windows_text_file(str(self.release_duration_sec))
        helpers.df_to_csv(self.zxyc_df, f = f'zxyc_from_disp__wx_{self.wx_enum}__hazard_type_{self.hazard_type}__rel_dur_{rel_dur_sec}_sec')

from py_lopa.classes.indoor_modeling_3 import Indoor_Modeling_3