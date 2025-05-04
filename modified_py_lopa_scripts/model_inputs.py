import sys
import copy
import math
from importlib.metadata import version
import pandas as pd

from py_lopa.data.debug_info import Debug_Info
from py_lopa.data.exception_enum import Exception_Enum
from py_lopa.calcs.consts import Consts
from py_lopa.classes.chems import Chems
import py_lopa.calcs.helpers as helpers
from py_lopa.classes.building import Building

consts = Consts()

cd = consts.CONSEQUENCE_DATA

class Model_Inputs:
    
    def __init__(self , inputs):

        self.inputs = copy.deepcopy(inputs)
        self.prep_keys()
        self.set_identification_string()
        self.clean_inputs()
        self.clean_attributes()
        self.enumerate_relief_phase()
        self.bldg_consequence_targs_dict()
        self.convert_to_internal_units()
        self.model_updates_for_catastrophic_failure()
        self.set_details_for_onsite_impact_evaluation()
        self.adjust_for_quick_test()
        self.prep_chems()
        self.get_version_info()
        self.model_updates_for_flammables()


    def prep_keys(self):

        # units of some input properties are different.  
        # conversions are in method 'convert_to_internal_units' 

        self.CHEM_MIX = 'chemical_mix' 
        self.MAX_HOLE_SZ_M = 'max_hole_size_in'
        self.PRESS_PA = 'pressure_psig'
        self.TEMP_K = 'temp_deg_c'
        self.MOL_FRACTIONS = 'mole_fracts'
        self.MASS_FRACTIONS = 'mass_fracts'
        self.COMPOSITION = 'composition'
        self.COMPOSITION_IS_IN_MOLES = 'comp_is_moles'
        self.RELEASE_ELEVATION_M = 'release_elevation_m'
        self.INHALATION = consts.INHALATION
        self.FLASH_FIRE = consts.FLASH_FIRE
        self.VAPOR_CLOUD_EXPLOSION = consts.VCE
        self.PV_BURST = consts.PV_BURST
        self.BLDG_NUMS = 'bldg_nums_low_med_high'
        self.BLDG_DISTS = 'bldg_dists_low_med_high'
        self.BLDG_HTS = 'bldg_hts_low_med_high'
        self.BLDG_ACH = 'bldg_ach_low_med_high'
        self.BLDG_NIGHT_OCC = 'bldg_nighttime_occupancy'
        self.BLDG_HVAC_SHUTOFF = 'bldg_hvac_shutoff_min_low_med_high'
        self.OFFSITE_DIST = 'offsite_dist_m'
        self.RELIEF_PHASE = 'relief_phase'
        self.RELIEF_CALC_AVAILABLE = 'relief_calc_available'
        self.RELIEF_CALC_RATE_KG_S = 'relief_calc_rate_lb_hr'
        self.CATASTROPHIC_VESSEL_FAILURE = 'catastrophic_vessel_failure'
        self.STORAGE_VOLUME_M3 = 'storage_volume_m3'
        self.STORAGE_MASS_KG = 'storage_mass_kg'
        self.RELEASE_ANGLE_DEGREES_FROM_HORIZONTAL = 'release_angle_degrees_from_horizontal'
        self.PERSONNEL_PER_10K_M2 = 'personnel_per_10k_m2'
        self.PERSONNEL_PER_10K_M2_NIGHT = 'personnel_per_10k_m2_night'
        self.RELEASE_INDOORS = 'release_indoors'
        self.AVAILABLE_POOL_AREA_M2 = 'available_pool_area_m2'
        self.RELEASE_LATITUDE_DEGREE_DECIMAL = 'release_latitude_degree_decimal'
        self.RELEASE_LONGITUDE_DEGREE_DECIMAL = 'release_longitude_degree_decimal'
        self.PERSONNEL_WORKING_AT_ELEVATION_NEAR_RELEASE_POINT = 'personnel_working_at_elevation_near_release_point'
        self.LOG_HANDLER = 'log_handler'
        self.KML_HANDLER = 'kml_handler'
        self.DISPLAY_DIAGNOSTIC_INFORMATION = 'display_diagnostic_information'
        self.USE_ONE_MET_CONDITION_WORST_CASE = 'use_one_met_condition_worst_case'
        self.ROOM_VOL_M3 = 'room_vol_m3'
        self.PRODUCTION_AREA_ACH = 'production_area_ach'
        self.ROOM_VENT_DIAMETER_IN = 'room_vent_diameter_in'
        self.ROOM_VENT_ELEVATION_M = 'room_vent_elevation_m'
        self.ROOM_VENT_RELEASE_ANGLE_DEGREES_FROM_HORIZONTAL = 'room_vent_release_angle_degrees_from_horizontal'
        self.EXIT_VAPOR_MASS_FRACTION = 'exit_vapor_mass_fraction'
        self.QUICK_TEST = 'quick_test'
        self.FULL_TEST = 'full_test'
        self.USE_MULTICOMPONENT_METHOD = 'use_multicomponent_method'
        self.USE_DOSE_AND_PROBIT = 'use_dose_and_probit'
        self.FLUID_SPEC = 'fluid_spec'
        self.DNV_BENCHMARK = 'dnv_benchmark'
        self.DESCRIPTION = 'description'
        self.DNV_BENCHMARK_SPECS = 'dnv_benchmark_specs'
        self.JET_FIRE_ANALYSIS = 'jet_fire_analysis'
        self.JET_FIRE_ONLY = 'jet_fire_only'
        self.MAX_RELEASE_DURATION_MIN = 'max_release_duration_min'
        self.GET_PHAST_DISCHARGE_ONLY = 'get_phast_discharge_only'

    def set_identification_string(self):
        # this id_string will be included in the interim results to help users see which run is in progress.
        # the units look as though they don't match up.  
        # however, this is using the property names associated with the input dict keys. 
        # the units haven't been converted from input units to internal units, so they will still have the units 
        # as displayed in the id_string.
        id_string = ''
        id_string += f'chems: {self.inputs[self.CHEM_MIX]}\n'
        id_string += f'press: {self.inputs[self.PRESS_PA]} psig\n'
        id_string += f'temp: {self.inputs[self.TEMP_K]} deg C\n'

        self.id_string = id_string

    def clean_inputs(self):
        
        relief_phase_key = self.RELIEF_PHASE
        chem_mix_key = self.CHEM_MIX
        comp_key = self.COMPOSITION
        bldg_dists_key = self.BLDG_DISTS
        bldg_nums_key = self.BLDG_NUMS
        bldg_hts_key = self.BLDG_HTS
        bldg_ach_key = self.BLDG_ACH
        bldg_hvac_key = self.BLDG_HVAC_SHUTOFF
        bldg_night_occ_key = self.BLDG_NIGHT_OCC
        mole_flag_key = self.COMPOSITION_IS_IN_MOLES
        inhal_flag_key = self.INHALATION
        ff_flag_key = self.FLASH_FIRE
        relief_calc_avail_flag_key = self.RELIEF_CALC_AVAILABLE
        catastrophic_vessel_failure_flag_key = self.CATASTROPHIC_VESSEL_FAILURE
        release_indoors_flag_key = self.RELEASE_INDOORS
        personnel_working_at_elevation_near_release_point_flag_key = self.PERSONNEL_WORKING_AT_ELEVATION_NEAR_RELEASE_POINT
        display_diagnostic_information_flag_key = self.DISPLAY_DIAGNOSTIC_INFORMATION
        use_one_met_condition_worst_case_flag_key = self.USE_ONE_MET_CONDITION_WORST_CASE
        quick_test_key = self.QUICK_TEST
        full_test_key = self.FULL_TEST
        use_multicomponent_method_key = self.USE_MULTICOMPONENT_METHOD
        use_dose_and_probit_key = self.USE_DOSE_AND_PROBIT
        vapor_cloud_explosion_key = self.VAPOR_CLOUD_EXPLOSION
        fluid_spec_key = self.FLUID_SPEC
        dnv_benchmark_key = self.DNV_BENCHMARK
        description_key = self.DESCRIPTION
        dnv_benchmark_specs_key = self.DNV_BENCHMARK_SPECS
        jet_fire_analysis_key = self.JET_FIRE_ANALYSIS
        jet_fire_only_key = self.JET_FIRE_ONLY
        pv_burst_flag_key = self.PV_BURST
        max_release_duration_min_key = self.MAX_RELEASE_DURATION_MIN
        get_phast_discharge_only_key = self.GET_PHAST_DISCHARGE_ONLY

        # since building data can be null, it is being excluded from this method of parsing different potential inputs.  
        # This now means that distances and heights must be received as numeric types.

        lists = [chem_mix_key, bldg_nums_key, bldg_dists_key, comp_key, bldg_hts_key, bldg_ach_key, bldg_hvac_key, bldg_night_occ_key]
        
        non_float_lists = [chem_mix_key, bldg_nums_key, bldg_night_occ_key]

        for key in lists:
            if isinstance(self.inputs[key], str):
                if '|' in self.inputs[key]:
                    self.inputs[key] = self.inputs[key].split("|")
                if self.inputs[key][0] == '[' and self.inputs[key][-1] == ']':
                    self.inputs[key] = self.inputs[key][1:-1]
                    self.inputs[key] = self.inputs[key].split(";")
            
            if not isinstance(self.inputs[key], list):
                self.inputs[key] = [self.inputs[key]]
        
        if len(self.inputs[chem_mix_key]) != len(self.inputs[comp_key]):
            raise Exception(Exception_Enum.MISMATCH_CHEMICAL_CONCENTRATION_COUNT)

        self.set_nulls_to_default()

        for key in lists:
            if key not in non_float_lists:
                for i in range(len(self.inputs[key])):
                    val = self.inputs[key][i]
                    if not pd.isna(val):
                        self.inputs[key][i] = float(val)

        flags = [mole_flag_key, inhal_flag_key, ff_flag_key, relief_calc_avail_flag_key, catastrophic_vessel_failure_flag_key, release_indoors_flag_key, personnel_working_at_elevation_near_release_point_flag_key, display_diagnostic_information_flag_key, use_one_met_condition_worst_case_flag_key, quick_test_key, full_test_key, use_multicomponent_method_key, use_dose_and_probit_key, vapor_cloud_explosion_key, dnv_benchmark_key, jet_fire_analysis_key, jet_fire_only_key, pv_burst_flag_key, get_phast_discharge_only_key]

        for key in flags:
            self.inputs[key] = self.parse_bool_from_input_flag(self.inputs[key])

        # components, building dists, building hts, bldg ach and bldg nighttime occupancy are lists of floats, handled separately.
        # since they are lists (and not floats), they are added to the non-float keys.
        nonfloat_keys = flags + lists +  [relief_phase_key, self.LOG_HANDLER, self.KML_HANDLER, fluid_spec_key, description_key, dnv_benchmark_specs_key]
        

        for key in self.inputs.keys():
            if key in nonfloat_keys:
                continue
            if self.inputs[key] is None:
                continue
            self.inputs[key] = float(self.inputs[key])

        self.inputs[self.RELEASE_ELEVATION_M] = max(0,self.inputs[self.RELEASE_ELEVATION_M])

        if self.inputs[self.VAPOR_CLOUD_EXPLOSION] or self.inputs[self.PV_BURST]:

            # explosions during the day are assumed to be much worse than at night
            self.inputs[self.USE_ONE_MET_CONDITION_WORST_CASE] = False
            self.inputs[self.FLASH_FIRE] = True
            
            # when vce is triggered, isolate to only the flammable release.  also, only the full hour of release time will be modeled (not controlled here)
            self.inputs[self.INHALATION] = False

            if self.inputs[self.PERSONNEL_PER_10K_M2_NIGHT] is None:
                self.inputs[self.PERSONNEL_PER_10K_M2_NIGHT] = self.inputs[self.PERSONNEL_PER_10K_M2]

            # the indoor model is strictly for handling exfiltration.  vce hazards inside are to be handled as outdoor models.  
            # congested volumes should be increased by one blast strength class for indoor releases, but no other conditions.
            self.inputs[self.RELEASE_INDOORS] = False

            # PV Burst model requires a catastrophic vessel rupture.  it will use the PWS isentropic pathway to availabe energy for the blast
            if self.inputs[self.PV_BURST]:
                self.inputs[self.CATASTROPHIC_VESSEL_FAILURE] = True

        if self.inputs[self.USE_ONE_MET_CONDITION_WORST_CASE]:
            self.inputs[self.BLDG_NIGHT_OCC] = 3 * [cd.OCC_UNOCCUPIED]

        self.VALID_HAZARDS = {
            cd.HAZARD_TYPE_FLASH_FIRE: self.inputs[self.FLASH_FIRE],
            cd.HAZARD_TYPE_INHALATION: self.inputs[self.INHALATION],
            cd.HAZARD_TYPE_VCE: self.inputs[self.VAPOR_CLOUD_EXPLOSION],
            cd.HAZARD_TYPE_PV_BURST: self.inputs[self.PV_BURST],
        }


    def is_iterable(self, val):
        ans = False
        try:
            iterator = iter(val)
            ans = True
        except:
            pass
        else:
            pass

        return ans

    def set_nulls_to_default(self):
        for key in self.inputs.keys():
            self.inputs[key] = self.set_vals_to_null(self.inputs[key])

        if pd.isna(self.inputs[self.RELIEF_CALC_RATE_KG_S]):
            rel = self.inputs[self.RELIEF_CALC_AVAILABLE]
            rel = self.parse_bool_from_input_flag(rel)
            if rel:
                raise Exception(Exception_Enum.RELIEF_RATE_REQUIRED_IF_RELIEF_CALC_AVAIALABLE)
            self.inputs[self.RELIEF_CALC_RATE_KG_S] = -1

        for i in range(len(self.inputs[self.BLDG_ACH])):
            val = self.inputs[self.BLDG_ACH][i]
            if pd.isna(val):
                self.inputs[self.BLDG_ACH][i] = 10
        
        for i in range(len(self.inputs[self.BLDG_HVAC_SHUTOFF])):
            val = self.inputs[self.BLDG_HVAC_SHUTOFF][i]
            if pd.isna(val):
                self.inputs[self.BLDG_HVAC_SHUTOFF][i] = 60

        if pd.isna(self.inputs[self.RELEASE_LATITUDE_DEGREE_DECIMAL]):
            self.inputs[self.RELEASE_LATITUDE_DEGREE_DECIMAL] = -1e6
        
        if pd.isna(self.inputs[self.RELEASE_LONGITUDE_DEGREE_DECIMAL]):
            self.inputs[self.RELEASE_LONGITUDE_DEGREE_DECIMAL] = -1e6

        for i in range(len(self.inputs[self.BLDG_DISTS])):
            val = self.inputs[self.BLDG_DISTS][i]
            if pd.isna(val):
                self.set_bldg_data_to_null(idx=i)
        
        if self.inputs[self.MAX_HOLE_SZ_M] is None:
            cat = self.inputs[self.CATASTROPHIC_VESSEL_FAILURE]
            cat = self.parse_bool_from_input_flag(cat)
            if not cat:
                raise Exception(Exception_Enum.HOLE_SIZE_REQUIRED_IF_FAILURE_NOT_CATASTROPHIC)
            self.inputs[self.MAX_HOLE_SZ_M] = -1

        if not isinstance(self.inputs[self.AVAILABLE_POOL_AREA_M2], str):
            if self.inputs[self.AVAILABLE_POOL_AREA_M2] is None:
                self.inputs[self.AVAILABLE_POOL_AREA_M2] = Consts.POOL_AREA_DEFAULT_M2
            else:
                if self.inputs[self.AVAILABLE_POOL_AREA_M2] < 0:
                    self.inputs[self.AVAILABLE_POOL_AREA_M2] = Consts.POOL_AREA_DEFAULT_M2

        if self.inputs[self.STORAGE_VOLUME_M3] is None and self.inputs[self.STORAGE_MASS_KG] is None:
            raise Exception(Exception_Enum.STORAGE_BASIS_NOT_DEFINED)
        

    def parse_bool_from_input_flag(self, inp_val):
        if isinstance(inp_val, bool):
            return inp_val
        if pd.isna(inp_val):
            return False
        inp_val = str(inp_val)
        c = inp_val[0].lower()
        return (c == 't') or (c == '1') or (c == 'y')

    def clean_attributes(self):
        #the keys of the input dictionary were stored as 
        #attributes of the "Model_Inputs" object.

        #This method converts the references to keys to instead refer
        #their values.

        #that is, instead of being model_input.inputs['key'], the input will be "model_input.key"

        keys = [
            'CHEM_MIX',
            'MAX_HOLE_SZ_M',
            'PRESS_PA',
            'TEMP_K',
            'COMPOSITION',
            'COMPOSITION_IS_IN_MOLES',
            'RELEASE_ELEVATION_M',
            'INHALATION',
            'FLASH_FIRE',
            'BLDG_NUMS',
            'BLDG_DISTS',
            'BLDG_HTS',
            'BLDG_ACH',
            'BLDG_HVAC_SHUTOFF',
            'OFFSITE_DIST',
            'RELIEF_PHASE',
            'RELIEF_CALC_AVAILABLE',
            'RELIEF_CALC_RATE_KG_S',
            'CATASTROPHIC_VESSEL_FAILURE',
            'STORAGE_VOLUME_M3',
            'STORAGE_MASS_KG',
            'RELEASE_ANGLE_DEGREES_FROM_HORIZONTAL',
            'PERSONNEL_PER_10K_M2',
            'PERSONNEL_PER_10K_M2_NIGHT',
            'RELEASE_INDOORS',
            'AVAILABLE_POOL_AREA_M2',
            'RELEASE_LATITUDE_DEGREE_DECIMAL',
            'RELEASE_LONGITUDE_DEGREE_DECIMAL',
            'PERSONNEL_WORKING_AT_ELEVATION_NEAR_RELEASE_POINT',
            'LOG_HANDLER',
            'KML_HANDLER',
            'DISPLAY_DIAGNOSTIC_INFORMATION',
            'USE_ONE_MET_CONDITION_WORST_CASE',
            'BLDG_NIGHT_OCC',
            'ROOM_VOL_M3',
            'PRODUCTION_AREA_ACH',
            'ROOM_VENT_DIAMETER_IN',
            'ROOM_VENT_ELEVATION_M',
            'ROOM_VENT_RELEASE_ANGLE_DEGREES_FROM_HORIZONTAL',
            'EXIT_VAPOR_MASS_FRACTION',
            'QUICK_TEST',
            'FULL_TEST',
            'USE_MULTICOMPONENT_METHOD',
            'USE_DOSE_AND_PROBIT',
            'VAPOR_CLOUD_EXPLOSION',
            'PV_BURST',
            'FLUID_SPEC',
            'DNV_BENCHMARK',
            'DESCRIPTION',
            'DNV_BENCHMARK_SPECS',
            'JET_FIRE_ANALYSIS',
            'JET_FIRE_ONLY',
            'MAX_RELEASE_DURATION_MIN',
            'GET_PHAST_DISCHARGE_ONLY',
        ]

        for key in keys:
            input_key = getattr(self, key)
            if input_key in self.inputs.keys():
                setattr(self, key, self.inputs[input_key])

        keys = [
            'MOL_FRACTIONS',
            'MASS_FRACTIONS'
        ]

        for key in keys:
            delattr(self, key)
        
        # delattr(self, 'inputs')

    def bldg_consequence_targs_dict(self):
        self.CONSEQUENCE_DISTS = {}
        self.CONSEQUENCE_DISTS[cd.KEYS_DIST_TO_OFFSITE_M] = self.OFFSITE_DIST

        # ONSITE "catastrophic" impact is 10xERPG-3 to 2000 ft onsite
        # OFFSITE "catastrophic" impact is 10xERPG-3 reaching 740 ft into community
        # the maximum distance considered by the application should be the greater of 
        # these two values.
        max_dist_onsite = 2000 / 3.28084
        max_dist_offsite = self.CONSEQUENCE_DISTS[cd.KEYS_DIST_TO_OFFSITE_M] + 740 / 3.28084
        max_dist = max(max_dist_onsite, max_dist_offsite)
        min_dist = min(max_dist_onsite, max_dist_offsite)
        
        
        max_height_m = 0
        for i in range(len(self.BLDG_DISTS)):
            val = self.BLDG_DISTS[i]
            if pd.isna(val):
                continue
            max_dist = max(val, max_dist)
            min_dist = min(val, min_dist)
            max_height_m = max(val, max_height_m)
            
        # the maximum height to evaluate is the height of the tallest roof plus a range to confirm higher val.
        # as of 5/5/22, this elev is 6 m.  Stored in Consts.ELEVATION_RANGE_FOR_CONC_EVAL_M.
        # the additional 6 m tests for peak concentrations in the vertical vicinity of a roof.
        # to account for localized wind effects, elevations within 6 m (+/-) of a target height
        # are evaluated.  The worst case concentration is used.

        max_height_m += Consts.ELEVATION_RANGE_FOR_CONC_EVAL_M

        self.bldgs = []

        for i in range(len(self.BLDG_NUMS)):
            
            # NOTE - bldg object will take hvac_shutoff time as an input in minutes.  
            # convert and store hvac_shutoff time in units of seconds.
            # the property is hvac_shutoff_sec.
            
            # for buildings that don't have a listed distance (dist set to None, etc.)
            # all data items for that building have been set to pd.NA.
            # this will check if the building number is pd.NA.  if so, then day occupancy
            # will also be set to pd.NA. 

            bldg_is_utilized = not pd.isna(self.BLDG_DISTS[i])
            d_o = cd.KEYS_BLDG_ALL_OCC[i]
            if not bldg_is_utilized:
                d_o = pd.NA

            n_o = self.BLDG_NIGHT_OCC[i]
            if pd.isna(n_o):
                n_o = cd.KEYS_BLDG_UNOCCUPIED
            if not isinstance(n_o, str):
                night_occ_str = cd.KEYS_BLDG_UNOCCUPIED
                if n_o >= cd.OCC_LOW and n_o < cd.OCC_MEDIUM:
                    night_occ_str = cd.KEYS_BLDG_LOW_OCC
                if n_o >= cd.OCC_MEDIUM and n_o < cd.OCC_HIGH:
                    night_occ_str = cd.KEYS_BLDG_MED_OCC
                if n_o >= cd.OCC_HIGH:
                    night_occ_str = cd.KEYS_BLDG_HIGH_OCC
                n_o = night_occ_str

            bldg = Building(
                occupancy= d_o,
                num = self.BLDG_NUMS[i],
                dist_m=self.BLDG_DISTS[i],
                ht_m=self.BLDG_HTS[i],
                ach = self.BLDG_ACH[i],
                hvac_shutoff_min=self.BLDG_HVAC_SHUTOFF[i],
                nighttime_occupancy = n_o,
                use_one_met_condition_worst_case = self.USE_ONE_MET_CONDITION_WORST_CASE,
                use_dose_and_probit=self.USE_DOSE_AND_PROBIT)
            
            self.bldgs.append(copy.deepcopy(bldg))
        


        self.CONSEQUENCE_DISTS[cd.KEYS_MAX_DIST_M] = max_dist
        self.CONSEQUENCE_DISTS[cd.KEYS_MIN_DIST_M] = min_dist

        self.CONSEQUENCE_DISTS[cd.KEYS_MAX_HEIGHT_M] = max_height_m

        delattr(self,'BLDG_NUMS')
        delattr(self,'BLDG_DISTS')
        delattr(self,'BLDG_HTS')
        delattr(self,'BLDG_ACH')
        delattr(self,'BLDG_HVAC_SHUTOFF')
        delattr(self,'BLDG_NIGHT_OCC')

    def enumerate_relief_phase(self):
        
        self.RELIEF_PHASE = str(self.RELIEF_PHASE)
        self.RELIEF_PHASE = self.RELIEF_PHASE.lower()

        rel_phase = consts.RELIEF_PHASE.VAPOR
        if self.RELIEF_PHASE[0] == "t" or self.RELIEF_PHASE[0] == "2":
            rel_phase = consts.RELIEF_PHASE.TWO_PHASE
        if self.RELIEF_PHASE[:3] == "let":
            rel_phase = consts.RELIEF_PHASE.LET_MODEL_DECIDE
        if self.RELIEF_PHASE[:3] == "liq":
            rel_phase = consts.RELIEF_PHASE.LIQUID
        
        self.RELIEF_PHASE = rel_phase

    def convert_to_internal_units(self):
        
        # this looks strange.  however, the input units are different than the property names.
        # rather than attempt to change them while importing and casting data types, converting
        # units after all property designations and other data is cleaned up.

        self.MAX_HOLE_SZ_M = self.MAX_HOLE_SZ_M / 3.28084 / 12
        self.PRESS_PA = (self.PRESS_PA + 14.6959) * 101325 / 14.6959
        self.TEMP_K = self.TEMP_K + 273.15
        self.RELIEF_CALC_RATE_KG_S = self.RELIEF_CALC_RATE_KG_S / Consts.LB_PER_KG / 3600
        self.RELEASE_ANGLE_DEGREES_FROM_HORIZONTAL *= math.pi / 180

        
    def model_updates_for_catastrophic_failure(self):
        # if a catastrophic leak is being modeled, set a default diameter for the model.
        # this will be updated later by the model.

        if self.CATASTROPHIC_VESSEL_FAILURE:
            self.MAX_HOLE_SZ_M = 20 / 12 / 3.28084
            self.RELIEF_CALC_AVAILABLE = False
            self.RELIEF_PHASE = consts.RELIEF_PHASE.LET_MODEL_DECIDE
            self.RELEASE_ANGLE_DEGREES_FROM_HORIZONTAL = 0

    def set_details_for_onsite_impact_evaluation(self):
        # the population per 10,000 m2 is provided as an input from the user
        # the impact area of the release to on-site workers will be calculated.
        # the required areas for serious, critical, & catastrophic categorization of
        # on-site impact are calculated here

        self.ONSITE_CATEGORY_IMPACT_AREA_M2 = {x: 10000 * cd.ONSITE_IMPACTED_PERSONNEL[x] / self.PERSONNEL_PER_10K_M2 for x in cd.CAT_ALL_CATEGORIES}
        self.ONSITE_CATEGORY_IMPACT_AREA_M2_NIGHT = {x: 0 for x in cd.CAT_ALL_CATEGORIES}

        if not self.USE_ONE_MET_CONDITION_WORST_CASE:
            self.ONSITE_CATEGORY_IMPACT_AREA_M2_NIGHT = {x: 10000 * cd.ONSITE_IMPACTED_PERSONNEL[x] / self.PERSONNEL_PER_10K_M2_NIGHT for x in cd.CAT_ALL_CATEGORIES}

        self.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS = copy.deepcopy(Consts.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS)

        if self.PERSONNEL_WORKING_AT_ELEVATION_NEAR_RELEASE_POINT:
            # if there are workers within the 6-meter (20 ft) vertical offset, this will model their potential impact

            e = self.RELEASE_ELEVATION_M
            t = Consts.ELEVATION_RANGE_FOR_CONC_EVAL_M
            elevs_around_release = [max(0, e-t), max(0, e-t/2), e, e+t/2, e + t]
            self.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS.extend(elevs_around_release)
        
        if self.RELEASE_ELEVATION_M < max(self.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS):
            if self.RELEASE_ELEVATION_M not in self.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS:
                self.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS.append(self.RELEASE_ELEVATION_M)

        self.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS = list(set(self.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS))

    def adjust_for_quick_test(self):
        # any tests that are specifically for benchmarking against the upstream model (DNV PYWPS), 
        # set it as a quick test in the target table model.  It will also be handled for analysis at the
        # completion of the target table run

        if self.DNV_BENCHMARK:
            self.QUICK_TEST = True
        
        # for models specified as being a quick test, the number of models and API calls will be greatly reduced
        if self.QUICK_TEST:
            # a maximum of one hazard should be considered.
            if self.INHALATION and self.FLASH_FIRE:
                self.FLASH_FIRE = False
            
            # since we have an option to only use one weather condition, we will activate that.  it only works for a nighttime condition.
            self.USE_ONE_MET_CONDITION_WORST_CASE = True
            
            # analyzing impact to buildings increases the number of elevations to analyze.  this removes the checked building properties.
            # debug xxx apple - testing timings for full set of elevations using dists and footprints calc
            bldg:Building
            for bldg in self.bldgs:
                bldg.dist_m = None
                bldg.ht_m = None

    def prep_chems(self):
        self.chems = Chems(mi=self)

        if not self.chems.inhal_credible:
            self.INHALATION = False

        if not self.chems.ff_credible:
            self.FLASH_FIRE = False

    def get_version_info(self):

        self.VERSION_PY_LOPA = f'{Debug_Info.current_version} from Debug Info'
        try:
            self.VERSION_PY_LOPA = version('py_lopa')
        except:
            pass
        
        self.VERSION_PY_PWS = 'N/A'
        try:
            self.VERSION_PY_PWS = version('pypws')
        except:
            pass
    
    def set_bldg_data_to_null(self, idx):

        self.inputs[self.BLDG_NUMS][idx] = pd.NA
        self.inputs[self.BLDG_DISTS][idx] = pd.NA
        self.inputs[self.BLDG_HTS][idx] = pd.NA
        self.inputs[self.BLDG_ACH][idx] = pd.NA
        self.inputs[self.BLDG_HVAC_SHUTOFF][idx] = pd.NA
        self.inputs[self.BLDG_NIGHT_OCC][idx] = pd.NA

    def set_vals_to_null(self, val):
        
        if isinstance(val, list):
            for i in range(len(val)):
                val[i] = self.set_vals_to_null(val[i])
        
        if val == '':
            val = pd.NA
        if isinstance(val, str):
            l = val.lower()
            if l == 'none' or l == 'null' or l == 'n/a':
                val = pd.NA
        
        return val
    
    def model_updates_for_flammables(self):
        if self.JET_FIRE_ONLY:
            self.JET_FIRE_ANALYSIS = True