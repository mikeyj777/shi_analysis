import copy
import json
import datetime
from datetime import datetime as dt

from pypws.enums import ResultCode

from py_lopa.model_work.model_controller import Model_Controller
from py_lopa.calcs import helpers
from py_lopa.calcs.consts import Consts

cd = Consts().CONSEQUENCE_DATA

class Model_Interface:
    
    def __init__(self):
        self.inputs = {}
        self.interim_data = []
        self.t0 = dt.now(datetime.UTC)
        self.time_stamp_utc = helpers.time_stamp_format(self.t0)
        self.data = None
        self.debug_data = None
        self.dispersion_plots_json = None
        self.dispersion_plots_pillow = None
        self.images_byte_strings = None
        self.mi = None
        self.chems = None
        self.dispersion_dicts = None
        self.timings = None

        # maintaining the below as properties to aid in header analysis
        self.vce_data = []
        self.phast_disp = {}
        self.phast_discharge = {}
        self.material = None
        self.flashresult = None


    def set_inputs_from_csv(self, infile):

        self.inputs = helpers.get_model_inputs_from_csv(infile)

    def set_inputs_from_web_form(self):
        pass

    def set_inputs_from_json(self, path_to_json_file = None, json_data = None):
        if path_to_json_file is not None:
            with open(path_to_json_file) as f:
                d = json.load(f)
        else:
            d = json.loads(json_data)
        
        main_key = 'PrimaryInputs'
        chems_key = 'ChemicalComponents'
        bldgs_key = 'BuildingInfo'
        if main_key not in d:
            main_key = 'AssesmentDetails'
            chems_key = 'InputChemicals'
            bldgs_key = 'InputNearbyBuildings'

        prim = d[main_key]
        use_multicomponent_method = True
        if 'UseMulticomponentMethod' in prim:
            use_multicomponent_method = prim['UseMulticomponentMethod']
        
        vapor_cloud_explosion = False
        if 'VceModel' in prim:
            vapor_cloud_explosion = prim['VceModel']


        release_latitude_degree_decimal = prim['ApproxLatitude']
        release_longitude_degree_decimal = prim['ApproxLongitude']
        inhalation = prim['StudyForInhalation']
        flash_fire = prim['StudyForFire']
        relief_calc_available = prim['ReliefCalculation']
        comp_is_moles = prim['MolarComposition']

        storage_press_psig = prim['StoragePressure']
        storage_temp_deg_c = prim['StorageTemperature']

        exit_press_psig = prim['ExitPressure']
        exit_temp_deg_c = prim['ExitTemperature']

        pressure_psig = storage_press_psig
        temp_deg_c = storage_temp_deg_c

        if relief_calc_available:
            pressure_psig = exit_press_psig
            temp_deg_c = exit_temp_deg_c

        relief_calc_rate_lb_hr = prim['ReliefRate']
        catastrophic_vessel_failure = prim['CatastrophicRelease']
        max_hole_size_in = prim['DischargeDiameter']
        release_elevation_m = prim['ReleaseElevation']

        rel_angle_input = prim['ReleaseOrientation']

        release_angle_degrees_from_horizontal = float(rel_angle_input) * 15

        rel_phase_input = prim['RelievingPhase']

        rel_phase_translation = {
            0:cd.PHASE_LET_MODEL_DECIDE,
            1:cd.PHASE_VAPOR,
            2:cd.PHASE_LIQUID,
            3:cd.PHASE_TWO_PHASE
        }


        relief_phase = rel_phase_translation[rel_phase_input]
        storage_volume_m3 = prim['MaximumVolume']
        storage_mass_kg = prim['MaximumMass']
        basis = prim['Basis']

        if basis == 'Mass':
            storage_volume_m3 = None
        else:
            storage_mass_kg = None
            
        available_pool_area_m2 = prim['AvailablePoolArea']
        release_indoors = prim['IndoorRelease']
        offsite_dist_m = prim['DistNonProcPlantPop']
        personnel_per_10k_m2 = prim['ProcAreaPopDensity']
        

        # display_diagnostic_information = prim['IncludeDiagnostics']
        # most times, on the desktop version, will want the diagnostic display.
        display_diagnostic_information = True
        personnel_working_at_elevation_near_release_point = prim['PersonnelAtReleaseElevation']
        use_one_met_condition_worst_case_input = prim['WeatherCondition']

        # for cases where VCE is credible, study the elevations around the release point to determine flammable clouds near the source.
        if vapor_cloud_explosion:
            personnel_working_at_elevation_near_release_point = True

        use_one_met_condit_translation = {
            'night': True,
            'both': False
        }

        use_one_met_condition_worst_case = use_one_met_condit_translation[use_one_met_condition_worst_case_input]

        personnel_per_10k_m2_night = prim['NightProcAreaPopDensity']

        chemical_mix = []
        composition = []

        chems = d[chems_key]

        for chem in chems:
            chemical_mix.append(chem['CASNumber'])
            composition.append(chem['Contribution'])

        bldg_nums_low_med_high = [0,0,0]
        bldg_dists_low_med_high = [0,0,0]
        bldg_hts_low_med_high = [0,0,0]
        bldg_ach_low_med_high = [0,0,0]
        bldg_nighttime_occupancy = [0,0,0]
        bldg_hvac_shutoff_min_low_med_high = [0,0,0]
        
        bldg_data = d[bldgs_key]
        
        night_occ_translation = {
            0: cd.KEYS_BLDG_LOW_OCC,
            1: cd.KEYS_BLDG_MED_OCC,
            2: cd.KEYS_BLDG_HIGH_OCC,
            3: cd.KEYS_BLDG_UNOCCUPIED,
            None: None
        }

        for bldg in bldg_data:
            day_occ = bldg['OccupancyLevel']
            bldg_nums_low_med_high[day_occ] = bldg['BuildingNumber']
            bldg_dists_low_med_high[day_occ] = bldg['DistanceFromRelease']
            bldg_hts_low_med_high[day_occ] = bldg['HVACHeight']
            bldg_ach_low_med_high[day_occ] = bldg['AirChanges']
            bldg_hvac_shutoff_min_low_med_high[day_occ] = bldg['AirHandlerShutOff']
            
            if not use_one_met_condition_worst_case:
                night_occ = bldg['NighttimeOccupancyLevel']
                bldg_nighttime_occupancy[day_occ] = night_occ_translation[night_occ]

        room_vol_m3 = None
        production_area_ach = None
        room_vent_diameter_in = None
        room_vent_elevation_m = None
        room_vent_release_angle_degrees_from_horizontal = None

        if release_indoors:
            if 'ProductionRoomVol' in prim:
                room_vol_m3 = prim['ProductionRoomVol']

            
            if 'ProductionRoomACH' in prim:
                production_area_ach = prim['ProductionRoomACH']
            
            
            if 'ProductionRoomVentDiameter' in prim:
                room_vent_diameter_in = prim['ProductionRoomVentDiameter']

            
            if 'ProductionRoomVentElevation' in prim:
                room_vent_elevation_m = prim['ProductionRoomVentElevation']
            
            
            if 'ProductionRoomVentOrientation' in prim:
                room_vent_angle_input = prim['ProductionRoomVentOrientation']
                room_vent_release_angle_degrees_from_horizontal = float(room_vent_angle_input) * 15
        
        exit_vapor_mass_fraction = None
        if 'ExitVaporMassFraction' in prim:
            exit_vapor_mass_fraction = prim['ExitVaporMassFraction']

        self.set_inputs_as_arguments(
            chemical_mix=chemical_mix,
            composition=composition,
            comp_is_moles=comp_is_moles,
            max_hole_size_in=max_hole_size_in,
            pressure_psig=pressure_psig,
            temp_deg_c=temp_deg_c,
            release_elevation_m=release_elevation_m,
            flash_fire=flash_fire,
            inhalation=inhalation,
            use_one_met_condition_worst_case=use_one_met_condition_worst_case,
            bldg_nums_low_med_high=bldg_nums_low_med_high,
            bldg_dists_low_med_high=bldg_dists_low_med_high,
            bldg_hts_low_med_high=bldg_hts_low_med_high,
            bldg_ach_low_med_high=bldg_ach_low_med_high,
            bldg_nighttime_occupancy=bldg_nighttime_occupancy,
            bldg_hvac_shutoff_min_low_med_high=bldg_hvac_shutoff_min_low_med_high,
            offsite_dist_m=offsite_dist_m,
            relief_phase=relief_phase,
            relief_calc_available=relief_calc_available,
            relief_calc_rate_lb_hr=relief_calc_rate_lb_hr,
            catastrophic_vessel_failure=catastrophic_vessel_failure,
            storage_volume_m3=storage_volume_m3,
            storage_mass_kg=storage_mass_kg,
            release_angle_degrees_from_horizontal=release_angle_degrees_from_horizontal,
            personnel_per_10k_m2=personnel_per_10k_m2,
            personnel_per_10k_m2_night=personnel_per_10k_m2_night,
            release_indoors=release_indoors,
            available_pool_area_m2=available_pool_area_m2,
            release_latitude_degree_decimal=release_latitude_degree_decimal,
            release_longitude_degree_decimal=release_longitude_degree_decimal,
            personnel_working_at_elevation_near_release_point=personnel_working_at_elevation_near_release_point,
            display_diagnostic_information=display_diagnostic_information,
            room_vol_m3=room_vol_m3,
            production_area_ach=production_area_ach,
            room_vent_diameter_in=room_vent_diameter_in,
            room_vent_elevation_m=room_vent_elevation_m,
            room_vent_release_angle_degrees_from_horizontal=room_vent_release_angle_degrees_from_horizontal,
            exit_vapor_mass_fraction = exit_vapor_mass_fraction,
            use_multicomponent_method = use_multicomponent_method,
            vapor_cloud_explosion=vapor_cloud_explosion,
        )

    def set_inputs_from_dict(self, input_dict):
        self.inputs = input_dict

    def set_inputs_as_arguments(self, 
        chemical_mix = ['ammonia'],
        composition = [1],
        comp_is_moles = True,
        max_hole_size_in = -1,
        pressure_psig = 50,
        temp_deg_c = -33,
        release_elevation_m = 10,
        flash_fire = True,
        inhalation = True,
        use_one_met_condition_worst_case=True,
        bldg_nums_low_med_high = ['18', '215', '54D'],
        bldg_dists_low_med_high = [10, 20, 30],
        bldg_hts_low_med_high = [0, 0, 0],
        bldg_ach_low_med_high = [10,10, 10],
        bldg_nighttime_occupancy = [cd.KEYS_BLDG_UNOCCUPIED, cd.KEYS_BLDG_UNOCCUPIED, cd.KEYS_BLDG_UNOCCUPIED], # xxx debug apple - for this analysis, we want to minimize processing time.  restricting building heights will ensure that a minimal number of elevations are evaluated.
        bldg_hvac_shutoff_min_low_med_high = [60, 60, 60],
        offsite_dist_m = 200,
        relief_phase = 'let model decide',
        relief_calc_available = False,
        relief_calc_rate_lb_hr = -1,
        catastrophic_vessel_failure = False,
        storage_volume_m3 = None,
        storage_mass_kg = None,
        release_angle_degrees_from_horizontal = 0,
        personnel_per_10k_m2 = 3,
        personnel_per_10k_m2_night=0,
        release_indoors = False,
        available_pool_area_m2 = Consts.POOL_AREA_DEFAULT_M2,
        release_latitude_degree_decimal = 36.522605,
        release_longitude_degree_decimal = -82.538582,
        personnel_working_at_elevation_near_release_point=False,
        log_handler=print,
        kml_handler=print,
        display_diagnostic_information=True,
        room_vol_m3 = 5000,
        production_area_ach = 3,
        room_vent_diameter_in = 50,
        room_vent_elevation_m = 10,
        room_vent_release_angle_degrees_from_horizontal = 0,
        exit_vapor_mass_fraction = None,
        full_test = False,
        quick_test = False,
        use_multicomponent_method = True,
        use_dose_and_probit = False,
        vapor_cloud_explosion = False,
        dnv_benchmark = False,
        description = None,
        dnv_benchmark_specs = None,
        jet_fire_analysis = False,
        jet_fire_only = False,
        pv_burst = False,
        max_release_duration_min = None,
        get_phast_discharge_only = False,
    ):

        self.inputs['chemical_mix'] = chemical_mix
        self.inputs['composition'] = composition
        self.inputs['comp_is_moles'] = comp_is_moles
        self.inputs['max_hole_size_in'] = max_hole_size_in
        self.inputs['pressure_psig'] = pressure_psig
        self.inputs['temp_deg_c'] = temp_deg_c
        self.inputs['release_elevation_m'] = release_elevation_m
        self.inputs['inhalation'] = inhalation
        self.inputs['flash_fire'] = flash_fire
        self.inputs['use_one_met_condition_worst_case'] = use_one_met_condition_worst_case
        self.inputs['bldg_nums_low_med_high'] = bldg_nums_low_med_high
        self.inputs['bldg_dists_low_med_high'] = bldg_dists_low_med_high
        self.inputs['bldg_hts_low_med_high'] = bldg_hts_low_med_high
        self.inputs['bldg_ach_low_med_high'] = bldg_ach_low_med_high
        self.inputs['bldg_nighttime_occupancy'] = bldg_nighttime_occupancy
        self.inputs['bldg_hvac_shutoff_min_low_med_high'] = bldg_hvac_shutoff_min_low_med_high
        self.inputs['offsite_dist_m'] = offsite_dist_m
        self.inputs['relief_phase'] = relief_phase
        self.inputs['relief_calc_available'] = relief_calc_available
        self.inputs['relief_calc_rate_lb_hr'] = relief_calc_rate_lb_hr
        self.inputs['catastrophic_vessel_failure'] = catastrophic_vessel_failure
        self.inputs['storage_volume_m3'] = storage_volume_m3
        self.inputs['storage_mass_kg'] = storage_mass_kg
        self.inputs['release_angle_degrees_from_horizontal'] = release_angle_degrees_from_horizontal
        self.inputs['personnel_per_10k_m2'] = personnel_per_10k_m2
        self.inputs['personnel_per_10k_m2_night'] = personnel_per_10k_m2_night
        self.inputs['release_indoors'] = release_indoors
        self.inputs['available_pool_area_m2'] = available_pool_area_m2
        self.inputs['release_latitude_degree_decimal'] = release_latitude_degree_decimal
        self.inputs['release_longitude_degree_decimal'] = release_longitude_degree_decimal
        self.inputs['personnel_working_at_elevation_near_release_point'] = personnel_working_at_elevation_near_release_point
        self.inputs['log_handler'] = log_handler
        self.inputs['kml_handler'] = kml_handler
        self.inputs['display_diagnostic_information'] = display_diagnostic_information
        self.inputs['room_vol_m3'] = room_vol_m3
        self.inputs['production_area_ach'] = production_area_ach
        self.inputs['room_vent_diameter_in'] = room_vent_diameter_in
        self.inputs['room_vent_elevation_m'] = room_vent_elevation_m
        self.inputs['room_vent_release_angle_degrees_from_horizontal'] = room_vent_release_angle_degrees_from_horizontal
        self.inputs['exit_vapor_mass_fraction'] = exit_vapor_mass_fraction
        self.inputs['full_test'] = full_test
        self.inputs['quick_test'] = quick_test
        self.inputs['use_multicomponent_method'] = use_multicomponent_method
        self.inputs['use_dose_and_probit'] = use_dose_and_probit
        self.inputs['vapor_cloud_explosion'] = vapor_cloud_explosion
        self.inputs['dnv_benchmark'] = dnv_benchmark
        self.inputs['description'] = description
        self.inputs['dnv_benchmark_specs'] = dnv_benchmark_specs
        self.inputs['jet_fire_analysis'] = jet_fire_analysis
        self.inputs['jet_fire_only'] = jet_fire_only
        self.inputs['pv_burst'] = pv_burst
        self.inputs['max_release_duration_min'] = max_release_duration_min
        self.inputs['get_phast_discharge_only'] = get_phast_discharge_only

    def run(self):
        log_handler = self.inputs['log_handler']
        self.mc = Model_Controller(inputs=self.inputs)
        self.mc.phast_discharge = self.phast_discharge
        self.mc.phast_disp = self.phast_disp
        self.mc.material = self.material
        self.mc.flashresult = self.flashresult
        self.mc.vce_data = self.vce_data

        self.resultCode = self.mc.run()

        # storing these objects at a higher level even if result code doesn't return success.
        # others are handled differently depending on result code.

        self.phast_discharge = self.mc.phast_discharge
        self.phast_disp = self.mc.phast_disp
        self.material = self.mc.material
        self.flashresult = self.mc.flashresult
        self.vce_data = self.mc.vce_data

        if self.resultCode == ResultCode.SUCCESS:
            if 'get_phast_discharge_only' in self.inputs:
                if self.inputs['get_phast_discharge_only']:
                    return ResultCode.SUCCESS
            if 'vapor_cloud_explosion' in self.inputs:
                if self.inputs['vapor_cloud_explosion']:
                    return ResultCode.SUCCESS
            self.data = self.mc.data
            self.output_list = self.mc.output_list
            self.debug_data = self.mc.debug_data
            self.dispersion_plots_json = self.mc.dispersion_plots
            self.dispersion_dicts = self.mc.dispersion_dicts
            self.timings = self.mc.timings
            try:
                self.dispersion_plots_pillow = helpers.get_pil_images_from_json_list(self.dispersion_plots_json)
                self.images_byte_strings = helpers.get_list_of_image_byte_strings_from_pil_image_list(self.dispersion_plots_pillow)
            except:
                pass
            self.mi = self.mc.mi
            self.chems = self.mc.chems
            del_t = dt.now(datetime.UTC) - self.t0
            runtime = helpers.get_timedelta_in_h_m_s(del_t.seconds)            

            log_handler(f'time to run: ' + runtime + '.')
        elif self.resultCode == ResultCode.NO_DISCHARGE_RECORDS_ERROR:
            self.data = self.mc.data
            log_handler(f'no discharge or dispersion emissions from at least one modeled duration')
        else:
            self.data = self.mc.data
            log_handler(f'an unexpected error has occurred.  please consult admin.')
    
        return self.resultCode
        

    def get_data(self):
        return self.data
