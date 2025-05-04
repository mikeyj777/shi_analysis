import pickle
import math

import copy
import pandas as pd
import numpy as np

from py_lopa.calcs.consts import Consts, Wx_Enum
from py_lopa.calcs import helpers
from py_lopa.calcs.flattening import Flattening
from py_lopa.classes.building import Building


cd = Consts().CONSEQUENCE_DATA

class Conseq_Assess:
    
    def __init__(self, chems = None, mi = None):
        self.flatten = None
        self.dose_probit_data_for_kml_zxy = None
        self.chems = chems
        self.mi = mi
        self.use_one_met_condition_worst_case=True
        if self.mi is not None:
            self.use_one_met_condition_worst_case = self.mi.USE_ONE_MET_CONDITION_WORST_CASE
        res = {
            Wx_Enum.NIGHT: {},
            Wx_Enum.DAY: {}
        }

        if self.use_one_met_condition_worst_case:
            res = {cd.WX_WORST_CASE: {}}
        self.onsite_offsite_conseq_dict = {
            cd.CLASS_ONSITE: copy.deepcopy(res),
            cd.CLASS_OFFSITE: copy.deepcopy(res)
        }

    def calc_probability_of_exposure(self, prob, return_floats = False):

        targ_dict = cd.PROBABILITY_OF_EXPOSURE_MAXIMUM_CLOSED_ENDED
        if return_floats:
            targ_dict = cd.PROBABILITY_OF_EXPOSURE_MAXIMUM_CLOSED_ENDED_FLOAT_VALUES

        for max_in_cat in targ_dict.keys():
            if prob <= max_in_cat:
                return targ_dict[max_in_cat]


    def assessment(self, dispersion, mass_released_kg, wx_enum, hazard_type, duration_sec):
        self.consequence_list = []
        self.dose_probit_data_for_kml_zxy = None
        self.duration_sec = duration_sec
        self.mass_released_kg = mass_released_kg
        self.wx_enum = wx_enum
        self.dispersion = dispersion
        self.analysis_df = dispersion.analysis_df
        self.conc_pfls = dispersion.conc_profiles
        self.haz_cat_conc_footprints_df = dispersion.haz_cat_conc_footprints_df
        # self.mi = dispersion.mi
        self.impact_areas_m2_df = {}
        self.hazard_type = hazard_type

        self.prepare_results_dicts()

        if not dispersion.has_vapor:
            return

        self.get_onsite_result()
        # self.get_building_result()
        # self.get_offsite_result()
        # self.get_lfl_distance()
    
    def prepare_results_dicts(self):

        self.onsite_offsite_conseq_dict[cd.CLASS_ONSITE][self.wx_enum][self.hazard_type] = {
            cd.CAT_TITLE: pd.NA,
            cd.IMPACT_DISTANCE_M: pd.NA,
            cd.IMPACT_AREA_M2: pd.NA,
            cd.REQUIRED_AREA_M2: pd.NA,
            cd.PROBABILITY_OF_EXPOSURE: pd.NA
        }

        if self.hazard_type == cd.HAZARD_TYPE_FLASH_FIRE:
            self.onsite_offsite_conseq_dict[cd.CLASS_ONSITE][self.wx_enum][self.hazard_type][cd.DISTANCE_TO_LFL_M] = pd.NA

        self.onsite_offsite_conseq_dict[cd.CLASS_OFFSITE][self.wx_enum][self.hazard_type] = pd.NA

        bldg:Building
        for bldg in self.mi.bldgs:
            bldg.building_conseq_result_dict[self.wx_enum][self.hazard_type] = pd.NA

    def get_onsite_result(self):
        
        req_area_m2_dict = copy.deepcopy(self.mi.ONSITE_CATEGORY_IMPACT_AREA_M2)
        if not self.use_one_met_condition_worst_case and self.wx_enum == Wx_Enum.NIGHT:
            req_area_m2_dict = copy.deepcopy(self.mi.ONSITE_CATEGORY_IMPACT_AREA_M2_NIGHT)
        
        flatten = Flattening(haz_cat_conc_footprints_df=self.haz_cat_conc_footprints_df, analysis_df=self.analysis_df)
        self.flatten =  flatten
        self.impact_areas_m2_df = flatten.get_areas()
        

        df = self.impact_areas_m2_df
        ser_crit_cat = [cd.CAT_SERIOUS, cd.CAT_CRITICAL, cd.CAT_CATASTROPHIC]
        critical_and_catastrophic = [cd.CAT_CRITICAL, cd.CAT_CATASTROPHIC]
        min_mod_ser = [cd.CAT_MINOR, cd.CAT_MODERATE, cd.CAT_SERIOUS]

        for haz_type in cd.HAZARD_ALL_TYPES:
            curr_haz_df = df[df[cd.KEYS_TARG_AND_TYPE_FLAM_OR_INHAL] == haz_type]
            if len(curr_haz_df) == 0:
                continue
            conseq_set = False
            
            # cat score is a ranking of the consequence category.  catastrophic is 5.  minor is 1.  
            # This is used for scoring which result is worse from an overall risk perspective.  
            # This allows for consequence to figure in POE when reporting what is worse.
            # for example, a Moderate with POE 0 will be returned as the governing target over SER POE -2.
            # it 

            self.onsite_offsite_conseq_dict[cd.CLASS_ONSITE][self.wx_enum][haz_type] = {
                cd.CAT_TITLE: pd.NA,
                cd.IMPACT_DISTANCE_M: pd.NA,
                cd.IMPACT_AREA_M2: pd.NA,
                cd.REQUIRED_AREA_M2: pd.NA,
                cd.PROBABILITY_OF_EXPOSURE: pd.NA
            }

            conseq_list = []

            for i in range(len(cd.CAT_ALL_CATEGORIES)-1,-1,-1):
                cat_score = i + 1
                cat = cd.CAT_ALL_CATEGORIES[i]
                req_area_m2 = req_area_m2_dict[cat]
                impact_dist_m = self.get_max_dist_at_conc_onsite_wrapper(haz_type=haz_type, haz_class=cd.CLASS_ONSITE, conseq_cat=cat)
                curr_haz_cat_row = curr_haz_df[curr_haz_df[cd.CAT_TITLE] == cat]
                impact_area_m2 = helpers.get_data_from_pandas_series_element(curr_haz_cat_row[cd.AREA_M2])

                prob = impact_area_m2 / req_area_m2
                prob = helpers.float_value_rounded_to_x_sig_figs(prob, 2)
                log_poe_float = self.calc_probability_of_exposure(prob = prob, return_floats=True)
                cat_plus_poe = cat_score + log_poe_float
                poe_string = self.calc_probability_of_exposure(prob=prob)

                if impact_area_m2 > req_area_m2 and cat in critical_and_catastrophic:
                    conseq_list.append({
                        cd.PROBABILITY_OF_EXPOSURE: poe_string,
                        cd.CAT_SCORE: cat_score,
                        cd.CAT_TITLE: cat,
                        cd.IMPACT_DISTANCE_M: impact_dist_m,
                        cd.IMPACT_AREA_M2: impact_area_m2,
                        cd.REQUIRED_AREA_M2: req_area_m2,
                        cd.CAT_PLUS_POE: cat_plus_poe,
                        cd.LOG_POE: log_poe_float
                    })
                    conseq_set = True
                    

                if impact_area_m2 > 0 and cat in min_mod_ser:
                    conseq_list.append({
                        cd.PROBABILITY_OF_EXPOSURE: poe_string,
                        cd.CAT_SCORE: cat_score,
                        cd.CAT_TITLE: cat,
                        cd.IMPACT_DISTANCE_M: impact_dist_m,
                        cd.IMPACT_AREA_M2: impact_area_m2,
                        cd.REQUIRED_AREA_M2: req_area_m2,
                        cd.CAT_PLUS_POE: cat_plus_poe,
                        cd.LOG_POE: log_poe_float
                    })
                    conseq_set = True
                    
            if not conseq_set:
                return
            
            if self.mi.USE_DOSE_AND_PROBIT and haz_type == cd.HAZARD_TYPE_INHALATION:
                governing_targ_dict = self.determine_governing_consequence(conseq_list, check_for_rank_when_partial_order_not_included=False)
                
                conseq = governing_targ_dict[cd.CAT_TITLE]
                
                
                if conseq in ser_crit_cat:
                    dose_prob_conseq_list = self.get_onsite_result_dose_probit(req_area_m2_dict)
                    if dose_prob_conseq_list is not None:
                        # merging items from dose-response (new) list into the existing conc-based (old) list.  
                        # nothing on the merged list should have a rank (conseq + poe) greater than what is in the new list
                        dose_prob_conseq_df = pd.DataFrame(dose_prob_conseq_list)
                        max_cat_plus_poe_dose_resp = dose_prob_conseq_df[cd.CAT_PLUS_POE].max()
                        i = 0
                        conc_based_list_updated = False
                        while i < len(conseq_list):
                            row = conseq_list[i]
                            if row[cd.CAT_PLUS_POE] >= max_cat_plus_poe_dose_resp:
                                conseq_list.pop(i)
                                conc_based_list_updated = True
                                continue
                            i += 1
                        
                        if conc_based_list_updated:
                            conseq_list.extend(dose_prob_conseq_list)
                
            
            governing_targ_dict = self.determine_governing_consequence(conseq_list=conseq_list)

            self.onsite_offsite_conseq_dict[cd.CLASS_ONSITE][self.wx_enum][haz_type] = governing_targ_dict

            # if not conseq_set:
            #     # self.set_consequence_by_tier_designation()
            #     self.set_conseq_by_analysis_df(class_ = cd.CLASS_ONSITE, haz_type = haz_type, ignore_higher_level_impacts=False)

            self.consequence_list = conseq_list
    
    def determine_governing_consequence(self, conseq_list, check_for_rank_when_partial_order_not_included = True):

        df = pd.DataFrame(conseq_list)
        # calculate raw probability of exposure
        
        # partial orders can create scenarios where a potentially higher-risk target can be masked.
        # for example, serious with poe -1.5 will report lower risk than moderate with poe 0.  in case the partial order can't be met, the serious target should win.
        # the user will need to be notified of this potential
        df['log_poe_no_partials'] = np.ceil(df['log_poe'])

        # get the 'cat plus poe' score.  This is showing that, once credit is taken for POE, which result is the highest on a risk scale?
        df[cd.CAT_PLUS_POE] = df[cd.CAT_SCORE] + df['log_poe']

        # sort on 'cat plus poe', and then by the raw cat score.  This will get the highest risk value to the top, and will allow any 
        df = df.sort_values([cd.CAT_PLUS_POE, cd.CAT_SCORE], ascending=[False, False])

        worst_case_partial_order_included_series = df.iloc[0]
        rank_for_worst_case_when_partial_order_included = helpers.get_data_from_pandas_series_element(worst_case_partial_order_included_series[cd.CAT_SCORE])
        poe_text_for_output = helpers.get_data_from_pandas_series_element(worst_case_partial_order_included_series[cd.PROBABILITY_OF_EXPOSURE])

        df_sorted_without_partial_order = copy.deepcopy(df)
        df_sorted_without_partial_order['cat_plus_poe_no_partials'] = df_sorted_without_partial_order[cd.CAT_SCORE] + df_sorted_without_partial_order['log_poe_no_partials']
        df_sorted_without_partial_order = df_sorted_without_partial_order.sort_values(['cat_plus_poe_no_partials', cd.CAT_SCORE], ascending=[False, False])

        # check for governing when partial order not included
        if check_for_rank_when_partial_order_not_included:
                
            worst_case_when_partial_order_NOT_included = df_sorted_without_partial_order.iloc[0]
            rank_for_worst_case_when_partial_order_NOT_included = helpers.get_data_from_pandas_series_element(worst_case_when_partial_order_NOT_included[cd.CAT_SCORE])
            
            if rank_for_worst_case_when_partial_order_NOT_included > rank_for_worst_case_when_partial_order_included:
                category_with_partial_order = helpers.get_data_from_pandas_series_element(worst_case_partial_order_included_series[cd.CAT_TITLE])
                poe_text_for_output = helpers.get_data_from_pandas_series_element(worst_case_partial_order_included_series[cd.PROBABILITY_OF_EXPOSURE])
                poe_text_no_partial_order = helpers.get_data_from_pandas_series_element(worst_case_when_partial_order_NOT_included[cd.PROBABILITY_OF_EXPOSURE])
                category_no_partial_order = helpers.get_data_from_pandas_series_element(worst_case_when_partial_order_NOT_included[cd.CAT_TITLE])
                poe_text_for_output = f'{category_with_partial_order} with POE {poe_text_for_output} also "{category_no_partial_order} with POE {poe_text_no_partial_order}" higher target may apply.'
            
            worst_case_partial_order_included_series[cd.PROBABILITY_OF_EXPOSURE] = poe_text_for_output

        output = worst_case_partial_order_included_series[[cd.CAT_TITLE, cd.IMPACT_DISTANCE_M, cd.IMPACT_AREA_M2, cd.REQUIRED_AREA_M2, cd.PROBABILITY_OF_EXPOSURE, cd.CAT_PLUS_POE]]

        output = output.to_dict()

        return output

    def get_onsite_result_dose_probit(self, req_area_m2_dict):

        ser_crit_cat = [cd.CAT_SERIOUS, cd.CAT_CRITICAL, cd.CAT_CATASTROPHIC]
        check_slot = not pd.isna(self.chems.ta.slot_ppm_n_min)
        check_probit = not pd.isna(self.chems.ta.probit_a)
        if not check_slot and not check_probit:
            return

        # a value below zero here is used to indicate that the slot or probit may not have been analyzed.
        
        slot_exceeded_distance = -1
        prob_1_pct_distance = -1
        
        # if slot or probit will be analyzed the minimum distance is set to zero.  this indicates a case where the 
        # dose accumulated is under the slot or prob threshold, even when fleeing from the source of an event.

        if check_slot:
            slot = self.chems.ta.slot_ppm_n_min
            slot_exceeded_distance = 0

        if check_probit:
            prob_1_pct_distance = 0
        
        ans = self.chems.ta.on_site_dose_probit_analysis(self.conc_pfls)
        
        df:pd.DataFrame = ans['conc_and_tox_profile']

        row_slot_threshold = None
        row_prob_threshold = None

        for _, row in df.iterrows():
            dw_dist_m = helpers.get_data_from_pandas_series_element(row['dw_dist_m'])
            if check_slot:
                dose_for_slot = helpers.get_data_from_pandas_series_element(row['slot_dose_ppm_n_min'])
                if dose_for_slot >= slot:
                    slot_exceeded_distance = dw_dist_m
                    row_slot_threshold = row
            if check_probit:
                prob = helpers.get_data_from_pandas_series_element(row['probablity_of_severe_injury'])
                if prob >= Consts.THRESHOLD_PROBABILITY_SEVERE_INJURY:
                    prob_1_pct_distance = dw_dist_m
                    row_prob_threshold = row
            
        
        # consider a scenario where the "required area" is 31,416 m2 to impact 1 person.  
        # That indicates that one person is expected on average to occupy 31,416 m2.  This is a circle with radius 100 m.
        # Consider an individual fleeing from a release.  The individual is accumulating dose as he flees (this assumes he can't hold his breath)
        # if he is far enough away from the release, he will not accumulate sufficient dose to have a severe injury.
        # the zone inside which an individual will accumulate sufficient dose to exceed the SLOT dose is considered the impact area.

        # critical and catastrophic require impact area to be greater than required area.
        # serious can take credit for POE
        slot_is_governing = True
        row = row_slot_threshold
        impact_distance_m = slot_exceeded_distance

        if check_slot:
            if check_probit and slot_exceeded_distance > prob_1_pct_distance:
                slot_is_governing = False
        else:
            slot_is_governing = False
        
        
        if not slot_is_governing:
            row = row_prob_threshold
            impact_distance_m = prob_1_pct_distance
        
        impact_area_m2 = 0
        if impact_distance_m > 0:
            time_min = row['time_min']
            dw_dist_m = impact_distance_m
            conc_ppm = row['conc_ppm']

            from py_lopa.calcs.pasquill_gifford_dispersion_calcs import Pasquill_Gifford_Dispersion_Calcs
            pg = Pasquill_Gifford_Dispersion_Calcs(ta = self.chems.ta)
            impact_area_m2 = pg.get_total_impacted_area_m2(wx_enum=self.wx_enum, x = dw_dist_m, conc_ppm=conc_ppm, time_min=time_min)
            self.dose_probit_data_for_kml_zxy = pg.data_for_kml_zxy

        conseq_list = []
        cat = cd.CAT_SERIOUS
        for i in range(len(ser_crit_cat)-1,-1,-1):
            cat_score = i + 3 # this will be used to judge who is worse when accounting for POE
            cat = ser_crit_cat[i]
            req_area_m2 = req_area_m2_dict[cat]
            if req_area_m2 <= 0:
                continue
            
            # get the area of the impacted zone.  This is a gaussian distribution between the source and the conc at the threshold beyond which people can
            # safely escape
            
            if impact_area_m2 >= req_area_m2:
                break


        if cat == cd.CAT_SERIOUS:
            if impact_distance_m == 0:
                impact_area_m2 = 0.01 * req_area_m2
    
        prob = min(1, impact_area_m2 / req_area_m2)
        prob = helpers.float_value_rounded_to_x_sig_figs(prob, 2)
        log_poe_float = self.calc_probability_of_exposure(prob = prob, return_floats=True)
        cat_plus_poe = cat_score + log_poe_float
        poe_string = self.calc_probability_of_exposure(prob = prob)

        conseq_list.append({
            cd.PROBABILITY_OF_EXPOSURE: poe_string,
            cd.CAT_SCORE: cat_score,
            cd.CAT_TITLE: cat,
            cd.IMPACT_DISTANCE_M: impact_distance_m,
            cd.IMPACT_AREA_M2: impact_area_m2,
            cd.REQUIRED_AREA_M2: req_area_m2,
            cd.CAT_PLUS_POE: cat_plus_poe,
            cd.LOG_POE: log_poe_float
        })

        if len(conseq_list) == 0:
            return

        return conseq_list




    def set_consequence_by_tier_designation(self):
        release_indoors = self.mi.release_indoors
        
        # self.onsite_offsite_conseq_dict[cd.CLASS_ONSITE][cd.HAZARD_TYPE_FLASH_FIRE][cd.RELEASE_TIER_DESIGNATION] = pd.NA

    def get_lfl_distance(self):

        if self.hazard_type != cd.HAZARD_TYPE_FLASH_FIRE:
            return

        min_ht_m = 0
        max_ht_m = Consts.ELEVATION_RANGE_FOR_CONC_EVAL_M

        if self.mi is not None:
            max_ht_m = max(self.mi.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS)
        
        lfl_volf = self.chems.flam_loc

        lfl_dist_m = self.get_min_dist_at_conc(lfl_volf, min_ht_m, max_ht_m)

        self.onsite_offsite_conseq_dict[cd.CLASS_ONSITE][self.wx_enum][cd.HAZARD_TYPE_FLASH_FIRE][cd.DISTANCE_TO_LFL_M] = lfl_dist_m


    def get_building_result(self):
        
        bldg:Building
        for bldg in self.mi.bldgs:
            if pd.isna(bldg.dist_m):
                continue
            bldg.prep_targ_dicts(self.use_one_met_condition_worst_case, self.wx_enum, self.hazard_type)
            bldg.set_conseq_targs(self.hazard_type, self.analysis_df)
            bldg.set_conseq_results(self.conc_pfls, self.use_one_met_condition_worst_case, self.wx_enum, self.hazard_type, self.duration_sec, self.chems)
                                
    def get_offsite_result(self):
        self.set_conseq_by_analysis_df(class_= cd.CLASS_OFFSITE, haz_type=self.hazard_type)

    def set_conseq_by_analysis_df(self, class_, haz_type, ignore_higher_level_impacts = False):
        df = self.analysis_df

        curr_haz_df = df[df[cd.KEYS_TARG_AND_TYPE_FLAM_OR_INHAL] == haz_type]
        if len(curr_haz_df) == 0:
            return
        curr_haz_class = curr_haz_df[curr_haz_df[cd.CLASS_TITLE] == class_]
        cats_to_check = cd.CAT_ALL_CATEGORIES

        if ignore_higher_level_impacts:
            cats_to_check = cats_to_check[0:len(cats_to_check)-len(cd.CAT_HIGHER_LEVEL_IMPACT)]
        
        min_ht_m = 0
        max_ht_m = Consts.ELEVATION_RANGE_FOR_CONC_EVAL_M
        if self.mi is not None:
            max_ht_m = max(self.mi.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS)

        onsite_consequence_set = False
        offsite_consequence_set = False
        for i in range(len(cats_to_check)-1,-1,-1):
            cat = cats_to_check[i]
            row = curr_haz_class[curr_haz_class[cd.CAT_TITLE] == cat]
            targ_conc_volf = helpers.get_data_from_pandas_series_element(row[cd.KEYS_TARG_AND_TYPE_CONC_VOLF])
            if pd.isna(targ_conc_volf):
                continue
            if not onsite_consequence_set:
                min_dist_m_to_conc = self.get_min_dist_at_conc(targ_conc_volf, min_ht_m, max_ht_m)
                if min_dist_m_to_conc is not pd.NA:
                    if class_ == cd.CLASS_ONSITE:
                        if min_dist_m_to_conc > 0:
                            onsite_consequence_set = True
                            self.onsite_offsite_conseq_dict[cd.CLASS_ONSITE][self.wx_enum][haz_type] = {
                                cd.CAT_TITLE: cat,
                                cd.IMPACT_DISTANCE_M: min_dist_m_to_conc,
                                cd.IMPACT_AREA_M2: pd.NA,
                                cd.REQUIRED_AREA_M2: pd.NA,
                                cd.PROBABILITY_OF_EXPOSURE: pd.NA
                            }

            if not offsite_consequence_set:
                max_dist_m_at_conc = self.get_max_dist_at_conc(targ_conc_volf, min_ht_m, max_ht_m)
                if class_ == cd.CLASS_OFFSITE:
                    targ_dist_m = helpers.get_data_from_pandas_series_element(row[cd.KEYS_TARG_AND_TYPE_DIST_M])
                    if pd.isna(targ_dist_m):
                        targ_dist_m = self.mi.CONSEQUENCE_DISTS[cd.KEYS_DIST_TO_OFFSITE_M]
                    targ_dist_m *= (1 - Consts.TOLERANCE_CONSEQUENCE_RESULT_PCT)
                    if not pd.isna(max_dist_m_at_conc):
                        if max_dist_m_at_conc >= targ_dist_m:
                            self.onsite_offsite_conseq_dict[cd.CLASS_OFFSITE][self.wx_enum][haz_type] = cat
                            offsite_consequence_set = True
            
            if offsite_consequence_set and onsite_consequence_set:
                break
    
    def get_min_dist_at_conc(self, targ_conc_volf, min_ht_m, max_ht_m):
        flat_proj = Flattening(conc_pfls=self.conc_pfls)
        min_dist_m = flat_proj.calc_min_dist_at_conc(targ_conc_volf, min_ht_m, max_ht_m)

        return min_dist_m

    def get_max_dist_at_conc(self, targ_conc_volf, min_ht_m, max_ht_m):
        # Flattening will analyze the max conc of the cloud to determine the extent of impact to 
        # minor/moderate if applicable for personnel out in the plant or inside production bldgs
        # flatten studies multiple elevations and projects the cloud extent each elevation to the ground,
        # this is effectively the "shadow" of the plume.
        # for example, if the required concentration is ERPG-3, flatten will calculate the distance of the "shadow" 
        # where concentration is ERPG-3 or greater.  
        
        flat_proj = Flattening(conc_pfls=self.conc_pfls)
        max_dist_m = flat_proj.calc_max_dist_at_conc(targ_conc_volf, min_ht_m, max_ht_m)

        return max_dist_m
    
    def get_max_dist_at_conc_onsite_wrapper(self, haz_type, haz_class, conseq_cat):

        df = self.analysis_df

        curr_haz_df = df[df[cd.KEYS_TARG_AND_TYPE_FLAM_OR_INHAL] == haz_type]
        if len(curr_haz_df) == 0:
            return pd.NA
        curr_haz_class = curr_haz_df[curr_haz_df[cd.CLASS_TITLE] == haz_class]

        row = curr_haz_class[curr_haz_class[cd.CAT_TITLE] == conseq_cat]
        targ_conc_volf = helpers.get_data_from_pandas_series_element(row[cd.KEYS_TARG_AND_TYPE_CONC_VOLF])

        min_ht_m = 0
        max_ht_m = Consts.ELEVATION_RANGE_FOR_CONC_EVAL_M
        if self.mi is not None:
            max_ht_m = max(self.mi.ONSITE_EVAL_CONC_FOOTPRINT_ELEVATIONS)

        max_dist = self.get_max_dist_at_conc(targ_conc_volf=targ_conc_volf, min_ht_m=min_ht_m, max_ht_m=max_ht_m)

        return max_dist

