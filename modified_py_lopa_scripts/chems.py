import numpy as np
import pandas as pd
import sys
import copy

from pypws.enums import Phase
from pypws.calculations import VesselLeakCalculation
from pypws.entities import FlashResult

from py_lopa.data.exception_enum import Exception_Enum
from py_lopa.calcs.consts import Consts
from py_lopa.calcs import helpers, thermo_pio
from py_lopa.classes.toxicological_analysis import Toxicological_Analysis
from py_lopa.phast_io import phast_prep
from py_lopa.data.tables import Tables
from py_lopa.calcs.get_phys_props import Get_Phys_Props

consts = Consts()

ra = consts.RELEASE_STREAM_ATTRIBS
cd = consts.CONSEQUENCE_DATA

class Chems:

    def __init__(self, mi = None):
        
        self.cheminfo = helpers.get_dataframe_from_csv(Tables().CHEM_INFO, encoding='cp1252')
        
        if mi is None:
            return

        self.chemmix_chemname = []
        self.material_ids = []
        self.tables = None
        self.ff_credible = False
        self.inhal_credible = False
        self.inhalation_locs = {}
        self.flam_loc = 0
        self.flam_loc_fracts = {}
        self.vapor_pressure_data = {}
        self.chem_of_greatest_hazard_vapor_volf_dose_and_probit_consts = {}

        self.chemmix_casno = mi.CHEM_MIX
        self.composition = mi.COMPOSITION
        self.comp_in_moles = mi.COMPOSITION_IS_IN_MOLES
        self.mi = mi
        self.ta = Toxicological_Analysis()
        self.ta.offsite_dist = self.mi.OFFSITE_DIST
        self.ta.catastrophic_release = self.mi.CATASTROPHIC_VESSEL_FAILURE
        self.shi_analysis_data = {}

        self.chems_to_cas_nos()
        # self.handle_chemicals_incompatible_with_pws_models()
        self.cas_nos_to_chems()
        self.check_if_chems_are_hazardous()

        self.mws = helpers.get_mws(chem_mix=self.chemmix_casno, cheminfo=self.cheminfo)
        self.material_ids = helpers.get_material_ids(self.chemmix_casno, self.cheminfo)
        self.normalize_to_mole_fraction()
        self.join_chems_lists_in_dataframe()

    def check_if_chems_are_hazardous(self):
        if not self.mi.FLASH_FIRE and not self.mi.INHALATION:
            raise Exception(Exception_Enum.NO_TARGET_STUDY_SELECTED)
        df = self.cheminfo

        # check all materials.
        # ff is considered credible if at least 1 chem has lel in cheminfo table
        # inhalation is considered credible if at least 1 chem has loc-3 in cheminfo table.
        # once you find a qualifying chemical, stop searching.

        for i in range(len(self.chemmix_casno)):
            
            row = df[df['cas_no'] == self.chemmix_casno[i]]
            
            if self.mi.FLASH_FIRE:
                lel = row['lel'].values[0]
                if not pd.isna(lel): 
                    if lel > 0 and lel < 1e6:
                        self.ff_credible = True
                        #if inhalation already found to be credible, stop looking
                        if self.inhal_credible:
                            break

            if self.mi.INHALATION:
                loc3 = row['loc_3'].values[0]
                if not pd.isna(loc3):
                    if loc3 > 0 and loc3 < 1e6:
                        self.inhal_credible = True
                        #if flash fire already found to be credible, stop looking
                        if self.ff_credible:
                            break

        if (self.mi.FLASH_FIRE and not self.ff_credible):
            self.mi.FLASH_FIRE = False
            self.mi.LOG_HANDLER('Flash Fire will not be analyzed as database currently has no associated flammability limits for the components of the releasing stream')
            
        if (self.mi.INHALATION and not self.inhal_credible):
            self.mi.INHALATION = False
            self.mi.LOG_HANDLER('Toxicity will not be analyzed as database currently has no associated toxic limits for the components of the releasing stream')

        # if not self.ff_credible and not self.inhal_credible:
        #     raise Exception(Exception_Enum.NO_HAZARDOUS_LIMITS_AVAILABLE)
        
    def chems_to_cas_nos(self):
        for c in range(len(self.chemmix_casno)):
            self.chemmix_casno[c] = helpers.get_cas_no(self.chemmix_casno[c], self.cheminfo)
    
    def cas_nos_to_chems(self):
        for c in range(len(self.chemmix_casno)):
            chemnm = helpers.get_chem_name(self.chemmix_casno[c], self.cheminfo)
            self.chemmix_chemname.append(chemnm)

    def handle_chemicals_incompatible_with_pws_models(self):
        # 9/12/24 - the PWS multicomponent model currently has issues with dispersing material mixtures containing air and water.
        # the plan until a fix is pushed to PWS is as follows:
        # 1. any components with air have the air component swapped out for nitrogen
        # 2. if water is in the mixture, the mixture model will be specified as pseudocomponent.

        # to replace air with nitrogen:
        # 1.  check for presence of air. 
            # the component of the mixture that is comprised of air will be added to any existing nitrogen in the mixture
        # 2.  check if nitrogen is in the mixture.  if so, store the index of nitrogen in the mixture

        # 9/20/24 - updates have been pushed to production.  results for materials with water/air should be accurate.
        # the call to this method has been commented out.

        idx_nitrogen = None
        idx_air = None
        air_component = None
        for i in range(len(self.chemmix_casno)):
            cas_no = self.chemmix_casno[i]
            if cas_no == cd.CAS_NOS_OF_CONCERN[cd.AIR]:
                air_component = self.composition[i]
                idx_air = i
            if cas_no == cd.CAS_NOS_OF_CONCERN[cd.WATER]:
                self.mi.USE_MULTICOMPONENT_METHOD = False
            if cas_no == cd.CAS_NOS_OF_CONCERN[cd.NITROGEN]:
                idx_nitrogen = i
        
        if air_component is not None:
            if idx_nitrogen is not None:
                self.composition[idx_nitrogen] += air_component
            else:
                self.chemmix_casno.append(cd.CAS_NOS_OF_CONCERN[cd.NITROGEN])
                self.composition.append(air_component)
        
            del self.chemmix_casno[idx_air]
            del self.composition[idx_air]

            self.mi.CHEM_MIX = self.chemmix_casno
            self.mi.COMPOSITION = self.composition

    def normalize_to_mole_fraction(self):
        if not self.comp_in_moles:
            for i in range(len(self.composition)):
                c = self.chemmix_casno[i]
                mw = self.mws[i]
                self.composition[i] /= mw

        tot = sum(self.composition)
        for i in range(len(self.composition)):
            self.composition[i] /= tot
    
    def join_chems_lists_in_dataframe(self):

        d = {
                ra.CAS_NO: self.chemmix_casno,
                ra.CHEM_NAME: self.chemmix_chemname,
                ra.MOLE_FRACT: self.composition,
                ra.MW: self.mws,
                ra.MATERIAL_COMPONENT_ID: self.material_ids
            }

        d_cols = d.keys()
        
        df = pd.DataFrame(d, columns = d_cols)

        self.releasing_stream_chem_data = copy.deepcopy(df)
        
        # since chemical list and composition may be updated in the "releasing_stream_chem_data" dataframe, set some potentially 
        # conflicting attributes to 'pd.NA'.  These attributes should be deleted, but it 
        # appears to require special handling.  In the case that this causes a future issue, simply setting the attribs to 'pd.NA'.

        attrib_to_set_to_none = ['chemmix_casno', 'chemmix_chemname', 'composition', 'composition_mass_fract', 'mws', 'phast_ids']

        for attrib in attrib_to_set_to_none:
            setattr(self, attrib, pd.NA)

    def final_system_checks(self, flashresult:FlashResult, state, material):
        
        tol_c = consts.TOLERANCE_DEGREES_C_FROM_RELIEF_PHASE

        mi = self.mi

        # state calcs are normally pt flash based.
        # previously, these were pulled from mi.
        # however, using multicomponent (and even potentially some pc) configurations
        # with two-phase results requires a liquid fraction be present in the flash calcs.
        # if a fraction is provided as a model input, it will be used in the flash calc.
        # most times it will not have been provided, so it is wrapped as None.

        lf = None
        if mi.EXIT_VAPOR_MASS_FRACTION is not None:
            lf = 1 - mi.EXIT_VAPOR_MASS_FRACTION

        # the 'mole_fraction'' attribute is the composition as entered. 
        # if entered as a vapor phase and the flash calc returns a
        # liquid, pad vessel with n2

        #get flash point and nbp before any potential padding
        # if self.mi.FLASH_FIRE:
        #     self.fp_K = thermo_pio.flash_point_deg_K(chems = self)
        # self.nb_K = thermo_pio.nbp_deg_K(chems = self)

        mi_rp = mi.RELIEF_PHASE
        rp_enum = consts.RELIEF_PHASE
        

        # if a relief calc is available, the reported relief exit conditions should have a
        # phase similar to what is reported from Phast.  If it is more than a few
        # degrees from a phase change, exit program.  If close enough, tweak temperature parameter.
        

        #Test if conditions are vastly different than those calculated by Phast
        if mi_rp == rp_enum.VAPOR:
            if flashresult.fluid_phase != Phase.VAPOUR:
                # with an available relief calc, the conditions should within a few degrees of appropriate phase
                if mi.RELIEF_CALC_AVAILABLE:
                    if abs(flashresult.dew_point_temperature - (mi.TEMP_K)) > tol_c:
                        raise Exception(Exception_Enum.INVALID_RELIEF_VAPOR_SPECIFIED_LIQUID_CALCULATED)
                    mi.TEMP_K = flashresult.dew_point_temperature + 0.01
                else:
                    #without a relief calc, pad the vessel with n2 until flash result is no longer subcooled.
                    mi.TEMP_K = self.pad_with_n2_to_flash_result_dew_point_temp_return_final_temp(flashresult = flashresult, mi = mi, state = state)
                    self.fix_mole_fractions()
                
                state = phast_prep.prep_state(
                    press_pa=mi.PRESS_PA, 
                    temp_K=mi.TEMP_K, 
                    lf=lf, 
                    use_multicomponent_modeling=mi.USE_MULTICOMPONENT_METHOD
                )
                material = phast_prep.prep_material(self)
                flashresult = phast_prep.flash_calc(state, material)

                while flashresult.fluid_phase != Phase.VAPOUR:
                    mi.TEMP_K += 1
                    state = phast_prep.prep_state(
                        press_pa=mi.PRESS_PA, 
                        temp_K=mi.TEMP_K, 
                        lf=lf, 
                        use_multicomponent_modeling=mi.USE_MULTICOMPONENT_METHOD
                    )
                    flashresult = phast_prep.flash_calc(state, material)

        if mi_rp == rp_enum.LIQUID:
            if flashresult.fluid_phase != Phase.LIQUID:
                if abs(mi.TEMP_K - flashresult.bubble_point_temperature) > tol_c:
                    raise Exception(Exception_Enum.INVALID_RELIEF_LIQUID_SPECIFIED_VAPOR_CALCULATED)
                mi.TEMP_K = flashresult.bubble_point_temperature - 1
                state = phast_prep.prep_state(
                    press_pa=mi.PRESS_PA, 
                    temp_K=mi.TEMP_K, 
                    lf=lf, 
                    use_multicomponent_modeling=mi.USE_MULTICOMPONENT_METHOD
                )
                flashresult = phast_prep.flash_calc(state, material)
                
                while flashresult.fluid_phase != Phase.LIQUID:
                    mi.TEMP_K -= 1
                    state = phast_prep.prep_state(
                        press_pa=mi.PRESS_PA, 
                        temp_K=mi.TEMP_K, 
                        lf=lf,
                        use_multicomponent_modeling=mi.USE_MULTICOMPONENT_METHOD
                    )
                    flashresult = phast_prep.flash_calc(state, material)
        
        if mi_rp == rp_enum.TWO_PHASE:
            if flashresult.fluid_phase != Phase.TWO_PHASE:
                self.mi.LOG_HANDLER(f'WARNING - Two-phase selected as source phase.  However, the source phase is calculated as {flashresult.fluid_phase.name}. It will be modeled as relieving from that phase.')


        self.calc_vapor_pressure_data(mi)

        ans = { 'state': state, 
                'flashresult': flashresult, 
                'model_inputs': mi,
                'material': material}

        return ans

    def fix_mole_fractions(self):

        df_release = self.releasing_stream_chem_data

        molfs = df_release[ra.MOLE_FRACT].to_list()

        tot_moles = sum(molfs)
        for i in range(len(molfs)):
            molfs[i] /= tot_moles

        self.releasing_stream_chem_data[ra.MOLE_FRACT] = molfs
        

    
    def calc_vapor_pressure_data(self, mi):
        vp_args = {
            'chems': self,
            'mi': mi
        }
        self.vapor_pressure_data = thermo_pio.vpress_pa_and_vapor_phase_comp_and_component_vapor_pressures(args=vp_args)

    def get_vapor_comp(self, data_dict):
        if sum(data_dict['ys']) > 0:
            return data_dict['ys']
        k_times_zi = data_dict['k_times_zi']
        return helpers.normalize_fractions(k_times_zi)


    def get_data_for_shi_analysis(self, mi, vessel_leak_calc):
        gpp_storage_condits = Get_Phys_Props(chem_mix=mi.CHEM_MIX, molar_basis=True, chem_composition=mi.COMPOSITION, temp_C=mi.TEMP_K - 273.15, press_psig=mi.PRESS_PA * 14.6959 / 101325 - 14.6959)
        gpp_storage_condits.run()
        data_dict_storage = gpp_storage_condits.data_dict
        ave_nbp_deg_c = data_dict_storage['ave_nbp_C']
        overall_ave_mw = data_dict_storage['ave_mw']

        fin_state = vessel_leak_calc.discharge_records[0].final_state
        disch_temp_c = fin_state.temperature - 273.15
        gpp_discharge_condits = Get_Phys_Props(chem_mix=mi.CHEM_MIX, molar_basis=True, chem_composition=mi.COMPOSITION, temp_C=disch_temp_c, press_psig = 0)
        gpp_discharge_condits.run()
        data_dict_discharge = gpp_discharge_condits.data_dict

        ans = {
            'ave_nbp_deg_c': ave_nbp_deg_c,
            'overall_ave_mw': overall_ave_mw,
            'discharge': {
                'temp_deg_c': disch_temp_c,
                'vapor_mol_composition': self.get_vapor_comp(data_dict_discharge),
                'overall_vapor_mole_fraction': data_dict_discharge['overall_vapor_mole_fraction'],
                'total_vapor_moles': self.get_vapor_moles(mi, data_dict_discharge),
                'component_vapor_moles': [self.get_vapor_moles(mi, data_dict_discharge, i) for i in range(len(mi.CHEM_MIX))],
                'total_moles': self.get_overall_moles(mi, data_dict_discharge),
                'component_total_moles': [self.get_overall_moles(mi, data_dict_discharge, i) for i in range(len(mi.CHEM_MIX))],
            },
            'storage': {
                'temp_deg_c': mi.TEMP_K - 273.15,
                'vapor_mol_composition': self.get_vapor_comp(data_dict_storage),
                'overall_vapor_mole_fraction': data_dict_storage['overall_vapor_mole_fraction'],
                'total_vapor_moles': self.get_vapor_moles(mi, data_dict_storage),
                'component_vapor_moles': [self.get_vapor_moles(mi, data_dict_storage, i) for i in range(len(mi.CHEM_MIX))],
                'total_moles': self.get_overall_moles(mi, data_dict_storage),
                'component_total_moles': [self.get_overall_moles(mi, data_dict_storage, i) for i in range(len(mi.CHEM_MIX))],
            }
        }

        return ans
    
    def get_vapor_moles(self, mi, data_dict, component_idx = None):
        overall_ave_mw = data_dict['ave_mw']
        mass_released_kg = mi.STORAGE_MASS_KG
        moles_released_kmol = mass_released_kg / overall_ave_mw
        overall_vapor_mole_fraction = data_dict['overall_vapor_mole_fraction']
        vapor_moles_released = moles_released_kmol * overall_vapor_mole_fraction
        if component_idx is not None:
            vapor_mol_composition = data_dict['vapor_mol_composition']
            vap_fract_component = vapor_mol_composition[component_idx]
            vapor_moles_released *= vap_fract_component
        return vapor_moles_released

    def get_overall_moles(self, mi, data_dict, component_idx = None):
        overall_ave_mw = data_dict['ave_mw']
        mass_released_kg = mi.STORAGE_MASS_KG
        moles_released_kmol = mass_released_kg / overall_ave_mw
        if component_idx is not None:
            molf = mi.COMPOSITION[component_idx]
            moles_released_kmol *= molf
        return moles_released_kmol



    def set_limits_of_concern(self, mi, flash_results, vessel_leak_calc:VesselLeakCalculation):
        df = self.cheminfo

        df_release = self.releasing_stream_chem_data

        data = thermo_pio.get_vapor_phase_comp_and_flash_calc_from_discharge_result(vessel_leak_calc, self, mi)
        self.vapor_molfs = data['vap_molfs']
        self.flash_data = data['flash_data']

        self.shi_analysis_data = self.get_data_for_shi_analysis(mi, vessel_leak_calc)

        # 'shi' is the Substance Hazard Index.  It is the vapor-phase conc / loc.  
        # the greater the value, the more severe the impact from that chemical

        # *** xxx debug apple
        # for this implementation, the shi calc will be run for both storage and discharge.
        # this will allow for analysis of shi versus consequence.
        # storage is run first and all the normal calculations continue.  the storage values are written to the 
        # appropriate dictionary.  then, the discharge is run.  it also stores to the appropriate dictionary.
        # discharge conditions also overwrites the needed calcs for the rest of the 
        # model to run so that the dispersion and results are based on the discharge condition.

        mixture_cas_nos = df_release[ra.CAS_NO].to_list()
        for condition in ['storage', 'discharge']:
            worst_shi = -1
            worst_shi_idx = -1
            targ_loc = 0
            lel_contribs = 0
            molfs = self.shi_analysis_data[condition]['vapor_mol_composition']
            for i in range(len(df_release)):
                c = mixture_cas_nos[i]
                molf = molfs[i]
                row = df[df['cas_no'] == c]
                if mi.INHALATION:
                    loc3 = row['loc_3'].values[0]
                    shi = -1
                    if not pd.isna(loc3) and loc3 < 1e6:
                        shi = molf / loc3
                    if shi > worst_shi and molf > 0:
                        worst_shi = shi
                        worst_shi_idx = i
                        self.shi_analysis_data[condition]['shi_tox'] = worst_shi
                        # the target concentration will be worst-case toxic conc, but corrected 
                        # for its overall concentration in the vapor phase. 
                        # that is, if chemical X is lethal at 100 ppm 
                        # and the vapor phase is 10% chemical X, then the 
                        # concentration of the vapor phase is toxic at 100 ppm / 10% = 1000 ppm

                        targ_loc = loc3 / molf
                
                if mi.FLASH_FIRE:
                    #determine average lel of vapor phase based on le chatlier's mixing rule
                    lel = 1000000000000 #large value to negate contribution of non-flammable to average LFL
                    lel_val = row['lel'].values[0]
                    if not pd.isna(lel_val):
                        if lel_val > 0 and lel_val < 1e6:
                            lel = lel_val 
                    lel_contribs += molf / lel

            # set this such that flammable models can utilize the building infilt analysis in site toxicological analysis
            self.ta.volf = 1
            

            if mi.INHALATION:

                if targ_loc <= 0 or targ_loc >= 1e6:
                    mi.LOG_HANDLER('inhalation assessment selected, but releasing stream is not estimated to be acutely toxic given the input conditions.')
                    mi.INHALATION = False
                else:
                    self.inhalation_locs = {}
                    self.inhalation_locs[3] = targ_loc
                    cas_worst_shi = mixture_cas_nos[worst_shi_idx]
                    row_worst_shi = df[df['cas_no'] == cas_worst_shi]
                    molf_worst_shi = self.vapor_molfs[worst_shi_idx]

                    # populate inhalation loc dict with LOC-1, -2, and -3
                    # the loc-3 value will be available for all chemicals with at least 
                    # one toxic limit.  any missing values will be 1/7th of the higher value.

                    
                    for i in range(2,0,-1):
                        loc = row_worst_shi['loc_' + str(i)].values[0]
                        if pd.isna(loc) or loc >= 1e6:
                            loc = self.inhalation_locs[i + 1] / 7
                            
                        self.inhalation_locs[i] = loc / molf_worst_shi

                        #stored as vol fract for use in phast model
                        self.inhalation_locs[i] /= 1e6
                    
                    self.inhalation_locs[3] /= 1e6

                    self.chem_of_greatest_hazard_vapor_volf_dose_and_probit_consts = {
                        ra.CAS_NO: mixture_cas_nos[worst_shi_idx],
                        cd.CONC_VOLF_TITLE: molf_worst_shi
                    }

                    volf_for_dose_probit = 1
                    if mi.USE_DOSE_AND_PROBIT:
                        volf_for_dose_probit = molf_worst_shi

                    self.ta.get_dose_and_probit_constants(cas_no = mixture_cas_nos[worst_shi_idx], volf = volf_for_dose_probit, cheminfo = self.cheminfo)

            if mi.FLASH_FIRE:
                
                # ave LFL = 1/lel_contribs.  If lel_contribs is less than 1e-6,
                # then the LFL > 1,000,000 ppm, which is not credible.

                if lel_contribs > 1e-6:
                    self.flam_loc = 1 / lel_contribs

                    self.shi_analysis_data[condition]['lfl'] = self.flam_loc

                    #stored as vol fract for use in phast model
                    self.flam_loc /= 1e6
                    
                    ids = [1, 2, 3]
                    fracts = [0.25, 0.50, 1.00]

                    self.flam_loc_fracts = {i : f * self.flam_loc for (i, f) in zip (ids, fracts)}
                    
                else:
                    mi.LOG_HANDLER('Flammable assessment selected, but releasing stream is not estimated to be flammable given the input conditions.')
                    mi.FLASH_FIRE = False
                    mi.VAPOR_CLOUD_EXPLOSION = False
        
        mi.VALID_HAZARDS[cd.HAZARD_TYPE_INHALATION] = mi.INHALATION
        if not mi.FLASH_FIRE:
            mi.VALID_HAZARDS[cd.HAZARD_TYPE_FLASH_FIRE] = False
            mi.VALID_HAZARDS[cd.HAZARD_TYPE_VCE] = False

    def pad_with_n2_to_flash_result_dew_point_temp_return_final_temp(self, flashresult, mi, state):

        cas_nos = self.releasing_stream_chem_data[ra.CAS_NO].to_list()
        if cd.NITROGEN_CAS_NO not in cas_nos:
            self.append_n2_to_relieving_stream_data()
        liquid_comp = thermo_pio.liquid_phase_comp_iterate_to_flashresult_at_bubble_point_temp(chems = self, state = state)
        self.releasing_stream_chem_data[ra.MOLE_FRACT] = liquid_comp
        material = phast_prep.prep_material(self)
        flashresult = phast_prep.flash_calc(state, material)

        # if after padding with n2, the flash calc still shows as liquid,
        # increment temperature until it is vapor. 

        lf = None
        if mi.EXIT_VAPOR_MASS_FRACTION is not None:
            lf = 1 - mi.EXIT_VAPOR_MASS_FRACTION

        if flashresult.fluid_phase != Phase.VAPOUR:
            if flashresult.dew_point_temperature - mi.TEMP_K > 5:
                raise Exception(Exception_Enum.ERROR_CANT_MODEL_VAPOR_RELIEF_FOR_GIVEN_COMPOSITION)

            mi.TEMP_K = flashresult.dew_point_temperature
            state = phast_prep.prep_state(
                press_pa=mi.PRESS_PA, 
                temp_K=mi.TEMP_K, 
                lf=lf, 
                use_multicomponent_modeling=mi.USE_MULTICOMPONENT_METHOD
            )
            flashresult = phast_prep.flash_calc(state, material)
            while flashresult.fluid_phase != Phase.VAPOUR:
                mi.TEMP_K += 1
                state = phast_prep.prep_state(
                    press_pa=mi.PRESS_PA, 
                    temp_K=mi.TEMP_K, 
                    lf=lf, 
                    use_multicomponent_modeling=mi.USE_MULTICOMPONENT_METHOD
                )
                flashresult = phast_prep.flash_calc(state, material)

        return mi.TEMP_K
    
    def append_n2_to_relieving_stream_data(self):

        df = self.cheminfo
        row = df[df['cas_no'] == cd.NITROGEN_CAS_NO]
        
        mws = self.releasing_stream_chem_data[ra.MW].to_list()

        mw_n2 = row['mw'].values[0]
        mat_id_n2 = row['mat_comp_id'].values[0]
        mws.append(mw_n2)
        
        # create a row to be added to the releasing stream chemical properties dataframe
        # for N2.  the "mol_fraction" for N2 is the initial guess on its liquid mol fract.  this value will be
        # optimized such that the system temp is bubble point at system pressure.

        df_n2 = pd.DataFrame({
            ra.CAS_NO: [cd.NITROGEN_CAS_NO],
            ra.CHEM_NAME: ['nitrogen'],
            ra.MOLE_FRACT: [0.001],
            ra.MW: [mw_n2],
            ra.MATERIAL_COMPONENT_ID: [mat_id_n2],
        })

        self.releasing_stream_chem_data = pd.concat([self.releasing_stream_chem_data, df_n2], ignore_index=True)