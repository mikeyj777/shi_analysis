import sys
import numpy as np
import pandas as pd

# import thermo
# from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL, IGMIX
# from thermo.interaction_parameters import IPDB
# from thermo.volume import VolumeLiquidMixture, VolumeLiquid

from py_lopa.calcs import helpers, dippr_eqns, thermo_pio
from py_lopa.calcs.consts import Consts
from py_lopa.data.tables import Tables


# pull chemical mixture data and return physical property data
# final data stored in Get_Phys_Props.data in JSON format.

class Get_Phys_Props:
    
    def __init__(self, chem_mix, molar_basis, chem_composition, temp_C, press_psig):
        self.data = None
        self.vp_data = None
        self.temp_C = temp_C
        self.press_psig = press_psig
        self.cheminfo = helpers.get_dataframe_from_csv(Tables().CHEM_INFO, encoding='cp1252')
        self.liq_dens_df = helpers.get_dataframe_from_csv(Tables().LIQ_DENSITY_DATA)
        self.initialize_properties()
        self.parse_inputs(chem_mix, chem_composition, molar_basis, temp_C, press_psig)

    def run(self):
        self.set_mol_and_mass_comps()
        self.calc_ave_mw()
        self.calc_nbp()
        self.get_liq_and_vap_compositions_ideal()
        if not self.good_data:
            self.make_output_dict_and_json_string()
            return
        self.get_vpress_data()
        self.get_liq_densities()
        self.calc_liquid_vol_comp()
        self.calc_liquid_level_and_vap_density()
        self.calc_ave_mw_vap()
        self.make_output_dict_and_json_string()

    def initialize_properties(self):
        self.good_data = True
        self.mass_composition = []
        self.molar_composition = []
        self.nbp_Cs = []
        self.liquid_mass_composition = []
        self.liquid_molar_composition = []
        self.liquid_vol_composition = []
        self.vapor_mass_composition = []
        self.vapor_molar_composition = []
        self.mws = []
        self.component_vapor_densities_lb_gal = []
        self.component_vapor_densities_kg_m3 = []
        self.liquid_level = None
        self.ave_mw = None
        self.nbp_K = None
        self.nbp_C = None
        self.liq_dens_data = {
            'densities_lb_gal': None,
            'densities_kg_m3':  None,
            'ave_dens_kg_m3': None,
            'ave_dens_lb_gal': None,
        }

    def parse_inputs(self, chem_mix, chem_composition, molar_basis, temp_C, press_psig):
        if isinstance(chem_mix, str):
            chem_mix = chem_mix.split('|')
        self.chem_mix = helpers.get_cas_nos(chem_mix, self.cheminfo)
        if isinstance(chem_composition, str):
            chem_composition = chem_composition.split('|')
        self.chem_composition = np.asarray(chem_composition)
        self.chem_composition = self.chem_composition.astype(np.float64)
        self.chem_composition = helpers.normalize_fractions(self.chem_composition)
        self.chem_composition = self.chem_composition.tolist()
        self.temp_K = temp_C + 273.15
        self.press_Pa = (press_psig + 14.6959) * 101325 / 14.6959
        self.parse_mol_basis_input(molar_basis)
        self.mws = helpers.get_mws(self.chem_mix, self.cheminfo)
    
    def parse_mol_basis_input(self, m):
        if type(m) is str:
            m = m[0].lower()
            if m == 't' or m == 'y':
                m = True
            else:
                m = False
        if type(m) is float:
            m = int(m)
        if type(m) is int:
            if m == 1:
                m = True
            else:
                m = False
        self.molar_basis = m

    def set_mol_and_mass_comps(self):
        if self.molar_basis:
            self.molar_composition = self.chem_composition
            self.mass_composition = helpers.mass_fracts_from_mol_fracts(self.chem_composition, self.mws)
        else:
            self.molar_composition = helpers.mol_fracts_from_mass_fracts(self.chem_composition, self.mws)
            self.mass_composition = self.chem_composition

    def calc_ave_mw(self):
        mws = np.asarray(self.mws)
        zs = np.asarray(self.molar_composition)
        self.ave_mw = zs.dot(mws.T)

    def calc_ave_mw_vap(self):
        mws = np.array(self.mws)
        ys = np.array(self.vapor_molar_composition)
        self.ave_mw_vap = ys.dot(mws.T)

    def calc_nbp(self):
        self.nbp_K = thermo_pio.nbp_deg_K(chems = None, x0 = 298.15, mixture_cas_nos = self.chem_mix, mixture_molfs = self.molar_composition, cheminfo = self.cheminfo)
        self.nbp_C = self.nbp_K - 273.15
        self.nbp_Cs = []
        for i in range(len(self.chem_mix)):
            comp_nbp_K = thermo_pio.nbp_deg_K(chems = None, x0 = 298.15, mixture_cas_nos = [self.chem_mix[i]], mixture_molfs = [1], cheminfo = self.cheminfo)
            comp_nbp_C = comp_nbp_K - 273.15
            self.nbp_Cs.append(comp_nbp_C)

    def get_liq_and_vap_compositions_ideal(self):
        data = thermo_pio.ideal_flash_calc_get_vf_xs_ys_mol_basis(chem_mix=self.chem_mix, mws= self.mws, overall_molfs=self.molar_composition, temp_K= self.temp_K, press_Pa=self.press_Pa, cheminfo=self.cheminfo)
        if sum(data['xs']) + sum(data['ys']) == 0:
            self.good_data = False
            return
        self.k_times_zi = data['k_times_zi']
        self.vp_data = data['vpress_data']
        self.vapor_molar_composition = data['ys']
        self.liquid_molar_composition = data['xs']
        self.vapor_mass_composition = helpers.mass_fracts_from_mol_fracts(self.vapor_molar_composition, self.mws)
        self.liquid_mass_composition = helpers.mass_fracts_from_mol_fracts(self.liquid_molar_composition, self.mws)
        self.overall_vapor_mol_percent = data['vf']
        if self.overall_vapor_mol_percent is None:
            self.good_data = False
            return
        self.overall_liquid_mol_percent = 1 - data['vf']
        self.overall_vapor_mass_fraction = data['overall_vapor_mass_fraction']
        if self.overall_vapor_mass_fraction is None:
            self.good_data = False
            return
        self.overall_liquid_mass_percent = 1 - data['overall_vapor_mass_fraction']

    def get_vpress_data(self):
        psat = np.array(self.vp_data)
        xs = np.array(self.molar_composition)
        partial_pressures_pa = np.multiply(xs, psat)
        v_press_pa_total = partial_pressures_pa.sum()
        v_press_psia_total = v_press_pa_total * 14.6959 / 101325
        self.v_press_psig_total = v_press_psia_total - 14.6959
        partial_vps_psig = partial_pressures_pa * 14.6959 / 101325 - 14.6959
        self.comp_vps_psig = partial_vps_psig


    def get_liq_densities(self):
        df = self.liq_dens_df
        mix = self.chem_mix
        tot_vol_m3 = 0
        densities_kg_m3 = []
        densities_lb_gal = []
        tot_vol_ok = True
        for i in range(len(mix)):
            row = df[df['cas_no'] == mix[i]]
            comp_dens_kg_m3 = 0
            comp_dens_lb_gal = 0
            if len(row) == 0:
                continue
            eqn_str = str(int(row['equation_id'].values[0]))
            coeffs = []
            for j in range(5):
                coeffs.append(row['liqdenscoeff' + str(j)].values[0])
                if pd.isna(coeffs[j]):
                    coeffs[j] = 0
            method_to_call = getattr(dippr_eqns, 'eqn_' + eqn_str)
            comp_dens_kmol_m3 = method_to_call(coeffs, self.temp_K)
            comp_dens_kmol_m3 = max(comp_dens_kmol_m3, 0)
            comp_dens_kg_m3 = comp_dens_kmol_m3 * self.mws[i]
            comp_dens_lb_gal = comp_dens_kg_m3 / Consts.GAL_PER_M3 * Consts.LB_PER_KG
            densities_lb_gal.append(comp_dens_lb_gal)
            densities_kg_m3.append(comp_dens_kg_m3)
            if comp_dens_kmol_m3 <= 0:
                continue
            tot_vol_m3 += self.liquid_mass_composition[i] / comp_dens_kg_m3
            
        
        ave_dens_kg_m3 = 0
        ave_dens_lb_gal = 0
        if tot_vol_m3 > 0  and tot_vol_ok:
            ave_dens_kg_m3 = 1 / tot_vol_m3
            ave_dens_lb_gal = ave_dens_kg_m3 / Consts.GAL_PER_M3 * Consts.LB_PER_KG
        
        self.liq_dens_data = {
                'ave_dens_kg_m3': ave_dens_kg_m3,
                'ave_dens_lb_gal': ave_dens_lb_gal,
                'densities_kg_m3': densities_kg_m3,
                'densities_lb_gal': densities_lb_gal
            }
    
    def calc_liquid_vol_comp(self):
        self.liquid_vol_composition = []
        ldd =  self.liq_dens_data['densities_lb_gal']
        for i in range(len(self.liquid_mass_composition)):
            liq_dens = ldd[i]
            vol_comp = 0
            if liq_dens > 0:
                vol_comp = self.liquid_mass_composition[i] / liq_dens
            self.liquid_vol_composition.append(vol_comp)
        tot_vol = sum(self.liquid_vol_composition)
        if tot_vol == 0:
            return
        self.liquid_vol_composition = [x / tot_vol for x in self.liquid_vol_composition]

    def calc_liquid_level_and_vap_density(self):
        liq_mass = self.overall_liquid_mass_percent
        vap_mass = 1 - liq_mass
        vap_vol = 0
        r = Consts.R_PA_M3_KMOL_DEGK
        self.component_vapor_densities_kg_m3 = []
        self.component_vapor_densities_lb_gal = []
        for i in range(len(self.chem_mix)):
            comp_vap_mass = vap_mass * self.vapor_mass_composition[i]
            comp_vap_mol = comp_vap_mass / self.mws[i]
            comp_vap_vol = comp_vap_mol * r * self.temp_K / self.press_Pa
            rho_kg_m3 = 0
            rho_lb_kg = 0
            if comp_vap_mol > 0:
                rho_kg_m3 = comp_vap_mass / comp_vap_vol
                rho_lb_kg = rho_kg_m3 * Consts.LB_PER_KG / Consts.GAL_PER_M3
            self.component_vapor_densities_kg_m3.append(rho_kg_m3)
            self.component_vapor_densities_lb_gal.append(rho_lb_kg)
            vap_vol += comp_vap_vol
        self.vap_dens_kg_m3 = 0
        if vap_vol > 0:
            self.vap_dens_kg_m3 = vap_mass / vap_vol

        self.liquid_level = 0
        liq_ave_dens_kg_m3 =  self.liq_dens_data['ave_dens_kg_m3']
        liq_vol = 0
        if liq_ave_dens_kg_m3 > 0:
            liq_vol = liq_mass / liq_ave_dens_kg_m3
            self.liquid_level = liq_vol / (liq_vol + vap_vol)
        self.overall_density_kg_m3 = 0
        self.overall_density_lb_gal = 0
        tot_vol = liq_vol + vap_vol
        if tot_vol > 0:
            self.overall_density_kg_m3 = 1 / tot_vol
            self.overall_density_lb_gal = self.overall_density_kg_m3 / Consts.GAL_PER_M3 * Consts.LB_PER_KG

    def make_output_dict_and_json_string(self):
        self.data_dict = {
            'temp_C': self.temp_C,
            'press_psig': self.press_psig,
            'overall_mass_percent': self.mass_composition,
            'overall_mole_percent': self.molar_composition,
            'component_nbp_C': self.nbp_Cs,
            'liquid_mass_composition': self.liquid_mass_composition,
            'liquid_molar_composition': self.liquid_molar_composition,
            'liquid_vol_composition': self.liquid_vol_composition,
            'vapor_mass_composition': self.vapor_mass_composition,
            'vapor_mol_composition': self.vapor_molar_composition,
            'mws': self.mws,
            'component_liq_densities_lb_gal': self.liq_dens_data['densities_lb_gal'],
            'component_liq_densities_kg_m3': self.liq_dens_data['densities_kg_m3'],
            'component_vap_densities_lb_gal': self.component_vapor_densities_lb_gal,
            'component_vap_densities_kg_m3': self.component_vapor_densities_kg_m3,
            'liquid_level': self.liquid_level,
            'ave_mw': self.ave_mw,
            'ave_nbp_C': self.nbp_C,
            'liq_ave_dens_kg_m3': self.liq_dens_data['ave_dens_kg_m3'],
            'liq_ave_dens_lb_gal': self.liq_dens_data['ave_dens_lb_gal'],
            'overall_density_kg_m3': self.overall_density_kg_m3,
            'overall_density_lb_gal':  self.overall_density_lb_gal,
            'total_vapor_pressure_psig_at_system_temp': self.v_press_psig_total,
            'component_vapor_pressures_psig_at_system_temp':  self.comp_vps_psig,
            'overall_vapor_mass_fraction': self.overall_vapor_mass_fraction,
            'overall_vapor_mole_fraction': self.overall_vapor_mol_percent,
            'ave_mw_vap': self.ave_mw_vap,
            'k_times_zi': self.k_times_zi,
        }
    
        self.data = helpers.dict_to_json(self.data_dict)