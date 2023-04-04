#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name: dFBA
date: 02.04.23
author: Snorre Sulheim

Holds the class used to run dFBA in python. With leakage. Code copied from outcompete prioject
"""

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from pathlib import Path
import scipy.stats as st
import seaborn as sns
import cobra
import leakage

from scipy.integrate import solve_ivp
from collections import defaultdict
import utils
import time
from itertools import chain
from optlang.symbolics import Zero
from random import shuffle
from collections.abc import MutableMapping

# import traceback

# Default values
DEFAULT_VALUES = {
 'Km': 10, # mmol/L (COMETS default is 0.01 mmol/cm^3)
 'Vmax': 20, #mmol/gDW/h
 'hill': 1,
 'time_step': 0.1, #h
 'max_cycles': 1000,
 'volume': 1, #1 cm^3
}

MEDIUM_EPSILON = 1e-4 # mMol
DIFF_METHOD = "solve_ivp"
GROWTH_EPSILON = 1e-8
EXCHANGE_FLUX_EPSILON = 1e-6
RATIO_CONVERGENCE_LIMIT = 0.01
RATIO_CUTOFF = 0.9

class dFBA(object):
    def __init__(self, iterations = 100, dt = 1, basis_medium = "M9", folder = "./dFBA", method = "pFBA", 
                 pfba_fraction = 1.0, dilution_interval = None, dilution_factor = 0.1, run_until_convergence = True, save_files = False, discard_padding_on_save = True, store_exchanges_flag = False, 
                     which_exchanges = "all", store_exchanges_rate = 1):
       self.models = {}
       self.method = method
       self.medium = Medium(iterations, dt, basis_medium)
       self.dt = dt # Hours
       self.iterations = iterations
       self.folder = Path(folder)
       self.folder.mkdir(exist_ok = True)
       self.pfba_fraction = pfba_fraction
       self.dilution_counter = 0
       self.dilution_interval = dilution_interval
       self.dilution_factor = dilution_factor
       self.run_until_convergence = run_until_convergence
       self.store_exchanges_flag = store_exchanges_flag
       self.which_exchanges = which_exchanges
       self.store_exchanges_rate = store_exchanges_rate



       self.dilution_species_ratio = None
       self.stop_flag = False
       self.save_files = False
       self.discard_padding = True
       


    def add_model(self, model_name, cbmodel, initial_biomass, add_leakage = False, leakage_dict = None):
            model = Model(cbmodel, initial_biomass, self.iterations, self.dt, method = self.method, 
                          pfba_fraction = self.pfba_fraction, add_leakage = add_leakage, leakage_dict = leakage_dict,
                          store_exchanges_flag = self.store_exchanges_flag, which_exchanges = self.which_exchanges,
                          store_exchanges_rate = self.store_exchanges_rate, medium = self.medium)
            
            self.models[model_name] = model

    def add_models(self, model_name_dict):
        for model_name, model_info in model_name_dict.items():
            cbmodel, initial_biomass = model_info
            model = Model(cbmodel, initial_biomass, self.iterations, self.dt, method = self.method, 
                          pfba_fraction = self.pfba_fraction, medium = self.medium)
            self.models[model_name] = model


    def initiate_medium_and_constraints(self):
        self.medium.get_medium_components_from_models(self.models)
        # for model_name, model in self.models.items():
        self.medium.store_medium(0)

    # @profile
    def run(self):
        start_time = time.time()
        self.initiate_medium_and_constraints()
        model_names = list(self.models.keys())
        for i in range(1, self.iterations):
            if self.stop_flag:
                break

            for model_name in model_names:
                self.simulate_model(model_name, i)
            shuffle(model_names)
            self.medium.store_medium(i)
            if self._all_models_infeasible():
                print("Simulation ended as all models are infeasible")
                break

            # Dilution
            if self.dilution_interval:
                self.dilution_counter += self.dt
                if self.dilution_counter >= self.dilution_interval:
                    # print("Dilute!")
                    self.dilute()
                    self.dilution_counter = 0 
    
        self.number_of_iterations = i
        # self.final_ratios = self._get_biomass_ratios()

        print("dFBA simulation took {0} seconds".format(int(time.time()-start_time)))
        # Store results
        self.store_results()

    def store_results(self, save_files = None, remove_zeros = True):
        if not save_files:
            save_files = self.save_files

        # Biomass 
        biomass_df = pd.DataFrame()
        biomass_df["Timepoint"] = np.arange(0, self.iterations)
        biomass_df["Time"] = biomass_df["Timepoint"]*self.dt

        for model_name, model in self.models.items():
            biomass_df[model_name] = model.biomass_array

        if remove_zeros:
            zero_idx = (biomass_df[[m for m in self.models.keys()]]==0).all(axis = 1)
            biomass_df = biomass_df.loc[~zero_idx,:]

        self.biomass_df = biomass_df

        # Media
        self.concentrations_df = self.medium.get_concentrations_df()

        # Exchanges
        if self.store_exchanges_flag:
            self.model_exchanges = dict()
            for model_name, model in self.models.items():
                self.model_exchanges[model_name] = model.get_exchange_fluxes()
        else:
            self.model_exchanges = None


        if self.discard_padding:
            # Remove the empty rows in biomass and concentration arrays
            idx = self.biomass_df.Timepoint < self.number_of_iterations
            self.biomass_df = self.biomass_df.loc[idx,:]

        if save_files:
            biomass_fn = self.folder / "biomass.csv"
            biomass_df.to_csv(biomass_fn)
            concentrations_fn = self.folder / "media.csv"
            self.concentrations_df.to_csv(concentrations_fn)

            for model_name, df in self.model_exchanges.items():
                ex_fn = self.folder / "exchanges_{0}.csv".format(model_name)
                df.to_csv(ex_fn)

    def _all_models_infeasible(self):
        for model_name, model in self.models.items():
            if model.status == "optimal":
                return False
        return True

    def dilute(self):
        # Diluted medium
        for m_id, conc in self.medium.concentrations.items():
            if self.medium.static_concentrations[m_id]:
                # Do nothing
                continue

            try:
                initial_conc = self.medium.initial_conditions[m_id]
            except KeyError:
                initial_conc = 0

            self.medium.concentrations[m_id] = conc*self.dilution_factor + initial_conc*(1-self.dilution_factor)
        # Dilute biomass
        for model_name, model in self.models.items():
            model.current_biomass = model.current_biomass * self.dilution_factor

        # Test convergence with respect to previous dilution
        convergence = self._test_dilution_convergence()
        if convergence and self.run_until_convergence:
            self.stop_flag = True
        elif max([v for k,v in self.dilution_species_ratio.items()])>RATIO_CUTOFF:
            self.stop_flag = True

    def _test_dilution_convergence(self):
        """
        This function tests whether the strain abundances have converged
        """
        prev_ratios = self.dilution_species_ratio
        self.dilution_species_ratio = self._get_biomass_ratios()

        if not prev_ratios:
            return False
        else:
            # Compare ratios
            ratio_diffs = []
            for key, value in prev_ratios.items():
                ratio_diffs.append(abs(value-self.dilution_species_ratio[key])/value)
            if all([x<RATIO_CONVERGENCE_LIMIT for x in ratio_diffs]):
                return True
            else:
                return False


    def get_final_biomass(self):
        temp = self._get_biomass_ratios()
        biomass_ratios = {"Rel. abundance {0}".format(key): value for (key, value) in temp.items()}
        biomass_ratios["Stable"] = int(self.is_coexistence())
        for model_name, model in self.models.items():
            key = "Final biomass {0}".format(model_name)
            biomass_ratios[key] = model.current_biomass
        return biomass_ratios


    def is_coexistence(self, n = None, minimum_biomass = None):
        """
        Check if the simulation reached a stable community

        Parameters
        ----------

        n:  int, optional
            Number of species required to be stable. Default = None, correspeonds to all species.
        minimum_biomass: float, optional
            Minimum biomass required for a model to be considered to be present. Default is the inoculum size.

        """
        if not n:
            n = len(self.models)

        if minimum_biomass:
            min_biomass_dict = {model_name: minimum_biomass for (model_name, _) in self.models.items()}
        else:
            min_biomass_dict = {model_name: model.initial_biomass for (model_name, model) in self.models.items()}


        biomass_ratios = self._get_biomass_ratios()
        # biomass = {model_name: model.current_biomass for (model_name, model) in self.models.items()}
        species_present = 0
        for model_name, model in self.models.items():
            if biomass_ratios[model_name] > 1 - RATIO_CUTOFF:
                if model.current_biomass > min_biomass_dict[model_name]:
                    species_present += 1

        stable_cultivation = species_present == n
        biomass_ratios["Stable"] = int(stable_cultivation)
        return stable_cultivation

                



    def _get_biomass_ratios(self):
        biomass_ratios = {}
        total_biomass = 0
        for model_name, model in self.models.items():
            biomass_ratios[model_name] = model.current_biomass
            total_biomass += model.current_biomass

        biomass_ratios = {key:value/total_biomass for (key, value) in biomass_ratios.items()}
        return biomass_ratios

    # @profile
    def simulate_model(self, model_name, i):
        model = self.models[model_name]
        #self._update_constraints(model)
        model.update_bounds_from_medium(self.medium)
        try:
            fluxes = solve_model(model, reactions = model.exchange_reactions + [model.objective])
        except cobra.exceptions.Infeasible:
            model.status = "infeasible"
        else:
            if fluxes[model.objective_id] < GROWTH_EPSILON:
                model.status == "no growth"
            elif model.cbmodel.solver.status != "optimal":
                model.status = "infeasible"
            else:
                model.status = "optimal"

        if model.status in ["optimal", "no growth"]:
            if model.status == "optimal":
                model.update_biomass(fluxes[model.objective_id], i)
            else:
                model.biomass_array[i] = model.current_biomass
            if self.store_exchanges_flag and (i%self.store_exchanges_rate == 0):
                model.store_exchanges(fluxes, i)
            self.medium.update_medium(fluxes, model, model.current_biomass)
        else:
            print("Model {0} is not feasible at timepoint {1}".format(model_name, i))



    # def _update_constraints(self, model):
    #     for m_id, conc in self.medium.concentrations.items():
    #         r_ex = model.met_ex_dict[m_id]
    #         lb_monod = monod(conc, model.km_dict[r_ex], model.Vmax_dict[r_ex])
    #         r = model.cbmodel.reactions.get_by_id(r_ex)
    #         r.lower_bound = -lb_monod

class Medium(object):
    def __init__(self, iterations, dt, basis_medium = "M9", store_concentrations = "all"):
        self.concentrations = defaultdict(float)
        self.static_concentrations = defaultdict(lambda: False)
        self.basis_medium = basis_medium
        self.iterations = iterations
        self.dt = dt
        self.store_concentrations = store_concentrations
        self.set_basis_medium()


    def set_basis_medium(self):
        if self.basis_medium == "M9":
            # utils.set_FBA_M9_basis(model.cbmodel)
            self.set_M9_basis_medium()
        else:
            raise NotImplementedError

    def set_store_concentrations(self, store_concentrations):
        self.store_concentrations = store_concentrations

    def get_medium_components_from_models(self, models):
        all_mets = []
        for _, model in models.items():
            ex_mets = list(model.met_ex_dict.keys())
            all_mets += ex_mets
        all_mets = list(set(all_mets))
        all_mets.sort()
        self.metabolites = all_mets
        self.met_conc_array = np.zeros((len(all_mets), self.iterations))

    def define_static_M9_medium(self):
        pass


    def define_initial_conditions(self, met_conc_dict):
        self.initial_conditions = met_conc_dict
        for m_id, conc in met_conc_dict.items():
            self.concentrations[m_id] = conc

    # @profile
    def update_medium(self, fluxes, model, model_biomass):
        # fluxes = fluxes.to_dict()
        mu = fluxes[model.objective_id]
        for r_id, flux in fluxes.items():
            if (r_id == model.objective_id):
                continue
            print(model.ex_met_dict)
            try:
                m_id = model.ex_met_dict[r_id]
            except KeyError:
                print('No {0}'.format(m_id))
                continue
            if flux!= 0:
                if self.static_concentrations[m_id]:
                    # Don't update medium
                    continue
                if DIFF_METHOD == "Euler":
                    delta_concentration = flux*model_biomass*self.dt
                    self.concentrations[m_id] = max(0, delta_concentration + self.concentrations[m_id])
                else:
                    dS_dt_fun = lambda t, S: model_biomass*np.exp(mu*t)*flux
                    sol = solve_ivp(dS_dt_fun, [0, self.dt], [0], rtol = 1e-5)
                    self.concentrations[m_id] = max(0, sol.y[0][-1] + self.concentrations[m_id])
    # @profile
    def store_medium(self, timestep):
        self.met_conc_array[:, timestep] = [self.concentrations[m_id] for m_id in self.metabolites]

    def get_concentrations_df(self, remove_zeros = True):
        df = pd.DataFrame(self.met_conc_array.T, columns = self.metabolites)
        df["Timepoint"] = np.arange(0, self.iterations)
        df["Time"] = df["Timepoint"]*self.dt

        if remove_zeros:
            df = df.loc[~(df.loc[:, self.metabolites]==0).all(axis =1),:]
            df = df.loc[:, ~(df == 0).all(axis = 0)]
        return df

    def set_M9_basis_medium(self):
        data_folder = utils.get_data_folder()
        m9_basis_fn = data_folder / "medium"/ "M9_basis_aerobe.csv"
        M9_df = pd.read_csv(m9_basis_fn, header = 0)
        self.basis_mets = list(M9_df["ID"])
        for _, row in M9_df.iterrows():
            m_id = row["ID"]
            self.concentrations[m_id] = 1000
            self.static_concentrations[m_id] = True



class Model(object):
    """
    Notes: 
     - monod equations are currently acting on exchange reactions, but should rather act on transporters
    """
    def __init__(self, cbmodel, initial_biomass, iterations, dt, method = "pFBA", pfba_fraction = 1.0,
                 add_leakage = False, leakage_dict = None, store_exchanges_flag = True, which_exchanges = "all",
                 store_exchanges_rate = 1, medium = None, unconstrain_basis_mets = True):
        self.cbmodel = cbmodel
        self.method = method
        self.pfba_fraction = pfba_fraction
        self.km_dict = defaultdict(lambda: DEFAULT_VALUES['Km'])
        self.Vmax_dict = defaultdict(lambda: DEFAULT_VALUES['Vmax'])
        self.hill_dict = defaultdict(lambda: DEFAULT_VALUES['hill'])
        self.leakage_dict = leakage_dict
        self.add_leakage = add_leakage
        self._lb = {}
        # self.store_original_lower_bounds()
        self.current_biomass = initial_biomass
        self.initial_biomass = initial_biomass
        self.objective_id = utils.find_objective(cbmodel)
        self.objective = self.cbmodel.reactions.get_by_id(self.objective_id)
        self.biomass_array = np.zeros(iterations)*np.nan
        self.biomass_array[0] = initial_biomass
        self.iterations = iterations
        self.status = 'optimal'
        self.medium = medium
        self.unconstrain_basis_mets = unconstrain_basis_mets
        self.dt = dt
        self.store_exchanges_flag = store_exchanges_flag
        self.store_exchanges_rate = store_exchanges_rate
        self.which_exchanges = which_exchanges
        self.prep_model()

    def prep_model(self):
        # First, close all uptake and open all secretion
        for r in self.cbmodel.exchanges:
            r.lower_bound = 0
            r.upper_bound = 1000

        # Add pFBA objective if chosen
        if self.method == "pFBA":
            self.pFBA_model = make_pfba_model(self.cbmodel, self.objective_id)
        elif self.method == "random_exchange":
            self.random_model, self.random_objective = make_random_model(self.cbmodel, self.objective_id)
            print(self.cbmodel.id, self.random_objective)
        elif self.method == 'FBA_with_leakage':
            self.cbmodel = make_leaky_model(self.cbmodel, self.objective_id)

        # if self.add_leakage:
        #     if self.leakage_dict is None:
        #         self.leakage_dict = DEFAULT_leakage_DICT
        #     self.leakage_dict["carbon_sources"] = get_leakage_carbon_sources(self.leakage_dict["carbon_sources"], 
        #                                                                  self.cbmodel)
        #     cd = self.leakage_dict
        #     self.cbmodel, self.leakage_constraint_dict = leakage.add_leakage_constraint(self.cbmodel, self.objective_id, 
        #                                     cd['carbon_sources'], cd['wc'], cd['we'], cd['wr'], cd['pmax'], 
        #                                     randomize_we = cd["randomize_we"])
        

        # Map dictionary mapping metabolites to reactions
        self.make_metabolite_reaction_mapping()

        # make list of exchange reactions
        self.exchange_reactions = self.cbmodel.exchanges
        self.exchange_reaction_ids = [r.id for r in self.exchange_reactions]
        self.exchange_metabolites = [self.ex_met_dict[r_id] for r_id in self.exchange_reaction_ids]

        # Make array for exchanges
        if self.store_exchanges_flag:
            exch_array_size = (self.iterations // self.store_exchanges_rate)
            if isinstance(self.which_exchanges, list):
                self.which_exchange_reactions = []
                for m_id in self.which_exchanges:
                    try:
                        r_id = self.met_ex_dict[m_id]
                    except KeyError:
                        pass
                    else:
                        self.which_exchange_reactions.append(r_id)
                
            # elif isinstance(self.which_exchanges, (bool, type(None))) and not self.which_exchanges:
            #     print("Selected no exchanges to store")
            #     self.store_exchanges_flag = False

            elif self.which_exchanges.lower() == "all":
                self.which_exchange_reactions = []
                if self.medium is not None:
                    remove_mets = [x for x in self.medium.basis_mets if x not in ["o2_e", 'pi_e', 'nh4_e']]
                else:
                    remove_mets = []
                for r_id in self.exchange_reaction_ids:
                    m_id = self.ex_met_dict[r_id]
                    if not m_id in remove_mets:
                        self.which_exchange_reactions.append(r_id)


                # Remove M9 default compounds except o2, nh4_e and pi_e
                self.exchange_flux_array = np.zeros((exch_array_size, len(self.which_exchange_reactions)))*np.nan
            else:
                print("Could interpret {0} to store exchanges")
                raise ValueError
            self.exchange_flux_array = np.zeros((exch_array_size, len(self.which_exchange_reactions)))*np.nan
            self.exchange_flux_timepoints = np.zeros(exch_array_size)*np.nan

        # Set very high vmax for basic mets in the medium to avoid that those (in particular oxygen) constrains the model
        if self.unconstrain_basis_mets:
            if self.medium is not None:
                for m_id in self.medium.basis_mets:
                    self.set_Vmax(m_id, 1000)
            else:
                print("Model has no medium")
                raise ValueError

    # def store_original_lower_bounds(self):
    #     for r in self.cbmodel.exchanges:
    #         self._lb[r.id] = r.lower_bound

    # @profile
    def update_biomass(self, growth_rate, i):
        #print(solution.fluxes)
        #print(self.objective)

        if DIFF_METHOD == "Euler":
            self.current_biomass = self.current_biomass + growth_rate*self.current_biomass*self.dt
            self.biomass_array[i] = self.current_biomass
        else:   
            mu = growth_rate
            dX_dt_fun = lambda t, X: mu*X[0]
            sol = solve_ivp(dX_dt_fun, [0, self.dt], [self.current_biomass], rtol = 1e-5)
            self.current_biomass = sol.y[0][-1]
            self.biomass_array[i] = sol.y[0][-1]


    

    def set_km(self, m_id, km):
        """
        Set km for a specific metabolite in mM
        """
        try:
            r_id = self.met_ex_dict[m_id]
        except KeyError:
            pass
            #print("Can't set Km for {0}; not in model {1}".format(m_id, self.cbmodel.id))
        else:
            self.km_dict[r_id] = km

    def set_Vmax(self, m_id, Vmax):
        """
        Set Vmax in mmol/gDW/h
        """
        try:
            r_id = self.met_ex_dict[m_id]
        except KeyError:
            pass
            #print("Can't set Vmax for {0}; not in model {1}".format(m_id, self.cbmodel.id))
        else:
            self.Vmax_dict[r_id] = Vmax
    
    def set_hill(self, m_id, hill):
        try:
            r_id = self.met_ex_dict[m_id]
        except KeyError:
            pass
            #print("Can't set hill coefficient for {0}; not in model {1}".format(m_id, self.cbmodel.id))
        else:
            self.hill_dict[r_id] = hill

    # @property
    # def exchange_reactions(self):
    #     return list(self.ex_met_dict.keys())
    
    
    def make_metabolite_reaction_mapping(self):
        self.met_ex_dict = {}
        self.ex_met_dict = {}
        for r in self.cbmodel.exchanges:
            (m, i) = r.metabolites.popitem()
            self.met_ex_dict[m.id] = r.id
            self.ex_met_dict[r.id] = m.id
    
    # @profile
    def update_bounds_from_medium(self, medium):
        if self.add_leakage:
            # update_mets = []
            for m_id, x in medium.concentrations.items():
                try:
                    r_id = self.met_ex_dict[m_id]
                except KeyError:
                    pass
                else:
                    r = self.cbmodel.reactions.get_by_id(r_id)
                    if m_id in self.leakage_dict["carbon_sources"]:
                        if x < MEDIUM_EPSILON:
                            r.lower_bound = 0
                            self.leakage_dict['wc'][m_id] = 1000
                        else:  
                            r.lower_bound = -1000
                            self.leakage_dict['wc'][m_id] = self.leakage_dict['wc0'][m_id]*(1 + (self.leakage_dict['Km']/x))
                    else:
                        if np.isinf(x):
                            r.lower_bound = -1000
                        else:
                            v_max = monod(x, self.km_dict[r_id], self.Vmax_dict[r_id])
                            r.lower_bound = -v_max

                # print(self.leakage_dict['wc'])
            self.leakage_constraint_dict = leakage.update_carbon_constraints(self.cbmodel, self.leakage_constraint_dict, self.leakage_dict["wc"], self.leakage_dict["carbon_sources"])
            # print(self.leakage_dict["wc"], self.leakage_dict["wc0"])
            self.cbmodel.constraints.leakage.set_linear_coefficients(coefficients = self.leakage_constraint_dict)
            

        else:

            for m_id, x in medium.concentrations.items():
                try:
                    r_id = self.met_ex_dict[m_id]
                except KeyError:
                    pass
                else:

                    r = self.cbmodel.reactions.get_by_id(r_id)
                    if np.isinf(x):
                        r.lower_bound = -1000
                    else:
                        v_max = monod(x, self.km_dict[r_id], self.Vmax_dict[r_id])
                        r.lower_bound = -v_max

    def store_exchanges(self, fluxes, i):
        if self.store_exchanges_flag:
            if i%self.store_exchanges_rate == 0:
                j = i//self.store_exchanges_rate
                self.exchange_flux_array[j, :] = [fluxes[r_id] for r_id in self.which_exchange_reactions]
                self.exchange_flux_timepoints[j]=i

    def get_exchange_fluxes(self):
        all_nan_timepoints = np.isnan(self.exchange_flux_array).all(axis = 1)
        header = [self.ex_met_dict[r_id] for r_id in self.which_exchange_reactions]
        df = pd.DataFrame(self.exchange_flux_array[~all_nan_timepoints, :], columns = header)
        df["Timepoint"] = self.exchange_flux_timepoints[~all_nan_timepoints]
        df.fillna(0, inplace = True)
        
        if self.which_exchanges == "all":
            all_zero_mets = (df.abs() < EXCHANGE_FLUX_EPSILON).all(axis = 0)
            df = df.loc[:, ~all_zero_mets]

        return df  


def get_leakage_carbon_sources(carbon_sources, cbmodel):
    if isinstance(carbon_sources, list):
        if np.all([c[-2:]=="_e" for c in carbon_sources]):
            return carbon_sources
    elif isinstance(carbon_sources, str):
        if carbon_sources[-2:] == "_e":
            return [carbon_sources]
        elif carbon_sources == "all":
            # All exchanges except those in M9 medium
            carbon_sources = []
            m9_mets = utils.get_m9_metabolite_ids()
            for r in cbmodel.exchanges:
                m, _ = r.metabolites.popitem()
                if not m.id in m9_mets:
                    carbon_sources.append(m.id)
            return carbon_sources
    # Something is wrong
    raise ValueError



def solve_model(model, reactions):
    max_growth_rate = model.cbmodel.slim_optimize()
    # print(model.cbmodel.summary())
    if model.method == "FBA":
        fluxes = get_fluxes(model.cbmodel, reactions)
    elif model.method == 'FBA_with_leakage':
        fluxes = _FBA_with_leakage(model.cbmodel)
    elif model.method == "pFBA":
        if not np.isnan(max_growth_rate):
            for r in model.cbmodel.exchanges:
                model.pFBA_model.reactions.get_by_id(r.id).bounds = r.bounds
            model.pFBA_model.reactions.get_by_id(model.objective_id).lower_bound = model.pfba_fraction * max_growth_rate
            model.pFBA_model.slim_optimize()
            fluxes = get_fluxes(model.pFBA_model, reactions)
        else:
            fluxes = get_fluxes(model.cbmodel, reactions)
    elif model.method == "random_exchange":
        if not np.isnan(max_growth_rate):
            for r in model.cbmodel.exchanges:
                model.random_model.reactions.get_by_id(r.id).bounds = r.bounds
            model.random_model.reactions.get_by_id(model.objective_id).lower_bound = model.pfba_fraction * max_growth_rate
            model.random_model.slim_optimize()
            fluxes = get_fluxes(model.random_model, reactions)
        else:
            fluxes = get_fluxes(model.cbmodel, reactions)
    else:
        raise NotImplementedError
    return fluxes

def make_pfba_model(model, objective_id):
    """
    It seems like it is much faster to pre-make a separate pFBA model rather than running pFBA at each timepoint
    """
    # Copy model
    pfba_model = model.copy()

    # estimate max_growth
    # max_growth = pfba_model.slim_optimize()
    # print(max_growth)
    # r_obj = pfba_model.reactions.get_by_id(objective_id)
    # print(r_obj)
    # pfba_model.reactions.get_by_id(objective_id).lower_bound = 0
    pfba_model.reactions.get_by_id(objective_id).objective_coefficient = 0
    reaction_variables = (
        (rxn.forward_variable, rxn.reverse_variable) for rxn in pfba_model.reactions
    )
    variables = chain(*reaction_variables)
    pfba_model.objective = pfba_model.problem.Objective(
        Zero, direction="min", sloppy=True, name="_pfba_objective")
    pfba_model.objective.set_linear_coefficients({v: 1.0 for v in variables})
    return pfba_model

def make_leaky_model(model, objective_id):
    leaky_model = model.copy()
    leak_reactions = []
    # leak_mets
    for m in leaky_model.metabolites:
        if (m.compartment == 'c') and (m.id[-2:]=='_c'):
            try:
                m_e = leaky_model.metabolites.get_by_id(m.id.replace('_c','_e'))
            except KeyError:
                # m_e = cobra.Metabolite(m.id.replace('_c','_e'))
                m_e = m.copy()
                m_e.id = m.id.replace('_c','_e')
                m_e.compartment = 'e'
                r_ex = leaky_model.add_boundary(m_e, type='exchange')
            else:
                r_ex = leaky_model.reactions.get_by_id('EX_{0}'.format(m_e))
            r_ex.upper_bound = 1000

            # Add leak reactions
            r = cobra.Reaction('LEAK_{0}'.format(m.id))
            r.add_metabolites({m:-1, m_e:1})
            r.bounds = (0, 0)
            leak_reactions.append(r)
            # leak_mets.append(m)
    leaky_model.add_reactions(leak_reactions)
    return leaky_model

def make_random_model(model, objective_id):
    # Copy model
    random_model = model.copy()
    random_objective = np.random.choice(model.exchanges)
    random_model.reactions.get_by_id(objective_id).objective_coefficient = 0
    random_model.reactions.get_by_id(random_objective.id).objective_coefficient = 1
    return random_model, random_objective.id

def get_fluxes(model, reactions):
    """
    Generates fast solution representation of the current solver state.
    From https://github.com/opencobra/cobrapy/issues/477#issuecomment-292715696
    """
    fluxes = {}
    if model.solver.status == 'optimal':
        var_primals = dict(zip(model.solver._get_variables_names(), model.solver._get_primal_values()))
        for rxn in reactions:
            fluxes[rxn.id] = var_primals[rxn.id] - var_primals[rxn.reverse_id]
    else:
        for rxn in reactions:
            fluxes[rxn.id] = 0
    return fluxes        


def monod(X, Km, Vmax, hill = 1):
    return Vmax*(X**hill/(Km+X**hill))


class ModelDict(MutableMapping):
    """A dictionary that applies an arbitrary key-altering
       function before accessing the keys"""

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        return self.store[key]

    def __setitem__(self, key, value):
        self.store[key] = value

    def __delitem__(self, key):
        del self.store[key]

    def __iter__(self):
        return iter(self.store)
    
    def __len__(self):
        return len(self.store)

    def __enter__(self):
        for model_name, model in self.store.items():
            model.__enter__()
        return self
    def __exit__(self, exc_type, exc_value, tb):
        if exc_type is not None:
            traceback.print_exception(exc_type, exc_value, tb)
        for model_name, model in self.store.items():
            model.__exit__(exc_type, exc_value, tb)



def _FBA_with_leakage(model, objective = 'BIOMASS_Ec_iJO1366_core_53p95M', ratio_of_optimum = 0.95, 
                     max_shadow_price = 1, min_shadow_price = 1e-6, slope = -1.38):
    # Simulate model
    solution = model.optimize()

    # Create a temporary instance of the model
    turnover_mets = get_turnover_mets(model, solution)
    shadow_prices = leakage.estimate_shadow_prices(model)
    with model as M:
        # Create a variable that scales the leakage
        scale = M.problem.Variable('leakage_scale')
        M.add_cons_vars([scale])
        constraints = []
    
        # # Add new reactions and mets for leakage
        # leak_reactions = []
        # leak_mets = []
        for i, m_id in enumerate(turnover_mets):
            # Get shadow price
            sp = -1*shadow_prices[m_id]
            if (sp > max_shadow_price) or np.isnan(sp) or (sp < min_shadow_price):
                # print(m_id, sp)
                continue
            m = M.metabolites.get_by_id(m_id)
            # Check if demand reaction exist
            # for r in m.reactions:
            #     if r.boundary:
            #         print(r)
            # Add excahnge
            # try:
            #     m_e = M.metabolites.get_by_id(m_id.replace('_c','_e'))
            # except KeyError:
            #     # m_e = cobra.Metabolite(m_id.replace('_c','_e'))
            #     m_e = m.copy()
            #     m_e.id = m_id.replace('_c','_e')
            #     m_e.compartment = 'e'
            #     M.add_boundary(m_e, type='exchange')
            # else:
            #     r_ex = model.reactions.get_by_id('EX_{0}'.format(m_e))
            #     r_ex.upper_bound = 1000
                
            # Add leak reaction
            r = M.reactions.get_by_id('LEAK_{0}'.format(m_id))
            r.bounds = (0, 1000)
            # if i > 90:
            #     break

            sp = shadow_prices[m.id]
            predicted_leakage = np.abs(sp)**(slope)
            constraint = M.problem.Constraint(
                r.flux_expression - scale*predicted_leakage,
                lb=0,
                ub=1)
            constraints.append(constraint)
            # print(m, r, sp)
        
        M.reactions.get_by_id(objective).objective_coefficient = 0
        M.reactions.get_by_id(objective).lower_bound = solution.objective_value*0.95
        M.add_cons_vars(constraints)
        obj = M.problem.Objective(M.solver.variables.leakage_scale, direction = 'max')
        M.objective = obj

        # Consider to don't get all fluxes, use get_fluxes
        new_solution = M.optimize(objective_sense=None)
        # print(M.summary())
        # print(new_solution.fluxes['BIOMASS_Ec_iJO1366_core_53p95M'])
    return new_solution.fluxes
    

def FBA_with_leakage(model, objective = 'BIOMASS_Ec_iJO1366_core_53p95M', ratio_of_optimum = 0.95, 
                     max_shadow_price = 1, min_shadow_price = 1e-6, slope = -1.38):
    # Simulate model
    solution = model.optimize()

    # Create a temporary instance of the model
    turnover_mets = get_turnover_mets(model, solution)
    shadow_prices = leakage.estimate_shadow_prices(model)
    with model as M:
        # Create a variable that scales the leakage
        scale = M.problem.Variable('leakage_scale')
        M.add_cons_vars([scale])
    
        # Add new reactions and mets for leakage
        leak_reactions = []
        leak_mets = []
        for i, m_id in enumerate(turnover_mets):
            # Get shadow price
            sp = -1*shadow_prices[m_id]
            if (sp > max_shadow_price) or np.isnan(sp) or (sp < min_shadow_price):
                # print(m_id, sp)
                continue
            m = M.metabolites.get_by_id(m_id)
            # Check if demand reaction exist
            # for r in m.reactions:
            #     if r.boundary:
            #         print(r)
            # Add excahnge
            try:
                m_e = M.metabolites.get_by_id(m_id.replace('_c','_e'))
            except KeyError:
                # m_e = cobra.Metabolite(m_id.replace('_c','_e'))
                m_e = m.copy()
                m_e.id = m_id.replace('_c','_e')
                m_e.compartment = 'e'
                M.add_boundary(m_e, type='exchange')
            else:
                r_ex = model.reactions.get_by_id('EX_{0}'.format(m_e))
                r_ex.upper_bound = 1000
                
            # Add leak reaction
            r = cobra.Reaction('LEAK_{0}'.format(m_id))
            r.add_metabolites({m:-1, m_e:1})
            r.bounds = (0, 1000)
            leak_reactions.append(r)
            leak_mets.append(m)
            # if i > 90:
            #     break
        M.add_reactions(leak_reactions)

        # Now add constraints on the leakage
        constraints = []
        for m, r in zip(leak_mets, leak_reactions):
            sp = shadow_prices[m.id]
            predicted_leakage = np.abs(sp)**(slope)
            constraint = M.problem.Constraint(
                r.flux_expression - scale*predicted_leakage,
                lb=0,
                ub=1)
            constraints.append(constraint)
            # print(m, r, sp)
        
        M.reactions.get_by_id(objective).objective_coefficient = 0
        M.reactions.get_by_id(objective).lower_bound = solution.objective_value*0.95
        M.add_cons_vars(constraints)
        obj = M.problem.Objective(M.solver.variables.leakage_scale, direction = 'max')
        M.objective = obj

        # Consider to don't get all fluxes, use get_fluxes
        new_solution = M.optimize(objective_sense=None)
        # print(M.summary())
        # print(new_solution.fluxes['BIOMASS_Ec_iJO1366_core_53p95M'])
    return new_solution.fluxes

def get_turnover_mets(model, solution, flux_epsilon = 1e-4):
    intracellular_metabolites = [m.id for m in model.metabolites if m.compartment =='c']
    S_matrix = cobra.util.create_stoichiometric_matrix(model, array_type='DataFrame')
    reactions_flux_bool = solution.fluxes.abs() > flux_epsilon
    
    # Remove non-intracellular metabolites and zero-flux reactions
    S_reduced = S_matrix.loc[intracellular_metabolites, reactions_flux_bool]

    # Remove metabolites that have all zero rows in the S_non_zero matrix
    nonzero_mets = (S_reduced != 0).any(axis = 1)
    nonzero_rxns = (S_reduced != 0).any(axis = 0)
    S = S_reduced.loc[nonzero_mets, nonzero_rxns]

    # Remove prot pool
    S = S.loc[S.index!='prot_pool',:]

    # Remove ACP metabolites (assumed to be to large since they are connected to a protein)
    S = S.loc[~S.index.str.contains('ACP'), :]

    # Remove mets that are already excreted
    already_exchanged_mets = get_already_exchanged_mets(model, solution)

    S = S.loc[~S.index.isin(already_exchanged_mets), :]
    # print(S.index, already_exchanged_mets)
    return S.index



def get_already_exchanged_mets(model, solution):
    # List already excreted metabolites
    exchanged_mets = []
    for r in model.boundary:
        flux = solution.fluxes[r.id]
        if flux != 0:
            exchanged_mets.append(list(r.metabolites.keys())[0].id)
    exchanged_mets = [m_id.replace('_e', '_c') for m_id in exchanged_mets]
    return exchanged_mets



if __name__ == '__main__':
    model = cobra.io.read_sbml_model('../models/e_coli/momentiJO1366.xml')
    model.solver = 'gurobi'
    # Initial conditions is 0.013 gDW/L (in total)
    model_name = 'Ecoli'
    model_name_dict = {model_name: [model, 0.006]}

    glucose_mM = utils.convert_gL_to_mM("C6H12O6", 4)
    D = dFBA(iterations = 3, dt = 0.1, method = "FBA_with_leakage", store_exchanges_flag = False)#(, pfba_fraction = 0.95)
    D.add_models(model_name_dict)
    # Set Km and vMax
    D.models[model_name].set_km("glc__D_e", 1)
    D.models[model_name].set_Vmax("glc__D_e", 10)



    D.medium.define_initial_conditions({"glc__D_e": glucose_mM})
    # D.medium.set_store_concentrations(["glc__D_e", "nh3_e"])
    D.run()
