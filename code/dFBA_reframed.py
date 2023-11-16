#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name: dFBA with reframed
date: 15.11.23
author: Snorre Sulheim

Holds the class used to run dFBA in python. With leakage. Code copied from outcompete project, but translated to 
use reframed library instead of cobrapy
"""
import reframed
import time
import numpy as np
from reframed.solvers.solution import Status
import utils
from collections import defaultdict
import pandas as pd
from pathlib import Path
from scipy.integrate import solve_ivp
from random import shuffle
import logging
import dotenv


REPO_PATH =  Path(dotenv.find_dotenv()).parent
TMP_FOLDER = REPO_PATH / 'tmp'
TMP_FOLDER.mkdir(exist_ok=True)
timestr = time.strftime("%Y%m%d_%H%M")
LOGFILE = TMP_FOLDER / f'dFBA_{timestr}.log'

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s %(levelname)s:%(name)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                     handlers=[logging.FileHandler(str(LOGFILE)),logging.StreamHandler()])
logging.captureWarnings(True)



# Default values
DEFAULT_VALUES = {
 'Km': 1, # mmol/L (COMETS default is 0.01 mmol/cm^3)
 'Vmax': 10, #mmol/gDW/h
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
MAX_UPTAKE = -50


# EXCLUDE_LEAKAGE = ['accoa', 'malcoa', 'lipidX', 'nadh', 'nadph', 'nadp', 'nad']

class dFBA(object):
    def __init__(self, iterations = 100, dt = 1, basis_medium = "M9", folder = "./dFBA", method = "pFBA", 
                 pfba_fraction = 1.0, dilution_interval = None, dilution_factor = 0.1, run_until_convergence = True, save_files = False, discard_padding_on_save = True, store_exchanges_flag = False, 
                     which_exchanges = "all", store_exchanges_rate = 1, leakage_params = {'slope': -3, 'intercept': -4}):
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
       self.leakage_params = leakage_params
       self.logger = logging.getLogger('dFBA')



       self.dilution_species_ratio = None
       self.stop_flag = False
       self.save_files = False
       self.discard_padding = True


    def add_model(self, model_name, cbmodel, initial_biomass):
            model = Model(cbmodel, initial_biomass, self.iterations, self.dt, method = self.method, 
                          pfba_fraction = self.pfba_fraction,
                          store_exchanges_flag = self.store_exchanges_flag, which_exchanges = self.which_exchanges,
                          store_exchanges_rate = self.store_exchanges_rate, medium = self.medium,
                          leakage_slope = self.leakage_slope)
            
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
            # print(i)
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

    # # @profile
    def simulate_model(self, model_name, i):
        model = self.models[model_name]
        #self._update_constraints(model)
        model.update_bounds_from_medium(self.medium)
        # try:
        fluxes = solve_model(model, reactions = model.exchange_reaction_ids + [model.objective_id], constraints = model.medium_constraints)
        # except cobra.exceptions.Infeasible:
        #     model.status = "infeasible"
        # else:
        if fluxes[model.objective_id] < GROWTH_EPSILON:
            model.status == "no growth"
        elif model.solution.status != Status.OPTIMAL:
            model.status = "infeasible"
        else:
            model.status = "optimal"
        # print(model.status, fluxes[model.objective_id])
        if model.status in ["optimal", "no growth"]:
            if model.status == "optimal":
                model.update_biomass(fluxes[model.objective_id], i)
            else:
                model.biomass_array[i] = model.current_biomass
            if self.store_exchanges_flag and (i%self.store_exchanges_rate == 0):
                model.store_exchanges(fluxes, i)
            self.medium.update_medium(fluxes, model, model.current_biomass)
        else:
            logger.info("Model {0} is not feasible at timepoint {1}".format(model_name, i))



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

    # # @profile
    def update_medium(self, fluxes, model, model_biomass):
        # fluxes = fluxes.to_dict()
        mu = fluxes[model.objective_id]
        for r_id, flux in fluxes.items():
            if (r_id == model.objective_id):
                continue
            # print(model.ex_met_dict)
            try:
                m_id = model.ex_met_dict[r_id]
            except KeyError:
                print('No met for {0}'.format(r_id))
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
    # # @profile
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
            m_id = 'M_'+row["ID"]
            self.concentrations[m_id] = 1000
            self.static_concentrations[m_id] = True



class Model(object):
    """
    Notes: 
     - monod equations are currently acting on exchange reactions, but should rather act on transporters
    """
    def __init__(self, cbmodel, initial_biomass, iterations, dt, method = "pFBA", pfba_fraction = 1.0,
                  store_exchanges_flag = True, which_exchanges = "all",
                 store_exchanges_rate = 1, medium = None, unconstrain_basis_mets = True, leakage_slope = -1.38):
        self.cbmodel = cbmodel
        self.method = method
        self.pfba_fraction = pfba_fraction
        self.km_dict = defaultdict(lambda: DEFAULT_VALUES['Km'])
        self.Vmax_dict = defaultdict(lambda: DEFAULT_VALUES['Vmax'])
        self.hill_dict = defaultdict(lambda: DEFAULT_VALUES['hill'])
        self._lb = {}
        # self.store_original_lower_bounds()
        self.current_biomass = initial_biomass
        self.initial_biomass = initial_biomass
        self.objective_id = model.biomass_reaction
        self.objective = self.cbmodel.reactions[self.objective_id]
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
        self.leakage_slope  = leakage_slope
        self.uptake_constraints = {}
        self.prep_model(medium)

    def prep_model(self, medium):
        # First, close all uptake and open all secretion
        for r_id in self.cbmodel.get_exchange_reactions():
            self.cbmodel.set_flux_bounds(r_id, 0, 1000)

        # Map dictionary mapping metabolites to reactions
        self.make_metabolite_reaction_mapping()

        # Get medium constraints
        if medium:
            self.update_bounds_from_medium(medium)
            
        # Add pFBA objective if chosen
        if self.method == "pFBA":
            pfba_solution, self.pFBA_solver, self.pFBA_constraints = make_pfba_model(self.cbmodel, self.objective_id, constraints= self.medium_constraints)
        
        # elif self.method == 'FBA_with_leakage':
        #     print("##### Leakage")
        #     self.cbmodel, self.leaky_reactions, self.leaky_mets, new_exchanges = make_leaky_model(self.cbmodel, self.objective_id)

        

        # make list of exchange reactions
        self.exchange_reaction_ids = self.cbmodel.get_exchange_reactions()
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

    # # @profile
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
        for r_id in self.cbmodel.get_exchange_reactions():
            r = self.cbmodel.reactions[r_id]
            m_id  = next(iter(r.stoichiometry))
            self.met_ex_dict[m_id] = r_id
            self.ex_met_dict[r_id] = m_id
    
    # @profile
    def update_bounds_from_medium(self, medium):
        self.medium_constraints = {}
        for m_id, x in medium.concentrations.items():
            try:
                r_id = self.met_ex_dict[m_id]
            except KeyError:
                logging.debug(f'Could not find {m_id}')
            else:
                r = self.cbmodel.reactions[r_id]
                if np.isinf(x):
                    self.medium_constraints[r_id] = (MAX_UPTAKE, r.ub)
                else:
                    v_max = monod(x, self.km_dict[r_id], self.Vmax_dict[r_id])
                    self.medium_constraints[r_id] = (-v_max, r.ub)

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


def make_pfba_model(model, objective_id, fraction_of_optimum = 1, constraints = None):
    """
    It seems like it is much faster to pre-make a separate pFBA model rather than running pFBA at each timepoint
    """
    solver = solver_instance(model)
    objective = {objective_id: 1}

    if isinstance(constraints, dict):
        pre_solution = FBA(model, objective, minimize = True, constraints = constraints, solver=solver)
        if pre_solution.status != Status.OPTIMAL:
            return pre_solution

        if (fraction_of_optimum is None) or (fraction_of_optimum == 1):
            solver.add_constraint('obj', objective, '=', pre_solution.fobj, update=False)
        else:
            solver.add_constraint('obj', objective, '>', fraction_of_optimum * pre_solution.fobj, update=False)

    if not hasattr(solver, 'pFBA_flag'):
        solver.pFBA_flag = True
        for r_id in model.reactions:
            if model.reactions[r_id].reversible:
                pos, neg = r_id + '_p', r_id + '_n'
                solver.add_variable(pos, 0, inf, update=False)
                solver.add_variable(neg, 0, inf,  update=False)
        solver.update()
        for r_id in model.reactions:
            if model.reactions[r_id].reversible:
                pos, neg = r_id + '_p', r_id + '_n'
                solver.add_constraint('c' + pos, {r_id: -1, pos: 1}, '>', 0, update=False)
                solver.add_constraint('c' + neg, {r_id: 1, neg: 1}, '>', 0, update=False)
        solver.update()

    objective = dict()
    for r_id in reactions:
        if model.reactions[r_id].reversible:
            pos, neg = r_id + '_p', r_id + '_n'
            objective[pos] = 1
            objective[neg] = 1
        else:
            objective[r_id] = 1
    if pre_solution:
        solution = solver.solve(objective, minimize=True, constraints=constraints)
        solution.pre_solution = pre_solution
    else:
        solution = None

    return solution, solver, objective


def solve_model(model, reactions, constraints):
    fba_solution = reframed.FBA(model.cbmodel, constraints=constraints)
    if fba_solution.status == Status.OPTIMAL:
        max_growth_rate = fba_solution.values
    # print(model.cbmodel.summary())
    if model.method == "FBA":
        model.solution = fba_solution
        fluxes = get_fluxes(model, fba_solution, reactions)
    elif model.method == 'FBA_with_leakage':
        if not np.isnan(max_growth_rate):
            fluxes = _FBA_with_leakage(model.cbmodel, model.objective_id, reactions = reactions, ratio_of_optimum = model.pfba_fraction, slope = model.leakage_slope)
        else:
            fluxes = get_fluxes(model.cbmodel, reactions)

    elif model.method == "pFBA":
        if not np.isnan(max_growth_rate):
            for r in model.cbmodel.exchanges:
                model.pFBA_model.reactions.get_by_id(r.id).bounds = r.bounds
            model.pFBA_model.reactions.get_by_id(model.objective_id).lower_bound = model.pfba_fraction * max_growth_rate
            sol = model.pFBA_model.slim_optimize()
            fluxes = get_fluxes(model.pFBA_model, reactions)
        else:
            fluxes = get_fluxes(model.cbmodel, reactions)
    else:
        raise NotImplementedError
    return fluxes

def get_fluxes(model, solution, reactions):
    """
    Generates fast solution representation of the current solver state.
    From https://github.com/opencobra/cobrapy/issues/477#issuecomment-292715696
    """
    fluxes = {}
    if solution.status == Status.OPTIMAL:
        fluxes = {r_id: solution.values[r_id] for r_id in reactions}
    else:
        fluxes = {r_id:0 for r_id in reactions}
    return fluxes   


def FBA_with_leakage(model, constraints, fraction_of_optimum = 0.9, only_turnover_metabolites = True, min_turnover = 1e-6):
    # Make sure exchange reactions are open
    for r_id in model.get_exchange_reactions():
        model.reactions[r_id].ub = 1000

    # Get pFBA solution
    pfba_solution = reframed.pFBA(model, constraints=constraints,obj_frac=1)
    growth = pfba_solution.values[model.biomass_reaction]
    
    # Get the list of potential leaky metabolites
    # Only metabolites with a turnover in pFBA solution
    candidate_metabolites = []
    candidate_metabolites_exchanges = []
    if only_turnover_metabolites:
        for m_id, turnover in pfba_solution.get_metabolites_turnover(model).items():
            if m_id.endswith('_c') and turnover > min_turnover:
                r_ex_id = f'R_EX_{m_id[2:-2]}_e'
                if r_ex_id in model.reactions:
                    candidate_metabolites.append(m_id)
                    candidate_metabolites_exchanges.append(r_ex_id)
    else:
        for m_id in model.metabolites:
            if m_id.endswith('_c'):
                r_ex_id = f'R_EX_{m_id[:-2]}_e'
                if model.reactions.get(r_ex_id):
                    candidate_metabolites.append(m_id)
                    candidate_metabolites_exchanges.append(r_ex_id)

    # Now only consider the exchanges that can have positive fluxes
    fva_results = reframed.FVA(model, obj_frac=fraction_of_optimum, constraints=constraints, reactions = candidate_metabolites_exchanges)
    leak_mets = []
    leak_exchanges = []
    for m_id, r_ex_id in zip(candidate_metabolites, candidate_metabolites_exchanges):
        lb, ub = fva_results[r_ex_id]
        if ub > 0:
            leak_mets.append(m_id)
            leak_exchanges.append(r_ex_id)

    metabolite_values, log_leakage_rates = predict_log_leakage_rates_from_shadow_prices(model, constraints, leak_mets, leak_exchanges)

    # Now add leakage variables and constraints
    # selected_mets = ['M_pyr_c', 'M_cit_c']
    solver = reframed.solver_instance(model)
    log_vars = []
    for m_id, r_ex_id in zip(leak_mets, leak_exchanges):
        r_ex = model.reactions[r_ex_id]
        
        # This seems like too much but...
        r_ex.ub = 1000
        if r_ex.lb < 0:
            # This will be a challenge in dFBA.... but it doesn't work for the log-constraint if r_ex.lb can be < 0
            continue
        if not log_leakage_rates.get(r_ex_id):
            continue
        # log_x = problem.addVar(vtype=gp.GRB.CONTINUOUS, name = f'log_{r_ex_id}', lb=-gp.GRB.INFINITY)
        solver.add_variable(var_id= f'log_{r_ex_id}', lb=-10, ub=10)
        x = solver.problem.getVarByName(r_ex_id)
        # x.ub = 1000

        log_x = solver.problem.getVarByName(f'log_{r_ex_id}')
        log_vars.append(f'log_{r_ex_id}')

        # Constrain log_x to be np.log10(x)
        # The parameter options are important
        c = solver.problem.addGenConstrLogA(x, log_x, 10, name = f'logc_{r_ex_id}', options='FuncPieceError=1e-5 FuncPieces=-2')
        solver.constr_ids.append(f'logc_{r_ex_id}')
    solver.problem.update()
    solver.update()

    # Now predict FBA
    t0 = time.time()
    r_growth = model.biomass_reaction
    lin_obj = {}
    for log_name in log_vars:
        r_ex_id = log_name[4:]
        # This is the linear part of (x-b)^2 = a^2 - 2xb + b^2
        lin_obj[log_name] = -2*log_leakage_rates[r_ex_id]
        
    # This is the quadratic part
    quad_obj = {}
    for key, _ in lin_obj.items():
        quad_obj.update({(key,key):1})

    constraints_ = {r_growth:growth*fraction_of_optimum} # Consider to alow a range of values {r_growth: (lb, ub)}
    constraints_.update(constraints)

    # Another key parameter
    solver.problem.params.BarHomogeneous = 1
    solution = solver.solve(lin_obj, quadratic=quad_obj, minimize=True, constraints=constraints_, get_values=True)
    return solution



def predict_log_leakage_rates_from_shadow_prices(model, constraints, leak_mets, leak_exchanges, 
                                                 slope = -3, intercept = -4, min_metabolite_value = 1e-7):
    solution = reframed.FBA(model, constraints=constraints, shadow_prices=True)
    predicted_log_leakage_rates = {}
    predicted_metabolite_values = {}
    for m_id, r_ex_id in zip(leak_mets, leak_exchanges):
        metabolite_value = - solution.shadow_prices[m_id]
        if np.isfinite(metabolite_value) and metabolite_value > 1e-7:
            # Use trendline from fit 
            lograte = intercept + slope*np.log10(metabolite_value)
            predicted_log_leakage_rates[r_ex_id] = lograte
            predicted_metabolite_values[r_ex_id] = metabolite_value
    return predicted_metabolite_values, predicted_log_leakage_rates


def monod(X, Km, Vmax, hill = 1):
    return Vmax*(X**hill/(Km+X**hill))



if __name__ == '__main__':
    if 0:
        # Test FBA with leakage
        model = reframed.load_cbmodel('../models/e_coli/iJO1366.xml')
        model.solver = 'gurobi'
        for r_id in model.get_exchange_reactions():
            model.reactions[r_id].ub = 1000

        constraints = {"R_EX_glc__D_e": -10}
        solution = FBA_with_leakage(model, constraints)
        print(solution.show_values('R_EX'))
        print(solution.show_values('R_BIOMASS'))

    if 1:
        # Test dFBA
        # model = cobra.io.read_sbml_model('../models/e_coli/momentiJO1366.xml')
        model = reframed.load_cbmodel('../models/e_coli/iJO1366.xml')
    
        model.solver = 'gurobi'
        model_name = 'E_coli'
        model_name_dict = {model_name: [model, 0.006]}

        glucose_mM = utils.convert_gL_to_mM("C6H12O6", 4)
        D = dFBA(iterations = 100, dt = 0.1, method = "FBA", store_exchanges_flag = False)#(, pfba_fraction = 0.95)
        D.add_models(model_name_dict)
        # Set Km and vMax
        D.models[model_name].set_km("glc__D_e", 1)
        D.models[model_name].set_Vmax("glc__D_e", 8)



        D.medium.define_initial_conditions({"M_glc__D_e": glucose_mM})
        # D.medium.set_store_concentrations(["glc__D_e", "nh3_e"])
        D.run()
        print(D.biomass_df)



