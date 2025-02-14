import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from pathlib import Path
import scipy.stats as st
import seaborn as sns
import reframed
import cobra

from scipy.interpolate import CubicSpline, UnivariateSpline



def get_concentrations(data_folder, organism):
    exometabolites_folder = Path(data_folder)
    
    # Filenames
    fn_exometabolites = exometabolites_folder / "{0}_exometabolites.csv".format(organism)
    fn_exometabolites_std = exometabolites_folder / "{0}_exometabolites_std.csv".format(organism)
    
    # Read files as dataframes
    df_exometabolites = pd.read_csv(fn_exometabolites, index_col=0)
    df_exometabolites_std = pd.read_csv(fn_exometabolites_std, index_col=0)

    # Drop columns that ends with 'MS' (They are not reported in paczia's figures)
    drop_cols = [x for x in df_exometabolites.columns if x.endswith(' MS')]
    df_exometabolites.drop(columns = drop_cols, inplace= True)
    df_exometabolites_std.drop(columns = drop_cols, inplace= True)
    
    # met_abbreviations = [x for x in df_exometabolites.columns]#.replace(" MS", "")
    # df_exometabolites.columns = met_abbreviations
    # df_exometabolites_std.columns = met_abbreviations
    
    return df_exometabolites, df_exometabolites_std



    #
def estimate_shadow_price_for_met(model, m, solution, delta = 0.1, existing_flux = None):
    with model as M:
        try:
            r = M.reactions.get_by_id('DM_{0}'.format(m.id))
        except KeyError:
            r = M.add_boundary(m, type = 'demand')
        old_lb = r.lower_bound
        if old_lb != 0:
            print("existing DM with constraints for ", m.id)
            r.bounds = (old_lb + delta, 1000)
        else:
            r.bounds = (delta, 1000)
        if existing_flux:
            r_id, flux = existing_flux
            M.reactions.get_by_id(r_id).bounds = (flux, flux)
        sp = (M.slim_optimize()-solution.objective_value)/delta
    return sp
    
def estimate_shadow_prices(model, intracellular_only = True, delta = 0.1, metabolites = []):
    intracellular_only = True
    wt_growth_rate = model.slim_optimize()
    shadow_prices = {}
    if not len(metabolites):
        metabolites = [m.id for m in model.metabolites]
    for i, m_id in enumerate(metabolites):
        # print(i, m_id)
        m = model.metabolites.get_by_id(m_id)
        if intracellular_only:
            if m.compartment != 'c':
                continue
        if m.id[:5]=='prot_':
            continue
        with model:
            try:
                r = model.reactions.get_by_id('DM_{0}'.format(m.id))
            except KeyError:
                r = model.add_boundary(m, type = 'demand')
            old_lb = r.lower_bound
            r.bounds = (old_lb + delta, 1000)
            shadow_prices[m.id] = (model.slim_optimize()-wt_growth_rate)/delta


def estimate_shadow_prices_reframed(model, constraints, intracellular_only = True, delta = 0.01, metabolites = []):
    intracellular_only = True
    solution = reframed.FBA(model, constraints=constraints)
    # Check if metabolite is already being secreted

    wt_growth_rate = solution.fobj
    shadow_prices = {}
    if not len(metabolites):
        metabolites = [m for m in model.metabolites]

    temp = model.copy()
    for i, m_id in enumerate(metabolites):
        try:
            r = temp.reactions[f'DM_{m_id}']
        except KeyError:
            temp.add_reaction_from_str(f'DM_{m_id}: {m_id} -->  [0, 0]')

    for i, m_id in enumerate(metabolites):
        predicted_flux, r_secretion = get_predicted_metabolite_secretion(model, solution, m_id)
        # print(i, m_id)
        m = model.metabolites[m_id]
        if intracellular_only:
            if m.compartment != 'c':
                continue
        if m.name.startswith('M_prot_'):
            continue
        sp_constraints = {f'DM_{m_id}': (delta, 10)}
        sp_constraints.update(constraints)
        if predicted_flux != 0:
            sp_constraints.update({r_secretion:(predicted_flux, predicted_flux)})

        sp_solution = reframed.FBA(temp, constraints=sp_constraints)
        if isinstance(sp_solution.fobj, float):
            shadow_prices[m_id] = (sp_solution.fobj-wt_growth_rate)/delta
        else:
            shadow_prices[m_id] = np.nan
        
    return shadow_prices


def get_predicted_metabolite_secretion_reframed(model, solution, m_id):
    mc_name = model.metabolites[m_id].name.split(" ")[0]
    predicted_flux = 0
    r_id = None
    for r_id in model.get_exchange_reactions():
        flux = solution.values[r_id]
        r = model.reactions[r_id]
        if  np.abs(flux) >1e-3:
            # print(r_id, r.name, flux, r.lb, r.ub)
            me_id = r.get_substrates()
            if len(me_id):
                me_name = model.metabolites[me_id[0]].name.split(' ')[0]
                if me_name == mc_name:
                    predicted_flux = flux
                    break
    return predicted_flux, r_id
        
    return shadow_prices