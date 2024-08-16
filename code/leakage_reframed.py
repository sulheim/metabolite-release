import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from pathlib import Path
import scipy.stats as st
import seaborn as sns
import reframed

def estimate_shadow_prices(model, constraints, intracellular_only = True, delta = 0.01, metabolites = []):
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


def get_predicted_metabolite_secretion(model, solution, m_id):
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
                
    
