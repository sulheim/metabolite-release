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
        # print(i, m_id)
        m = model.metabolites[m_id]
        if intracellular_only:
            if m.compartment != 'c':
                continue
        if m.name.startswith('M_prot_'):
            continue
        sp_constraints = {f'DM_{m_id}': (delta, 10)}
        sp_constraints.update(constraints)
        sp_solution = reframed.FBA(temp, constraints=sp_constraints)
        if isinstance(sp_solution.fobj, float):
            shadow_prices[m_id] = (sp_solution.fobj-wt_growth_rate)/delta
        else:
            shadow_prices[m_id] = np.nan
        
    return shadow_prices