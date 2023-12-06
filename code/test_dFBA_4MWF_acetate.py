#!/usr/bin/env python
# coding: utf-8


import reframed
from pathlib import Path
import pandas as pd
import utils
from dFBA_reframed import dFBA
import numpy as np

def set_auxotrophy_dict(auxotrophy_dict):
    constraints_dict = {}
    for species_abbr, species_dict in auxotrophy_dict.items():
        species_constraint_dict = {}
        if species_dict.get('vitamins'):
            for met_id in species_dict['vitamins']:
                r_id = f'R_EX_{met_id}_e'
                species_constraint_dict[r_id] = (-0.001, 0)
        if species_dict.get('amino acids'):
            for met_id in species_dict['amino acids']:
                r_id = f'R_EX_{met_id}_e'
                species_constraint_dict[r_id] = (-0.2, 0)
        constraints_dict[species_abbr] = species_constraint_dict
    return constraints_dict

if __name__ == '__main__':

    growth_rates_fn = '/Users/ssulheim/Library/CloudStorage/OneDrive-SINTEF/UNIL/MWF/experimental work/230525_growth_phenotyping/fitted_growth_parameters.csv'
    growth_rates_df = pd.read_csv(growth_rates_fn, index_col=0)

    cs = 'Acetate'
    growth_rate_dict = {}
    for species_abbr in growth_rates_df.Species.unique():
        idx = (growth_rates_df['Species']==species_abbr) & (growth_rates_df['Carbon source']==cs)
        growth_rate_dict[species_abbr] = growth_rates_df.loc[idx, 'r'].values[0]


    model_folder = Path('/Users/ssulheim/git/mwf_gems/models/gf_models/carveme_TFA')
    model1_fn = model_folder / 'polished_At.xml'
    model2_fn = model_folder / 'polished_Ct.xml'
    model3_fn = model_folder / 'polished_Oa.xml'
    model4_fn = model_folder / 'polished_Ml.xml'

    model_name_dict = {}
    yield_dict = {}
    vmax_dict = {}
    for fn, name in zip([model1_fn, model2_fn, model3_fn, model4_fn],['At', 'Ct', 'Oa', 'Ml']):
        model = reframed.load_cbmodel(fn)
        model.reactions.R_EX_h_e.lb = -1000
        model.reactions.R_EX_o2_e.lb = -1000
        model_name_dict[name] = [model, 0.01]
        model.solver = 'gurobi'
        fba_result = reframed.FBA(model, constraints={'R_EX_glc__D_e':0, 'R_EX_cit_e':0, 'R_EX_ac_e':-10})
        yield_dict[name] = fba_result.fobj/10
        vmax_dict[name] = growth_rate_dict[name]/yield_dict[name]


    auxotrophy_dict = {
        'Ml': {'amino acids': ['pro__L', 'cys__L'], 'vitamins':['thm', 'btn']},
        'Oa': {'vitamins': ['thm']}
        }
    auxotrophy_constraints = set_auxotrophy_dict(auxotrophy_dict)


    N_runs = 200
    for j in range(N_runs):
        cs_formula = model.metabolites['M_ac_c'].metadata['FORMULA']
        ac_mM = utils.convert_gL_to_mM(cs_formula, 2)
        dt = 0.25
        total_time = 58
        iterations = int(total_time//dt)
        dilution_interval = 20
        D = dFBA(iterations = iterations, dt = dt, method = "FBA_with_leakage", store_exchanges_flag = False, fraction_of_optimum=0.9, 
                 dilution_interval=dilution_interval, folder = 'dFBA_191123')
        D.medium.define_initial_conditions({"M_ac_e": ac_mM})
        D.add_models(model_name_dict, auxotrophy_constraints=auxotrophy_constraints)

        # Set Km and vMax
        noise = np.random.normal(0, 1, len(model_name_dict))
        for i, name in enumerate(model_name_dict.keys()):
            D.models[name].set_km("M_ac_e", 1)
            D.models[name].set_Vmax("M_ac_e", max(0, vmax_dict[name]+noise[i]))

        # D.medium.set_store_concentrations(["glc__D_e", "nh3_e"])
        D.run()
        # D.store_results(True)
        D.store_results(True, tag = j+1)


