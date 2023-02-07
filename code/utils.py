#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name: utils.py
date: 15.01.21
author: Snorre Sulheim

A script for testing the genome-scale metabolic models (GEMs) used in the SEP-outcompete project
"""

import cobra
from pathlib import Path
import pandas as pd
from dotenv import find_dotenv

# File names (relative to root folder)
MODEL_INFO_FN = r"models\list_of_models.csv"
KM_DEFAULT= 0.01 #mmol/L
VMAX_DEFAULT = 10 # mmol/gDW/h
M9_AEROBE_FN = "M9_basis_aerobe.csv"

def get_results_folder():
    return get_project_root() / "results"

def get_original_models_folder():
    return get_project_root() / "models" / "original_models"

def get_project_root():
    denv_path = find_dotenv()
    if len(denv_path):
        root_folder = Path(denv_path).parent
    else:
        root_folder = Path(__file__).parent.parent.parent
    return root_folder

def get_data_folder():
    return get_project_root() / "data"

def get_log_folder():
    return get_project_root() / "results" / "log"

def get_curated_models_folder():
    return get_project_root() / 'models' / 'curated_models'

def get_yaml_model(strain_name = None, strain_abb = None):
    if strain_name:
        model_info = get_model_info(strain_name = strain_name)
    else:
        model_info = get_model_info(strain_abb = strain_abb)

    model_fn = model_info["Abbreviation"]
    return cobra.io.load_yaml_model(str(get_curated_models_folder() / "{0}.yml".format(model_fn)))

def get_model_info(strain_name = None, strain_abb = None):
    """
    Returns info about:
    a) One strain if strain_name or strain_abb are provided (as dict)
    b) All strains if not (as dataframe)
    """
    model_info = pd.read_csv(get_project_root() / MODEL_INFO_FN, sep = ",", header = 0, skipinitialspace=True)
    
    if strain_name and strain_name.lower() in list(model_info["Name"].str.lower()):
        row_idx = model_info["Name"].str.lower() == strain_name.lower()
        return model_info.loc[row_idx, :].iloc[0].to_dict()

    elif strain_abb and strain_abb.lower() in list(model_info["Abbreviation"].str.lower()):
        row_idx = model_info["Abbreviation"].str.lower() == strain_abb.lower()
        return model_info.loc[row_idx, :].iloc[0].to_dict()

    else:
        print("Model name/abb {0}/{1} not found in list".format(strain_name, strain_abb))
        return model_info

def set_FBA_medium(model, carbon_sources = "glc__D_e", rate = 4):
    """
    carbon_sources can be a string or a list of strings (met_IDs).
    Assumes a m9 medium as a basis.
    rate can be a float/int or a list of floats / ints equal to the length of carbon_sources
    in g/L


    """
    set_FBA_M9_basis(model)

    # Make sure carbon sources is a list
    if isinstance(carbon_sources, str):
        carbon_sources = [carbon_sources]

    if isinstance(rate, list):
        if len(rate)!= len(carbon_sources):
            raise ValueError
    else:
        rate = [rate]*len(carbon_sources)

    for m_id, ra in zip(carbon_sources, rate):
        try:
            m = model.metabolites.get_by_id(m_id)
        except KeyError:
            continue
        r = get_exchange_reaction(m)
        if r:
            r.lower_bound = -abs(ra)


def set_FBA_M9_basis(model):
    """
    Set the lower bound to -1000 for exchange reactions providing the basis (minerals, salts, ammonium, phosphate) for M9 minimal medium
    """
    m9_basis_fn = get_data_folder() / "medium" / M9_AEROBE_FN
    M9_df = pd.read_csv(m9_basis_fn, header = 0)
    
    # First, close all uptake and open all secretion
    for r in model.exchanges:
        r.lower_bound = 0
        r.upper_bound = 1000

    # Then set M9 medium
    for _, row in M9_df.iterrows():
        m_id = row["ID"]
        try:
            m = model.metabolites.get_by_id(m_id)
        except KeyError:
            continue
        r = get_exchange_reaction(m)
        if r:
            #print(m.id, r.id, r.reaction)
            r.lower_bound = -1000

def get_M9_basis_exchange_reactions(model):
    m9_basis_fn = get_data_folder() / "medium" / M9_AEROBE_FN
    M9_df = pd.read_csv(m9_basis_fn, header = 0)
    exchange_reactions = []
    # Then set M9 medium
    for _, row in M9_df.iterrows():
        m_id = row["ID"]
        try:
            m = model.metabolites.get_by_id(m_id)
        except KeyError:
            continue
        r = get_exchange_reaction(m)
        exchange_reactions.append(r.id)
    return exchange_reactions

def get_m9_metabolite_ids():
    m9_basis_fn = get_data_folder() / "medium" / M9_AEROBE_FN
    M9_df = pd.read_csv(m9_basis_fn, header = 0)
    return list(M9_df.ID)

def get_exchange_reaction(metabolite):
    for r in metabolite.reactions:
        if len(r.metabolites) == 1:
            return r
    # No exchange reaction, raise error
    print("No exchange reaction for ", metabolite.id)
    return None

def MM(conc, km = KM_DEFAULT, Vmax = VMAX_DEFAULT):
    """
    Michelis Menten kinetics
    """
    return Vmax*conc/(km+conc)


def convert_gL_to_mmol(formula_string, gL, space_width = 1):
    """
    Convert concentration in g/L to mmol for a given metabolite formula and a given grid space width.

    Args:
        formula_string (str): Chemical formula
        gL (float, int): The medium concentration of this medium compound in gram/L
        space_width (float, int): The space width used in the respective comets simulation in cm 
    """
    formula = cobra.core.formula.Formula(formula_string)
    mW = formula.weight
    volume = (space_width*1e-1)**3 # Volume in dm^3 == L
    # Calculate amount in mmol
    mol = gL*volume/mW
    mmol = mol*1e3
    return mmol

def convert_gL_to_mM(formula_string, gL):
    """
    Convert concentration in g/L to mM for a given metabolite formula and a given grid space width.

    Args:
        formula_string (str): Chemical formula
        gL (float, int): The medium concentration of this medium compound in gram/L
    """
    formula = cobra.core.formula.Formula(formula_string)
    mW = formula.weight
    # Calculate amount in mmol
    M = gL/mW
    mM = M*1e3
    return mM

def get_biomass_reaction(model):
    for r in model.reactions:
        if "biomass" in r.id.lower():
            return r
        elif "biomass" in r.name.lower():
            return r
    return None

def convert_gL_to_mmol_cobra(metabolite, gL, space_width):
    """
    DEPRECATED
    Convert concentration in g/L to mmol for a given metabolite and a given grid space width. '

    Args:
        metabolite (cobra.Metabolite): The metabolite that you want to calculate the medium amount for. Has to feature a chemical formula
        gL (float, int): The medium concentration of this medium compound in gram/L
        space_width (float, int): The space width used in the respective comets simulation in cm 
    """
    mW = metabolite.formula_weight # in g/mol
    volume = (space_width*1e-1)**3 # Volume in dm^3 == L
    # Calculate amount in mmol
    mol = gL*volume/mW
    mmol = mol*1e3
    return mmol

def get_exchange_metabolites(model, reaction_list):
    mets = []
    with model:
        for r_id in reaction_list:
            r = model.reactions.get_by_id(r_id)
            (m, i) = r.metabolites.popitem() 
            mets.append(m.id)
    return mets

def get_selected_carbon_sources(carbon_sources_string):

    select_few_mets = False
    if carbon_sources_string is None:
        fn = get_data_folder() / "selected_carbon_sources.csv"
    elif (isinstance(carbon_sources_string, str)):
        if carbon_sources_string.split(".")[-1] in ["csv", "tsv", "xlsx", "txt"]:
            fn = get_data_folder() / carbon_sources_string
        else:
            select_few_mets = True
    df = pd.read_csv(fn, index_col = 0)

    if select_few_mets:
        carbon_source_ids = [x.strip() for x in carbon_sources_string.split(",")]
        df = df.loc[df.index.isin(carbon_source_ids), :]
    return df

def get_strain_matrix():
    fn = get_data_folder() / "strain_met_matrix.csv"
    df = pd.read_csv(fn, index_col = 0)
    df.set_index("ID", inplace = True)
    return df

def get_models(model_names = None):
    """
    Get the cobra models for the model names given (list). If model_names = None, all models are loaded. 
    """
    model_folder = get_curated_models_folder()
    models = {}
    if not model_names:
        for model_fn in model_folder.glob("*.yml"):
            model_name = model_fn.stem
            model = cobra.io.load_yaml_model(model_fn)
            model.id = model_name
            models[model_name] = model
    else:
        for model_name in model_names:
            model_fn = model_folder / "{0}.yml".format(model_name)
            model = cobra.io.load_yaml_model(model_fn)
            model.id = model_name
            models[model_name] = model
    return models

def find_objective(model):
    objectives =  []
    for r in model.reactions:
        if r.objective_coefficient == 1:
            objectives.append(r.id)

    if len(objectives) == 1:
        return objectives[0]
    else:
        return objectives

# def is_transport_reaction(r):
#     base_mets = list(set([m.id.rsplit("_")[0] for m in r.metabolites.keys()]))
#     remove_mets = ['h2o', 'h', 'atp', 'adp', 'na1', 'k', 'amp', 'pi', 'pep', 'pyr']
#     other_mets = [x for x in base_mets if x not in remove_mets]
#     if len(other_mets) <= 1:
#         return True
#     else:
#         return False

def make_formula_sheet():
    xml_fn = "../../models/curated_models/E_coli.xml"
    model = cobra.io.read_sbml_model(xml_fn)
    id_formula_dict = {}
    for m in model.metabolites:
        if m.compartment == "e":
            id_formula_dict[m.id] = m.formula
    s = pd.Series(id_formula_dict)
    s.name = "Formula"
    fn = get_data_folder() / "formula.csv"
    s.to_csv(fn)

def get_formula(met_id):
    fn = get_data_folder() / "formula.csv"
    s = pd.read_csv(fn, index_col = 0, squeeze = True).to_dict()
    return s[met_id]

def get_formula_dict():
    fn = get_data_folder() / "formula.csv"
    s = pd.read_csv(fn, index_col = 0, squeeze = True).to_dict()
    return s    

if __name__ == '__main__':
    # print(get_model_info("bacillus subtilis"))
    if 0:
        model = cobra.io.load_yaml_model(r"../../models/curated_models/E_coli.yml")
        m = model.metabolites.get_by_id("glu__L_e")
        convert_gL_to_mmol(m, 4, 1)

    if 0:
        get_selected_carbon_sources()
    if 0:
        get_strain_matrix()
    if 0:
        print(get_project_root())
    if 1:
        make_formula_sheet()
        get_formula("glc__D_e")