import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from pathlib import Path
import scipy.stats as st
import seaborn as sns

def get_leakage(data_folder, organism, time):
    exometabolites_folder = Path(data_folder)
    
    # Filenames
    fn_exometabolites = exometabolites_folder / "{0}_exometabolites.csv".format(organism)
    fn_exometabolites_std = exometabolites_folder / "{0}_exometabolites_std.csv".format(organism)
    fn_OD = exometabolites_folder / "{0}_OD.csv".format(organism)
    
    # Read files as dataframes
    df_exometabolites = pd.read_csv(fn_exometabolites, index_col=0)
    df_exometabolites_std = pd.read_csv(fn_exometabolites_std, index_col=0)
    df_OD = pd.read_csv(fn_OD, index_col=0)
    
    met_abbreviations = [x.replace(" MS", "") for x in df_exometabolites.columns]
    leakage_list = []
    leakage_uncertainty_list = []
    for met_abbrv in met_abbreviations:
        leakage, std = estimate_leakage_rate_for_met(time, met_abbrv, df_OD, df_exometabolites, df_exometabolites_std)
        leakage_list.append(leakage)
        leakage_uncertainty_list.append(std)
    #print(leakage_list)
    leakage_df = pd.DataFrame({"Metabolite": met_abbreviations, "Leakage (uM/h/OD)":leakage_list,"Leakage std": leakage_uncertainty_list})#, columns = ["Time", , "Leakage std"])
    return leakage_df
    
    

def estimate_leakage_rate_for_met(time, met_abbrv, df_OD, df_exometabolites, df_exometabolites_std):
    """
    Estimate leakage by dividing the change in concentration of two timepoints by the OD of the mean.
    Time = 6 or 7 seems to be mid-exponential phase
    """
    OD_mean = df_OD.loc[time, "OD mean"]
    OD_std = df_OD.loc[time, "OD std"]
    
    # met_conc
    met_conc_1 = df_exometabolites.loc[time-1, met_abbrv]
    met_conc_2 = df_exometabolites.loc[time+1, met_abbrv]
    met_conc_1_std = df_exometabolites_std.loc[time-1, met_abbrv]
    met_conc_2_std = df_exometabolites_std.loc[time+1, met_abbrv]
    
    # leakage rate
    delta_time = 2 # Divide by delta time = 2 hours
    leakage_per_h = (met_conc_2 - met_conc_1) / delta_time 
    leakage_rate_per_h_per_OD = leakage_per_h / OD_mean
    
    if leakage_per_h != 0:
        # Error propagation
        # Assume no covariance; https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulae
        leakage_per_h_std = np.sqrt(met_conc_1_std**2 + met_conc_2_std**2)/delta_time
        #print(leakage_per_h_std, leakage_per_h, OD_std, OD_mean)
        rel_leakage_std = np.sqrt((leakage_per_h_std/leakage_per_h)**2+(OD_std/OD_mean)**2)
        leakage_std = rel_leakage_std*leakage_rate_per_h_per_OD
    else:
        leakage_std = 0
    return leakage_rate_per_h_per_OD, leakage_std
    
    
    #
def estimate_shadow_prices(model, intracellular_only = True, epsilon = 0.1):
    intracellular_only = True
    wt_growth_rate = model.slim_optimize()
    shadow_prices = {}
    for m in model.metabolites:
        if intracellular_only:
            if m.compartment != 'c':
                continue
        with model:
            try:
                r = model.reactions.get_by_id('DM_{0}'.format(m.id))
            except KeyError:
                r = model.add_boundary(m, type = 'demand')
            old_lb = r.lower_bound
            r.bounds = (old_lb + epsilon, 1000)
            shadow_prices[m.id] = (model.slim_optimize()-wt_growth_rate)/epsilon
        
    return shadow_prices