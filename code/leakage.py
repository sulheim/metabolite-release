import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from pathlib import Path
import scipy.stats as st
import seaborn as sns
# import leakage

from scipy.interpolate import CubicSpline, UnivariateSpline

"""
E. coli
OD to gDW -  multiple values used
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3663833/ [0.396,0.515]
- https://www.sciencedirect.com/science/article/pii/S1096717618303537?via%3Dihub 0.32 # 
"""

b_subtilis_N_per_OD =  2*1e8 #cells per mL at OD660 = 1 https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=4&id=105286
b_subtilis_weight = 2.2*1e-13 # https://bionumbers.hms.harvard.edu/bionumber.aspx?id=115203&ver=1&trm=bacillus+weight&org=

gDW_per_OD = {
    'e_coli': 0.32, # gDW/L https://www.sciencedirect.com/science/article/pii/S1096717618303537?via%3Dihub
    'yeast' : np.mean([0.644, 0.848]),
    'b_licheniformis': 0.32, # Assume same weight as E.coli (but I think it is bigger -> maybe higher number?), calculation is worng #b_subtilis_weight*b_subtilis_N_per_OD*1000
    'c_glutamicum': 0.32 #Assume same as E. coli
}

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

def get_leakage(data_folder, organism, time, unit = '/gDW', method = 'one-way-diff', only_significant_changes = True):
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
        # leakage, std = estimate_leakage_rate_for_met(time, met_abbrv, df_OD, df_exometabolites, df_exometabolites_std)
        leakage = estimate_leakage_rate_for_met_no_std(time, met_abbrv, df_OD, df_exometabolites, df_exometabolites_std,
                                                       organism = organism, unit = unit, method = method, 
                                                       only_significant_changes=only_significant_changes)
        leakage_list.append(leakage)
        # leakage_uncertainty_list.append(std)
    #print(leakage_list)
    if unit == 'OD':
        name = "Leakage (mM/OD/h)"
    else:
        name = "Leakage (mmol/gDW/h)"
    leakage_df = pd.DataFrame({"Metabolite": met_abbreviations, name:leakage_list})#,"Leakage std": leakage_uncertainty_list})#, columns = ["Time", , "Leakage std"])
    return leakage_df
    
    
def estimate_leakage_rate_for_met_no_std(time, met_abbrv, df_OD, df_exometabolites, df_exometabolites_std, organism, method = 'spline', 
                                            unit = 'OD', only_significant_changes = True):
    # print("Method: ", method)     
    if method == 'spline':
        w = 1/df_OD['OD std']
        w[np.isnan(w)] = np.min(w)
        cs_od = UnivariateSpline(df_OD.index, df_OD['OD mean'], k = 2, s=0.2, ext =3, w=w*1e-1)
        mean_std = np.nanmean(df_exometabolites_std[met_abbrv])
        met_conc =  df_exometabolites[met_abbrv]
        met_conc[np.isnan(met_conc)]=0
        cs_met = UnivariateSpline(df_exometabolites.index, met_conc, s = mean_std*2, k = 1)
        leakage_per_h = cs_met.derivatives(time)[1]
        leakage_rate_per_h_per_OD = leakage_per_h / cs_od(time)

    elif method == 'two-ways-diff':
        OD_mean = df_OD.loc[time, "OD mean"]
        
        # met_conc
        met_conc_1 = df_exometabolites.loc[time-1, met_abbrv]
        met_conc_2 = df_exometabolites.loc[time+1, met_abbrv]    

        # leakage rate
        delta_time = 2 # Divide by delta time = 2 hours
        leakage_per_h = (met_conc_2 - met_conc_1) / delta_time 
        leakage_rate_per_h_per_OD = leakage_per_h / OD_mean
    else:
        # One-way diff
        t1 = np.floor(time)
        t2 = np.ceil(time)
        if t2==t1:
            t2 = t1+1

        OD_mean = df_OD.loc[[t1, t2], "OD mean"].mean()
        # OD = df_OD.loc[t1, 'OD mean'] + (df_OD.loc[t2, 'OD mean']-df_OD.loc[t1, 'OD mean'])*(time-t1)
        delta_time = 1 # Divide by delta time = 2 hours
        # met_conc
        met_conc_1 = df_exometabolites.loc[t1, met_abbrv]
        met_conc_2 = df_exometabolites.loc[t2, met_abbrv]
        met_conc_1_std = df_exometabolites_std.loc[t1, met_abbrv]
        met_conc_2_std = df_exometabolites_std.loc[t2, met_abbrv]

        # Test stats(Welch's t-test)
        n_obs = 3
        t_score, p_val = st.ttest_ind_from_stats(met_conc_1, met_conc_1_std, n_obs, met_conc_2, met_conc_2_std, n_obs, equal_var = False)
        if only_significant_changes and (p_val > 0.05):
            leakage_rate_per_h_per_OD = 0
        else:
            # leakage rate
            leakage_per_h = (met_conc_2 - met_conc_1) / delta_time 
            leakage_rate_per_h_per_OD = leakage_per_h / OD_mean
            
    # The /1000 is to convert from uM to mM
    if unit == 'OD':
        return leakage_rate_per_h_per_OD/1000
    else:
        # OD to gDW
        return leakage_rate_per_h_per_OD/gDW_per_OD[organism]/1000

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

def get_glucose_uptake_rate(data_folder, organism, time, unit = 'gDW', method = 'spline'):
    """
    Unit can be gDW or OD
    """
    exometabolites_folder = Path(data_folder)
    
    # Filenames
    fn_glucose = exometabolites_folder / "{0}_glucose.csv".format(organism)
    fn_OD = exometabolites_folder / "{0}_OD.csv".format(organism)

    # Read files as dataframes
    df_glucose = pd.read_csv(fn_glucose, index_col=0)
    df_OD = pd.read_csv(fn_OD, index_col=0)
    df_glucose.loc[0,:]=[20,0] # From paper

    if method == 'spline':
        cs_od = UnivariateSpline(df_OD.index, df_OD['OD mean'] , s=1)
        cs_glc = UnivariateSpline(df_glucose.index, df_glucose['Glucose mean'] , s=0.5)
        glc_per_h = cs_glc.derivatives(time)[1]/cs_od(time)
        glc_per_h_per_od = glc_per_h/cs_od(time)
    elif method == 'two-ways-diff':
        glc_per_h = (df_glucose.loc[time+1, 'Glucose mean']-df_glucose.loc[time-1, 'Glucose mean'])/2
        glc_per_h_per_od = glc_per_h/df_OD.loc[time, "OD mean"]
    else:
        # One-way diff
        # One-way diff
        t1 = np.floor(time)
        t2 = np.ceil(time)
        if t2==t1:
            t2 = t1+1

        OD_mean = df_OD.loc[[t1, t2], "OD mean"].mean()
        glc_per_h = (df_glucose.loc[t2, 'Glucose mean']-df_glucose.loc[t1, 'Glucose mean'])/1.0
        glc_per_h_per_od = glc_per_h/OD_mean

    
    Mw_glc = 180.156
    glc_mM_per_h_per_od = glc_per_h_per_od/Mw_glc*1000 # Convert from g/L to mM
    
    if unit == 'OD':
        return glc_mM_per_h_per_od
    else:
        # OD to gDW

        return glc_mM_per_h_per_od/gDW_per_OD[organism]

# get growth rate
def get_growth_rate(data_folder, organism, time, method = 'spline'):
    exometabolites_folder = Path(data_folder)
    fn_OD = exometabolites_folder / "{0}_OD.csv".format(organism)
    df_OD = pd.read_csv(fn_OD, index_col=0)
    if method == 'spline':
        cs = UnivariateSpline(df_OD.index, df_OD['OD mean'] , s=1)
        x = np.linspace(0, df_OD.index.max(), 100)
        OD = cs(x)
        ln_OD = np.log(OD)
        idx = np.argmin(np.abs(x-time))
        mu = (ln_OD[idx+1]-ln_OD[idx])/(x[idx+1]-x[idx])

    else:
        df_OD['ln(OD)'] = np.log(df_OD['OD mean'])
        mu = (df_OD.loc[time+1,'ln(OD)'] - df_OD.loc[time-1,'ln(OD)'])/2
    return mu


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
        
    return shadow_prices