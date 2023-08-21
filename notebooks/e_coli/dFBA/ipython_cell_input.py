import cobra

model = cobra.io.read_sbml_model('../../../models/e_coli/iJO1366.xml')

# model = cobra.io.read_sbml_model('../../../models/e_coli/iJR904.xml')
# model = cobra.io.read_sbml_model('../../models/e_coli/e_coli_core.xml')

model.solver = 'gurobi'
# Initial conditions is 0.013 gDW/L (in total)
model_name = 'Ecoli'
model_name_dict = {model_name: [model, 0.006]}

glucose_mM = utils.convert_gL_to_mM("C6H12O6", 4)
D = FBA_leak.dFBA(iterations = 3, dt = 0.2, method = "FBA_with_leakage", store_exchanges_flag = False)#(, pfba_fraction = 0.95)
D.add_models(model_name_dict)

# Set Km and vMax
D.models[model_name].set_km("glc__D_e", 1)
D.models[model_name].set_Vmax("glc__D_e", 10)



D.medium.define_initial_conditions({"glc__D_e": glucose_mM})
# D.medium.set_store_concentrations(["glc__D_e", "nh3_e"])
D.run()
