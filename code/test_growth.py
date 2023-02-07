import cobra

def check_fva_exchange_reactions(model):
    
if __name__ == '__main__':
    model = cobra.io.read_sbml_model('momentiJO1366.xml')
    solution = model.optimize()
    print(model.summary())