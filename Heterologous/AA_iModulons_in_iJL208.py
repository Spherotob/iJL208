"""
- Integrate AA iModulons into M. florum model iJL208
- Test if AA or other auxotrophies can be lifted
"""

# %%
import cobra
from os.path import abspath, join

import pandas as pd


# %% load models
# M. florum
model = cobra.io.load_json_model(join(abspath("..\Tutorials"), "final_iJL208_28_07_2020.json"))
model_h = model.copy()
model_tot = model.copy() # carries all AA iModulons
# E. coli iML1515
iML1515 = cobra.io.load_json_model("iML1515.json")

# %% test M. florum model
sol = model.optimize()
print(model.summary())


# %% load AA iModulon data and genes
g_data = pd.read_excel("AA_iModulon_summary.xls")
# get iModulon names
iMod_names = g_data.loc[:, "iModulon Name"]


# %% integrate iModulon genes in model
iMod_models = {}
iMod_rxn_data = {}
iMod_rxns_total = []
for iMod in iMod_names:
    rxns_data = {}
    with model_h:
        # get iModulon data
        iMod_data = g_data.loc[g_data["iModulon Name"] == iMod].squeeze() # squeeze into Series
        # get genes
        iMod_bnumbers = iMod_data["B-numbers"].split(";")
        # get reaction of genes
        iMod_rxns = []
        iMod_rxns_unique = []
        for b in iMod_bnumbers:
            try:
                rxns = iML1515.genes.get_by_id(b).reactions
                iMod_rxns.append(rxns)        
                for r in rxns:
                    if r not in iMod_rxns_unique: iMod_rxns_unique.append(r)    
                    if r not in iMod_rxns_total: iMod_rxns_total.append(r)
            except:
                print("Gene " + b + " not in iML1515")
        # print reactions  
        for r in iMod_rxns_unique:
            print(r)
        
        # integrate reaction into M. florum model
        model_h.add_reactions(iMod_rxns_unique)
        # integrate reactions into all-iModulons model
        model_tot.add_reactions(iMod_rxns_unique)
        # get novel compounds
        novel_cmp = []
        for r in iMod_rxns_unique:
            for m in r.metabolites:
                try:
                    model.metabolites.get_by_id(m.id)
                except:
                    novel_cmp.append(m.id)
                    
        # save model
        iMod_models[iMod] = model_h.copy()
        
        # test iModulon reactions
        iMod_rxns_flux = {}
        active_flux_count = 0
        for r in iMod_rxns_unique:
            model_h.objective = r.id
            sol = model_h.optimize()
    #         print(model_h.metabolites.gln__L_c.shadow_price)
            if sol.status=="optimal":
                iMod_rxns_flux[r.id] = sol.objective_value
            else:
                iMod_rxns_flux[r.id] = None 
                
            # active flux count
            if iMod_rxns_flux[r.id] is not None and iMod_rxns_flux[r.id] != 0:
                active_flux_count += 1
                
        rxns_data["fluxes"] = iMod_rxns_flux
        rxns_data["active_reactions"] = active_flux_count
        
        iMod_rxn_data[iMod] = rxns_data
        
# %% investigate iJL208 with all AA iModulons
iMod_rxns_flux = {}
active_flux_count = 0
for r in iMod_rxns_total:

    with model_tot:
        model_tot.objective = r.id
        sol = model_tot.optimize()

        if sol.status=="optimal":
            iMod_rxns_flux[r.id] = sol.objective_value
        else:
            iMod_rxns_flux[r.id] = None 
            
        # active flux count
        if iMod_rxns_flux[r.id] is not None and iMod_rxns_flux[r.id] != 0:
            active_flux_count += 1
