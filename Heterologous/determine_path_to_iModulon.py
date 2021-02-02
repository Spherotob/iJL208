"""
Determine metabolic path from native M. florum to iModulon AA synthesis
"""

import cobra
from optlang.symbolics import Zero

from os.path import abspath, join
from itertools import chain

import pandas as pd
import numpy as np


# %% load models
# M. florum
model = cobra.io.load_json_model(join(abspath("..\Tutorials"), "final_iJL208_28_07_2020.json"))
model_original = model.copy()
model_a = model.copy() # augmented with E. coli reactions/genes
# E. coli iML1515
iML1515 = cobra.io.load_json_model("iML1515.json")

# enforce growth
# model_a.reactions.get_by_id("BIOMASS_step3_c").lower_bound = 0.1
# lift secretion bounds
model_a.reactions.get_by_id("EX_ac_e").bounds = (0, 1000)
model_a.reactions.get_by_id("EX_lac__L_e").bounds = (0, 1000)

# enable hugher carbon uptake fluxes
# model_a.reactions.EX_sucr_e.lower_bound = -10
# model_a.reactions.EX_rib__D_e.lower_bound = -10
# model_a.reactions.EX_fru_e.lower_bound = -10
# relax NOX2 bounds
# model_a.reactions.NOX2.upper_bound = 1000

# %% user parameter
filename = "AA_iModulon_in_iJL208_minimal_yield.xlsx"

bio = "BIOMASS_step3_c"     # biomass formation equation ID

add_transport_reactions = True # add transport/exchange reactions for non-native metabolites introduced

id_prefix_ecoli = "_ecoli"     # denotes reactions from iML1515 in iJL208
id_prefix_iMod = "_iMod"     # denotes iModulon reactions
# mode for determining set of heterologous reaction towards iModulons
mode = "min_yield" # "max_yield": Maximize AA yield/flux; "min_yield": Guarantee a minimal AA yield
bio_min_flux = 0.05 # minimal growth rate relative to maximum growth(~0.22 1/h)


# amino acids
AA_IDs = ["gly", "arg__L", "ala__L", "leu__L", "thr__L", "arg__L", "asn__L", "glu__L",
      "asp__L", "gln__L", "his__L", "ile__L", "lys__L", "met__L", "phe__L", "ser__L",
      "trp__L", "tyr__L", "val__L", "cys__L"]

# %% restrict carbon uptake reactions
# determine flux distribution
bound_factor = 1 # bounds on exchange flux based on WT flux distribution times this factor
sol_wt = model.optimize()
for ex_rxn in model.exchanges:
    if len(ex_rxn.reactants) == 1:
        ex_met = ex_rxn.reactants[0]
    else:
        continue
    if "C" in list(ex_met.elements.keys()) and sol_wt.fluxes.loc[ex_rxn.id] < 0:
        model_a.reactions.get_by_id(ex_rxn.id).lower_bound = sol_wt.fluxes.loc[ex_rxn.id]*bound_factor
        print(ex_rxn.id + ": " + str(model_a.reactions.get_by_id(ex_rxn.id).lower_bound))
# set bounds manually
model_a.reactions.DM_gly_c.lower_bound = sol_wt.fluxes.loc["DM_gly_c"]*bound_factor
model_a.reactions.EX_1pyr5c.upper_bound = sol_wt.fluxes.loc["EX_1pyr5c"]*bound_factor


# %% calculate wild-type growth
sol = model_a.optimize()
growth_WT = sol.objective_value
print(model_a.summary())
# save model
model_p = model_a.copy()

# %% augment iJL208 with iML1515 reactions (gene associated)
ecoli_rxns = []
def create_rxn_from_object(rxnID, rxnObj):
    # copy reaction
    rxn_copy = cobra.Reaction(id=rxnID,
                              name=rxnObj.name,
                              subsystem=rxnObj.subsystem,
                              lower_bound=rxnObj.lower_bound,
                              upper_bound=rxnObj.upper_bound)
    # add metabolites
    rxn_copy.add_metabolites(rxnObj.metabolites)
    # parse gene associations
    rxn_copy.gene_reaction_rule = rxnObj.gene_reaction_rule   
    return rxn_copy

iJL208_rxns = [rxn.id for rxn in model.reactions]
for rxn in iML1515.reactions:
    # check if reaction is associated with genes
    if len(rxn.genes)==0:
        continue
        # pass
    # reaction already in iJL208?
    if rxn.id in iJL208_rxns:
        continue
    # do not include transmembrane protein/reactions
    comp = []
    for met in rxn.metabolites.keys(): comp.append(met.compartment)
    if len(np.unique(comp)) > 1:
        continue
    
    
    
    # copy reaction
    rxn_copy = create_rxn_from_object(rxn.id + id_prefix_ecoli, rxn)
    # add reaction to iJL208
    model_a.add_reaction(rxn_copy)
    # save reaction
    ecoli_rxns.append(rxn_copy.id)

# test augmented reactions
augment_fluxes = pd.Series(dtype="float64")
for r in model_a.reactions:
    if id_prefix_ecoli in r.id:
        with model_a:
            model_a.objective = r.id
            sol = model_a.slim_optimize()
            augment_fluxes = augment_fluxes.append(pd.Series([sol], [r.id]))


# %% load AA iModulon data and genes
g_data = pd.read_excel("AA_iModulon_summary.xls")
# get iModulon names
iMod_names = g_data.loc[:, "iModulon Name"]

# define functions
def add_iModulon_to_model(model, g_data, iMod):
       
    # get iModulon data
    iMod_data = g_data.loc[g_data["iModulon Name"] == iMod].squeeze() # squeeze into Series
    # get genes
    iMod_bnumbers = iMod_data["B-numbers"].split(";")
    # get reaction of genes
    iMod_rxns_unique = []
    for b in iMod_bnumbers:
        try:
            rxns = iML1515.genes.get_by_id(b).reactions       
            for r in rxns:
                if r not in iMod_rxns_unique: iMod_rxns_unique.append(r)    
        except:
            print("Gene " + b + " not in iML1515")
            
    # mark iModulon reactions in M. florum model
    iMod_rxns = []
    for iMod_rxn in iMod_rxns_unique:
        
        # add exchange reaction (secretion only) for novel, non-native metabolites
        if add_transport_reactions:
            for met in iMod_rxn.metabolites:
                try:
                    model_original.metabolites.get_by_id(met.id)
                    # metabolite already in original iJL208, do nothing
                except:
                    # metabolite not in model, add secretion reaction
                    scrt_rxn_id = "EX_" + met.id
                    scrt_rxn = cobra.Reaction(id=scrt_rxn_id,
                                              lower_bound=0,
                                              upper_bound=1000)
                    # check if metabolite was introduced by other heterologous reactions
                    try:
                        scrt_met = model.metabolites.get_by_id(met.id)
                    except:                        
                        scrt_met = cobra.Metabolite(id=met.id,
                                                formula=met.formula,
                                                name=met.name,
                                                charge=met.charge,
                                                compartment=met.compartment)
                    
                    scrt_rxn.add_metabolites({scrt_met: -1})
                    model.add_reaction(scrt_rxn)
                        
        
        try:
            rxn = model.reactions.get_by_id(iMod_rxn.id + id_prefix_ecoli)
            # switch reactions to savely change id
            rxn_id_change = create_rxn_from_object(iMod_rxn.id + id_prefix_iMod, rxn)         
            model.add_reaction(rxn_id_change)
            model.reactions.get_by_id(iMod_rxn.id + id_prefix_ecoli).delete()
        
        except:
            print("iModulon reaction not in model: " + iMod_rxn.id)
            rxn_id_change = create_rxn_from_object(iMod_rxn.id + id_prefix_iMod, iMod_rxn)
            model.add_reaction(rxn_id_change)

        iMod_rxns.append(rxn_id_change.id)
        
        
        
        
    return iMod_rxns
     
    
def disable_amino_acid_uptake(model, amino_acid):
    # disable amino acid uptake (or respective di-peptide)
    for ex_rxn in model.exchanges:
        if amino_acid in ex_rxn.id:
            ex_rxn.bounds = (0, 0)
            
        # consider misspelling of ala__L in di-peptide
        if amino_acid == "ala__L":
            if "ala_L" in ex_rxn.id:
                ex_rxn.bounds = (0, 0)
                
    return model
                
    

# %% integrate iModulon genes in model
iMod_models = {}
iMod_rxn_data = {}
iMod_AA_association = {}
iMod_AA_production = {}
iMod_AA_solutions = {}
for iMod in iMod_names:
    rxns_data = {}
    with model_a:
        
        # add iModulon to model
        iMod_rxns = add_iModulon_to_model(model_a, g_data, iMod)
      
        # determine which AA can be synthesized by iModulon reactions
        AA_synthesis = []
        for iMod_rxn_id in iMod_rxns:
            iMod_rxn = model_a.reactions.get_by_id(iMod_rxn_id)
            for met, stoic in iMod_rxn.metabolites.items():
                if met.id[:-2] in AA_IDs and met.id not in AA_synthesis and met.id[-2:]=="_c":
                    if stoic < 0 and iMod_rxn.lower_bound < 0:
                        AA_synthesis.append(met.id)
                    elif stoic > 0 and iMod_rxn.upper_bound > 0:
                        AA_synthesis.append(met.id)                    
        iMod_AA_association[iMod] = AA_synthesis
        
        # calculate path towards iModulon AA synthesis/overprodction
        iMod_AA_production[iMod] = {}
        iMod_AA_solutions[iMod] = {}
        
        for AA_syn in AA_synthesis:   

            # protect model
            with model_a:

                # disable all E. coli reactions with AA participation to avoid production outside iModulon
                for aa_rxn in model_a.metabolites.get_by_id(AA_syn).reactions:
                    if id_prefix_ecoli in aa_rxn.id:
                        model_a.reactions.get_by_id(aa_rxn.id).bounds = (0, 0)
                        pass
                        
                # disable amino acid uptake (or respective di-peptide)     
                model_a = disable_amino_acid_uptake(model_a, AA_syn[:-2])
                # # disable amino acid uptake (or respective di-peptide)
                # for ex_rxn in model_a.exchanges:
                #     if AA_syn[:-2] in ex_rxn.id:
                #         ex_rxn.bounds = (0, 0)
                #         pass
                #     # consider misspelling of ala__L in di-peptide
                #     if AA_syn[:-2] == "ala__L":
                #         if "ala_L" in ex_rxn.id:
                #             ex_rxn.bounds = (0, 0)
                #             pass
                  
                           
                  
                # create amino acid exchange
                aa_met = model_a.metabolites.get_by_id(AA_syn)
                ex_aa_id = "EX_" + AA_syn
                ex_aa = cobra.Reaction(id=ex_aa_id,
                                       lower_bound=0,
                                       upper_bound=1000)
                ex_aa.add_metabolites({aa_met: -1})
                model_a.add_reaction(ex_aa)
                
                
                ################
                if mode == "max_yield":               
                    # maximize AA yield/rate, pFBA
                    model_a.objective = ex_aa_id
                    sol = model_a.optimize()
                    if sol.status == "optimal":
                        # calculate pFBA solution
                        sol = cobra.flux_analysis.pfba(model_a)  
                    else:
                        sol = []
                    
                                       
                elif mode =="min_yield":
                    # only guarantee a minimal AA flux/yield
                    
                    # # get e coli reactions
                    # ecoli_rxn_obj = []
                    # for rxn in model_a.reactions:
                    #     if id_prefix_ecoli in rxn.id:
                    #         ecoli_rxn_obj.append(rxn)
                    # # determine maximum growth
                    # # block E. coli reactions
                    # with model_a:
                    #     # get e coli reactions
                    #     ecoli_rxn_obj = []
                    #     for rxn in model_a.reactions:
                    #         if id_prefix_ecoli in rxn.id:
                    #             model_a.reactions.get_by_id(rxn.id).bounds = (0, 0)
                                           
                    #     max_growth = model_a.slim_optimize(error_value=-1)
                    # print(max_growth)
                        
                    if growth_WT != -1:
                        # enforce flux on Biomass  
                        with model_a:
                            model_a.reactions.get_by_id(bio).lower_bound = growth_WT*bio_min_flux
                            model_a.reactions.get_by_id(bio).upper_bound = growth_WT*bio_min_flux
                            model_a.reactions.get_by_id(ex_aa_id).lower_bound = 1e-2
                            # maximize AA synthesis
                            model_a.objective = ex_aa_id
                            model_a.objective_direction = "min"
                            sol = model_a.slim_optimize(error_value=-1)
                            if sol != -1:
                                
                                # conduct pfba 
                                # sol = cobra.flux_analysis.pfba(model_a)
                                
                                # conduct pfba manually, only with E. coli reactions
                                # get e coli reactions
                                ecoli_rxn_obj = []
                                for rxn in model_a.reactions:
                                    if id_prefix_ecoli in rxn.id:
                                        ecoli_rxn_obj.append(rxn)
                                        
                                reaction_variables = ((rxn.forward_variable,
                                            rxn.reverse_variable) for rxn in ecoli_rxn_obj)        
                                variables = chain(*reaction_variables)   
                                
                                model_a.objective = model_a.problem.Objective(
                                    Zero, direction="min", sloppy=False, name="_pfba_objective"
                                    )
                                model_a.objective.set_linear_coefficients({v: 1.0 for v in variables})
                                
                                sol = model_a.optimize()
                            else:
                                sol = []
                        
                    else:
                        sol = []
                            
                            
                            
                
                if sol: 
                    # save AA synthesis rate
                    iMod_AA_production[iMod][AA_syn] = sol[ex_aa_id]
                    # save flux solution
                    iMod_AA_solutions[iMod][AA_syn] = sol
                    # # determine necessary E. coli heterologous reactions
                    # for f in sol.fluxes.index:
                    #     if f.endswith(id_prefix_ecoli) and sol.fluxes[f] != 0:
                    #         iMod_AA_solutions[iMod][AA_syn].append(f[:-len(id_prefix_ecoli)])
                            
                else:
                    iMod_AA_production[iMod][AA_syn] = None
        
        
# %% Analyze solutions
iMod_res = pd.DataFrame(columns=["iModulon",
                                 "Amino acid",
                                 "E. coli connect reactions", 
                                 "E. coli connect reaction fluxes",
                                 "E.coli connect genes", "E. coli connect subsystems",
                                 "Amino acid synthesis rate [mmol/gDW/h]",
                                 "Maximum growth rate [1/h]",
                                 "Amino acid co-synthesis",
                                 "Amino acid co-synthesis growth [1/h]"])

# loop through solutions
c = 0
for iMod in iMod_AA_solutions.items():
    for AA in iMod[1].items():
        # save AA synthesis rate
        iMod_res.loc[c, "Amino acid synthesis rate [mmol/gDW/h]"] \
            = iMod_AA_production[iMod[0]][AA[0]]
        
        # load solution
        sol = AA[1]
        # make new line
        iMod_res.loc[c, "iModulon"] = iMod[0]
        iMod_res.loc[c, "Amino acid"] = AA[0]
        # determine necessary E. coli heterologous reactions
        ecoli_rxns = []
        ecoli_subsystems = []
        ecoli_genes = []
        ecoli_rxns_flux = []
        for f in sol.fluxes.index:
            if f.endswith(id_prefix_ecoli) and abs(sol.fluxes[f]) > 0:
                ecoli_id = f[:-len(id_prefix_ecoli)]
                # save reaction ID
                ecoli_rxns.append(ecoli_id)
                # save flux
                ecoli_rxns_flux.append(str(round(sol.fluxes[f], 4)))
                # save subsystem
                ecoli_subsystems.append(iML1515.reactions.get_by_id(ecoli_id).subsystem)
                # save related genes
                for g in iML1515.reactions.get_by_id(ecoli_id).genes:
                    ecoli_genes.append(g.id)
                
                
        iMod_res.loc[c, "E. coli connect reactions"] = ", ".join(ecoli_rxns)
        iMod_res.loc[c, "E. coli connect reaction fluxes"] = ", ".join(ecoli_rxns_flux)
        iMod_res.loc[c, "E. coli connect subsystems"] = ", ".join(list(np.unique(ecoli_subsystems)))
        iMod_res.loc[c, "E.coli connect genes"] = ", ".join(list(np.unique(ecoli_genes)))
        
        c += 1
        


# %% check essentiality of E. coli reactions and producce simulation data for minimal reaction sets

writer = pd.ExcelWriter(filename, engine="xlsxwriter")

obj_val = []
essential_reactions = {}
iMod_info_save = {}
max_growth_solutions = {}
for i, iMod_AA in iMod_res.iterrows():
    
    # get heterologous E. coli reactions to connect iModulon
    rxns_connect = iMod_AA["E. coli connect reactions"].split(", ")
    # get AA 
    AA_syn = iMod_AA["Amino acid"]
    
    with model_p:
        # add iModulon
        iMod_rxns_unique = add_iModulon_to_model(model_p, g_data, iMod_AA["iModulon"])
        
        # disable amino acid uptake (or respective di-peptide)     
        model_p = disable_amino_acid_uptake(model_p, AA_syn[:-2])
        # disable amino acid uptake (or respective di-peptide)
        # for ex_rxn in model_p.exchanges:
        #     if AA_syn[:-2] in ex_rxn.id:
        #         ex_rxn.bounds = (0, 0)
        #         pass
        #     # consider misspelling of ala__L in di-peptide
        #     if AA_syn[:-2] == "ala__L":
        #         if "ala_L" in ex_rxn.id:
        #             ex_rxn.bounds = (0, 0)
        #             pass
         
        # create amino acid exchange
        aa_met = model_p.metabolites.get_by_id(AA_syn)
        ex_aa_id = "EX_" + AA_syn
        ex_aa = cobra.Reaction(id=ex_aa_id,
                               lower_bound=0,
                               upper_bound=1000)
        ex_aa.add_metabolites({aa_met: -1})
        model_p.add_reaction(ex_aa)
        
        # add heterologous E. coli reactions
        ecoli_rxns = []
        if rxns_connect[0]:          
            for rxnID in rxns_connect:
                # copy iML1515 reaction
                rxn_ecoli = iML1515.reactions.get_by_id(rxnID)
                rxn_copy = create_rxn_from_object(rxnID + id_prefix_ecoli, rxn_ecoli)
                # add to model
                model_p.add_reaction(rxn_copy)
                ecoli_rxns.append(rxn_copy.id)

            
        
        # test essentiality of E. coli reactions
        model_p.reactions.get_by_id(bio).lower_bound = growth_WT*bio_min_flux   # enforce biomass formation
        model_p.objective = ex_aa_id
        is_essential = []
        dict_key = iMod_AA["iModulon"] + "; " + AA_syn
        essential_reactions[dict_key] = []
        for ecoli_rxn in ecoli_rxns:
            # save bounds
            bounds_save = model_p.reactions.get_by_id(ecoli_rxn).bounds
            # delete reactions
            model_p.reactions.get_by_id(ecoli_rxn).bounds = (0, 0)
            # calculate AA sythesis rate
            sol_AA = model_p.slim_optimize(error_value=-1)
            if sol_AA <= 0:
                # reaction is essential for AA synthesis
                is_essential.append(True)
                essential_reactions[dict_key].append(ecoli_rxn)
            else:
                is_essential.append(False)
                continue

            # restore reaction bounds
            model_p.reactions.get_by_id(ecoli_rxn).bounds = bounds_save
            
        model_p.reactions.get_by_id(bio).lower_bound = 0   # relax biomass formation  constraints gain 
            
        
        # get maximum flux towards AA
        sol_max_aa = model_p.optimize()
        
        # calculate maximum growth rate
        model_p.objective = bio
        model_p.reactions.get_by_id(bio).lower_bound = 0
        model_p.reactions.get_by_id(bio).upper_bound = 1000
        sol_max_mu = model_p.optimize()
      
        # save model
        model_id = (iMod_AA["iModulon"] + "; " + AA_syn).replace("/", "_")
        cobra.io.save_json_model(model_p, "Models/" + model_id + ".json")
        # max_growth_solutions[(iMod_AA["iModulon"] + "; " + AA_syn).replace("/", "_")]\
        #     = model_p.copy()
        
        # test if auxotrophy of other AAs are lifted at the same time
        aa_cosynth = {}
        for amino_acid in AA_IDs:
            if amino_acid != AA_syn[:-2] and amino_acid != "cys__L":
                with model_p:
                    # disable amino acid uptake
                    model_p = disable_amino_acid_uptake(model_p, amino_acid)
                    # optimize for growth
                    sol_aa_syn = model_p.slim_optimize(error_value=-1)
                    # is growth possible?
                    if sol_aa_syn > 1e-2:
                        aa_cosynth[amino_acid] = str(round(sol_aa_syn, 3))
        
        
      
    # make extra sheet with additional information
    iMod_info = pd.DataFrame(columns=["Origin",
                                      "Reaction ID",
                                      "Reaction name",
                                      "Genes",
                                      "Reaction flux (max yield) [mmol/gDW/h]",
                                      "Reaction flux (max growth) [mmol/gDW/h]",
                                      "Reaction formula",
                                     ])
    c = 0
        
    # update table 
    ecoli_rxns = []
    ecoli_subsystems = []
    ecoli_genes_total = []
    ecoli_rxns_flux = []
    for r in essential_reactions[dict_key]:
        ecoli_id = r[:-len(id_prefix_ecoli)]
        # save reaction ID
        ecoli_rxns.append(ecoli_id)
        # save flux
        ecoli_rxns_flux.append(str(sol_max_aa.fluxes[r]))
        # save subsystem
        ecoli_subsystem = iML1515.reactions.get_by_id(ecoli_id).subsystem
        ecoli_subsystems.append(ecoli_subsystem)
        
        # save related genes
        ecoli_genes = []
        for g in iML1515.reactions.get_by_id(ecoli_id).genes:
            ecoli_genes.append(g.id)
            ecoli_genes_total.append(g.id)
            
        # fill extra sheet
        # heterologous E. coli reactions
        iMod_info.loc[c, "Origin"] = "E. coli heterologous"
        iMod_info.loc[c, "Reaction ID"] = ecoli_id
        iMod_info.loc[c, "Reaction name"] = iML1515.reactions.get_by_id(ecoli_id).name
        iMod_info.loc[c, "Genes"] = ", ".join(list(np.unique(ecoli_genes)))
        iMod_info.loc[c, "Reaction flux (max yield) [mmol/gDW/h]"] = round(sol_max_aa.fluxes[r], 4)
        iMod_info.loc[c, "Reaction flux (max growth) [mmol/gDW/h]"] = round(sol_max_mu.fluxes[r], 3)
        iMod_info.loc[c, "Reaction formula"] = iML1515.reactions.get_by_id(ecoli_id).build_reaction_string()
        c += 1
     
        
     
    # update summary  
    iMod_AA.loc["E. coli connect reactions"] = ", ".join(ecoli_rxns)
    iMod_AA.loc["E. coli connect reaction fluxes"] = ", ".join(ecoli_rxns_flux)
    iMod_AA.loc["E. coli connect subsystems"] = ", ".join(list(np.unique(ecoli_subsystems)))
    iMod_AA.loc["E.coli connect genes"] = ", ".join(list(np.unique(ecoli_genes_total)))
    iMod_AA.loc["Amino acid synthesis rate [mmol/gDW/h]"] = round(sol_max_aa.objective_value, 4)
    iMod_AA.loc["Maximum growth rate [1/h]"] = round(sol_max_mu.objective_value, 3)
    
    # add co-synthesis of AAs
    if len(aa_cosynth) > 0:
        iMod_AA.loc["Amino acid co-synthesis"] = ", ".join(list(aa_cosynth.keys()))
        iMod_AA.loc["Amino acid co-synthesis growth [1/h]"] = ", ".join(list(aa_cosynth.values()))
        
     
    # fill extra sheet with iModulon reaction data
    for r in iMod_rxns_unique:
        # get ID
        ecoli_id = r[:-len(id_prefix_iMod)]
        # save related genes
        ecoli_genes = []
        for g in iML1515.reactions.get_by_id(ecoli_id).genes:
            ecoli_genes.append(g.id)
        
        iMod_info.loc[c, "Origin"] = "iModulon"
        iMod_info.loc[c, "Reaction ID"] = ecoli_id
        iMod_info.loc[c, "Reaction name"] = iML1515.reactions.get_by_id(ecoli_id).name
        iMod_info.loc[c, "Genes"] = ", ".join(list(np.unique(ecoli_genes)))
        iMod_info.loc[c, "Reaction flux (max yield) [mmol/gDW/h]"] = round(sol_max_aa.fluxes[r], 4)
        iMod_info.loc[c, "Reaction flux (max growth) [mmol/gDW/h]"] = round(sol_max_mu.fluxes[r], 3)
        iMod_info.loc[c, "Reaction formula"] = iML1515.reactions.get_by_id(ecoli_id).build_reaction_string()
        c += 1
        
      
    iMod_info_save[(iMod_AA["iModulon"] + "; " + AA_syn).replace("/", "_")] = iMod_info
                
    
    
    


# %% save data frame
iMod_res.to_excel(writer,
                  sheet_name="Summary", 
                  index=False, )
for key, table in iMod_info_save.items():
    table.to_excel(writer,
                  sheet_name=key, 
                  index=False)  
writer.save()

