'''
new:
Volume implementation
simple diffusion algorithm

Adding more commentary
'''
print('program started')
print('loading modules')
import numpy
#import cProfile  as profile
import itertools
import numpy as np
import os
import re
import gurobipy
import libsbml
from copy import deepcopy
import scipy.sparse as sparse
#from openpyxl import load_workbook
import matplotlib.pyplot as plt
#from multiple_cellsV2.py import *
INF = float('inf')
NAN = float('nan')
print('loading modules completed')


####################################################################################
def run_model(environment,sbml,exchange_reaction_list):
    '''Runs some tests on the recon_2.2 model that is 
       stored in the folder d.models'''
    
    v, f_opt,exchange_list=max_fluxes(sbml,environment,exchange_reaction_list)
    return exchange_list,f_opt
def get_formula(sID, sbml):
    '''
    Get the formula of species with ID sID
    '''
    return get_notes_field(sID, 'FORMULA', sbml)

def get_notes_field(eID, name, sbml):
    '''
    Gets the notes field.
    '''
    element = sbml.getModel().getElementBySId(eID)
    notes = element.getNotesString()
    f = re.search(name + ':([^<]+)', notes)
    return f.group(1).strip() if f is not None else ''

def get_source_reactions(sbml):
    '''Determine source and sink reactions'''
    model = sbml.getModel()
    
    rID_list = []

    # strip out format used in recon 2.1
    species = model.getSpecies('M_carbon_e')
    if species:
        species.setBoundaryCondition(True)

    for reaction in model.getListOfReactions():
        nS, nP = 0, 0
        #outs=[]
        for reactant in reaction.getListOfReactants():
            sID = reactant.getSpecies()
            #print('sID: ,',sID)
            try:
                if not model.getSpecies(sID).getBoundaryCondition():
                    nS += 1
            except:
                print('reaction: ', reaction.getId())
                print('sID reactant: ',sID)
                print(model.getSpecies(sID).getBoundaryCondition())
        for product in reaction.getListOfProducts():
            sID = product.getSpecies()
            db1=model.getSpecies(sID)
            try:
                model.getSpecies(sID).getBoundaryCondition()
            except:
                print('reaction: ', reaction.getId())
                print('sID product: ',sID)
                print(model.getSpecies(sID).getBoundaryCondition())
            if not model.getSpecies(sID).getBoundaryCondition():
                nP += 1
        if (nS == 0) or (nP == 0):
            rID_list.append(reaction.getId())
    return rID_list


def max_fluxes(sbml,environment,exchange_reaction_list):
    '''
    Written to mimic neilswainston matlab function maxFluxes
    '''
    #get the variables to maximize (objective) and their weight (objective_coeff).
    objective,objective_coeff=read_objective()
    #dictionary containing the maximum flux for the nutrients in the extracellular matrix.
    #The 'evironment' variables are passed down from the tissue level model.
    carbon_sources={
                    'EX_glc(e)':environment[1],
                    'R_EX_lac_L_LPAREN_e_RPAREN_':environment[2],
                     'EX_o2(e)':environment[0]
                    }
    #A list of the freely available nutrients. These are not moddelled in the tissue model.
    media=[
        'EX_ca2(e)',
        'EX_cl(e)',
        'EX_fe2(e)',
        'EX_fe3(e)',
        'EX_h(e)',
        'EX_h2o(e)',
        'EX_k(e)',
        'EX_na1(e)',
        'EX_nh4(e)',
        'EX_so4(e)',
        'EX_pi(e)']
    model = sbml.getModel()
    listOfNodes=[]
    listOfEdges=[]
    
    
    #create list of cell uptake and secretion reactions
    
    normoxic=False
    #create a list of the exchange reactions that are allowed.
    allowed_exchanges=read_allowed_exchanges()
    #And block every exchange reaction that is not on the allowed list.
    block_illegal_exchanges(sbml,exchange_reaction_list,allowed_exchanges)
    #Run the model. f_opt is the value reached by the objective function. v is a list of all the reactions and their respective fluxes.
    v, f_opt = max_flux(sbml, carbon_sources, objective, objective_coeff, normoxic, media)
    print('objective value: ',f_opt)
    print('----')
##    print ('%s (%s):\t%g' % (carbon_source,
##                            'normoxic' if normoxic else 'anaerobic',
##                            f_opt))
    #write results to files and print some
    # legend=open('legend.txt', 'w')
    interaction=open('exchange.txt','w')

    #create list of fluxes that need to be returned to the tissue model
    exchange_list=['oxygen','glucose','lactate']
    #Iterate through the list. The important reactions, the non-zero exchange reactions and the reactions to be returned to the tissue model are printed here
    for j, reaction in enumerate(model.getListOfReactions()):
        #Get the reaction ID of the current item of the list.
        ID=reaction.getId()
##        legend.write(str(j)+" "+ID)
##        legend.write("\n")
        if 'R_DM_atp_c_'==ID:
            result_atp=v[j]
            print('result_atp1: ', result_atp)
            exchange_reaction=ID+': '+str(v[j])
##            interaction.write(exchange_reaction)
##            interaction.write("\n")
        elif 'R_biomass_reaction2' in ID:
            result_biomass2=v[j]
            print('result_biomass: ',result_biomass2)
            exchange_reaction=ID+': '+str(v[j])
##            interaction.write(exchange_reaction)
##            interaction.write("\n")
        elif 'R_DM_custombiomass'==ID:
            result_biomass=v[j]
            print('result_biomass: ',result_biomass)
            exchange_reaction=ID+': '+str(v[j])
##            interaction.write(exchange_reaction)
##            interaction.write("\n")
        elif 'R_EX_lac_D_LPAREN_e_RPAREN_'==ID:
            
            result_lacD=v[j]
            print('result_lac (D): ',result_lacD)
            exchange_reaction=ID+': '+str(v[j])
##            interaction.write(exchange_reaction)
##            interaction.write("\n")
        elif 'R_L_LACt2r'==ID:
            result_lacL=v[j]
            print('lactate(L): ',result_lacL)
            exchange_reaction=ID+': '+str(v[j])
##            interaction.write(exchange_reaction)
##            interaction.write("\n")
        elif 'R_EX_lac_L_LPAREN_e_RPAREN_'==ID:
            result_lacL=v[j]
            exchange_list[2]=result_lacL
            print('result_lac (L): ',result_lacL)
            exchange_reaction=ID+': '+str(v[j])
##            interaction.write(exchange_reaction)
##            interaction.write("\n")
        elif 'R_HEX1'==ID:
            result_lacL=v[j]
            print('R_HEX1: ',result_lacL)
            exchange_reaction=ID+': '+str(v[j])
##            interaction.write(exchange_reaction)
##            interaction.write("\n")
##        elif 'R_GLCt4'==ID:
##            result_lacL=v[j]
##            print('C0_glucose: ',result_lacL)
##            exchange_reaction=ID+': '+str(v[j])
##            interaction.write(exchange_reaction)
##            interaction.write("\n") 
        elif 'R_EX_glc_LPAREN_e_RPAREN_' in ID:
            exchange_list[1]=v[j]
            result_glc=v[j]
            print('result_glc: ',result_glc)
            exchange_reaction=ID+': '+str(v[j])
##            interaction.write(exchange_reaction)
##            interaction.write("\n")
        elif 'R_EX_o2_LPAREN_e_RPAREN_'==ID:
            exchange_list[0]=v[j]
            print("EX_oxygen: ",v[j])
            exchange_reaction=ID+': '+str(v[j])
##            interaction.write(exchange_reaction)
##            interaction.write("\n")
        elif ID in exchange_list and v[j]!=0:
            print(ID,":",v[j])
            exchange_reaction=ID+': '+str(v[j])
##            interaction.write(exchange_reaction)
##            interaction.write("\n")
        elif '_EX_' in ID and v[j]!=0:
            print(ID,": ",v[j])
            exchange_reaction=ID+': '+str(v[j])
##            interaction.write(exchange_reaction)
##            interaction.write("\n")
        elif '_DM_' in ID and v[j]!=0:
            print(ID,": ",v[j])
            exchange_reaction=ID+': '+str(v[j])
##            interaction.write(exchange_reaction)
##            interaction.write("\n")
         
        # if abs(v[j])>0: #vanaf welke waarde de resulaten worden geprint (1e-10)
        #      rID = reaction.getId()
        #      txt, txt_formula, innodes, outnodes=mydisplay_reaction_and_formula_Nodes(rID, sbml, v[j]>0)
        #      reactionnode=format_for_ID_SBML(rID)
        #      if reactionnode not in listOfNodes:
        #         listOfNodes.append(reactionnode)
        #      for n in innodes+outnodes:
        #          if not n in listOfNodes: listOfNodes.append(n)
        #      for n in innodes:
        #          if not (n, reactionnode, v[j]) in listOfEdges:
        #             listOfEdges.append((n, reactionnode, str(v[j])))
        #      for n in outnodes:
        #          if not (reactionnode, n) in listOfEdges:
        #             listOfEdges.append((reactionnode, n, str(v[j])))
        #  #print ('%s\t%s\t%g\t[%s,\t%s]' % ('reaction', reactionnode, abs(v[j]), txt, txt_formula))
    # legend.close()
    interaction.close()
    # nodef=open('tumor.txt', 'w')
    # nodef.write("\n".join(listOfNodes))
    # nodef.write("\n")
    # nodef.write("\n".join([";".join(r) for r in listOfEdges]))
    # nodef.close()
    if f_opt=='inf' or f_opt==0:
        print('simulation failed, results unreliable.')
    return v, f_opt,exchange_list

def genTxtAndTxtFormulaAndNodes(sbml, l=[]):
    txt = ''
    txt_formula = ''
    nodesl=[]
    for index, reactant in enumerate(l):
        if index > 0:
            txt += '+ '
            txt_formula += '+ '
        stoich = reactant.getStoichiometry()
        if stoich == int(stoich):
            stoich = int(stoich)
        if stoich != 1:
            txt += str(stoich) + ' '
            txt_formula += str(stoich) + ' '
        sID = reactant.getSpecies()
        txt += sID
        formula = get_formula(sID, sbml)
        nodesl.append(sID)
        if not formula:
            formula = '?'
        txt_formula += formula
        len_txt = len(txt)
        len_txt_formula = len(txt_formula)
        if len_txt > len_txt_formula:
            txt_formula = txt_formula.ljust(len_txt)
        elif len_txt < len_txt_formula:
            txt = txt.ljust(len_txt_formula)
        txt += ' '
        txt_formula += ' '
    return txt, txt_formula, nodesl
def mydisplay_reaction_and_formula_Nodes(rID, sbml, forward=True):
    '''
    Display reaction, with formulae below each reactant.
    '''
    model = sbml.getModel()
    reaction = model.getReaction(rID)

    itxt = '--> '
    itxt_formula = '--> '

    rtxt, rtxt_formula, rnodes = genTxtAndTxtFormulaAndNodes(sbml, reaction.getListOfReactants())
    ptxt, ptxt_formula, pnodes = genTxtAndTxtFormulaAndNodes(sbml, reaction.getListOfProducts())
    if forward:
        txt = rtxt + itxt + ptxt
        txt_formula = rtxt_formula + itxt_formula + ptxt_formula
        innodes = rnodes
        outnodes = pnodes
    else:
        txt = ptxt + itxt + rtxt
        txt_formula = ptxt_formula + itxt_formula + rtxt_formula
        outnodes = rnodes
        innodes = pnodes
    return txt, txt_formula, innodes, outnodes


def max_flux(sbml, carbon_sources, objective, objective_coeff, normoxic, media):
    '''
    Written to mimic neilswainston matlab function maxFlux
    '''
    model=sbml.getModel()
    #set_infinite_bounds(sbml)
    # block import reactions
    block_all_imports(sbml)
    # define carbon sources/import possibilities
    for carbon_source in carbon_sources:
        cs_amount=carbon_sources[carbon_source]
        set_import_bounds(sbml, carbon_source, cs_amount)
    # define media (unlimited import reactions)
    set_import_bounds(sbml, media, INF)
    if normoxic:
        set_import_bounds(sbml, 'EX_o2(e)', INF) #allow unlimited oxygen
##    set_import_bounds(sbml, 'R_h2o(e)', INF) #allow unlimited water
    # specify objective and set it to maximise
    change_objective(sbml, objective,objective_coeff)

    # avoid infinities
    # set a limit to the maximum value of the objective function
    obj_max = 1e19+1
    change_rxn_bounds(sbml, objective, obj_max, 'u')
    #run the model
    v, f_opt = optimize_cobra_model(sbml)
    #if the value of the objective function is high enough, round to infinity
    if f_opt > 0.9 * obj_max:
        f_opt = INF
    return v, f_opt


def read_sbml(filename):
    '''
    Read an SBML file from specified path.

    '''
    reader = libsbml.SBMLReader()
    sbml = reader.readSBMLFromFile(filename)
    return sbml


def block_all_imports(sbml):
    '''
    Written to mimic neilswainston matlab function blockAllImports
    '''
    model = sbml.getModel()
    #cycle through the list of reactions.
    for rID in get_source_reactions(sbml):
        reaction = model.getReaction(rID)
        nR, nP = 0, 0
        #nR will become the number of reactants that not boundary conditions
        for reactant in reaction.getListOfReactants():
            sID = reactant.getSpecies()
            if not model.getSpecies(sID).getBoundaryCondition():
                nR += 1
        #nP will become the number of products that are not boundary conditions        
        for product in reaction.getListOfProducts():
            sID = product.getSpecies()
            if not model.getSpecies(sID).getBoundaryCondition():
                nP += 1
        #The kineticLaw contains the maximum and minimum flux of a reactions
        kineticLaw = reaction.getKineticLaw()
        if (nR == 1) and (nP == 0):
            kineticLaw.getParameter('LOWER_BOUND').setValue(0)
        if (nR == 0) and (nP == 1):
            kineticLaw.getParameter('UPPER_BOUND').setValue(0)


def change_rxn_bounds(sbml, rxn_name_list, value, bound_type='b'):
    '''
    Written to mimic the matlab function changeRxnBounds from
    http://opencobra.sf.net/
    '''
    # convert single entries to lists
    if isinstance(rxn_name_list, str):
        rxn_name_list = [rxn_name_list]
    if isinstance(value, (int, float, complex)):
        value = [value] * len(rxn_name_list)
    if isinstance(bound_type, str):
        bound_type = [bound_type] * len(rxn_name_list)
    for index, rID in enumerate(rxn_name_list):
        reaction = get_reaction_by_id(sbml, rID)
        if not reaction:
            print ('reaction ',rID,' not found' % rID)
        else:
            kineticLaw = reaction.getKineticLaw()
            if bound_type[index] in ['l', 'b']:
                kineticLaw.getParameter('LOWER_BOUND').setValue(value[index])
            if bound_type[index] in ['u', 'b']:
                kineticLaw.getParameter('UPPER_BOUND').setValue(value[index])


def change_objective(sbml, rxn_name_list, objective_coeff=1):
    '''
    Written to mimic the matlab function changeObjective from
    http://opencobra.sf.net/
    '''
    model = sbml.getModel()
    for reaction in model.getListOfReactions():
        kineticLaw = reaction.getKineticLaw()
        kineticLaw.getParameter('OBJECTIVE_COEFFICIENT').setValue(0)
    # convert single entries to lists
    if isinstance(rxn_name_list, str):
        rxn_name_list = [rxn_name_list]
    if isinstance(objective_coeff, (int, float, complex)):
        objective_coeff = [objective_coeff] * len(rxn_name_list)
    for index, rID in enumerate(rxn_name_list):
        reaction = get_reaction_by_id(sbml, rID)
        if not reaction:
            print ('reaction %s not found' % rID)
        else:
            kineticLaw = reaction.getKineticLaw()
            #print(objective_coeff[index])
            kineticLaw.getParameter('OBJECTIVE_COEFFICIENT').setValue(
                objective_coeff[index])

def format_for_SBML_ID(txt):
    '''
    Written to mimic the matlab function formatForSBMLID from
    http://opencobra.sf.net/
    '''
    txt = 'R_' + txt
    for symbol, replacement in [
            ('-', '_DASH_'),
            ('/', '_FSLASH_'),
            ('\\', '_BSLASH_'),
            ('(', '_LPAREN_'),
            (')', '_RPAREN_'),
            ('[', '_LSQBKT_'),
            (']', '_RSQBKT_'),
            (',', '_COMMA_'),
            ('.', '_PERIOD_'),
            ('\'', '_APOS_'),
            ('&', '&amp'),
            ('<', '&lt'),
            ('>', '&gt'),
            ('"', '&quot')]:
        txt = txt.replace(symbol, replacement)
    return txt

def format_for_ID_SBML(txt):
    '''
    Written to mimic the matlab function formatForSBMLID from
    http://opencobra.sf.net/
    '''
    for symbol, replacement in [
            ('-', '_DASH_'),
            ('/', '_FSLASH_'),
            ('\\', '_BSLASH_'),
            ('(', '_LPAREN_'),
            (')', '_RPAREN_'),
            ('[', '_LSQBKT_'),
            (']', '_RSQBKT_'),
            (',', '_COMMA_'),
            ('.', '_PERIOD_'),
            ('\'', '_APOS_'),
            ('&', '&amp'),
            ('<', '&lt'),
            ('>', '&gt'),
            ('"', '&quot')]:
        txt = txt.replace(replacement, symbol)
    return txt


def optimize_cobra_model(sbml):
    '''
    Replicate Cobra command optimizeCbModel(model,[],'one').
    '''
    bound = INF
    cobra = convert_sbml_to_cobra(sbml, bound)

    N, L, U = cobra['S'], list(cobra['lb']), list(cobra['ub'])
    f, b = list(cobra['c']), list(cobra['b'])
    #run the model
    v_sol, f_opt, _ = easy_lp(f, N, b, L, U, sbml,one=False)
    return v_sol, f_opt


def convert_sbml_to_cobra(sbml, bound=INF):
    '''
    Get Cobra matrices from SBML model.
    '''
    model = sbml.getModel()
    #create empty (sparse) matrix
    S = sparse.lil_matrix((model.getNumSpecies(), model.getNumReactions()))
    #lower bound, upper bound, weight in objective function,b,reversible (true/false), speciesID
    lb, ub, c, b, rev, sIDs = [], [], [], [], [], []
    
    for species in model.getListOfSpecies():   
        b.append(0.)
    sIDs = [species.getId() for species in model.getListOfSpecies()]
    for j, reaction in enumerate(model.getListOfReactions()):
        kinetic_law = reaction.getKineticLaw()
        rxn_lb = kinetic_law.getParameter('LOWER_BOUND').getValue()
        rxn_ub = kinetic_law.getParameter('UPPER_BOUND').getValue()
        rxn_c = kinetic_law.getParameter('OBJECTIVE_COEFFICIENT').getValue()                 
        rxn_rev = reaction.getReversible()
        sIDs.append(species.getId())
        #heuristics + create sparse matrix S

        h_c_reactant=False
        h_m_product=False

        
        for reactant in reaction.getListOfReactants():
            sID = reactant.getSpecies()

            if 'M_h_i' in sID:
                h_c_reactant=True

            s = reactant.getStoichiometry()
            if not model.getSpecies(sID).getBoundaryCondition():
                i = sIDs.index(sID)
                S[i, j] = S[i, j] - s
        for product in reaction.getListOfProducts():
            sID = product.getSpecies()
            if 'M_atp_' in sID:
                atp_product=True
            if 'M_h_m' in sID:
                h_m_product=True
            #no carbon fixation
            s = product.getStoichiometry()
            if not model.getSpecies(sID).getBoundaryCondition():
                i = sIDs.index(sID)
                S[i, j] = S[i, j] + s
            
        
##        #remove the artificial constrains (that are only there to prevent overloading inferior solvers)
##        if rxn_lb==-1000:
##            rxn_lb=-INF
##        if rxn_ub==1000:
##            rxn_ub=INF
            
            
        if rxn_lb < -bound:
            rxn_lb = -bound
        if rxn_ub > bound:
            rxn_ub = bound
        if rxn_lb < 0:
            rxn_rev = True

        #heuristic proton block
        if h_c_reactant==True and h_m_product==True:
            if rxn_lb!=0:
                print('heuristic proton block: ',reaction.getId())
                rxn_lb=0
        # carbon fixation
##        if co2_product==True and co2_reactant==False and rxn_lb!=0:
##            rxn_lb=0
##            print('carbon fixation filter: ',reaction.getId())
        #annoyance filter
        if 'M_datp_c' in reaction.getListOfProducts():
            rxn_lb=0
            rxn_ub=0
        
        lb.append(rxn_lb)
        ub.append(rxn_ub)
        c.append(rxn_c)
        rev.append(rxn_rev)
    lb, ub, c, b = numpy.array(lb), numpy.array(ub), numpy.array(c), \
        numpy.array(b)
    rev = numpy.array(rev)
    cobra = {'S': S, 'lb': lb, 'ub': ub, 'c': c, 'b': b, 'rev': rev}
    return cobra


def easy_lp(f, a, b, vlb, vub, sbml, one=False):
    '''
    Optimize lp using Gurobi.
    '''
    # create gurobi model
    lp = gurobipy.Model()
    lp.Params.OutputFlag = 0
    lp.Params.FeasibilityTol = 1e-9  # as per Cobra
    lp.Params.OptimalityTol = 1e-9  # as per Cobra
    rows, cols = a.shape
    # add variables (the reactions) to model
    for j in range(cols):
        LB = vlb[j]
        if LB == -INF:
            LB = -gurobipy.GRB.INFINITY
        UB = vub[j]
        if UB == INF:
            UB = gurobipy.GRB.INFINITY
##        if f[j]!=0:
##            print(f[j])
        lp.addVar(lb=LB, ub=UB, obj=f[j])
    lp.update()
    lpvars = lp.getVars()
    # iterate over the rows of S adding each row into the model
    # thus adding each species to the model
    S = a.tocsr()
    for i in range(rows):
        start = S.indptr[i]
        end = S.indptr[i + 1]
        variables = [lpvars[j] for j in S.indices[start:end]]
        coeff = S.data[start:end]
        expr = gurobipy.LinExpr(coeff, variables)

        lp.addConstr(lhs=expr, sense=gurobipy.GRB.EQUAL, rhs=b[i]) #lijkt te vaak plaats te vinden, veel onaffe constraint in lp-file
    lp.update()

    add_constrains_from_file(sbml,lp,lpvars)
##    Add_solubiltiy_constraint(sbml,lpvars,lp)
    block_all_DM(sbml)
    lp.ModelSense = -1 #set: gurobi wil maximize objective function
##    lp.write("C:\\Users\\Jelte\\MEP\\recon.lp")
##    lp.write("C:\\Users\\Jelte\\MEP\\recon.mps")
    #run the model, for real now.
    lp.optimize()

    u = numpy.empty(len(f))
    v = numpy.empty(len(f))
    u[:] = NAN
    v[:] = NAN
    f_opt = NAN
    conv = False
    if lp.Status == gurobipy.GRB.OPTIMAL:
        f_opt = lp.ObjVal
        conv = True
        v = [var.x for var in lp.getVars()]


    if f_opt == -0.0:
        f_opt = 0.0

    return v, f_opt, conv


def get_reaction_by_id(sbml, rID):
    '''Gets the reaction by id.'''
    model = sbml.getModel()
    reaction = model.getReaction(rID)
    if not reaction:
        # try cobra replacements
        rID = format_for_SBML_ID(rID)
        reaction = model.getReaction(rID)
    if not reaction:
        # try removing trailing underscore
        if rID[-1] == '_':
            rID = rID[:-1]
        reaction = model.getReaction(rID)
    if not reaction:
        # try adding '_in'
        reaction = model.getReaction(rID + '_in')
    if not reaction:
        # try known alternatives
        rID_map = {
            'R_DM_atp_c': 'R_HKt',  # alternative ATPase
            # alternative C10:0
            'R_EX_HC02175_LPAREN_e_RPAREN': 'R_EX_dca_LPAREN_e_RPAREN_',
            # alternative C12:0
            'R_EX_HC02176_LPAREN_e_RPAREN': 'R_EX_ddca_LPAREN_e_RPAREN_',
            # alternative C22:0
            'R_EX_docosac': 'R_EX_docosac_LPAREN_e_RPAREN_',
        }
        if rID in rID_map:
            rID = rID_map[rID]
            reaction = get_reaction_by_id(sbml, rID)
    return reaction


def set_import_bounds(sbml, rxn_name_list, value):
    '''Sets the import bounds.'''
    model = sbml.getModel()
    # convert single entries to lists
    if isinstance(rxn_name_list, str):
        rxn_name_list = [rxn_name_list]
    if isinstance(value, (int, float, complex)): #long doesn't exsist in 3.6 anymore
        value = [value] * len(rxn_name_list)
    for index, rID in enumerate(rxn_name_list):
        reaction = get_reaction_by_id(sbml, rID)
        if not reaction:
            print ('reaction',rID,' not found' % rID)
        else:
            nR, nP = 0, 0
            for reactant in reaction.getListOfReactants():
                sID = reactant.getSpecies()
                if not model.getSpecies(sID).getBoundaryCondition():
                    nR += 1
            for product in reaction.getListOfProducts():
                sID = product.getSpecies()
                if not model.getSpecies(sID).getBoundaryCondition():
                    nP += 1
            kineticLaw = reaction.getKineticLaw()
            val = abs(value[index])
            if (nR == 0) and (nP == 1):
                kineticLaw.getParameter('UPPER_BOUND').setValue(val)
            elif (nR == 1) and (nP == 0):
                kineticLaw.getParameter('LOWER_BOUND').setValue(-val)
            else:
                print ('reaction %s not import' % rID)


def set_infinite_bounds(sbml):
    '''
    Set default bounds to INF, rather than 1000 (say)
    '''
    model = sbml.getModel()
    for reaction in model.getListOfReactions():
        kineticLaw = reaction.getKineticLaw()
        param = kineticLaw.getParameter('LOWER_BOUND')
        if param.getValue() < -100:
            param.setValue(-INF)
        param = kineticLaw.getParameter('UPPER_BOUND')
        if param.getValue() > 100:
            param.setValue(INF)

def ID_to_Var(sbml,lpvars,reacID):
    model=sbml.getModel()
    variable=None
    for R, reaction in enumerate(model.getListOfReactions()):
        if reaction.getId()==reacID:
            variable=lpvars[R]
            return variable
    print('error: ReactionID not in list: ',reacID)

def ID_to_Var2(sbml,lpvars,reacID):
    model=sbml.getModel()
    variable_list=[]
    for R, reaction in enumerate(model.getListOfReactions()):
        if reacID in reaction.getId():
            variable=lpvars[R]
            variable_list.append(variable)
    if len(variable_list)==0:
        print('error: ReactionID not in list: ',reacID)
    else:      
        return variable_list
    
    
def create_translation_list(sbml):
    model=sbml.getModel()
    LoR=model.getListOfReactions()
    namelist=[]
    for reaction in LoR:
        namelist.append(reaction.getId())
    legend=open('legend2.txt', 'w')
    
    for R, item in enumerate(namelist):
        regel_nummer=str(R)
        regel=regel_nummer+" "+item
        legend.write(regel)
        legend.write("\n")
    legend.close()

def reaction_finder(sbml):
    print('search reactions')
    model=sbml.getModel()
    
    for reaction in model.getListOfReactions():
        pID_match=False
        rID_match=False
        for reactant in reaction.getListOfReactants():
            rID=reactant.getSpecies()
            if rID=='M_glu_L_c' or rID=='M_akg_c':
                rID_match=True
        for product in reaction.getListOfProducts():
            pID=product.getSpecies()
            if pID=='M_glu_L_c' or pID=='M_akg_c':
                pID_match=True
        if pID_match==True and rID_match==True:
            print('match: ',reaction.getId())             
    return

def subsystems():
    #read
    print('loading file')
    with open("C:\\Users\\Jelte\\MEP\models\\recon2.2_rep6_fork1.xml","r") as file:
        bestand=file.readlines()
    print('Done')
    #sanitize
    print('sanitize and convert to list')
    regel_lijst=[]
    for line in bestand:
        new_line=line.replace("\n",'')
        regel_lijst.append(new_line)
    print('length list: ', len(regel_lijst))
    print('Done')
    new_line=None
    bestand=None
    file=None
    print('parsing file')
    zone=None
    inReaction=False
    lijst={}
    for line in regel_lijst:
        if '<listOfReactions>' in line:
            zone='reactions'
        elif zone=='reactions' and '<reaction' in line:
            found_subsystem=False
            inReaction=True
            #get Id
            start=line.find('id="')+4
            end=line.find('"',start)
            identifier=line[start:end]
            #print(identifier)

        elif inReaction==True and zone=='reactions' and 'SUBSYSTEM:' in line:
            start=line.find('SUBSYSTEM:')+11
            end=line.find('<',start)
            subsystem=line[start:end]
            found_subsystem=True
            
        elif inReaction==True and '</reaction>' in line:
            lijst[identifier]=subsystem
            #print(identifier,': ',subsystem)
            if found_subsystem==False:
                subsystem='other'
                
            inReaction=False
        else:
            pass
        
    #create overview of exsisting subsystems
    unieke_lijst=[]
    for value in lijst.values():
        if value not in unieke_lijst:
            unieke_lijst.append(value)
    subsystem_list=open('subsystem_list.txt','w')
    unieke_lijst.sort()
    for item in unieke_lijst:
        subsystem_list.write(item+'\n')
    subsystem_list.close()
    print('Done, script completed')
    return lijst

def add_constrains_from_file(sbml,lp,lpvars):
    '''
add constraints from file. Keep in mind that this function add constraints and does not overwrite any other constraints.
/ob     use to define reactions and their respective weights for the objective function
/block  block all reactions following
/ub     set the upperbound for a reaction ie. R_EX_glc=9 sets the maximum flux for R_EX_glc to 9
/lb     set the lowerbound for a reaction
/c      set a new constraint.
        /c
        Reaction_1=2
        Reaction_2=5
        90
        /c=
        generates the contraint Reaction_1*2+Reaction_2*5=90
        another example
        /c
        Reaction_1=-0.5
        Reaction_2=5
        10
        /c<
        generates the constraint Reaction_1*-0.5+Reaction_2<10
    '''
    print('aplying constraints from file')
    #read
    with open("constraints.txt","r") as file:
        bestand=file.readlines()

    #sanitize
    regel_lijst=[]
    for line in bestand:
        new_line=line.replace("\n",'')
        regel_lijst.append(new_line)
    #print(regel_lijst)
    new_line=None
    bestand=None
    file=None

    constr_nr=0
    soort=None
    for line_nr, line in enumerate(regel_lijst):
        #ignore commented lines
        if line[0:1]=='#':
            pass
        elif line[0:2]=="/c":
            #new constraint starts
            
            
            start_line_nr=line_nr
            if soort=='constrain' and variable_list!=[]:
                #deze \c sluit een constraint
                #save constraint
                expr = gurobipy.LinExpr(coeff_list, variable_list)
                if len(line)<3 or line[2:3]=="=":
                    lp.addConstr(lhs=expr, sense=gurobipy.GRB.EQUAL, rhs=const)
                elif line[2:3]==">":
                    lp.addConstr(lhs=expr, sense=gurobipy.GRB.GREATER_EQUAL, rhs=const)
                elif line[2:3]=="<":
                    lp.addConstr(lhs=expr, sense=gurobipy.GRB.LESS_EQUAL, rhs=const)
                else:
                    print('add constraint: sign unclear, constraint not implemented')
                #simulated by print commando's
                print(variable_list)
                print(coeff_list)
            variable_list=[]
            coeff_list=[]
            const=0
            constr_nr=constr_nr+1
            soort='constrain'
        elif line[0:3]=="/ub":
            soort='upperbound'
        elif line[0:3]=="/lb":
            soort='lowerbound'
        elif line[0:3]=="/ob":
            soort='objective'
        elif line[0:6]=='/block':
            soort='block'
        
        elif soort=='constrain' and "=" in line:
            #analyse new line
            var_half=True
            variable=''
            coeff=''
            for colomn_nr in range(len(line)):
                if line[colomn_nr:colomn_nr+1]=="=":
                    var_half=False
                elif var_half==True:
                    variable=variable+line[colomn_nr:colomn_nr+1]
                elif var_half==False:
                    coeff=coeff+line[colomn_nr:colomn_nr+1]
            #convert coeff string to number
            try:
                coeff=float(coeff)
            except:
                print('conversion to float failed')
            #save found variable and coefficient to the list
            Gvariable=ID_to_Var(sbml,lpvars,variable)
            if Gvariable!=None:
                variable_list.append(Gvariable)
                coeff_list.append(coeff)
        elif soort=='constrain' and "=" not in line:
            coeff=''
            for colomn_nr in range(len(line)):
                coeff=coeff+line[colomn_nr:colomn_nr+1]
            try:
                const=float(coeff)
            except:
                print('conversion to number failed, constant ignored')
        elif soort=='upperbound':
            var_half=True
            variable=''
            coeff=''
            for colomn_nr in range(len(line)):
                if line[colomn_nr:colomn_nr+1]=="=":
                    var_half=False
                elif var_half==True:
                    variable=variable+line[colomn_nr:colomn_nr+1]
                elif var_half==False:
                    coeff=coeff+line[colomn_nr:colomn_nr+1]
            try:
                coeff=float(coeff)
            except:
                print('conversion to float failed')
            Gvariable=ID_to_Var(sbml, lpvars, variable)
            if Gvariable!=None:
                lp.addConstr(lhs=Gvariable,sense=gurobipy.GRB.LESS_EQUAL,rhs=coeff)
                #print('upper bound set for ', variable, ' to ', coeff)

        elif soort=='lowerbound':
            var_half=True
            variable=''
            coeff=''
            for colomn_nr in range(len(line)):
                if line[colomn_nr:colomn_nr+1]=="=":
                    var_half=False
                elif var_half==True:
                    variable=variable+line[colomn_nr:colomn_nr+1]
                elif var_half==False:
                    coeff=coeff+line[colomn_nr:colomn_nr+1]
            try:
                coeff=float(coeff)
            except:
                print('conversion to float failed')
            Gvariable=ID_to_Var(sbml, lpvars, variable)
            if Gvariable!=None:
                lp.addConstr(lhs=Gvariable,sense=gurobipy.GRB.GREATER_EQUAL,rhs=coeff)
                print('lower bound set for ', variable, ' to ', coeff)
        elif soort=='block':
            variable=''
            for colomn_nr in range(len(line)):
                variable=variable+line[colomn_nr:colomn_nr+1]
            Gvariable_list=ID_to_Var2(sbml, lpvars, variable)
            for Gvariable in Gvariable_list:
                if Gvariable!=None:
                    lp.addConstr(lhs=Gvariable,sense=gurobipy.GRB.EQUAL,rhs=0)
                    #print('block: ',variable)
                else:
                    print('block failed: ',variable)
    lp.update()
    print('processing contraints file completed')
############ end constraints from file

def load_subsystem_file():
    with open("allowed_subsystem_list.txt","r") as file:
        bestand=file.readlines()

    #sanitize
    regel_lijst=[]
    for line in bestand:
        new_line=line.replace("\n",'')
        regel_lijst.append(new_line)
    #print(regel_lijst)
    new_line=None
    bestand=None
    file=None
    return regel_lijst

def compat_recon1(sbml):
    '''
    for some reason more important than expected. Originaly needed for recon1, now for recon2.2_rep6_fork1 as well. The reason for this is unclear.
    '''
    model=sbml.getModel()
    #check for the first few reaction if the kineticLaws are present in the model
    try:
        for k, reaction in enumerate(model.getListOfReactions()) and k<10:
            kinetic_law = reaction.getKineticLaw()
            rxn_lb = kinetic_law.getParameter('LOWER_BOUND').getValue()
        
    except:
        for k, reaction in enumerate(model.getListOfReactions()):
            #print(reaction.getId())
            name=reaction.getId()

            rev=reaction.getReversible()
            if reaction.getReversible()==True:
                kineticLaw=reaction.createKineticLaw()
                kl=kineticLaw.createParameter()
                kl.setId('LOWER_BOUND')
                kl.setValue(-INF)
                k2=kineticLaw.createParameter()
                k2.setId('UPPER_BOUND')
                k2.setValue(INF)
                k3=kineticLaw.createParameter()
                k3.setId('OBJECTIVE_COEFFICIENT')
                k3.setValue(0)
            elif reaction.getReversible()==False:
                kineticLaw=reaction.createKineticLaw()
                kl=kineticLaw.createParameter()
                kl.setId('LOWER_BOUND')
                kl.setValue(0)
                k2=kineticLaw.createParameter()
                k2.setId('UPPER_BOUND')
                k2.setValue(INF)
                k3=kineticLaw.createParameter()
                k3.setId('OBJECTIVE_COEFFICIENT')
                k3.setValue(0)


def read_objective():
    objective=[]
    objective_coeff=[]
    #read
    with open("constraints.txt","r") as file:
        bestand=file.readlines()

    #sanitize
    regel_lijst=[]
    for line in bestand:
        new_line=line.replace("\n",'')
        regel_lijst.append(new_line)
    new_line=None
    bestand=None
    file=None

    constr_nr=0
    soort=None
    for line_nr, line in enumerate(regel_lijst):
        if line[0:1]=='#':
            pass
        elif line[0:3]=="/ob" and soort!='objective':
            soort='objective'
        elif "/" in line:
            soort=None
        elif soort=='objective' and "/" not in line:
            name_half=True
            name=''
            coeff=''
            for colomn_nr in range(len(line)):
                if line[colomn_nr:colomn_nr+1]=="=":
                    name_half=False
                elif name_half==True:
                    name=name+line[colomn_nr:colomn_nr+1]
                elif name_half==False:
                    coeff=coeff+line[colomn_nr:colomn_nr+1]
            try:
                coeff=float(coeff)
            except:
                print('conversion to float failed')
            objective.append(str(name))
            objective_coeff.append(coeff)
    return objective,objective_coeff



def check_number(number):
    try:
        float(number)
        return True
    except:
        return False
def block_all_DM(sbml):
    '''blocks all demand reactions unless they are in the objective function'''
    model=sbml.getModel()
    for reaction in model.getListOfReactions():
        KL=reaction.getKineticLaw()
        if '_DM_' in reaction.getId() and KL.getParameter('OBJECTIVE_COEFFICIENT').getValue()!=0:
            KL.getParameter('LOWER_BOUND').setValue(0)
            KL.getParameter('UPPER_BOUND').setValue(0)

def remove_subsystems(sbml):
    #remove subsystems
    model=sbml.getModel()
    lijst=subsystems()
    system_list=load_subsystem_file()
    for reaction in model.getListOfReactions():
        if (lijst[reaction.getId()] not in system_list) and (reaction.getId() not in system_list):
            KL=reaction.getKineticLaw()
            KL.getParameter('LOWER_BOUND').setValue(0)
            KL.getParameter('UPPER_BOUND').setValue(0)


def test(sbml):
    #content to be added depending on what is required. Currently
    return

    
def read_allowed_exchanges():
    with open("allowed_exchanges.txt","r") as file:
        bestand=file.readlines()
    #sanitize
    regel_lijst=[]
    for line in bestand:
        new_line=line.replace("\n",'')
        regel_lijst.append(new_line)
    #print(regel_lijst)
    new_line=None
    bestand=None
    file=None
    return regel_lijst

def block_illegal_exchanges(sbml,exchange_list,allowed_exchanges):
    model=sbml.getModel()
    for reaction in model.getListOfReactions():
        if reaction.getId() in exchange_list and reaction.getId() not in allowed_exchanges:
            KL=reaction.getKineticLaw()
            KL.getParameter('LOWER_BOUND').setValue(0)
            KL.getParameter('UPPER_BOUND').setValue(0)
            
def preload_model():
    tests_path = os.path.dirname(__file__)
    model_path = os.path.join(tests_path, 'models')
    model_path = os.path.normpath(model_path)
    #name='recon_2.2_orig'
    #name="Recon2.v02_mod2"
    #name='recon2.2_rep6_fork1'
    name='recon2.2_minimal11'
    print (name)
    model_filename = os.path.join(model_path, name + '.xml')
    print(model_filename)
    sbml = read_sbml(model_filename)
    compat_recon1(sbml)
    #remove_subsystems(sbml)
    create_translation_list(sbml)
    return sbml

def generate_list_of_exchange_reactions(sbml):
    model=sbml.getModel()
    print('generating list of exchange reactions')
    exchange_reaction_list=[]
    for reaction in model.getListOfReactions():
        comp_c=False
        comp_e=False
        for reactant in  reaction.getListOfReactants():
            ID=reactant.getSpecies()
            species=model.getSpecies(ID)
            if species.getCompartment()=='c':
               comp_c=True
            if species.getCompartment()=='e':
               comp_e=True
        for product in  reaction.getListOfProducts():
            ID=product.getSpecies()
            species=model.getSpecies(ID)
            if species.getCompartment()=='c':
               comp_c=True
            if species.getCompartment()=='e':
               comp_e=True
        if comp_c==True and comp_e==True:
            exchange_reaction_list.append(reaction.getId())
    print('exchange list completed')
    return exchange_reaction_list

def create_volume_list(NoC,shape,normalize):
    Volume=np.ones(NoC+1)
    totalVolume=NoC
    if shape=='cilinder' or shape=='cylinder' or shape=='cillinder':
        #cilinder (ignore height)
        i=0
        totalVolume=0
        while i<NoC+1:
            Volume[[i]]=pi*((i+2)**2)-((i+1)**2)
            totalVolume=Volume[[i]]+totalVolume
            i=i+1
    elif shape=='sphere':
        #sphere, from in->out
        i=0
        while i<NoC+1:
            Volume[[i]]=(4/3)*pi*(((i+2)^3)-(i+1^3))
            totalVolume=Volume[[i]]+totalVolume
            i=i+1
    elif shape=='slab' or 'linear' or 'costom':
        pass
    #normalize
    if normalize==True:
        i=0
        while i<NoC+1:
            Volume[[i]]=Volume[[i]]/(totalVolume/NoC)
            i=i+1
    return Volume
def totalResistance(cellNumber,Rm,Rd,volume):
    #calculates the total restistance between the vessel and the cell
    Resistance=0
    i=0
    while i<=cellNumber:
        if i==0:
            Resistance+=(Rm/volume[i])
        else:
            if volume[i-1]>volume[i]:
                Resistance+=(Rd/volume[i-1])
            else:
                Resistance+=(Rd/volume[i])
        i+=1
    return Resistance
###########################################################################################
#Tissue Scale model
#Settings
title='Real values' #The title for the graphs
NoC=5 #Number of Cells/compartments (layers)
location=0
shape='slab'
normalize=False
Rd_ox=50000/1
Rd_glc=149253.7313/1
Rd_lac=1
Rm_ox= 100000
Rm_glc=100000
Rm_lac=1
I_max_ox=5000
I_max_glc=600000


I_max_lac=0
d=1
deadcells=[]
#deadcells=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52]
#constants
pi=3.1415926535897950288


#create starting matrices

MaxImatrix=np.zeros((NoC+1,3))
Imatrix=np.zeros((NoC+2,3))
matrix=np.zeros((NoC+1,3))
matrix[0,0]=I_max_ox
matrix[0,1]=I_max_glc
matrix[0,2]=I_max_lac
MaxImatrix[0,0]=I_max_ox
MaxImatrix[0,1]=I_max_glc
MaxImatrix[0,2]=I_max_lac
#write compartments size

normalize=False
volume=create_volume_list(NoC,shape,normalize)
print('volume vector: ',volume)
#user input
#volume=[1,1]
#append ones uptil correct length is reached.
while len(volume)<NoC+1:
    volume.append(1)
    
#pre-load model
sbml=preload_model()
exchange_reaction_list=generate_list_of_exchange_reactions(sbml)
#prepare for iteration
#environment=np.dtype('float64')
no_change=False
debug=0
stop=False
round_counter=0
#starting outerloop
print('starting outer loop')
while no_change==False and stop==False:
    #save current configuration
    matrix_old=deepcopy(matrix)
    #prepare for iteration
    i=0
    progress=True
    deadCell=False
    SumOx=0
    SumGlc=0
    SumLac=0
    #LostOx=0
    print('round: ',round_counter)
    print('--------------------------------------')
    #start iteration: each cell once
    while NoC>i and progress==True and deadCell==False:
        print('----------------------------------------')
        #create environment, the maximum fluxes available for the cell-model
        #resistance
        RtotOx=totalResistance(i,Rm_ox,Rd_ox,volume)
        RtotGlc=totalResistance(i,Rm_glc,Rd_glc,volume)
        ## RtotLac=Rtotglc=totalResistance(i,Rm_lac,Rd_lac,volume)
        #create input vector
        environment=[0,0,0]
        environment[0]=(matrix[i,0]*Rm_ox/(RtotOx))/volume[i]
        environment[1]=(matrix[i,1]*Rm_glc/(RtotGlc))/volume[i]
        environment[2]=(matrix[i+1,2])/volume[i]
        #save input vector
        MaxImatrix[i,0]=environment[0]
        MaxImatrix[i,1]=environment[1]
        MaxImatrix[i+1,2]=environment[2]
        print('environment: ',environment)
        #single cell model
        #dead cell
        if i in deadcells:
            exchange_list=np.zeros(3)
            f_opt=-1
            print('dead cell')
        else: #not dead cell
            print('run model. Cellnumber= ',i)
            exchange_list,f_opt=run_model(environment,sbml,exchange_reaction_list)
        if round_counter==2 and i==0:
            debug=True
        #check for dead cells
        if f_opt==0 or f_opt=='inf':
            deadCell=True
            exchange_list[0]=0
            exchange_list[1]=0
            exchange_list[2]=0
        #
        
        #volume correction
        exchange_list[0]=(exchange_list[0])*volume[[i]]
        exchange_list[1]=(exchange_list[1])*volume[[i]]
        exchange_list[2]=(exchange_list[2])*volume[[i]]
        #use results to calculate totals
        SumOx=(-1*exchange_list[0])+SumOx
        SumGlc=(-1*exchange_list[1])+SumGlc
        SumLac=(-1*exchange_list[2])+SumLac
        #save exchange in matrix
        Imatrix[i,0]=exchange_list[0]
        Imatrix[i,1]=exchange_list[1]
        Imatrix[i,2]=exchange_list[2]
        #debug
##        if environment[1]!=matrix[i,1]/Volume[[i]]:
##            stop=True
##            print('environment and matrix diverted')
##            raise
        
        #check if results are possible and correct rounding issues
        difference=matrix[i,0]+exchange_list[0]
        if difference<0:
            if difference<-1e-7-Dox:
                stop=True
                print('oxygen used highe than oxygen availabe')

        difference=matrix[i,1]+exchange_list[1]
        if difference<0:
            if difference<-1e-7-Dglc:
                stop=True
                print('glucose used higher than glucose available')

            
            #print('glucose used higher than glucose available')
        difference=matrix[i+1,2]+exchange_list[2]
        if difference<0:
            if difference<-1e-7:
                stop=True
                print('lactate used higher than lactate available')
            else:
                exchange_list[2]=-1*matrix[i+1,2]
                print('small lactate deviation corrected')
        difference=0      
        #debug
##        if round_counter==2 and i==38:
##            pass
        #apply modifications
        matrix[i+1,0]=environment[0]+exchange_list[0]
        matrix[i+1,1]=environment[1]+exchange_list[1]
        matrix[i,2]=environment[2]+exchange_list[2]
        #check if rerun is required
        print('SumOx: ',SumOx)
        
        print('totalOx: ', SumOx)
        #Debug: check for lower than zero values
        if matrix[i+1,0]<0:
            print('Oxygen flux under zero', matrix[i+1,0])
            matrix[i+1,0]=0
            stop=True
            raise ValueError('Oxygen flux under zero')
        if matrix[i+1,1]<0:
            print('Glycose flux onder zero ',matrix[i+1,1])
            difference=environment[1]+exchange_list[1]
            print(difference)
            raise ValueError('glucosoe flux under zero')
                       
            stop=True
        if matrix[i+1,2]<0:
            print('Lactate flux under zero')
            matrix[i+1,2]=0
            raise ValueError('Lactate flux under zero')
            stop=True
        #check whether a naturally dead cell was simulated
        if exchange_list[0]==0 and exchange_list[1]==0 and exchange_list[2]==0 and i not in deadcells:
            progress=False
            #move to next cell
        i=i+1
    #test for changes
    print(matrix)
    print(matrix_old)
    if np.array_equal(matrix,matrix_old):
        no_change=True
        diffs=matrix[:,:]-matrix_old[:,:]
        print(diffs)
    else:
        diffs=matrix[:,:]-matrix_old[:,:]
        print(diffs)                          
    print('--------------------------------------')
    #round counter
    round_counter=round_counter+1
    if round_counter>500:
        print('Maximum round count reached. Aborting simulation. Show most recent results')
        stop=True

#remove artefact for max fluxes plot
if deadCell==True and i<NoC:
    #Thus the artefact is present
    print('removing artefact')
    #cycle though the not simulated dead-cells
    while i<NoC:
        RtotOx=totalResistance(i,Rm_ox,Rd_ox,volume)
        RtotGlc=totalResistance(i,Rm_glc,Rd_glc,volume)
        ## RtotLac=Rtotglc=totalResistance(i,Rm_lac,Rd_lac,volume)
        #create input vector
        environment=[0,0,0]
        environment[0]=(matrix[i,0]*Rm_ox/(RtotOx))/volume[i]
        environment[1]=(matrix[i,1]*Rm_glc/(RtotGlc))/volume[i]
        environment[2]=(matrix[i+1,2])/volume[i]
        #save input vector
        MaxImatrix[i,0]=environment[0]
        MaxImatrix[i,1]=environment[1]
        MaxImatrix[i+1,2]=environment[2]
        #exchanges are zero, cells are dead
        matrix[i+1,0]=environment[0]
        matrix[i+1,1]=environment[1]
        matrix[i,2]=environment[2]
        i+=1
    pass
    

    
#wipe some variables
matrix_old=None
environment=None
exchange_list=None
#If needed print warning
if stop==True:
    print('Something went wrong during the simulation')
print('Simulation completed. Processing resuts for visualisation...')
#process matix for userfriendly output

LacEx=-SumLac
print('oxygen used: ',-SumOx)
print(matrix[:,0])
print('lactate used: ', SumLac)
#create x-axis
i=0
x_axis=[]

while i<NoC:
    x_axis.append(i)
    i=i+1
#x-axis completed
#plot fluxes through lumens
plt.plot(x_axis,MaxImatrix[0:NoC,0],'bD-')
plt.plot(x_axis,MaxImatrix[0:NoC,1],'gv-')
plt.plot(x_axis,MaxImatrix[1:NoC+1,2],'rx-')
plt.legend(['oxygen','glucose','lactate'])
plt.xlabel('Cell number')
plt.ylabel('Maximum possible flux through L')
text='tissue oxygen usage: '+ str(SumOx)
plt.figtext(0.6,0.9,text)
text='tissue lactate excretion: '+ str(LacEx)
plt.figtext(0.6,0.95,text)
plt.title(title)
plt.savefig('max_fluxes.eps')

plt.show()
plt.clf()
#create actual flux matrix
print(matrix)
AFmatrix=np.zeros((NoC+1,3))
x_axis.append(i*d)
#x_axis.append(0)
i=NoC-1
print(AFmatrix)
while i>-1:
##    print('i= ',i)
##    print(Imatrix[i,0])
##    print(AFmatrix[i+1,0])
    AFmatrix[i,0]=AFmatrix[i+1,0]-Imatrix[i,0]

    AFmatrix[i,1]=AFmatrix[i+1,1]-Imatrix[i,1]
    AFmatrix[i,2]=AFmatrix[i+1,2]+Imatrix[i,2]
    i=i-1
print(AFmatrix[:,0])
#plt.plot(x_axis,AFmatrix[:,0],'b-')
plt.plot(x_axis,AFmatrix[:,0],'bD-')
plt.plot(x_axis,AFmatrix[:,1],'gv-')
plt.plot(x_axis,AFmatrix[:,2],'rx-')
plt.legend(['oxygen','glucose','lactate'])
plt.xlabel('L number, distance')
plt.ylabel('Actual flux through L')
text='tissue oxygen usage: '+ str(SumOx)
plt.figtext(0.6,0.9,text)
text='tissue lactate excretion: '+ str(LacEx)
plt.figtext(0.6,0.95,text)
plt.title(title)
print(AFmatrix)
plt.savefig('actual_fluxes.eps')
plt.show()
print(title)
print('Done')

