import efmtool
import numpy as np

def get_solutions(efms, reactions):
    t = efms.transpose()
    res = []
    for efm in t:
        flux_v = {}
        for i in range(len(reactions)):
            rxn_id = reactions[i]
            v = efm[i]
            if v != 0:
                flux_v[rxn_id] = v
        res.append(flux_v)
    return res

def efm(model):
    metabolites = list(map(lambda x: x.id, model.metabolites))
    reactions = list(map(lambda x: x.id, model.reactions))
    
    s_array = []
    rev = []
    for r in reactions:
        r_array = []
        rxn = model.reactions.get_by_id(r)
        rxn_metabolites = dict(map(lambda x: (x[0].id, x[1]), rxn.metabolites.items()))
        for m in metabolites:
            v = 0
            if m in rxn_metabolites:
                v = rxn_metabolites[m]
            r_array.append(v)
        s_array.append(r_array)
        rev.append(1 if rxn.lower_bound < 0 and rxn.upper_bound > 0 else 0)
    S = np.array(s_array)
    S = S.transpose()
    
    efms = efmtool.calculate_efms(S, rev, reactions, metabolites)
    
    return get_solutions(efms, reactions)
    