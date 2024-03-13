import logging
import efmtool
import numpy as np

logger = logging.getLogger(__name__)


def get_solutions(efms, reactions):
    t = efms.transpose()
    res = []
    for fm in t:
        flux_v = {}
        for i in range(len(reactions)):
            rxn_id = reactions[i]
            v = fm[i]
            if v != 0:
                flux_v[rxn_id] = v
        res.append(flux_v)
    return res


def efm(model):
    metabolites = list(map(lambda x: x.id, model.metabolites))
    reactions = []
    for r in model.reactions:
        lb = r.lower_bound
        ub = r.upper_bound
        if lb == ub and lb == 0:
            logger.warning(f'ignore zero bound reaction {r} [{lb}, {ub}]')
        else:
            reactions.append(r.id)

    s_array = []
    rev = []
    flipped = set()
    for r in reactions:
        rxn = model.reactions.get_by_id(r)
        lb = rxn.lower_bound
        ub = rxn.upper_bound

        rxn_metabolites = dict(map(lambda x: (x[0].id, x[1]), rxn.metabolites.items()))

        if lb < 0 and ub == 0:  # flip if reverse reaction
            flipped.add(rxn.id)
            rxn_metabolites = {k: v * -1 for k, v in rxn_metabolites.items()}

        r_array = [rxn_metabolites.get(m, 0) for m in metabolites]

        """
        # old stuff
        for m in metabolites:
            v = 0
            if m in rxn_metabolites:
                v = rxn_metabolites[m]
            r_array.append(v)
        """
        s_array.append(r_array)
        rev.append(1 if lb < 0 < ub else 0)  # 1 if reversible else 0
    S = np.array(s_array)
    S = S.transpose()

    res = efmtool.calculate_efms(S, rev, reactions, metabolites)

    return get_solutions(res, reactions)
