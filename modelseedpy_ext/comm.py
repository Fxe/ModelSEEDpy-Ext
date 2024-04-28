from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Tuple, Union
from cobra.core import Model, Reaction, Metabolite


class Comm(Model):

    def __init__(self, id_or_model: Union[str, "Model", None] = None, name: Optional[str] = None):
        super().__init__(id_or_model, name)


class CommFactory:

    def __init__(self):
        self.models = {}

    def with_model(self, model, abundance, index):
        if index in self.models:
            raise ValueError(f'Invalid index {index}. Already taken')
        self.models[index] = (model, abundance)
        return self

    def build_comm_biomass(self):
        r_bio_sum = Reaction('bio1', 'bio_com', '', 0, 1000)
        r_bio_sum.add_metabolites({
            model_com.metabolites.cpd11416_cA: -0.4,
            model_com.metabolites.cpd11416_cR: -0.6
        })
        model_com.add_reactions([r_bio_sum])

        model_com.objective = 'bio1'
        model_com.reactions.bio1.annotation = {'sbo': 'SBO:0000629'}
        model_com.reactions.bio1_A.annotation = {'sbo': 'SBO:0000629'}
        model_com.reactions.bio1_R.annotation = {'sbo': 'SBO:0000629'}

    def build(self):
        model = Comm('model')

        m_pointer = {}

        ex_ids = set()
        for model in models.values():
            ex_ids |= {list(r.metabolites)[0].id for r in model.exchanges}

        for token, model in models.items():
            m_mapper = {}
            for m in model.metabolites:
                if m.id in ex_ids:
                    # print(m)
                    if m.id not in ex:
                        m_copy = Metabolite(m.id, m.formula, m.name, m.charge, m.compartment)
                        m_mapper[m.id] = m_copy
                else:
                    m_id_copy = m.id[:-1] + token
                    # print(m_id_copy, m.id)
                    m_copy = Metabolite(m_id_copy, m.formula, m.name[:-1] + token, m.charge, m.compartment[:-1] + token)
                    m_mapper[m.id] = m_copy
            m_pointer[token] = m_mapper

        for token in models:
            model.add_metabolites(list(m_pointer[token].values()))

        r_pointer = {}
        for token, model in models.items():
            r_mapper = {}
            for r in model.reactions:
                if not r.id.startswith('EX_'):
                    if r.id.startswith('bio'):
                        r_copy_id = r.id + '_' + token
                    else:
                        r_copy_id = r.id[:-1] + token
                    r_copy = Reaction(r_copy_id, r.name[:-1] + token, r.subsystem, r.lower_bound, r.upper_bound)
                    r_copy_m = dict(map(lambda x: (model_com.metabolites.get_by_id(m_pointer[token][x[0].id].id), x[1]),
                                        r.metabolites.items()))
                    r_copy.add_metabolites(r_copy_m)
                    r_mapper[r_copy.id] = r_copy
                else:
                    # print(r)
                    pass
            r_pointer[token] = r_mapper

        for token in models:
            model_com.add_reactions(list(r_pointer[token].values()))

        return model