class GSPComBuilder:
    """
    Assembles Community model for denitrification
    """

    def __init__(self, model_rhoda, model_acido, media_com):
        """

        :param model_rhoda: Model
        :param model_acido: Model
        :param media_com: Media
        """
        self.models = {model_rhoda: "R", model_acido: "A"}
        self.media_com = media_com

    def a(self, model_com):
        ## Replace to cytochromes
        model_com.reactions.rxn01806_cA.add_metabolites(
            {
                model_com.metabolites.cpd00004_cA: +1,
                model_com.metabolites.cpd00003_cA: -1,
                model_com.metabolites.cpd18072_cA: +1,
                model_com.metabolites.cpd18074_cA: -1,
            }
        )
        model_com.reactions.rxn01806_cR.add_metabolites(
            {
                model_com.metabolites.cpd00004_cR: +1,
                model_com.metabolites.cpd00003_cR: -1,
                model_com.metabolites.cpd18072_cR: +1,
                model_com.metabolites.cpd18074_cR: -1,
            }
        )

    def external_denitri(self, model):
        from cobra.core import Metabolite, Reaction

        print("external denitri")
        delete = {
            "rxn10577_cR",
            "rxn10577_cA",
            "rxn08966_cA",
            "rxn08966_eR",
            "rxn09008_cR",
            "rxn09008_cA",
            "rxn05890_cR",
            "rxn11932_cA",
            "rxn11932_cR",
            "rxn05627_cR",
            "rxn05627_cA",
            "rxn05625_cA",
        }
        model.remove_reactions(delete)
        model.add_metabolites([Metabolite("cpd00418_e0", "NO", "NO_e0")])
        rxnR = Reaction("rxn09004_cR", "rxn09004_cR", "", -1000, 1000)
        rxnR.add_metabolites(
            {
                model.metabolites.cpd00209_cR: 1,
                model.metabolites.cpd00075_cR: -1,
                model.metabolites.cpd00209_e0: -1,
                model.metabolites.cpd00075_e0: 1,
            }
        )
        rxnA = Reaction("rxn09004_cA", "rxn09004_cA", "", -1000, 1000)
        rxnA.add_metabolites(
            {
                model.metabolites.cpd00209_cA: 1,
                model.metabolites.cpd00075_cA: -1,
                model.metabolites.cpd00209_e0: -1,
                model.metabolites.cpd00075_e0: 1,
            }
        )
        ex = Reaction("EX_cpd00418_e0", "EX_cpd00418_e0", "", 0, 1000)
        ex.add_metabolites({model.metabolites.cpd00418_e0: -1})
        model.add_reactions([rxnR, rxnA, ex])
        move = {
            "cpd00528_cR": ["cpd00528_e0", None],
            "cpd00528_cA": ["cpd00528_e0", None],
            "cpd00659_cR": ["cpd00659_e0", None],
            "cpd00659_cA": ["cpd00659_e0", None],
            "cpd00418_cA": ["cpd00418_e0", None],
            "cpd00418_cR": ["cpd00418_e0", None],
            "cpd00075_cR": ["cpd00075_e0", {"rxn14428_cR"}],
            # 'cpd00075_cR': ['cpd00075_e0', {'rxn09001_cR', 'rxn14427_cR', 'rxn09003_cR', 'rxn14428_cR'}],
        }

        for src_cpd_id in move:
            dst_cpd_id = move[src_cpd_id][0]
            if dst_cpd_id not in model.metabolites:
                print("!")
                break
            dst_cpd = model.metabolites.get_by_id(dst_cpd_id)
            src_cpd = model.metabolites.get_by_id(src_cpd_id)
            target_rxn_ids = move[src_cpd_id][1]
            for rxn in src_cpd.reactions:
                if target_rxn_ids is None or rxn.id in target_rxn_ids:
                    rxn.add_metabolites(
                        {
                            src_cpd: -1 * rxn.metabolites[src_cpd],
                            dst_cpd: rxn.metabolites[src_cpd],
                        }
                    )
        return model

    def b(self, model_com):
        ## Replace to cytochromes
        model_com.reactions.rxn01806_cA.add_metabolites(
            {
                model_com.metabolites.cpd00004_cA: +1,
                model_com.metabolites.cpd00003_cA: -1,
                model_com.metabolites.cpd15560_cA: +1,
                model_com.metabolites.cpd15561_cA: -1,
            }
        )
        model_com.reactions.rxn01806_cR.add_metabolites(
            {
                model_com.metabolites.cpd00004_cR: +1,
                model_com.metabolites.cpd00003_cR: -1,
                model_com.metabolites.cpd15560_cR: +1,
                model_com.metabolites.cpd15561_cR: -1,
            }
        )

    def build(self):
        from modelseedpy import MSBuilder

        ex_ids = set(self.media_com.get_media_constraints().keys())
        ex_ids.add("cpd00075_e0")
        ex = {}
        m_pointer = {}

        from cobra.core import Model, Reaction, Metabolite

        model_com = Model("com")
        for model, token in self.models.items():
            metabolites = []
            for m in model.metabolites:
                if m.id in ex_ids:
                    # print(m)
                    if m.id not in ex:
                        m_copy = Metabolite(
                            m.id, m.formula, m.name, m.charge, m.compartment
                        )
                        m_pointer[m.id] = m.id
                        metabolites.append(m_copy)
                else:
                    m_id_copy = m.id[:-1] + token
                    # print(m_id_copy, m.id)
                    m_copy = Metabolite(
                        m_id_copy,
                        m.formula,
                        m.name[:-1] + token,
                        m.charge,
                        m.compartment[:-1] + token,
                    )
                    m_pointer[m.id] = m_copy.id
                    metabolites.append(m_copy)
            model_com.add_metabolites(metabolites)
            reactions = []
            for r in model.reactions:
                if not r.id.startswith("EX_"):
                    if r.id.startswith("bio"):
                        r_copy_id = r.id + "_" + token
                    else:
                        r_copy_id = r.id[:-1] + token
                    r_copy = Reaction(
                        r_copy_id,
                        r.name[:-1] + token,
                        r.subsystem,
                        r.lower_bound,
                        r.upper_bound,
                    )
                    r_copy_m = dict(
                        map(
                            lambda x: (
                                model_com.metabolites.get_by_id(m_pointer[x[0].id]),
                                x[1],
                            ),
                            r.metabolites.items(),
                        )
                    )
                    r_copy.add_metabolites(r_copy_m)
                    reactions.append(r_copy)
                    # print(r_copy.name)
                else:
                    # print(r)
                    pass
            model_com.add_reactions(reactions)
            MSBuilder.add_exchanges_to_model(model_com, 'e' + token)

        MSBuilder.add_exchanges_to_model(model_com, 'e0')
        print('Exchanges ', len(model_com.exchanges), 'e0')

        r_bio_sum = Reaction("bio1", "bio_com", "", 0, 1000)
        r_bio_sum.add_metabolites(
            {
                model_com.metabolites.cpd11416_cA: -0.4,
                model_com.metabolites.cpd11416_cR: -0.6,
            }
        )
        rxn09008_cR = Reaction("rxn09008_cR", "rxn09008_cR", "", -1000, 1000)
        rxn09008_cR.add_metabolites(
            {
                model_com.metabolites.cpd00418_cR: -1,
                model_com.metabolites.cpd00418_eA: 1,
            }
        )
        model_com.add_reactions([r_bio_sum, rxn09008_cR])
        model_com.objective = "bio1"
        model_com.reactions.bio1.annotation = {"sbo": "SBO:0000629"}
        model_com.reactions.bio1_A.annotation = {"sbo": "SBO:0000629"}
        model_com.reactions.bio1_R.annotation = {"sbo": "SBO:0000629"}
        for r in model_com.exchanges:
            r.lower_bound = 0
            r.upper_bound = 1000
        for k, (lb, ub) in self.media_com.get_media_constraints().items():
            if "EX_" + k in model_com.reactions:
                r = model_com.reactions.get_by_id("EX_" + k)
                r.lower_bound = lb
                r.upper_bound = ub
        model_com.reactions.EX_cpd00528_e0.upper_bound = 100
        model_com.reactions.EX_cpd00659_e0.upper_bound = 100
        model_com.reactions.rxn05488_cA.add_metabolites(
            {
                model_com.metabolites.cpd00067_e0: +1,
                model_com.metabolites.cpd00067_cA: -1,
            }
        )
        # rxn = Reaction('rxn05625_cR', 'rxn05625_cR', '', -1000, 1000)
        # rxn.add_metabolites({model_com.metabolites.cpd00075_cR: -1, model_com.metabolites.cpd00075_e0: 1})
        # model_com.add_reactions([rxn])

        # model_com.reactions.rxn11937_cR.lower_bound = 0 #-1000 # block nosZ
        return model_com

    def media_fix(self, model_com, media_com):
        mediacompounds = []
        for ex_id, v in model_com.medium.items():
            compound_id = ex_id[:-3][3:]
            print(
                list(
                    filter(
                        lambda x: x["id"] == compound_id + "_e0",
                        model_o["modelcompounds"],
                    )
                )[0]["compound_ref"]
            )
            mediacompounds.append(
                {
                    "compound_ref": "489/6/15/compounds/id/" + compound_id,
                    "concentration": 20,
                    "id": compound_id,
                    "inchikey": "",
                    "maxFlux": v,
                    "minFlux": -100,
                    "name": compound_id,
                    "smiles": "",
                }
            )
        data = media_com.get_data()
        data["mediacompounds"] = mediacompounds
        return data

    def fix_kbase_object_data(self, model_acido, model_o, model1, model2):
        """
        model1: model_com
        model2: kbase_model
        model_o: kbase_model (object)
        """
        model_o["genome_ref"] = model_acido.genome_ref
        model_o["template_ref"] = model_acido.template_ref

        for r in model1.reactions:
            if r.id not in model2.reactions:
                cpd_id = list(r.metabolites)[0]
                if len(r.metabolites) > 1:
                    print(cpd_id.id, r.metabolites)
                else:
                    modelReactionReagents = list(
                        map(
                            lambda x: {
                                "coefficient": x[1],
                                "modelcompound_ref": "~/modelcompounds/id/" + x[0].id,
                            },
                            r.metabolites.items(),
                        )
                    )
                    print(cpd_id.id, modelReactionReagents)

                    compound_id = cpd_id.id
                    mr = {
                        "aliases": [],
                        "dblinks": {},
                        "direction": ">",
                        "edits": {},
                        "gapfill_data": {},
                        "id": "OUT_" + compound_id + "e",
                        "maxforflux": 1000,
                        "maxrevflux": 0,
                        "modelReactionProteins": [],
                        "modelReactionReagents": modelReactionReagents,
                        "modelcompartment_ref": "~/modelcompartments/id/c0",
                        "name": "SINK OF " + compound_id,
                        "numerical_attributes": {},
                        "probability": 0,
                        "protons": 0,
                        "reaction_ref": "~/template/reactions/id/rxn00000_c",
                        "string_attributes": {},
                    }
                    model_o["modelreactions"].append(mr)
