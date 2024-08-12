from modelseedpy import MSBuilder, MSATPCorrection, MSMedia, MSGapfill


MEDIA_GENOME_SCALE = medium_gapfill = {
    'EX_cpd00029_e0': 20.0,
 'EX_cpd00013_e0': 100.0,
 'EX_cpd00001_e0': 100.0,
 'EX_cpd00218_e0': 0.0, #niacin
 'EX_cpd00220_e0': 100.0,
 'EX_cpd00305_e0': 100.0,
 'EX_cpd00393_e0': 100.0,
 'EX_cpd03424_e0': 100.0,
 'EX_cpd00443_e0': 0.0, #ABEE
 'EX_cpd00644_e0': 0.0002281,
 'EX_cpd00263_e0': 100.0,
 'EX_cpd00048_e0': 100.0,
 'EX_cpd00009_e0': 100.0,
 'EX_cpd00242_e0': 29.759425,
 'EX_cpd00205_e0': 1.3415688,
 'EX_cpd00063_e0': 100.0,
 'EX_cpd00971_e0': 34.9324073,
 'EX_cpd00099_e0': 100.0,
 'EX_cpd00254_e0': 100.0,
 'EX_cpd00030_e0': 100.0,
 'EX_cpd00058_e0': 100.0,
 'EX_cpd00034_e0': 100.0,
 'EX_cpd10515_e0': 100.0,
 'EX_cpd00149_e0': 100.0,
 'EX_cpd00244_e0': 100.0,
 'EX_cpd11574_e0': 100.0,
 'EX_cpd15574_e0': 100.0,
 'EX_cpd00067_e0': 100.0,
 'EX_cpd00209_e0': 10.0}


def load_genome_from_annotation_file(filename):
    from pandas import read_csv
    from modelseedpy.core.msgenome import MSGenome, MSFeature

    df = read_csv(filename, sep='\t', index_col=0)
    features = {}
    for feature_id, d in df.iterrows():
        gene_id = d['gene_id']
        feature = MSFeature(gene_id, '')
        for rast_feature in d['RAST'].split('; '):
            feature.add_ontology_term('RAST', rast_feature)
        if feature.id not in features:
            features[feature.id] = feature
        else:
            raise ValueError('duplicate gene ID')
    genome = MSGenome()
    genome.add_features(list(features.values()))

    return genome


class SynComStudy:

    def __init__(self, template):
        self.rast = None
        self.kbase = None
        self.genome_acido = None
        self.genome_rhoda = None
        self.template = template

    def build_isolates(self):
        model_acido = self.build_isolate('3H11', self.genome_acido)
        model_rhoda = self.build_isolate('R12', self.genome_rhoda)

        return model_acido, model_rhoda

    def build_isolate(self, model_id, genome):
        model = self.build_model(model_id, genome)
        model_gapfill = model.copy()

        atp_tests = self.atp_correction(model)
        gap_fill = self.gap_fill(model, atp_tests)

        rxn_new, ex_new = self._integrate_solution(model_gapfill, gap_fill)
        print(f'new reactions {len(rxn_new)} new exchanges {len(ex_new)}')
        model_gapfill.medium = {k: v for k, v in MEDIA_GENOME_SCALE.items() if k in model_gapfill.reactions}

        return model_gapfill

    @staticmethod
    def _add_atpm(model):
        from cobra.core import Reaction
        if 'ATPM_c0' not in model.reactions:
            atpm = Reaction(f'ATPM_c0', f'ATPM', 'ATPM', 0, 1000)
            atpm.add_metabolites({
                model.metabolites.cpd00001_c0: -1,
                model.metabolites.cpd00002_c0: -1,
                model.metabolites.cpd00008_c0: 1,
                model.metabolites.cpd00009_c0: 1,
                model.metabolites.cpd00067_c0: 1,
            })
            model.add_reactions([atpm])

    def build_model(self, model_id, genome):
        """
        build base model
        :param model_id:
        :param genome:
        :return:
        """
        b = MSBuilder(genome, self.template, model_id)
        model_base = b.build(model_id, annotate_with_rast=False)
        SynComStudy._add_atpm(model_base)
        return model_base

    def atp_correction(self, model):
        media = MSMedia.from_dict({
            'cpd00001': 1000,
            'cpd00067': 1000,
            'cpd00209': 1,
            'cpd00029': 1,
        })
        media.id = 'nitrate'
        media.name = 'nitrate'
        media.get_media_constraints()
        atp_correction = MSATPCorrection(model, self.template, [media],
                                         'c0', atp_hydrolysis_id='ATPM_c0',
                                         load_default_medias=False)
        tests = atp_correction.run_atp_correction()

        return tests

    def _get_solution(self, gf, sol, model):
        res = {}
        for rxn in gf.gfmodel.reactions:
            v = round(sol.fluxes[rxn.id], 9)
            if v != 0:
                if rxn.id not in model.reactions:
                    if rxn.id[:-1] in self.template.reactions:
                        # print(v, rxn.id, rxn.build_reaction_string(True))
                        lb = -1000
                        ub = 1000
                        if v > 0:
                            lb = 0
                        elif v < 0:
                            ub = 0
                        res[rxn.id[:-1]] = (lb, ub)
        return res

    def _integrate_solution(self, model, gap_fill_solution):
        added_reactions = []
        for rxn_id, (lb, ub) in gap_fill_solution.items():
            template_reaction = self.template.reactions.get_by_id(rxn_id)
            model_reaction = template_reaction.to_reaction(model)
            model_reaction.lower_bound = lb
            model_reaction.upper_bound = ub
            added_reactions.append(model_reaction)
        model.add_reactions(added_reactions)
        add_exchanges = MSBuilder.add_exchanges_to_model(model)

        return added_reactions, add_exchanges

    def gap_fill(self, model, tests, min_biomass=0.02):
        model.objective = 'bio1'
        gapfill = MSGapfill(model,
                            default_gapfill_templates=[self.template],
                            test_conditions=tests,
                            default_target='bio1')

        gapfill.gfmodel.reactions.bio1.lower_bound = min_biomass
        gapfill.gfmodel.medium = MEDIA_GENOME_SCALE

        gapfill_fba = gapfill.gfmodel.optimize()

        gapfill_solution = self._get_solution(gapfill, gapfill_fba, model)

        return gapfill_solution


def build_genome_annotation_file(genome, gene_id_remap, starts_with_filter='FW510-R12',
                                 filename_psort='psortb_r12_results.tsv',
                                 filename_out='annotation_rhoda.tsv'):
    """
    :param genome: rhoda or acido genome
    :param gene_id_remap: dict mapping feature_id to id from expression data
    :param starts_with_filter: FW510-R12 for rhoda GW101-3H11 for acido
    :param filename_psort: psort output file
    :param filename_out: output annotation_file
    :return:
    """
    psort_data = {}
    with open(filename_psort, 'r') as fh:
        h = fh.readline()
        l = fh.readline()
        header = {v: i for i, v in enumerate(h.split('\t'))}
        while l:
            _p = [x.strip() for x in l.split('\t')]
            psort_data[_p[header['SeqID']]] = [
                _p[header['Final_Localization']],
                _p[header['Final_Score']],
                _p[header['Final_Localization_Details']],
                _p[header['Secondary_Localization']],
            ]
            l = fh.readline()
    remap = {}
    for k, v in gene_id_remap.items():
        if k.startswith(starts_with_filter):
            if v not in remap:
                remap[v] = k
            else:
                print(k)
    with open(filename_out, 'w') as fh:
        d = [
            'feature_id', 'gene_id', 'RAST',
            'psort_loc', 'psort_loc_score', 'psort_loc_details', 'psort_loc_sec']
        fh.write('\t'.join(d) + '\n')
        for f in genome.features:
            d = [
                f.id,
                remap[f.id],
                '; '.join(f.ontology_terms.get('RAST', [])),
                psort_data[f.id][0],
                psort_data[f.id][1],
                psort_data[f.id][2],
                psort_data[f.id][3]
            ]
            fh.write('\t'.join(d) + '\n')


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


def comm(models):
    from cobra.core import Model, Reaction, Metabolite
    model_com = Model('com')
    ex_ids = set()
    for model in models.values():
        ex_ids |= {list(r.metabolites)[0].id for r in model.exchanges}

    ex = {}
    m_pointer = {}
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
        model_com.add_metabolites(list(m_pointer[token].values()))

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

    from modelseedpy import MSBuilder

    exchanges = MSBuilder.add_exchanges_to_model(model_com)

    return model_com

## NEW STUFF
from modelseedpy import MSBuilder, MSATPCorrection, MSMedia, MSGapfill


MEDIA_GENOME_SCALE = medium_gapfill = {
    'EX_cpd00029_e0': 20.0,
 'EX_cpd00013_e0': 100.0,
 'EX_cpd00001_e0': 100.0,
 'EX_cpd00218_e0': 0.0, #niacin
 'EX_cpd00220_e0': 100.0,
 'EX_cpd00305_e0': 100.0,
 'EX_cpd00393_e0': 100.0,
 'EX_cpd03424_e0': 100.0,
 'EX_cpd00443_e0': 0.0, #ABEE
 'EX_cpd00644_e0': 0.0002281,
 'EX_cpd00263_e0': 100.0,
 'EX_cpd00048_e0': 100.0,
 'EX_cpd00009_e0': 100.0,
 'EX_cpd00242_e0': 29.759425,
 'EX_cpd00205_e0': 1.3415688,
 'EX_cpd00063_e0': 100.0,
 'EX_cpd00971_e0': 34.9324073,
 'EX_cpd00099_e0': 100.0,
 'EX_cpd00254_e0': 100.0,
 'EX_cpd00030_e0': 100.0,
 'EX_cpd00058_e0': 100.0,
 'EX_cpd00034_e0': 100.0,
 'EX_cpd10515_e0': 100.0,
 'EX_cpd00149_e0': 100.0,
 'EX_cpd00244_e0': 100.0,
 'EX_cpd11574_e0': 100.0,
 'EX_cpd15574_e0': 100.0,
 'EX_cpd00067_e0': 100.0,
 'EX_cpd00209_e0': 10.0}


def load_genome_from_annotation_file(filename):
    from pandas import read_csv
    from modelseedpy.core.msgenome import MSGenome, MSFeature

    df = read_csv(filename, sep='\t', index_col=0)
    features = {}
    for feature_id, d in df.iterrows():
        gene_id = d['gene_id']
        feature = MSFeature(gene_id, '')
        for rast_feature in d['RAST'].split('; '):
            feature.add_ontology_term('RAST', rast_feature)
        if feature.id not in features:
            features[feature.id] = feature
        else:
            raise ValueError('duplicate gene ID')
    genome = MSGenome()
    genome.add_features(list(features.values()))

    return genome


class SynComStudy:

    def __init__(self, template):
        self.rast = None
        self.kbase = None
        self.genome_acido = None
        self.genome_rhoda = None
        self.template = template

    def build_isolates(self):
        model_acido = self.build_isolate('3H11', self.genome_acido)
        model_rhoda = self.build_isolate('R12', self.genome_rhoda)

        return model_acido, model_rhoda

    def build_isolate(self, model_id, genome):
        model = self.build_model(model_id, genome)
        model_gapfill = model.copy()

        atp_tests = self.atp_correction(model)
        gap_fill = self.gap_fill(model, atp_tests)

        rxn_new, ex_new = self._integrate_solution(model_gapfill, gap_fill)
        print(f'gapfill reactions {len(rxn_new)} new exchanges {len(ex_new)}')
        model_gapfill.medium = {k: v for k, v in MEDIA_GENOME_SCALE.items() if k in model_gapfill.reactions}

        return model_gapfill

    @staticmethod
    def _add_atpm(model):
        from cobra.core import Reaction
        if 'ATPM_c0' not in model.reactions:
            atpm = Reaction(f'ATPM_c0', f'ATPM', 'ATPM', 0, 1000)
            atpm.add_metabolites({
                model.metabolites.cpd00001_c0: -1,
                model.metabolites.cpd00002_c0: -1,
                model.metabolites.cpd00008_c0: 1,
                model.metabolites.cpd00009_c0: 1,
                model.metabolites.cpd00067_c0: 1,
            })
            model.add_reactions([atpm])

    def build_model(self, model_id, genome):
        """
        build base model
        :param model_id:
        :param genome:
        :return:
        """
        b = MSBuilder(genome, self.template, model_id)
        model_base = b.build(model_id, annotate_with_rast=False)
        SynComStudy._add_atpm(model_base)
        return model_base

    def atp_correction(self, model):
        media = MSMedia.from_dict({
            'cpd00001': 1000,
            'cpd00067': 1000,
            'cpd00209': 1,
            'cpd00029': 1,
        })
        media.id = 'nitrate'
        media.name = 'nitrate'
        media.get_media_constraints()
        atp_correction = MSATPCorrection(model, self.template, [media],
                                         'c0', atp_hydrolysis_id='ATPM_c0',
                                         load_default_medias=False)
        tests = atp_correction.run_atp_correction()

        return tests

    def _get_solution(self, gf, sol, model):
        res = {}
        for rxn in gf.gfmodel.reactions:
            v = round(sol.fluxes[rxn.id], 9)
            if v != 0:
                if rxn.id not in model.reactions:
                    if rxn.id[:-1] in self.template.reactions:
                        # print(v, rxn.id, rxn.build_reaction_string(True))
                        lb = -1000
                        ub = 1000
                        if v > 0:
                            lb = 0
                        elif v < 0:
                            ub = 0
                        res[rxn.id[:-1]] = (lb, ub)
        return res

    def _integrate_solution(self, model, gap_fill_solution):
        added_reactions = []
        for rxn_id, (lb, ub) in gap_fill_solution.items():
            template_reaction = self.template.reactions.get_by_id(rxn_id)
            model_reaction = template_reaction.to_reaction(model)
            model_reaction.lower_bound = lb
            model_reaction.upper_bound = ub
            added_reactions.append(model_reaction)
        model.add_reactions(added_reactions)
        add_exchanges = MSBuilder.add_exchanges_to_model(model)

        return added_reactions, add_exchanges

    def gap_fill(self, model, tests, min_biomass=0.02):
        model.objective = 'bio1'
        gapfill = MSGapfill(model,
                            default_gapfill_templates=[self.template],
                            test_conditions=tests,
                            default_target='bio1')

        gapfill.gfmodel.reactions.bio1.lower_bound = min_biomass
        gapfill.gfmodel.medium = MEDIA_GENOME_SCALE

        gapfill_fba = gapfill.gfmodel.optimize()

        gapfill_solution = self._get_solution(gapfill, gapfill_fba, model)

        return gapfill_solution


def build_genome_annotation_file(genome, gene_id_remap, starts_with_filter='FW510-R12',
                                 filename_psort='psortb_r12_results.tsv',
                                 filename_out='annotation_rhoda.tsv'):
    """
    :param genome: rhoda or acido genome
    :param gene_id_remap: dict mapping feature_id to id from expression data
    :param starts_with_filter: FW510-R12 for rhoda GW101-3H11 for acido
    :param filename_psort: psort output file
    :param filename_out: output annotation_file
    :return:
    """
    psort_data = {}
    with open(filename_psort, 'r') as fh:
        h = fh.readline()
        l = fh.readline()
        header = {v: i for i, v in enumerate(h.split('\t'))}
        while l:
            _p = [x.strip() for x in l.split('\t')]
            psort_data[_p[header['SeqID']]] = [
                _p[header['Final_Localization']],
                _p[header['Final_Score']],
                _p[header['Final_Localization_Details']],
                _p[header['Secondary_Localization']],
            ]
            l = fh.readline()
    remap = {}
    for k, v in gene_id_remap.items():
        if k.startswith(starts_with_filter):
            if v not in remap:
                remap[v] = k
            else:
                print(k)
    with open(filename_out, 'w') as fh:
        d = [
            'feature_id', 'gene_id', 'RAST',
            'psort_loc', 'psort_loc_score', 'psort_loc_details', 'psort_loc_sec']
        fh.write('\t'.join(d) + '\n')
        for f in genome.features:
            d = [
                f.id,
                remap[f.id],
                '; '.join(f.ontology_terms.get('RAST', [])),
                psort_data[f.id][0],
                psort_data[f.id][1],
                psort_data[f.id][2],
                psort_data[f.id][3]
            ]
            fh.write('\t'.join(d) + '\n')


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


def comm(models):
    from cobra.core import Model, Reaction, Metabolite
    model_com = Model('com')
    ex_ids = set()
    for model in models.values():
        ex_ids |= {list(r.metabolites)[0].id for r in model.exchanges}

    ex = {}
    m_pointer = {}
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
        model_com.add_metabolites(list(m_pointer[token].values()))

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

    from modelseedpy import MSBuilder

    exchanges = MSBuilder.add_exchanges_to_model(model_com)

    return model_com


import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


class CommPlots:

    def __init__(self, model):
        self.linestyle_3h11 = ':'
        self.linestyle_r12 = '-.'
        self.color_acetate = 'black'
        self.model = model

    def generate_solutions(self, cobra, rates, expected_growth_60R_40A):
        """
        generates growth solutions and max ATPM of R12 and 3H11
        :param cobra:
        :param rates:
        :param expected_growth_60R_40A:
        :return:
        """
        self.model.medium = {'EX_cpd00067_e0': 100.0,
                             'EX_cpd00058_e0': 100.0,
                             'EX_cpd00001_e0': 1000.0,
                             'EX_cpd00971_e0': 34.9324073,
                             'EX_cpd00013_e0': 1000.0,
                             'EX_cpd00244_e0': 100.0,
                             'EX_cpd00205_e0': 1.3415688,
                             'EX_cpd00009_e0': 100.0,
                             'EX_cpd11574_e0': 100.0,
                             'EX_cpd00305_e0': 100.0,
                             'EX_cpd00048_e0': 100.0,
                             'EX_cpd00209_e0': 12.0,
                             'EX_cpd00254_e0': 100.0,
                             'EX_cpd03424_e0': 100.0,
                             'EX_cpd15574_e0': 100.0,
                             'EX_cpd10515_e0': 100.0,
                             'EX_cpd00149_e0': 100.0,
                             'EX_cpd00029_e0': 20.0,
                             'EX_cpd00034_e0': 100.0,
                             'EX_cpd00063_e0': 100.0,
                             'EX_cpd00030_e0': 100.0,
                             'EX_cpd00099_e0': 100.0}
        solution_array = {}

        self.model.objective_direction = 'max'
        for i in range(len(rates['EX_cpd00029_e0'])):
            self.model.objective = 'bio1'

            uptake_ac = rates['EX_cpd00029_e0'][i]
            uptake_no3 = rates['EX_cpd00209_e0'][i] if rates['EX_cpd00209_e0'][i] < 0 else 0
            uptake_no2 = rates['EX_cpd00075_e0'][i] if rates['EX_cpd00075_e0'][i] < 0 else 0
            # uptake_ac = rates['EX_cpd00029_e0'][i]
            _medium_up = {
                'EX_cpd00029_e0': uptake_ac,
                'EX_cpd00209_e0': uptake_no3,
                'EX_cpd00075_e0': uptake_no2,
            }
            for ex_id in {'EX_cpd00029_e0', 'EX_cpd00209_e0', 'EX_cpd00075_e0'}:
                rxn_ex = self.model.reactions.get_by_id(ex_id)
                rxn_ex.lower_bound = _medium_up[ex_id]
                # print(i, rxn_ex.id, rxn_ex.lower_bound)
            sol = cobra.flux_analysis.pfba(self.model)
            solution_array[i] = {'solution': sol}
            extra_growth = sol.fluxes['bio1'] - expected_growth_60R_40A[i]

            print(i, 'GROWTH 0 ATPM, exceess', sol.fluxes['bio1'], extra_growth)

            if extra_growth > 0:
                self.model.reactions.bio1.lower_bound = expected_growth_60R_40A[i]
                self.model.objective = 'ATPM_cA'
                sol_atpm = cobra.flux_analysis.pfba(self.model)
                atpm_max_A = sol_atpm.fluxes['ATPM_cA']
                solution_array[i]['max_atpm_a'] = atpm_max_A
                self.model.objective = 'ATPM_cR'
                sol_atpm = cobra.flux_analysis.pfba(self.model)
                atpm_max_R = sol_atpm.fluxes['ATPM_cR']
                solution_array[i]['max_atpm_r'] = atpm_max_A
                print(i, 'MAX isolate ATPM', atpm_max_A, atpm_max_R)
            # solution_array[i] = {'solution': sol}

            self.model.reactions.bio1.lower_bound = 0
            self.model.reactions.ATPM_cA.lower_bound = 0
            self.model.objective = 'bio1'

        return solution_array

    def generate_growth_gap_solutions(self, cobra, expected_growth_60R_40A, rates, solution_array):
        """
        Generate ATPM solutions with minimum growth rate from expected growth data
        :param cobra:
        :param expected_growth_60R_40A:
        :param rates:
        :param solution_array:
        :return:
        """
        solution_exp = {}

        self.model.reactions.ATPM_cA.lower_bound = 0
        self.model.reactions.ATPM_cR.lower_bound = 0
        self.model.reactions.bio1.lower_bound = 0
        self.model.optimize()

        for i in range(len(expected_growth_60R_40A)):
            uptake_acetate = rates['EX_cpd00029_e0'][i]
            uptake_nitrate = rates['EX_cpd00209_e0'][i] if rates['EX_cpd00209_e0'][i] < 0 else 0
            uptake_nitrite = rates['EX_cpd00075_e0'][i] if rates['EX_cpd00075_e0'][i] < 0 else 0

            self.model.reactions.EX_cpd00029_e0.lower_bound = uptake_acetate
            self.model.reactions.EX_cpd00209_e0.lower_bound = uptake_nitrate
            self.model.reactions.EX_cpd00075_e0.lower_bound = uptake_nitrite

            self.model.reactions.ATPM_cA.lower_bound = solution_array[i]['max_atpm_a'] * 0.4
            self.model.reactions.bio1.lower_bound = expected_growth_60R_40A[i]
            self.model.objective = 'ATPM_cR'

            print('Acetate:', uptake_acetate,
                  'NO3', uptake_nitrate,
                  'NO2', uptake_nitrite,
                  'min community biomass', expected_growth_60R_40A[i], 'mATP 3H11',
                  self.model.reactions.ATPM_cA.lower_bound)

            solution = cobra.flux_analysis.pfba(self.model)

            solution_exp[i] = solution

            print('Found solution mATP R12', solution_exp[i].fluxes['ATPM_cR'], 'solution status', solution.status)

        return solution_exp

    @staticmethod
    def constants():
        atpm_max_cA_isolate = 9.67570638881224
        atpm_max_cR_isolate = 23.801193202545065
        time_steps = [23, 19, 11, 18, 24, 24, 23]
        expected_growth_60R_40A = [3.98504E-06, 4.40109E-05, 3.79833E-05, 0.00001748, 8.79E-06, 0.000008511,
                                   -1.46191E-06]
        rates = {
            # acetate ac
            'EX_cpd00029_e0': [-0.00796337, -0.077158311, -0.126205391, -0.047720639, -0.044917204, -0.0163845,
                               -0.018624539],
            # nitrate NO3
            'EX_cpd00209_e0': [-0.015765185, -0.260613674, -0.336007637, 0.000367934, -0.00050591, 0.000160971,
                               -0.00033594],
            # nitrite NO2
            'EX_cpd00075_e0': [-0.000300664, 0.145752713, 0.30446552, -0.112949487, -0.100448798, -0.06937439,
                               0.000601328],
            # Nitrous oxide N2O
            'EX_cpd00659_e0': [0.002927536, 0.007122807, -0.018121213, -3.70367E-05, 2.77775E-05, -0.000138889, 0]
        }

        return atpm_max_cA_isolate, atpm_max_cR_isolate, time_steps, expected_growth_60R_40A, rates

    @staticmethod
    def get_exp_syncom():
        exp_data_symcom = {
            'time': [0, 23, 42, 53, 71, 95, 119, 142],
            'biomass': [0.00008208, 0.000173736, 0.001009944, 0.00142776, 0.0017424, 0.00195336, 0.002157624, 0.002124],
            'acetate': [18.0313603, 17.8482028, 16.3821949, 14.9939356, 14.1349641, 13.0569512, 12.6637232, 12.2353588],
            'no3': [9.14450462, 8.78190537, 3.83024557, 0.13416156, 0.14078438, 0.12864254, 0.13250585, 0.12477923],
            'no2': [0.09977915, 0.09286388, 2.86216542, 6.21128614, 4.17819538, 1.76742424, 0.10243888, 0.11626942],
            'n2o': [0, 0.06733333, 0.20266667, 0.00333333, 0.00266667, 0.00333333, 0, 0]
        }
        df_exp = {
            'i': [],
            'acetate': [],
            'biomass': [],
            'no3': [],
            'no2': [],
            'n2o': []
        }
        for i in range(7):
            df_exp['i'].append(exp_data_symcom['time'][i])
            df_exp['biomass'].append(exp_data_symcom['biomass'][i])
            df_exp['acetate'].append(exp_data_symcom['acetate'][i])
            df_exp['no3'].append(exp_data_symcom['no3'][i])
            df_exp['no2'].append(exp_data_symcom['no2'][i])
            df_exp['n2o'].append(exp_data_symcom['n2o'][i])

        return df_exp

    @staticmethod
    def exp_atpm_pred(model, od_coeff, d):
        model.reactions.ATPM_c0.lower_bound = 0
        res = {x: [] for x in ['total_time', 'gdw', 'ac', 'no3', 'no2', 'n2o', 'atpm', 'atpm_per_gdw_per_h']}
        for i in range(len(d) - 1):
            d_start = d[i]
            d_end = d[i + 1]
            time = d_end['time'] - d_start['time']
            dt_acetate = (d_end['Acetate'] - d_start['Acetate']) / time
            dt_OD = (d_end['OD600'] - d_start['OD600']) / time
            dt_no3 = (d_end['NO3'] - d_start['NO3']) / time
            dt_no2 = (d_end['NO2'] - d_start['NO2']) / time
            dt_n2o = (d_end['N2O'] - d_start['N2O']) / time
            res['total_time'].append(time)
            res['gdw'].append(dt_OD * od_coeff)
            res['ac'].append(dt_acetate)
            res['no3'].append(dt_no3)
            res['no2'].append(dt_no2)
            res['n2o'].append(dt_n2o)
            max_atpm = None
            if dt_acetate < 0 and dt_no3 < 0:
                model.objective = 'ATPM_c0'
                model.reactions.bio1.lower_bound = dt_OD * od_coeff
                model.reactions.EX_cpd00209_e0.lower_bound = dt_no3
                model.reactions.EX_cpd00029_e0.lower_bound = dt_acetate
                try:
                    sol_max_atpm = model.optimize()
                    if sol_max_atpm.status == 'optimal':
                        max_atpm = sol_max_atpm.fluxes['ATPM_c0']
                except Exception as e:
                    pass

            res['atpm'].append(max_atpm)
            res['atpm_per_gdw_per_h'].append(max_atpm / time if max_atpm else None)

        return pd.DataFrame(res)

    def generate_total_acc_data(self, time_steps, solution_exp):
        df_array = {
            'i': [],
            'gdw_3h11': [],
            'gdw_r12': [],
            'gdw_total': [],
            'acetate': [],
            'NO3': [],
            'NO2': [],
            'N2O': [],
            'N2': [],
            'NO': [],
        }
        acc_r_growth_acc = 0.00005472
        acc_a_growth_acc = 0.00002736
        acc_acetate = 18.0313603
        acc_no3 = 9.14450462
        acc_no2 = 0.09977915
        acc_no = 0
        acc_n2o = 0
        acc_n2 = 0
        time = 0
        for t_index in range(len(time_steps)):
            sol_exp = solution_exp[t_index]
            growth_3H11 = sol_exp.fluxes['bio1_A']
            growth_R12 = sol_exp.fluxes['bio1_R']
            up_acetate = sol_exp.fluxes['EX_cpd00029_e0']
            up_no3 = sol_exp.fluxes['EX_cpd00209_e0']
            up_no2 = sol_exp.fluxes['EX_cpd00075_e0']
            up_n2 = sol_exp.fluxes['EX_cpd00528_e0']
            up_n2o = sol_exp.fluxes['EX_cpd00659_e0']
            up_no = sol_exp.fluxes['EX_cpd00418_e0']

            iterations = time_steps[t_index]
            for it in range(iterations):
                df_array['i'].append(time)
                df_array['gdw_3h11'].append(acc_a_growth_acc)
                df_array['gdw_r12'].append(acc_r_growth_acc)
                df_array['gdw_total'].append(acc_r_growth_acc + acc_a_growth_acc)
                df_array['acetate'].append(acc_acetate)
                df_array['NO3'].append(acc_no3)
                df_array['NO2'].append(acc_no2)
                df_array['N2O'].append(acc_n2o)
                df_array['N2'].append(acc_n2)
                df_array['NO'].append(acc_no)
                acc_r_growth_acc += growth_R12
                acc_a_growth_acc += growth_3H11
                acc_no3 += up_no3
                acc_no2 += up_no2
                acc_n2o += up_n2o
                acc_n2 += up_n2
                acc_no += up_no
                acc_acetate += up_acetate
                time += 1

        print('Biomass 3H11', df_array['gdw_3h11'][-1],
              'Biomass R12', df_array['gdw_r12'][-1],
              'Biomass Total', df_array['gdw_3h11'][-1] + df_array['gdw_r12'][-1],
              'Time', df_array['i'][-1],
              'Acetate', df_array['acetate'][-1],
              'NO3', df_array['NO3'][-1],
              'NO2', df_array['NO2'][-1],
              'N2O', df_array['N2O'][-1],
              'N2', df_array['N2'][-1])

        return df_array

    def plot_total_acc(self, df_array, df_exp):


        fig, ax = plt.subplots()
        ax2 = ax.twinx()
        fig.set_size_inches(25.7, 8.27)
        # sns.lineplot(data = df_array[0], x='i', y='value', hue='reaction', marker='o', sort = False, ax=ax, linestyle='-')
        sns.lineplot(data=df_array, x='i', y='gdw_3h11', ax=ax2, linestyle=self.linestyle_3h11, color='sienna')
        sns.lineplot(data=df_array, x='i', y='gdw_r12', ax=ax2, linestyle=self.linestyle_r12, color='sienna')
        sns.lineplot(data=df_array, x='i', y='gdw_total', ax=ax2, linestyle='--', color='sienna')
        sns.scatterplot(data=df_exp, x='i', y='biomass', ax=ax2, color='sienna')
        sns.lineplot(data=df_array, x='i', y='acetate', ax=ax, color=self.color_acetate)
        sns.scatterplot(data=df_exp, x='i', y='acetate', ax=ax, color=self.color_acetate)
        sns.lineplot(data=df_array, x='i', y='NO3', ax=ax, color='blue')
        sns.scatterplot(data=df_exp, x='i', y='no3', ax=ax, color='blue')
        sns.lineplot(data=df_array, x='i', y='NO2', ax=ax, color='green')
        sns.scatterplot(data=df_exp, x='i', y='no2', ax=ax, color='green')
        sns.lineplot(data=df_array, x='i', y='N2O', ax=ax, color='purple')
        sns.scatterplot(data=df_exp, x='i', y='n2o', ax=ax, color='purple')
        sns.lineplot(data=df_array, x='i', y='NO', ax=ax, color='orange', linestyle='--')
        sns.lineplot(data=df_array, x='i', y='N2', ax=ax, color='cyan')

        plt.legend(handles=[
            mpatches.Patch(color=self.color_acetate, label='Acetate'),
            mpatches.Patch(color='blue', label='NO3'),
            mpatches.Patch(color='green', label='NO2'),
            mpatches.Patch(color='orange', label='NO'),
            mpatches.Patch(color='purple', label='N2O'),
            mpatches.Patch(color='cyan', label='N2'),
            mpatches.Patch(label='Dot experimental values', facecolor=None, color='white'),
            mpatches.Patch(edgecolor='sienna', label='Biomass 3H11', linestyle=self.linestyle_3h11, facecolor='white'),
            mpatches.Patch(edgecolor='sienna', label='Biomass R12', linestyle=self.linestyle_r12, facecolor='white'),
            mpatches.Patch(edgecolor='sienna', label='Biomass Total', linestyle='--', facecolor='white'),
        ])
        ax.set_title('Totals')
        ax.set_xlabel('time')
        ax.set_ylabel('mM')
        ax2.set_ylabel('gDW')
        sns.despine(fig, ax)

    def generate_total_uptake_data(self, time_steps, solution_exp):
        df_array = {
            'i': [],
            'gdw_3h11': [],
            'gdw_r12': [],
            'atpm_3h11': [],
            'atpm_r12': [],
            'acetate': [],
            'NO3': [],
            'NO2': [],
            'N2O': [],
            'N2': [],
            'NO': [],
        }
        time = 0
        for t_index in range(len(time_steps)):
            sol_exp = solution_exp[t_index]
            growth_3H11 = sol_exp.fluxes['bio1_A']
            growth_R12 = sol_exp.fluxes['bio1_R']
            atpm_3H11 = sol_exp.fluxes['ATPM_cA']
            atpm_R12 = sol_exp.fluxes['ATPM_cR']
            print(atpm_3H11, atpm_R12)
            up_acetate = sol_exp.fluxes['EX_cpd00029_e0']
            up_no3 = sol_exp.fluxes['EX_cpd00209_e0']
            up_no2 = sol_exp.fluxes['EX_cpd00075_e0']
            up_n2 = sol_exp.fluxes['EX_cpd00528_e0']
            up_n2o = sol_exp.fluxes['EX_cpd00659_e0']
            up_no = sol_exp.fluxes['EX_cpd00418_e0']

            iterations = time_steps[t_index]
            for it in range(iterations):
                df_array['i'].append(time)
                df_array['gdw_3h11'].append(growth_3H11)
                df_array['gdw_r12'].append(growth_R12)
                df_array['atpm_3h11'].append(atpm_3H11)
                df_array['atpm_r12'].append(atpm_R12)
                df_array['acetate'].append(up_acetate)
                df_array['NO3'].append(up_no3)
                df_array['NO2'].append(up_no2)
                df_array['N2O'].append(up_n2o)
                df_array['N2'].append(up_n2)
                df_array['NO'].append(up_no)
                time += 1

        return df_array

    def plot_total_uptake(self, df_array):
        fig, ax = plt.subplots(3, 1, height_ratios=[3, 1, 1])
        fig.set_size_inches(25.7, 8.27)

        sns.lineplot(data=df_array, x='i', y='acetate', ax=ax[0], color=self.color_acetate)
        sns.lineplot(data=df_array, x='i', y='NO3', ax=ax[0], color='blue')
        sns.lineplot(data=df_array, x='i', y='NO2', ax=ax[0], color='green')
        sns.lineplot(data=df_array, x='i', y='N2O', ax=ax[0], color='purple')
        sns.lineplot(data=df_array, x='i', y='NO', ax=ax[0], color='orange', linestyle='--')
        sns.lineplot(data=df_array, x='i', y='N2', ax=ax[0], color='cyan')

        sns.lineplot(data=df_array, x='i', y='gdw_3h11', ax=ax[1], linestyle=self.linestyle_3h11,
                     color='sienna')
        sns.lineplot(data=df_array, x='i', y='gdw_r12', ax=ax[1], linestyle=self.linestyle_r12,
                     color='sienna')
        sns.lineplot(data=df_array, x='i', y='atpm_3h11', ax=ax[2], linestyle=self.linestyle_3h11,
                     color='dimgray')
        sns.lineplot(data=df_array, x='i', y='atpm_r12', ax=ax[2], linestyle=self.linestyle_r12,
                     color='dimgray')
        ax[0].set_xlabel('')
        ax[1].set_xlabel('')
        ax[2].set_xlabel('time')
        ax[0].set_ylabel('mM/gDW/h')
        ax[1].set_ylabel('Biomass 1/h')
        ax[2].set_ylabel('mATP mM/gDW/h')
        sns.despine(fig, ax[0])
        sns.despine(fig, ax[1])
        sns.despine(fig, ax[2])

    def generate_organism_uptake_data(self, time_steps, solution_exp, monitor=None):
        if monitor is None:
            monitor = {
                'ac': 'cpd00029_e0'
            }
        df_array = []
        for i in range(3):
            df_array.append({
                'i': [],
                'growth': [],
                'acetate': [],
                'no3': [],
                'no2': [],
                'no': [],
                'n2o': [],
                'n2': [],
                'leu': []
            })
        time = 0
        for t_index in range(len(time_steps)):
            sol_exp = solution_exp[t_index]
            growth_3H11 = sol_exp.fluxes['bio1_A']
            growth_R12 = sol_exp.fluxes['bio1_R']
            leu_total, leu_acido, leu_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00107_e0', sol_exp))
            ac_total, ac_acido, ac_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00029_e0', sol_exp))
            no3_total, no3_acido, no3_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00209_e0', sol_exp))
            no2_total, no2_acido, no2_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00075_e0', sol_exp))
            no_total, no_acido, no_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00418_e0', sol_exp))
            n2o_total, n2o_acido, n2o_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00659_e0', sol_exp))
            n2_total, n2_acido, n2_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00528_e0', sol_exp))

            print(ac_total, ac_acido, ac_rhoda)

            iterations = time_steps[t_index]
            for it in range(iterations):
                df_array[0]['i'].append(time)
                df_array[1]['i'].append(time)
                df_array[2]['i'].append(time)
                df_array[1]['growth'].append(sol_exp.fluxes['bio1_A'])
                df_array[2]['growth'].append(sol_exp.fluxes['bio1_R'])
                df_array[0]['acetate'].append(ac_total)
                df_array[1]['acetate'].append(ac_acido)
                df_array[2]['acetate'].append(ac_rhoda)

                df_array[0]['no3'].append(no3_total)
                df_array[1]['no3'].append(no3_acido)
                df_array[2]['no3'].append(no3_rhoda)

                df_array[0]['no2'].append(no2_total)
                df_array[1]['no2'].append(no2_acido)
                df_array[2]['no2'].append(no2_rhoda)

                df_array[0]['no'].append(no_total)
                df_array[1]['no'].append(no_acido)
                df_array[2]['no'].append(no_rhoda)

                df_array[0]['n2o'].append(n2o_total)
                df_array[1]['n2o'].append(n2o_acido)
                df_array[2]['n2o'].append(n2o_rhoda)

                df_array[0]['n2'].append(n2_total)
                df_array[1]['n2'].append(n2_acido)
                df_array[2]['n2'].append(n2_rhoda)

                df_array[0]['leu'].append(leu_total)
                df_array[1]['leu'].append(leu_acido)
                df_array[2]['leu'].append(leu_rhoda)
                time += 1

        return df_array

    def plot_organism_uptake_data(self, df_array, title='Organism Uptake/Secretion rates'):
        fig, ax_arr = plt.subplots(3, 1, height_ratios=[3, 1, 1])
        fig.set_size_inches(25.7, 8.27)
        ax = ax_arr[0]
        ax2 = ax_arr[1]
        ax3 = ax_arr[2]

        sns.lineplot(data=df_array[1], x='i', y='no3', ax=ax3, color='blue', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='no3', ax=ax3, color='blue', linestyle=self.linestyle_r12)

        sns.lineplot(data=df_array[1], x='i', y='no2', ax=ax3, color='green', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='no2', ax=ax3, color='green', linestyle=self.linestyle_r12)

        sns.lineplot(data=df_array[1], x='i', y='leu', ax=ax2, color='green', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='leu', ax=ax2, color='green', linestyle=self.linestyle_r12)

        sns.lineplot(data=df_array[1], x='i', y='acetate', ax=ax, color='black', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='acetate', ax=ax, color='black', linestyle=self.linestyle_r12)

        sns.lineplot(data=df_array[1], x='i', y='no', ax=ax, color='orange', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='no', ax=ax, color='orange', linestyle=self.linestyle_r12)

        sns.lineplot(data=df_array[1], x='i', y='n2o', ax=ax, color='purple', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='n2o', ax=ax, color='purple', linestyle=self.linestyle_r12)

        sns.lineplot(data=df_array[1], x='i', y='n2', ax=ax, color='cyan', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='n2', ax=ax, color='cyan', linestyle=self.linestyle_r12)
        ax.set_title(title)

        ax.legend(handles=[
            mpatches.Patch(edgecolor=self.color_acetate, label='Acetate (3H11)', facecolor='white',
                           linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor=self.color_acetate, label='Acetate (R12)', facecolor='white',
                           linestyle=self.linestyle_r12),
            mpatches.Patch(edgecolor='orange', label='NO (3H11)', facecolor='white', linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor='orange', label='NO (R12)', facecolor='white', linestyle=self.linestyle_r12),
            mpatches.Patch(edgecolor='purple', label='N2O (3H11)', facecolor='white', linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor='purple', label='N2O (R12)', facecolor='white', linestyle=self.linestyle_r12),
            mpatches.Patch(edgecolor='cyan', label='N2 (3H11)', facecolor='white', linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor='cyan', label='N2 (R12)', facecolor='white', linestyle=self.linestyle_r12),
        ])
        ax2.legend(handles=[
            mpatches.Patch(edgecolor='green', label='Leucine (3H11)', facecolor='white', linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor='green', label='Leucine (R12)', facecolor='white', linestyle=self.linestyle_r12),
        ])
        ax3.legend(handles=[
            mpatches.Patch(edgecolor='blue', label='NO3 (3H11)', facecolor='white', linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor='blue', label='NO3 (R12)', facecolor='white', linestyle=self.linestyle_r12),
            mpatches.Patch(edgecolor='green', label='NO2 (3H11)', facecolor='white', linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor='green', label='NO2 (R12)', facecolor='white', linestyle=self.linestyle_r12),
        ])

        ax.set_ylabel('mM/gDW/h')
        ax2.set_ylabel('mM/gDW/h')
        ax3.set_ylabel('mM/gDW/h')
        ax.set_xlabel('')
        ax2.set_xlabel('')
        ax3.set_xlabel('time')
        sns.despine(fig, ax)
        sns.despine(fig, ax2)
        sns.despine(fig, ax3)

    def generate_leucine_efflux_data(self, cobra, solution_atpm, expected_growth_60R_40A, rates, per_atpm_to_leucine):
        solution_leucine = {}
        for i in range(len(expected_growth_60R_40A)):
            solution_leucine[i] = None
            self.model.reactions.EX_cpd00029_e0.lower_bound = rates['EX_cpd00029_e0'][i]
            self.model.reactions.EX_cpd00209_e0.lower_bound = rates['EX_cpd00209_e0'][i] if \
            rates['EX_cpd00209_e0'][
                i] < 0 else 0
            self.model.reactions.EX_cpd00075_e0.lower_bound = rates['EX_cpd00075_e0'][i] if \
            rates['EX_cpd00075_e0'][
                i] < 0 else 0
            self.model.reactions.ATPM_cA.lower_bound = solution_atpm[i].fluxes['ATPM_cA'] * (1 - per_atpm_to_leucine)
            self.model.reactions.ATPM_cR.lower_bound = solution_atpm[i].fluxes['ATPM_cR'] * (1 - per_atpm_to_leucine)
            self.model.reactions.bio1.lower_bound = expected_growth_60R_40A[i]
            self.model.objective = 'LeuE_cA'
            self.model.objective_direction = 'min'
            try:
                solution_leucine[i] = cobra.flux_analysis.pfba(self.model)
                print(i, solution_leucine[i].fluxes['LeuE_cA'])
            except Exception as ex:
                pass
        self.model.objective_direction = 'max'
        self.model.objective = 'bio1'
        return solution_leucine

    @staticmethod
    def _get_cpd_acc(model, cpd_id, sol):
        cpd_acc = {}
        metabolite = model.metabolites.get_by_id(cpd_id)
        for r in metabolite.reactions:
            rxn_s = r.metabolites[metabolite]
            v = sol.fluxes[r.id]
            cmps = frozenset(r.compartments)
            if v != 0:
                if cmps not in cpd_acc:
                    cpd_acc[cmps] = 0
                cpd_acc[cmps] += v * rxn_s
        return cpd_acc

    @staticmethod
    def _cpd_acc_to_a_r_t(cpd_acc):
        """
        parse cpd acc to acido, rhoda, total
        """
        return cpd_acc.get(frozenset({'e0'}), 0), cpd_acc.get(frozenset({'e0', 'cA'}), 0), cpd_acc.get(frozenset({'e0', 'cR'}), 0)
