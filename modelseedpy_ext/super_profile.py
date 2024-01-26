from modelseedpy.core.msatpcorrection import MSATPCorrection
from modelseedpy.core.msgenome import normalize_role
from cobra.core import Reaction
from cobra.flux_analysis import pfba


SEED_AA_FLAVOR = {
            'pyr': 'cpd00020_c0',
            'leu': 'cpd00107_c0',
            'lys': 'cpd00039_c0',
            'his': 'cpd00119_c0',
            'ile': 'cpd00322_c0',
            'thr': 'cpd00161_c0',
            'trp': 'cpd00065_c0',
            'tyr': 'cpd00069_c0',
            'ser': 'cpd00054_c0',
            'met': 'cpd00060_c0',
            'cys': 'cpd00084_c0',
            'arg': 'cpd00051_c0',
            'asn': 'cpd00132_c0',
            'asp': 'cpd00041_c0',
            'ala': 'cpd00035_c0',
            'gln': 'cpd00053_c0',
            'glu': 'cpd00023_c0',
            'gly': 'cpd00033_c0',
            'val': 'cpd00156_c0',
            'pro': 'cpd00129_c0',
            'phe': 'cpd00066_c0',
        }

SEED_FLAVOR_2 = {
            'glc': '',
            'atp': 'cpd00002_c0',
            'adp': 'cpd00008_c0',
            'h2o': 'cpd00001_c0',
            'h_c': 'cpd00067_c0',
            'h_e': 'cpd00067_e0',
            'pi': 'cpd00009_c0',
            'nad': 'cpd00003_c0',
            'nadh': 'cpd00004_c0',
            'nadp': 'cpd00006_c0',
            'nadph': 'cpd00005_c0',
            'q8': 'cpd15560_c0',
            'q8h2': 'cpd15561_c0',
            'accoa': 'cpd00022_c0',
            'ac': 'cpd00029_c0',
            'coa': 'cpd00010_c0',
        }

BIGG_AA_FLAVOR = {
            'pyr': 'pyr_c',
            'leu': 'leu__L_c',
            'lys': 'lys__L_c',
            'his': 'his__L_c',
            'ile': 'ile__L_c',
            'thr': 'thr__L_c',
            'trp': 'trp__L_c',
            'tyr': 'tyr__L_c',
            'ser': 'ser__L_c',
            'met': 'met__L_c',
            'cys': 'cys__L_c',
            'arg': 'arg__L_c',
            'asn': 'asn__L_c',
            'asp': 'asp__L_c',
            'ala': 'ala__L_c',
            'gln': 'gln__L_c',
            'glu': 'glu__L_c',
            'gly': 'gly_c',
            'val': 'val__L_c',
            'pro': 'pro__L_c',
            'phe': 'phe__L_c',
        }

BIGG_FLAVOR_2 = {
            'glc': '',
            'atp': 'atp_c',
            'adp': 'adp_c',
            'h2o': 'h2o_c',
            'h': '',
            'pi': 'pi_c',
            'nad': 'nad_c',
            'nadh': 'nadh_c',
            'nadp': 'nadp_c',
            'nadph': 'nadph_c',
            'q8': 'q8_c',
            'q8h2': 'q8h2_c',
            'accoa': 'accoa_c',
            'ac': 'ac_c',
            'coa': 'coa_c',
        }


class Profiler:
    def __init__(self, master):
        self.master = master
        pass


class PathwayModuleProfiler:
    def __init__(self, master):
        self.master = master

    def profile_genome(self, genome_id):
        pass


class FunctionProfiler:
    def __init__(self, master):
        self.master = master

    def profile_genome(self, genome_id):
        genome = self.master.fn_get_genome(genome_id)
        for f in genome.features:
            for f_str in f.functions:
                nmz = normalize_role(f_str)
                if nmz not in nmz_to_role:
                    nmz_to_role[nmz] = set()
                    nmz_to_genome_gene[nmz] = {}
                nmz_to_role[nmz].add(f_str)
                if genome_id not in nmz_to_genome_gene[nmz]:
                    nmz_to_genome_gene[nmz][genome_id] = set()
                nmz_to_genome_gene[nmz][genome_id].add(f.id)


class AtpCoreProfiler:
    def __init__(self, master, medias, atp_hydrolysis_id, extracell="e0"):
        self.master = master
        self.atp_hydrolysis_id = atp_hydrolysis_id
        self.medias = medias
        self.comp = {"EX_cpd00001_e0": 1000, "EX_cpd00067_e0": 1000}
        self.extracell = extracell

    def profile_genome(self, genome_id):
        model = self.master.fn_get_model(genome_id)
        model.objective = self.atp_hydrolysis_id
        return self.test_core_model(model, self.medias, cmp=self.extracell)

    def test_media(self, model, media, exs, cmp):
        media_const = media.get_media_constraints(cmp)
        medium_ = {}
        for cpd_id in media_const:
            if cpd_id not in exs:
                return None
        for cpd_id in media_const:
            lb, ub = media_const[cpd_id]
            if cpd_id in exs:
                rxn_exchange = exs[cpd_id]
                medium_[rxn_exchange.id] = -1 * lb
        medium_.update(self.comp)
        model.medium = medium_
        # print(media.id, model.medium)
        sol = model.optimize()
        return sol

    def get_metabolite_exchanges(self, model):
        metabolite_exchanges = {}
        for rxn_exchange in model.exchanges:
            metabolites = rxn_exchange.metabolites
            if len(metabolites) == 1:
                metabolite_exchanges[list(metabolites)[0].id] = rxn_exchange
            else:
                print("ignore", rxn_exchange)
        return metabolite_exchanges

    def test_core_model(self, model, medias, v=False, cmp="e0"):
        res = {}
        media_out = {}
        metabolite_exchanges = self.get_metabolite_exchanges(model)
        for media_id in medias:
            # print(media_id)
            sol = self.test_media(model, medias[media_id], metabolite_exchanges, cmp)
            media_out[media_id] = None
            res[media_id] = 0
            if sol and sol.status == "optimal":
                ex_out = dict(
                    filter(
                        lambda x: x[0].startswith("EX_") and x[1] != 0,
                        sol.fluxes.to_dict().items(),
                    )
                )
                media_out[media_id] = ex_out
                # print(media_id, sol.objective_value, ex_out)
                atp_val = sol.objective_value
                if atp_val > 0:
                    for cpd_id in medias[media_id].get_media_constraints(cmp):
                        media_cpd_val = sol.fluxes[metabolite_exchanges[cpd_id].id]
                        if media_cpd_val == 0:
                            atp_val = 0
                            # print('reject', media_id)
                res[media_id] = atp_val
        return res, media_out


class ProfilerAA:

    def __init__(self, model):
        self.model = model
        self.amino_acids = BIGG_AA_FLAVOR
        self.others = BIGG_FLAVOR_2
        self.test_reactions = {}

    def build_synth_metabolite_test(self, metabolite):
        if metabolite.id in self.model.metabolites:
            m = self.model.metabolites.get_by_id(metabolite.id)
            rxn_test = Reaction(f'test_{m.id}', f'Test {m.id} [{m.name}]', 'TEST', 0, 0)
            rxn_test.add_metabolites({
                m: -1
            })
            return rxn_test
        else:
            return None

    def build_all_synth_test(self):
        _to_add = []
        for m in self.model.metabolites:
            rxn_test = self.build_synth_metabolite_test(m)
            if rxn_test.id not in self.model.reactions:
                self.test_reactions[m.id] = rxn_test
                _to_add.append(rxn_test)
        self.model.add_reactions(_to_add)

    def build_tests2(self):
        rxn_test_atp = Reaction(f'test_atp', f'Test ATP', 'TEST', 0, 0)
        rxn_test_atp.add_metabolites({
            self.model.metabolites.get_by_id(self.others['h2o']): -1,
            self.model.metabolites.get_by_id(self.others['atp']): -1,
            self.model.metabolites.get_by_id(self.others['adp']): 1,
            self.model.metabolites.get_by_id(self.others['pi']): 1,
            self.model.metabolites.get_by_id(self.others['h_c']): 1,
        })
        rxn_test_nad = Reaction(f'test_nad', f'Test NAD', 'TEST', 0, 0)
        rxn_test_nad.add_metabolites({
            self.model.metabolites.get_by_id(self.others['nad']): -1,
            self.model.metabolites.get_by_id(self.others['nadh']): 1,
            self.model.metabolites.get_by_id(self.others['h_e']): -1,
        })
        rxn_test_nadp = Reaction(f'test_nadp', f'Test NADP', 'TEST', 0, 0)
        rxn_test_nadp.add_metabolites({
            self.model.metabolites.get_by_id(self.others['nadp']): -1,
            self.model.metabolites.get_by_id(self.others['nadph']): 1,
            self.model.metabolites.get_by_id(self.others['h_e']): -1,
        })
        rxn_test_q8 = Reaction(f'test_q8', f'Test Q8', 'TEST', 0, 0)
        rxn_test_q8.add_metabolites({
            self.model.metabolites.get_by_id(self.others['q8']): -1,
            self.model.metabolites.get_by_id(self.others['q8h2']): 1,
            self.model.metabolites.get_by_id(self.others['h_e']): -2,
        })
        rxn_accoa_coa = Reaction('test_accoa_coa', 'test_accoa_coa', 'TEST', 0, 0)
        rxn_accoa_coa.add_metabolites({
            self.model.metabolites.get_by_id(self.others['h2o']): -1,
            self.model.metabolites.get_by_id(self.others['accoa']): -1,
            self.model.metabolites.get_by_id(self.others['coa']): 1,
            self.model.metabolites.get_by_id(self.others['ac']): 1,
            self.model.metabolites.get_by_id(self.others['h_c']): 1,
        })
        to_add = [rxn_test_atp, rxn_test_nad, rxn_test_nadp, rxn_test_q8, rxn_accoa_coa]
        _to_add = []
        for r in to_add:
            if r.id not in self.model.reactions:
                _to_add.append(r)
        self.model.add_reactions(_to_add)
        added = {}
        for r in to_add:
            if r.id in self.model.reactions:
                added[r.id] = self.model.reactions.get_by_id(r.id)
        return added

    def build_tests(self):
        for aa, cpd_id in self.amino_acids.items():
            cpd = self.model.metabolites.get_by_id(cpd_id)
            rxn_test = Reaction(f'test_{aa}', f'Test {cpd.id} [{cpd.name}]', 'TEST', 0, 0)
            rxn_test.add_metabolites({
                cpd: -1
            })

            if rxn_test.id not in self.model.reactions:
                self.test_reactions[aa] = rxn_test
        self.model.add_reactions(list(self.test_reactions.values()))

    def profile_genome(self, genome_id):
        result = {}
        for test_id, r in self.test_reactions.items():
            self.model.objective = r.id
            model_reaction = self.reactions.get_by_id(r.id)
            model_reaction.lower_bound = 0
            model_reaction.upper_bound = 1
            solution = pfba(self.model)
            obj = solution.fluxes[r.id]
            if obj > 0 and solution.status == 'optimal':
                result[test_id] = solution
            else:
                result[r.id] = None
            r.upper_bound = 0

        return result


class AtpCoreGapfillProfiler:
    def __init__(self, master, atp_hydrolysis_id):
        self.master = master
        self.atp_hydrolysis_id = atp_hydrolysis_id
        self.atp_basic_medias = None
        self.atp_medias = None
        self.res = {}
        self.template = None
        self.min_gap = {
            "Glc/O2": 5,
            "Etho/O2": 0.01,
            "Ac/O2": 1,
            "Pyr/O2": 3,
            "Glyc/O2": 2,
            "Fum/O2": 3,
            "Succ/O2": 2,
            "Akg/O2": 2,
            "LLac/O2": 2,
            "Dlac/O2": 2,
            "For/O2": 2,
            "For/NO3": 1.5,
            "Pyr/NO": 2.5,
            "Pyr/NO2": 2.5,
            "Pyr/NO3": 2.5,
            "Pyr/SO4": 2.5,
            "Glc/DMSO": 5,
            "Glc/TMAO": 5,
        }

    def profile_genome(self, genome_id):
        result = {}
        model = self.master.fn_get_model(genome_id)
        atp_method = MSATPCorrection(
            model,
            self.template,
            self.atp_basic_medias,
            atp_hydrolysis_id=self.atp_hydrolysis_id,
        )
        atp_method.evaluate_growth_media()
        for media in atp_method.media_gapfill_stats:
            result[media.name] = atp_method.media_gapfill_stats[media]
        for media_id in self.min_gap:
            if media_id in self.atp_medias:
                res_min_gap = atp_method.msgapfill.run_gapfilling(
                    self.atp_medias[media_id],
                    self.atp_hydrolysis_id,
                    minimum_obj=self.min_gap[media_id],
                    binary_check=True,
                )
                result[media_id] = res_min_gap
        return result

    def profile(self):
        pass


class GeProf:
    def __init__(self, genome_ids, fn_get_model, fn_get_genome, genome_order=None):
        self.genome_order = None
        self.genome_ids = genome_ids
        self.fn_get_model = fn_get_model
        self.fn_get_genome = fn_get_genome

    def profile(self):
        pass


class TestModel:

    def __init__(self, model):
        self._model = model
        self.test_reactions = None
        self.test = None

    def add_tests(self):
        amino_acids = """
        cpd00107_c0
        cpd00039_c0
        cpd00119_c0
        cpd00322_c0
        cpd00161_c0
        cpd00065_c0
        cpd00069_c0
        cpd00054_c0
        cpd00060_c0
        cpd00084_c0
        cpd00051_c0
        cpd00132_c0
        cpd00041_c0
        cpd00035_c0
        cpd00053_c0
        cpd00023_c0
        cpd00033_c0
        cpd00156_c0
        cpd00129_c0
        cpd00066_c0
        """
        self.test = [self._model.metabolites.get_by_id(m) for m in amino_acids.split()]
        from cobra.core import Reaction

        for t in self.test:
            rxn_test = Reaction(f'test_{t.id}', f'Test {t.name}', 'TEST', 0, 0)
            rxn_test.add_metabolites({
                t: -1
            })
            self.test_reactions[t] = rxn_test
        self._model.add_reactions(self.test_reactions)

        return self.test_reactions

    def run_tests(self):
        test_solutions = {}
        for m, r in self.test_reactions.items():
            self._model.objective = r.id
            model_rxn = self._model.reactions.get_by_id(r.id)
            model_rxn.upper_bound = 1
            solution = cobra.flux_analysis.pfba(self._model)
            obj = solution.fluxes[model_rxn.id]
            if obj > 0 and solution.status == 'optimal':
                test_solutions[model_rxn.id] = solution
            else:
                test_solutions[model_rxn.id] = None
            model_rxn.upper_bound = 0

        return test_solutions
