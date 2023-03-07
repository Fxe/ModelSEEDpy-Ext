from modelseedpy.core.msatpcorrection import MSATPCorrection
from modelseedpy.core.msgenome import normalize_role


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
    
    def __init__(self, master, medias, atp_hydrolysis_id, extracell='e0'):
        self.master = master
        self.atp_hydrolysis_id = atp_hydrolysis_id
        self.medias = medias
        self.comp = {'EX_cpd00001_e0': 1000, 'EX_cpd00067_e0': 1000}
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
        #print(media.id, model.medium)
        sol = model.optimize()
        return sol
    
    def get_metabolite_exchanges(self, model):
        metabolite_exchanges = {}
        for rxn_exchange in model.exchanges:
            metabolites = rxn_exchange.metabolites
            if len(metabolites) == 1:
                metabolite_exchanges[list(metabolites)[0].id] = rxn_exchange
            else:
                print('ignore', rxn_exchange)
        return metabolite_exchanges
    
    def test_core_model(self, model, medias, v=False, cmp='e0'):
        res = {}
        media_out = {}
        metabolite_exchanges = self.get_metabolite_exchanges(model)
        for media_id in medias:
            #print(media_id)
            sol = self.test_media(model, medias[media_id], metabolite_exchanges, cmp)
            media_out[media_id] = None
            res[media_id] = 0
            if sol and sol.status == 'optimal':
                ex_out = dict(filter(lambda x: x[0].startswith('EX_') and x[1] != 0, sol.fluxes.to_dict().items()))
                media_out[media_id] = ex_out
                #print(media_id, sol.objective_value, ex_out)
                atp_val = sol.objective_value
                if atp_val > 0:
                    for cpd_id in medias[media_id].get_media_constraints(cmp):
                        media_cpd_val = sol.fluxes[metabolite_exchanges[cpd_id].id]
                        if media_cpd_val == 0:
                            atp_val = 0
                            #print('reject', media_id)
                res[media_id] = atp_val
        return res, media_out


class AtpCoreGapfillProfiler:
    
    def __init__(self, master, atp_hydrolysis_id):
        self.master = master
        self.atp_hydrolysis_id = atp_hydrolysis_id
        self.atp_basic_medias = None
        self.atp_medias = None
        self.res = {}
        self.template = None
        self.min_gap = {
            'Glc/O2': 5,
            'Etho/O2': 0.01,
            'Ac/O2': 1,
            'Pyr/O2': 3,
            'Glyc/O2': 2,
            'Fum/O2': 3,
            'Succ/O2': 2,
            'Akg/O2': 2,
            'LLac/O2': 2,
            'Dlac/O2': 2,
            'For/O2': 2,
            'For/NO3': 1.5,
            'Pyr/NO': 2.5,
            'Pyr/NO2': 2.5,
            'Pyr/NO3': 2.5,
            'Pyr/SO4': 2.5,
            'Glc/DMSO': 5,
            'Glc/TMAO': 5,
        }
    
    def profile_genome(self, genome_id):
        result = {}
        model = self.master.fn_get_model(genome_id)
        atp_method = MSATPCorrection(model, self.template, self.atp_basic_medias, atp_hydrolysis_id=self.atp_hydrolysis_id)
        atp_method.evaluate_growth_media()
        for media in atp_method.media_gapfill_stats:
            result[media.name] = atp_method.media_gapfill_stats[media]
        for media_id in self.min_gap:
            if media_id in self.atp_medias:
                res_min_gap = atp_method.msgapfill.run_gapfilling(self.atp_medias[media_id], self.atp_hydrolysis_id, minimum_obj=self.min_gap[media_id], binary_check=True)
                result[media_id] = res_min_gap
        return result
    
    def profile(self):
        pass


class GeProf:
    
    def __init__(self, genome_ids, fn_get_model, fn_get_genome, genome_order = None):
        self.genome_order = None
        self.genome_ids = genome_ids
        self.fn_get_model = fn_get_model
        self.fn_get_genome = fn_get_genome
    
    def profile(self):
        pass
