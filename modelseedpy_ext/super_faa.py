from modelseedpy import MSGenome
from modelseedpy.core.msgenome import MSFeature
from modelseedpy import RastClient
from modelseedpy_ext.re.hash_seq import HashSeq


class SuperFaa:

    def __init__(self):
        self.genomes = {}

        self.genome_feature_hash = None
        self.super_faa = None

        self.blank_annotation = None
        self.annotated_genes = None
        self.done = None
        self.rast_batch_size = 10000

    def aaa(self):
        self.genome_feature_hash = {}
        self.super_faa = {}
        for genome_id in self.genomes:
            self.genome_feature_hash[genome_id] = {}
            genome = self.genomes[genome_id]
            for f in genome.features:
                if len(f.seq) > 0:
                    hseq = HashSeq(f.seq)
                    self.genome_feature_hash[genome_id][f.id] = hseq.hash_value
                    if hseq.hash_value not in self.super_faa:
                        self.super_faa[hseq.hash_value] = f.seq

    def write(self, filename):
        with open(filename, 'w') as fh:
            for h, seq in self.super_faa.items():
                fh.write(f'>{h}\n')
                fh.write(f'{seq}\n')

    def _rast_batch(self, features: list, rast):
        _g = MSGenome()
        _g.add_features(features)
        rast.annotate_genome(_g)
        for _a in _g.features:
            blank = True
            for _rast in _a.ontology_terms.get('RAST', []):
                if _rast not in self.annotated_genes:
                    self.annotated_genes[_rast] = set()
                blank = False
                self.annotated_genes[_rast].add(_a.id)
            if blank:
                self.blank_annotation.add(_a.id)
            self.done.add(_a.id)

    def rast(self):
        self.blank_annotation = set()
        self.annotated_genes = {}
        self.done = set()

        rast = RastClient()

        _f = []
        for h, seq in self.super_faa.items():
            if h not in self.done and h not in self.blank_annotation:
                _f.append(MSFeature(h, seq))
                if len(_f) >= self.rast_batch_size:
                    self._rast_batch(_f, rast)
                    _f = []

        if len(_f) > 0:
            self._rast_batch(_f, rast)

# junk

import math
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


class LinearNmz:

    def __init__(self, v_min, v_max, m=1000):
        self.v_min = v_min
        self.v_max = v_max
        self.m = m

    def f(self, x):
        if x < self.v_min:
            raise ValueError('-')
        if x > self.v_max:
            raise ValueError('+')
        return ((x - self.v_min) / (self.v_max - self.v_min)) * self.m

    def get_color_gradient(self, a, b):
        cmap_name = 'custom_gradient'
        cmap = LinearSegmentedColormap.from_list(cmap_name, [a, b], N=self.m)
        return cmap


class GrowMeta:

    def __init__(self,
                 big_file='/scratch/fliu/data/GROW/Genome-Resolved-Open-Watersheds-database-GROWdb/Geospatial Data/GROWdb_with_vars_20220715.csv'):
        self.sample_data = {}
        self.sample_lu = {}
        self.df_meta = pd.read_csv(big_file)
        for row_id, d in self.df_meta.iterrows():
            self.sample_data[d['SampleName']] = d

        self.df_lu = pd.read_csv('./data/LandUseCat.csv', index_col='ID')
        for row_id, d in self.df_lu.iterrows():
            self.sample_lu[row_id] = d

    def loc_key(self, k):
        if k not in self.sample_data:
            if k.startswith('KC:'):
                k = k.split('_', 1)[1]
            return k.replace('.', '-').split('_MT_')[0]
        else:
            return k

    def get_meta(self, sample_id):
        return self.sample_data[self.loc_key(sample_id)]

    def get_lu(self, sample_id):
        return self.sample_lu[self.loc_key(sample_id)]

    def get_meta_value(self, sample_id, prop):
        v = self.sample_data[self.loc_key(sample_id)][prop]
        if pd.isna(v):
            return None
        return v

    def get_lu_value(self, sample_id, prop):
        v = self.sample_lu[self.loc_key(sample_id)][prop]
        if pd.isna(v):
            return None
        return v


class CHenryData:

    def __init__(self, filename, models, exclude=('id', 'name')):
        if exclude is None:
            self.exclude_columns = set()
        else:
            self.exclude_columns = set(exclude)
        self.df = pd.read_csv(filename, index_col=0)
        self.model_ids = {x for x in self.df.columns if x not in self.exclude_columns}
        self.model_indv = {'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'}
        self.model_indv_data = {x: {y: {} for y in self.model_ids} for x in self.model_indv}
        self.compartments = set()
        self.features = set()
        self.feature_column = {}
        self.zero_rows = set()
        self.models = models
        self.feature_order = None
        self.model_ids_order = None
        self.bio_to_c_index = None

        for row_id, d in self.df.iterrows():
            _p = d['id'].split('_')
            cmp = _p[-1]
            self.compartments.add(cmp)
            total = 0
            if cmp.startswith('c') and _p[0].startswith('rxn'):
                self.features.add(_p[0])
            for m_id in self.model_ids:
                total += math.fabs(d[m_id])
            if total == 0:
                self.zero_rows.add(d['id'])
            if cmp.startswith('c') and _p[0].startswith('rxn'):
                for k in self.model_ids:
                    self.model_indv_data[cmp][k][_p[0]] = d[k]

    def get_model_organism_nmz_coef(self):
        model_organism_nmz_coef = {c_index: {} for c_index in self.models}
        for row_id, d in self.df.iterrows():
            rxn_id = d['id']
            if rxn_id in self.bio_to_c_index:
                # cmp = _p[-1]
                for model_id in self.model_ids:
                    model_organism_nmz_coef[self.bio_to_c_index[rxn_id]][model_id] = d[model_id]
                #    model_organism_relab_prot[model_id][cmp] = d[model_id]
        return model_organism_nmz_coef

    def write_model_nmz_file2(self, c_index, filename='/scratch/fliu/data/GROW/c1_nmz.tsv'):
        model_organism_nmz_coef = self.get_model_organism_nmz_coef()
        with open(filename, 'w') as fh:
            l = list(self.model_indv_data[c_index])
            fh.write('reaction_id\t' + '\t'.join(l) + '\n')
            for f in self.feature_order:
                data = [f]
                for model_id in self.model_ids_order:
                    coef = model_organism_nmz_coef[c_index][model_id]
                    # coef = 1
                    v = self.model_indv_data[c_index][model_id].get(f, 0)
                    data.append(v / coef)
                fh.write('\t'.join([str(x) for x in data]) + '\n')

    def write_model_nmz_file(self, model_organism_nmz_coef, filename='/scratch/fliu/data/GROW/c1_nmz.tsv'):
        with open(filename, 'w') as fh:
            c_index = 'c1'
            l = list(self.model_indv_data[c_index])
            fh.write('reaction_id\t' + '\t'.join(l) + '\n')
            for f in self.feature_order:
                data = [f]
                for model_id in self.model_ids_order:
                    coef = model_organism_nmz_coef[c_index][model_id]
                    # coef = 1
                    v = self.model_indv_data[c_index][model_id].get(f, 0)
                    data.append(v / coef)
                fh.write('\t'.join([str(x) for x in data]) + '\n')
