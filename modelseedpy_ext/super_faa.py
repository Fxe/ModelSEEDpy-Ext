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
                    if hseq.hash_value not in self.super_faa:
                        self.genome_feature_hash[genome_id][f.id] = hseq.hash_value
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
