from modelseedpy_ext.utils import progress
import pyarrow as pa
import pyarrow.parquet as pq
import os
from modelseedpy import MSGenome
import json


class PanGenomeAnalysis:

    def __init__(self):
        self.genomes = {}
        self.genome_to_clade = {}
        self.clade_to_genome = {}
        self.feature_to_genome = {}
        self.taxonomy = {}
        self.clade_motupan = {}

    def get_genome_faa(self, g_id, folder='/scratch/shared/CDM/salterns/data/All_MAGs_proteins'):
        genome = self.genomes.get(g_id)
        if genome is None:
            _genome_path = f'{folder}/{g_id}.faa'
            if os.path.exists(_genome_path):
                genome = MSGenome.from_fasta(_genome_path)
                self.genomes[g_id] = genome
                return genome
            else:
                return None
        return genome

    def load_motupan(self, motupan_folder='/scratch/shared/CDM/salterns/analysis/motupan/'):
        for f in os.listdir(motupan_folder):
            if f.endswith('.json'):
                clade_id = f.split('.motupan')[0]
                if clade_id not in self.clade_motupan:
                    with open(f'{motupan_folder}/{f}', 'r') as fh:
                        _motu_data = json.load(fh)
                        rnd_id = list(_motu_data)[0]
                        self.clade_motupan[clade_id] = _motu_data[rnd_id]
                else:
                    raise ValueError(f'duplicate clade id: {clade_id}')

    def get_clade_taxonomy(self, clade_id):
        counter = {}
        for genome_id in self.clade_to_genome[clade_id]:
            _tx = self.taxonomy.get(genome_id)
            if _tx:
                if _tx not in counter:
                    counter[_tx] = 0
                counter[_tx] += 1
        return counter

    def load_motupan_genomes(self):
        for clade_id in self.clade_motupan:
            self.clade_to_genome[clade_id] = set()
            _motu_data = self.clade_motupan[clade_id]
            clade_genomes = {x['name'] for x in _motu_data['genomes']}
            for genome_id in clade_genomes:
                genome = self.get_genome_faa(genome_id)
                if genome:
                    self.clade_to_genome[clade_id].add(genome_id)
                    self.genome_to_clade[genome_id] = clade_id

    def load_taxonomy(self):
        import pandas as pd
        _taxonomy = {}
        _taxonomy.update(
            pd.read_csv('/home/fliu/cliff_mags/gtdb_table_arc.csv', index_col=0).to_dict()['Classification'])
        _taxonomy.update(
            pd.read_csv('/home/fliu/cliff_mags/gtdb_table_bac.csv', index_col=0).to_dict()['Classification'])
        for genome_id in self.genomes:
            if genome_id + '__' in _taxonomy:
                self.taxonomy[genome_id] = _taxonomy[genome_id + '__']
            else:
                print('not found', genome_id)

    def load_all_genomes(self):
        with open('/home/fliu/cliff_mags/data/ani_library_derep_mags.txt', 'r') as fh:
            for l in fh.readlines():
                genome_id = l.strip().split('/')[-1][:-3]
                genome = self.get_genome_faa(genome_id)
                for f in genome.features:
                    if f.id not in self.feature_to_genome:
                        self.feature_to_genome[f.id] = genome_id
        with open('/home/fliu/cliff_mags/data/ani_library_rep_mags.txt', 'r') as fh:
            for l in fh.readlines():
                genome_id = l.strip().split('/')[-1][:-3]
                genome = self.get_genome_faa(genome_id)
                for f in genome.features:
                    if f.id not in self.feature_to_genome:
                        self.feature_to_genome[f.id] = genome_id
        with open('/home/fliu/cliff_mags/data/ani_library_ext_mags.txt', 'r') as fh:
            for l in fh.readlines():
                genome_id = l.strip().split('/')[-1][:-3]
                genome = self.get_genome_faa(genome_id)
                for f in genome.features:
                    if f.id not in self.feature_to_genome:
                        self.feature_to_genome[f.id] = genome_id


def read_ani():
    def read_ani_out(f):
        ani_score = {}
        with open(f, 'r') as fh:
            l = fh.readline()
            while (l):
                _p = l.split('\t')
                g1, g2, score, k1, k2 = _p
                if g1 not in ani_score:
                    ani_score[g1] = {}
                ani_score[g1][g2] = (score, k1, k2)
                l = fh.readline()
        return ani_score

    ani_score_gca = read_ani_out('./ani_gtdb_gca2.txt')
    ani_score_gcf = read_ani_out('./ani_gtdb_gcf2.txt')
    ani_score_self = read_ani_out('./ani_self2.txt')


def read_mmseqs(filename='./data/mmseqs/output/mmseqs_cluster.tsv'):
    feature_to_rep_seq_h = {}
    ct = 0
    with open(filename, 'r') as fh:
        l = fh.readline()
        a, b = [s.strip() for s in l.split('\t')]
        while l:
            if a:
                feature_to_rep_seq_h[b] = a
            l = fh.readline()
            if l:
                a, b = [s.strip() for s in l.split('\t')]
            ct += 1


def write(feature_to_rep_seq_h, output_file='/scratch/fliu/data/cliff/db/mmseqs_clusters.parquet'):
    col_seq_id = []
    col_cluster_id = []
    for s_id, c_id in progress(feature_to_rep_seq_h.items()):
        col_seq_id.append(bytes(bytearray.fromhex(s_id)))
        col_cluster_id.append(bytes(bytearray.fromhex(c_id)))
    table = pa.Table.from_arrays([pa.array(col_cluster_id), pa.array(col_seq_id)],
                                 names=['cluster_id', 'seq_id'])
    pq.write_table(table, output_file)
