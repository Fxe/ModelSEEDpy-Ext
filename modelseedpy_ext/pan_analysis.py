import logging
import os
import json
import subprocess
import pandas as pd
from pandas import DataFrame
from modelseedpy import RastClient, MSGenome
from modelseedpy.core.msgenome import normalize_role

logger = logging.getLogger(__name__)


class PanAnalysis:

    def __init__(self, analysis_id, genomes, kbase_mags):
        self.id = analysis_id
        self.genomes = genomes
        self.kbase_mags = kbase_mags
        self.kbase_mags_faa = set()
        if not os.path.exists(self.id):
            print(f'create folder {self.id}')
            os.mkdir(self.id)
        self.df_motu_output = None
        self.cluster_functions = None
        self.cluster_functions_uni = None
        self.cluster_functions_cft = None
        self.cluster_functions_type = None
        self.function_to_clusters = None
        self.completeness = {}

        self.seq_annotation = None
        self.bad_seqs = None
        self.genomes_objects = {}

    def index_kbase_mags(self, kbase_api):
        for object_id in self.kbase_mags:
            object_id = object_id + '.RAST' # dumb hack
            genome = kbase_api.get_from_ws(object_id, 132089)
            genome_faa_filename = f'{os.path.abspath(self.id)}/{genome.id}.faa'
            genome.to_fasta(genome_faa_filename)
            self.kbase_mags_faa.add(genome_faa_filename)
            self.genomes_objects[object_id] = genome

    def index_function_to_clusters(self):
        self.function_to_clusters = {}
        for cluster_id in self.cluster_functions:
            for gene_id in self.cluster_functions[cluster_id]:
                function_strings = self.cluster_functions[cluster_id][gene_id]
                if function_strings is not None:
                    for s in function_strings:
                        nmz_role = normalize_role(s)
                        if nmz_role not in self.function_to_clusters:
                            self.function_to_clusters[nmz_role] = set()
                        self.function_to_clusters[nmz_role].add(cluster_id)

    def run_motu(self) -> DataFrame:
        output_motu_file = f'{os.path.abspath(self.id)}/motu_output.txt'
        if not os.path.exists(output_motu_file):
            input_faa_file = f'{self.id}/motu_input.txt'
            if not os.path.exists(input_faa_file):
                logger.info(f'create mOTU input file {input_faa_file}')
                with open(input_faa_file, 'w') as fh:
                    for s in self.genomes:
                        fh.write(f'/home/fliu/mice/gtdb_prokka/{s}.faa\n')
                    for s in self.kbase_mags_faa:
                        fh.write(f'{s}\n')

            cmd = [
                "python3",
                "/home/fliu/python3/mOTUlizer/mOTUlizer/bin/mOTUpan.py",
                '--txt',
                '--faas', os.path.abspath(input_faa_file),
                '--output', output_motu_file
            ]

            logger.info(f'running mOTU {cmd}')
            ret = subprocess.run(cmd, capture_output=True)

        self.df_motu_output = pd.read_csv(output_motu_file, sep='\t', comment='#', index_col=0)

        return self.df_motu_output

    def get_genome_object(self, genome_id):
        if genome_id not in self.genomes_objects:
            f = f'/home/fliu/mice/gtdb_prokka/{genome_id}.faa'
            logger.debug(f'fetch {genome_id} {f}')

            genome = MSGenome.from_fasta(f, split=' # ')
            self.genomes_objects[genome_id] = genome
        return self.genomes_objects[genome_id]

    def process_df_motu_output(self):
        self.cluster_functions = {}
        self.cluster_functions_type = {}

        for row_id, d in self.df_motu_output.iterrows():
            self.cluster_functions_type[row_id] = d['type']
            functions = {}
            for genome_id in d['genomes'].split(';'):
                genome = self.get_genome_object(genome_id)
                for gene_id in d['genes'].split(';'):
                    if gene_id in genome.features:
                        feature = genome.features.get_by_id(gene_id)
                        if feature.seq in self.seq_annotation:
                            functions[gene_id] = set(self.seq_annotation[feature.seq])
                        elif feature.seq in self.bad_seqs:
                            functions[gene_id] = None
                        else:
                            print('error')
            self.cluster_functions[row_id] = functions

        self.cluster_functions_uni = {}
        self.cluster_functions_cft = {}
        for c_id in self.cluster_functions:
            gene_f = self.cluster_functions[c_id]
            all_f = set()
            for k in gene_f.values():
                if k is None:
                    all_f.add(None)
                else:
                    all_f |= k
            if len(all_f) == 1:
                self.cluster_functions_uni[c_id] = list(all_f)[0]
            else:
                self.cluster_functions_cft[c_id] = list(all_f)

    def run_rast(self):
        self.seq_annotation = {}
        self.bad_seqs = set()
        rast = RastClient()
        genome_ids = set(self.genomes) | set(self.genomes_objects)
        for genome_id in genome_ids:
            genome = self.get_genome_object(genome_id)
            output_rast_file = f'{os.path.abspath(self.id)}/{genome_id}.json'
            if os.path.exists(output_rast_file):
                logger.debug(f'file [{output_rast_file}] found, reading previous annotation')
                with open(output_rast_file, 'r') as fh:
                    annotation = json.load(fh)
                    for feature_id in annotation:
                        f = genome.features.get_by_id(feature_id)
                        for role in annotation[feature_id]:
                            f.add_ontology_term('RAST', role)
            else:
                logger.debug(f'file [{output_rast_file}] not found, running RAST')
                rast.annotate_genome(genome)
                logger.debug(f'file [{output_rast_file}] saving to local cache')
                with open(output_rast_file, 'w') as fh:
                    annotation = {}
                    for f in genome.features:
                        if 'RAST' in f.ontology_terms:
                            annotation[f.id] = sorted(list(f.ontology_terms['RAST']))
                    fh.write(json.dumps(annotation))

            for f in genome.features:
                if 'RAST' in f.ontology_terms:
                    if f.seq not in self.seq_annotation:
                        self.seq_annotation[f.seq] = set(f.ontology_terms['RAST'])
                else:
                    self.bad_seqs.add(f.seq)
