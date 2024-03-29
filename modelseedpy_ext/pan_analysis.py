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

    DEFAULT_COMPLETENESS = 95
    DEFAULT_CONTAMINATION = 0

    def __init__(self, analysis_id, genomes, kbase_mags):
        self.id = analysis_id
        self.genomes = genomes
        self.kbase_mags = kbase_mags
        self.kbase_mags_faa = {}
        if not os.path.exists(self.id):
            print(f"create folder {self.id}")
            os.mkdir(self.id)
        self.df_motu_output = None
        self.cluster_functions = None
        self.cluster_functions_uni = None
        self.cluster_functions_cft = None
        self.cluster_functions_type = None
        self.function_to_clusters = None
        self.contamination = {}
        self.completeness = {}
        self.motu_ret = None
        self.fn_get_genome = None
        self.motu_path = '/home/fliu/python3/mOTUlizer/mOTUlizer/bin/'
        self.temp_folder = None

        self.seq_annotation = None
        self.bad_seqs = None
        self.genomes_objects = {}

    def set_temp_folder(self, path):
        if not os.path.exists(path):
            logger.info(f"mkdirs path [{path}]")
            os.mkdir(path)
        if os.path.isdir(path):
            self.temp_folder = path
        else:
            raise ValueError(f'Path exists and not a folder: [{path}]')

    def index_kbase_mags(self, kbase_api):
        for object_id in self.kbase_mags:
            object_id = object_id + ".RAST"  # dumb hack
            genome = kbase_api.get_from_ws(object_id, 132089)
            genome_faa_filename = f"{os.path.abspath(self.id)}/{genome.id}.faa"
            if self.temp_folder:
                genome_faa_filename = f"{self.temp_folder}/{genome.id}.faa"
            if not os.path.exists(genome_faa_filename):
                logger.info(f"write FAA file [{genome_faa_filename}]")
                genome.to_fasta(genome_faa_filename)
            self.kbase_mags_faa[object_id] = genome_faa_filename
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

    @staticmethod
    def _get_genome_faa_file_path(o):
        f = '/scratch/gtdb_r207/faa_files/'
        if o.startswith('GB_GCA_'):
            p = o.split('GB_GCA_')
            genome_id = f'{o[3:]}_protein.faa'
            f += f'GCA/{p[1][0:3]}/{p[1][3:6]}/{p[1][6:9]}/{genome_id}'
        elif o.startswith('RS_GCF_'):
            p = o.split('RS_GCF_')
            genome_id = f'{o[3:]}_protein.faa'
            f += f'GCF/{p[1][0:3]}/{p[1][3:6]}/{p[1][6:9]}/{genome_id}'
        else:
            raise ValueError(f'Error id: [{o}]')
        return f

    def make_motu_input_file(self, motu_input="motu_input.txt"):
        output_folder = self.temp_folder if self.temp_folder else os.path.abspath(self.id)
        input_faa_file = f"{output_folder}/{motu_input}"
        genome_ids = set()

        lines = []
        for o in self.genomes:
            f = '/scratch/gtdb_r207/faa_files/'
            if o.startswith('GB_GCA_'):
                p = o.split('GB_GCA_')
                genome_id = f'{o[3:]}_protein.faa'
                f += f'GCA/{p[1][0:3]}/{p[1][3:6]}/{p[1][6:9]}/{genome_id}'
            elif o.startswith('RS_GCF_'):
                p = o.split('RS_GCF_')
                genome_id = f'{o[3:]}_protein.faa'
                f += f'GCF/{p[1][0:3]}/{p[1][3:6]}/{p[1][6:9]}/{genome_id}'
            else:
                raise ValueError(f'Error id: [{o}]')
            genome_ids.add(f.split('/')[-1])
            lines.append(f)
        for o in self.kbase_mags_faa:
            genome_ids.add(o.split('/')[-1])
            lines.append(o)
        for k in lines:
            if not os.path.exists(k) or not os.path.isfile(k):
                raise ValueError(f'Genome FAA missing: [{o}]')

        with open(input_faa_file, "w") as fh:
            for line in lines:
                fh.write(line + "\n")

        return genome_ids

    def make_motu_completeness_file(self, genome_ids, motu_completeness="motu_completeness_input.txt"):
        output_folder = self.temp_folder if self.temp_folder else os.path.abspath(self.id)
        input_completeness_file = f"{output_folder}/{motu_completeness}"
        with open(input_completeness_file, "w") as fh:
            fh.write("Bin Id\tCompleteness\tContamination\n")
            for genome_id in genome_ids:
                completeness = self.completeness.get(
                    genome_id, self.DEFAULT_COMPLETENESS
                )
                contamination = self.contamination.get(
                    genome_id, self.DEFAULT_CONTAMINATION
                )
                fh.write(f"{genome_id}\t{completeness}\t{contamination}\n")

    def run_motu_python(self, threads=None, max_it=20, method='default', clean_after=False, quiet=True):
        from mOTUlizer.classes.mOTU import mOTU
        faas = {s: self._get_genome_faa_file_path(s) for s in self.genomes}
        faas.update(self.kbase_mags_faa)

        gene_clusters_dict = {}

        checkm = { k: self.completeness.get(k, self.DEFAULT_COMPLETENESS) for k in faas.keys()}

        precluster = False
        output = f'{self.temp_folder}/mmseqs'

        if not threads:
            import multiprocessing
            threads = multiprocessing.cpu_count()

        motu = mOTU(
            name=self.id,
            faas=faas,
            gene_clusters_dict=gene_clusters_dict,
            genome_completion_dict=checkm,
            max_it=max_it,
            threads=threads,
            precluster=precluster,
            method=method,
            output=output,
            clean_after=clean_after,
            quiet=quiet
        )

        from io import StringIO
        io = StringIO(motu.pretty_pan_table())
        df = pd.read_csv(io, sep="\t", comment="#", index_col=0)

        return df

    def run_motu(
        self,
        force_rerun=False,
        motu_input="motu_input.txt",
        motu_output="motu_output.txt",
        motu_completeness="motu_completeness_input.txt",
        env=None
    ) -> DataFrame:
        output_folder = self.temp_folder if self.temp_folder else os.path.abspath(self.id)
        output_motu_file = f"{output_folder}/{motu_output}"
        if not os.path.exists(output_motu_file) or force_rerun:
            input_completeness_file = f"{output_folder}/{motu_completeness}"
            input_faa_file = f"{output_folder}/{motu_input}"
            if not os.path.exists(input_faa_file):
                raise ValueError(f'Input file does not exists: [{input_faa_file}]')
            if not os.path.exists(input_completeness_file):
                raise ValueError(f'Input completeness file does not exists: [{input_completeness_file}]')

            cmd = [
                "python3",
                f"{self.motu_path}mOTUpan.py",
                "--checkm",
                input_completeness_file,
                "--txt",
                "--faas",
                os.path.abspath(input_faa_file),
                "--output",
                output_motu_file,
            ]

            logger.info(f"running mOTU {cmd}")
            if env:
                self.motu_ret = subprocess.run(cmd, capture_output=True, env=env)
            else:
                self.motu_ret = subprocess.run(cmd, capture_output=True)

        self.df_motu_output = pd.read_csv(
            output_motu_file, sep="\t", comment="#", index_col=0
        )

        return self.df_motu_output

    def get_genome_object(self, genome_id):
        if self.fn_get_genome:
            return self.fn_get_genome(genome_id)
        if genome_id not in self.genomes_objects:
            f = f"/home/fliu/mice/gtdb_prokka/{genome_id}.faa"
            logger.debug(f"fetch {genome_id} {f}")

            genome = MSGenome.from_fasta(f, split=" # ")
            self.genomes_objects[genome_id] = genome
        return self.genomes_objects[genome_id]

    def process_df_motu_output(self):
        self.cluster_functions = {}
        self.cluster_functions_type = {}

        for row_id, d in self.df_motu_output.iterrows():
            self.cluster_functions_type[row_id] = d["type"]
            functions = {}
            for genome_id in d["genomes"].split(";"):
                genome = self.get_genome_object(genome_id)
                for gene_id in d["genes"].split(";"):
                    if gene_id in genome.features:
                        feature = genome.features.get_by_id(gene_id)
                        if feature.seq in self.seq_annotation:
                            functions[gene_id] = set(self.seq_annotation[feature.seq])
                        elif feature.seq in self.bad_seqs:
                            functions[gene_id] = None
                        else:
                            print("error")
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
            output_rast_file = f"{os.path.abspath(self.id)}/{genome_id}.json"
            if os.path.exists(output_rast_file):
                logger.debug(
                    f"file [{output_rast_file}] found, reading previous annotation"
                )
                with open(output_rast_file, "r") as fh:
                    annotation = json.load(fh)
                    for feature_id in annotation:
                        f = genome.features.get_by_id(feature_id)
                        for role in annotation[feature_id]:
                            f.add_ontology_term("RAST", role)
            else:
                logger.debug(f"file [{output_rast_file}] not found, running RAST")
                rast.annotate_genome(genome)
                logger.debug(f"file [{output_rast_file}] saving to local cache")
                with open(output_rast_file, "w") as fh:
                    annotation = {}
                    for f in genome.features:
                        if "RAST" in f.ontology_terms:
                            annotation[f.id] = sorted(list(f.ontology_terms["RAST"]))
                    fh.write(json.dumps(annotation))

            for f in genome.features:
                if "RAST" in f.ontology_terms:
                    if f.seq not in self.seq_annotation:
                        self.seq_annotation[f.seq] = set(f.ontology_terms["RAST"])
                else:
                    self.bad_seqs.add(f.seq)
