import logging
import os
import subprocess
import multiprocessing

logging.getLogger(__name__)


class MMseqs2:

    def __init__(self, genomes, path='.'):
        self.mmseqs_bin = 'mmseqs'
        self.path = path
        self.genomes = genomes

        self.feature_genome = {}
        for genome in self.genomes:
            for f in genome.features:
                if f.id not in self.feature_genome:
                    self.feature_genome[f.id] = genome
                else:
                    raise ValueError('genome feature_id not unique')

    def _build_master_faa(self, filepath=None, force=False):
        if filepath is None:
            filepath = self.path + '/master.faa'
        if force or not os.path.exists(filepath):
            with open(filepath, 'w') as fh:
                for genome in self.genomes:
                    for f in genome.features:
                        if len(f.seq) > 0:
                            fh.write(f'>{f.id}\n')
                            fh.write(f'{f.seq}\n')

        return filepath

    def cluster(self, master_faa_path=None, cluster_method='easy-cluster', threads=None, coverage=0.8, coverage_mode=0):
        """
        (1) a maximum E-value threshold (option -e [0,\infty[) computed according to the gap-corrected
        Karlin-Altschul statistics using the ALP library.

        (2) a minimum coverage (option -c [0,1], which is defined by the number of aligned residue pairs
        divided by either the maximum of the length of query/centre and target/non-centre
        sequences alnRes/max(qLen,tLen) (default mode, --cov-mode 0), or by the length of the
        target/non-centre sequence alnRes/tLen (--cov-mode 1), or by the length of the query/centre
        alnRes/qLen (--cov-mode 2). Read more about how coverage is computed in section How to set the
        right alignment coverage to cluster.

        (3) a minimum sequence identity (--min-seq-id [0,1]) with option --alignment-mode 3 defined as the
        number of identical aligned residues divided by the number of aligned columns including internal gap
        columns, or, by default, defined by a highly correlated measure, the equivalent similarity score of
        the local alignment (including gap penalties) divided by the maximum of the lengths of the two locally
        aligned sequence segments. The score per residue equivalent to a certain sequence identity is obtained
        by a linear regression using thousands of local alignments as training set.
        """
        if threads is None:
            threads = multiprocessing.cpu_count()

        mmseqs_bin = self.mmseqs_bin
        mmseqs_mode = cluster_method  # easy-cluster easy-linclust

        param_a = f'{self.path}/output/mmseqs'
        param_b = f'{self.path}/output'

        cmd = [self.mmseqs_bin, cluster_method,
               '--threads', str(threads),
               '--min-seq-id', '0.0',

               # coverage and coverage mode
               '--cov-mode', str(coverage_mode),
               '-c', str(coverage),

               master_faa_path, param_a, param_b]

        print(cmd)

        output = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )

        return output

    @staticmethod
    def read_cluster_file(cluster_filepath):
        clusters = {}
        with open(cluster_filepath, 'r') as fh:
            line = fh.readline()
            while line:
                rep, member = line.split()
                if rep not in clusters:
                    clusters[rep] = set()
                clusters[rep].add(member)
                line = fh.readline()

        return clusters

    def build_pangenome(self, pangenome_id: str, cluster_filepath, pangenome_name=None):
        from cobrakbase.core.kbasegenome.pangenome import KBasePangenome
        from cobrakbase.core.kbasegenome.ortholog_family import OrthologItem, OrthologFamily

        if pangenome_name is None:
            pangenome_name = pangenome_id

        clusters = self.read_cluster_file(cluster_filepath)

        res = KBasePangenome(pangenome_id, pangenome_name, '', [str(g.info) for g in self.genomes])

        orthologs = []
        for cluster_id, members in clusters.items():
            items = []
            rep = self.feature_genome[cluster_id].features.get_by_id(cluster_id)
            functions = {}
            for m in members:
                genome = self.feature_genome[cluster_id]
                feature = self.feature_genome[cluster_id].features.get_by_id(cluster_id)
                _f = '; '.join(feature.functions)
                if _f not in functions:
                    functions[_f] = 0
                functions[_f] += 1
                items.append(OrthologItem(m, 0, genome.info))
            selected_function = ''
            selected_function_count = 0
            for _f, _c in functions.items():
                if _c > selected_function_count:
                    selected_function_count = _c
                    selected_function = _f
            orthologs.append(OrthologFamily(cluster_id, items,
                                            function=selected_function,
                                            protein_translation=rep.protein_translation))

        res.add_orthologs(orthologs)

        return res
