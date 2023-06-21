import os
import subprocess


class MMseqs2:

    def __init__(self, genomes, path='.'):
        self.mmseqs_bin = 'mmseqs'
        self.path = path
        self.genomes = genomes

    def _build_master_faa(self, filepath, force=False):
        if force or not os.path.exists(filepath):
            with open(filepath, 'w') as fh:
                for genome in self.genomes:
                    for f in genome.features:
                        fh.write(f'>{f.id}\n')
                        fh.write(f'{f}\n')

    def cluster(self, cluster_method='easy-cluster', threads=None, coverage=0.8, coverage_mode=0):
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
            threads = 1

        output_folder = self.path
        mmseqs_bin = self.mmseqs_bin
        mmseqs_mode = cluster_method  # easy-cluster easy-linclust
        param_threads = '--threads ' + threads
        param_min_seq_id = '0.0'

        # coverage and coverage mode
        param_cov_mode = coverage_mode
        param_c = coverage

        param_master_faa = f'{output_folder}/master.faa'
        param_a = f'{output_folder}/master_set_mmseqs/mmseqs_'
        param_b = f'{output_folder}/master_set_mmseqs'

        cmd = [self.mmseqs_bin, cluster_method, param_threads, param_min_seq_id, param_cov_mode, param_c,
               param_master_faa, param_a, param_b]

        output = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )

        return output
