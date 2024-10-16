from modelseedpy_ext.ani.skani import _read_search_output_as_parquet, _read_triangle_output_as_parquet
from modelseedpy_ext.ani import GenomeReference
import pandas as pd
import polars as pl


class ANIDao:

    def __init__(self):
        pass

    def get_ani(self, genome_ids):
        pass


class ExpMethod:

    def __init__(self, all_vs_all):
        self.clusters = {}
        self.all_vs_all = all_vs_all

    def compute_new_clusters(self, candidates):
        visited = set()
        remaining_candidates = set(candidates)
        while len(remaining_candidates) > 0:
            remaining_candidates = self._pop_representative_circun(self.clusters,
                                                                   remaining_candidates,
                                                                   visited)
        return visited

    def _get_best_quality(self, candidates):
        return list(candidates)[0]

    def _pop_representative_circun(self, clusters, candidates, visited):
        # select the best candidate for new representative
        rep = self._get_best_quality(candidates)

        # setup singleton cluster
        self.clusters[rep] = {rep: 100}  # self ani assume 100

        # get entities in the candidate pool that has ANI range of ani_radius
        for o in self.all_vs_all.filter(
                (pl.col("Query_file1") == rep) |
                (pl.col("Query_file2") == rep)
        ).rows(named=True):
            if o['Query_file1'] == rep:
                clusters[rep][o['Query_file2']] = o['ANI']
            else:
                clusters[rep][o['Query_file1']] = o['ANI']

        visited |= set(self.clusters[rep])
        return candidates - visited


class NrExp:
    MOCK_PRE_COMP_DATA = {
        'fastgenomics': (
            '/scratch/fliu/data/ani/library_morgan.txt',
            '/scratch/fliu/data/ani/morgan_vs_gtdb_r220_rep.out',
            '/scratch/fliu/data/ani/skani/morgan.triangle.slow.out'),
        'alexey': (
            '/scratch/fliu/data/ani/library_alexey.txt',
            '/scratch/fliu/data/ani/alexey_vs_gtdb_r220_rep.out',
            '/scratch/fliu/data/ani/skani/alexey.triangle.slow.out'),
        'grow': (
            '/scratch/fliu/data/ani/library_grow_kbase.txt',
            '/scratch/fliu/data/ani/skani/grow_vs_gtdb_r220_rep.slow.out',
            '/scratch/fliu/data/ani/grow_kbase_vs_self.out'),
        'mucc': (
            '/scratch/fliu/data/ani/library_mucc_kbase.txt',
            '/scratch/fliu/data/ani/mucc_kbase_vs_gtdb_r220_rep.out',
            '/scratch/fliu/data/ani/mucc_kbase_vs_self.out'),
        'anme': (
            '/scratch/fliu/data/ani/library_anme_kbase.txt',
            '/scratch/fliu/data/ani/skani/anme_vs_gtdb_r220_rep.slow.out',
            '/scratch/fliu/data/ani/skani/anme.triangle.slow.out'),
    }

    def __init__(self, ani_radius=95):
        self.path_to_exp_library = '/scratch/fliu/data/ani/skani/exp'
        self.rep_exp_counter = 0
        self.rep_clusters = {}
        self.member_to_cluster = {}
        self.dao_ani = None
        self.ani_method = None
        self.rep_exp_libraries = {}
        self.ani_radius = ani_radius

    def run_ani(self, query_library, rep_library):

        if rep_library == 'gtdb_r220_rep' and query_library in self.MOCK_PRE_COMP_DATA:
            file_vs_rep = self.MOCK_PRE_COMP_DATA[query_library][1]
            print('load pre calc dataset: ', file_vs_rep)
            return pl.from_arrow(_read_search_output_as_parquet(file_vs_rep, sep='\t'))
        elif query_library == rep_library and query_library in self.MOCK_PRE_COMP_DATA:
            file_vs_self = self.MOCK_PRE_COMP_DATA[query_library][2]
            print('load pre calc dataset: ', file_vs_self)
            return pl.from_arrow(_read_triangle_output_as_parquet(file_vs_self, sep='\t'))
        else:
            rep_library_file = self.rep_exp_libraries[rep_library][:-11]
            genome_files = self.MOCK_PRE_COMP_DATA[query_library][0]
            print(f'compute new ANI {genome_files} {rep_library_file}')
            ani_matrix, output, cmd = self.ani_method.distance(genome_files,
                                                               rep_library_file,
                                                               f'{self.path_to_exp_library}/{query_library}_{rep_library}.out')
            if ani_matrix is not None:
                return pl.from_arrow(ani_matrix)
            return ani_matrix

    def expand(self, query_library):
        # library_members = self.dao_ani.get_members_from_library(query_library)
        library_members = set(pd.read_csv(self.MOCK_PRE_COMP_DATA[query_library][0], header=None).to_dict()[0].values())
        print('library members:', len(library_members))

        added = set()
        candidates = set(library_members)

        for exp_library in self.rep_exp_libraries:
            if len(candidates) == 0:
                print('STOP all members assigned to clusters')
                break
            print(exp_library)

            df_vs_rep = self.run_ani(query_library, exp_library)
            if df_vs_rep is not None:
                print(f'number of ani pairs ({exp_library}):', len(df_vs_rep))

                # select all members with ANI >= radius
                to_previous_representatives = df_vs_rep.filter(
                    pl.col("ANI") >= self.ani_radius).group_by(
                    "Query_file", maintain_order=True).max()

                for o in to_previous_representatives.rows(named=True):
                    rep = o['Ref_file']
                    genome_id = o['Query_file']
                    if genome_id not in self.member_to_cluster:

                        if rep not in self.rep_clusters:
                            self.rep_clusters[rep] = {}
                        self.rep_clusters[rep][genome_id] = o['ANI']
                        self.member_to_cluster[genome_id] = rep
                        added.add(o['Query_file'])
                print('members added to rep clusters:', len(added))
                candidates -= added
                print('members candidate to new rep clusters:', len(candidates))
            else:
                print(f'found no matches against {exp_library}')

        if len(candidates) > 0:
            print('Generate new candidate clusters')

            df_all_vs_all = self.run_ani(query_library, query_library)
            df_all_vs_all = df_all_vs_all.filter(pl.col(
                "Query_file1").is_in(candidates) & pl.col(
                "Query_file2").is_in(candidates)).filter(pl.col("ANI") >= self.ani_radius)

            exp_method = ExpMethod(df_all_vs_all)
            visited = exp_method.compute_new_clusters(candidates)

            print('visited nodes:', len(visited))
            print('new clusters', len(exp_method.clusters))

            self.catalog_new_representatives(exp_method, query_library)

            self.rep_exp_counter += 1

    def catalog_new_representatives(self, exp_method, query_library):
        if len(exp_method.clusters) > 0:
            genome_files = [GenomeReference.from_file(x) for x in exp_method.clusters.keys()]
            exp_library_file = f'{self.path_to_exp_library}/iteration_{self.rep_exp_counter}.txt'
            ani_library = self.ani_method.build_library(genome_files, exp_library_file)
            self.rep_exp_libraries[query_library] = ani_library

            for k, members in exp_method.clusters.items():
                if k not in self.member_to_cluster:
                    if k not in self.rep_clusters:
                        self.rep_clusters[k] = {}
                        for m in members:
                            if m not in self.member_to_cluster:
                                self.member_to_cluster[m] = k
                                self.rep_clusters[k][m] = members[m]
            # self.rep_clusters.update(exp_method.clusters)
