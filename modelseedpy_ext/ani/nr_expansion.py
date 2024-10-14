from modelseedpy_ext.ani.skani import _read_search_output_as_parquet, _read_triangle_output_as_parquet
from modelseedpy_ext.ani import GenomeReference
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
        self.clusters[rep] = {rep}

        # get entities in the candidate pool that has ANI range of ani_radius
        for o in self.all_vs_all.filter(
                (pl.col("Query_file1") == rep) |
                (pl.col("Query_file2") == rep)
        ).rows(named=True):
            if o['Query_file1'] == rep:
                clusters[rep].add(o['Query_file2'])
            else:
                clusters[rep].add(o['Query_file1'])

        visited |= set(self.clusters[rep])
        return candidates - visited


class NrExp:

    def __init__(self):
        self.path_to_exp_library = '/scratch/fliu/data/ani/skani/exp'
        self.rep_exp_counter = 0
        self.rep_clusters = {}
        self.dao_ani = None
        self.ani_method = None
        self.rep_exp_libraries = {}
        self.ani_radius = 95

    def expand(self, exp_id, genome_files, file_vs_rep, file_vs_self):

        df_vs_rep = pl.from_arrow(_read_search_output_as_parquet(file_vs_rep, sep='\t'))
        print('number of ani pairs:', len(df_vs_rep))

        to_previous_representatives = df_vs_rep.filter(
            pl.col("ANI") >= self.ani_radius).group_by(
            "Query_file", maintain_order=True).max()

        added = set()
        for o in to_previous_representatives.rows(named=True):
            rep = o['Ref_file'].split('/')[-1]
            genome_id = o['Query_file'].split('/')[-1]
            if rep not in self.rep_clusters:
                self.rep_clusters[rep] = {}
            self.rep_clusters[rep][genome_id] = o['ANI']
            added.add(o['Query_file'])
        print('members added to rep clusters:', len(added))

        for exp_iteration, exp_library in self.rep_exp_libraries.items():
            self.ani_method.distance(genome_files,
                                     exp_library[:-11],
                                     f'{self.path_to_exp_library}/{exp_id}_iteration{exp_iteration}.out')
            print(exp_library)

        candidates_for_representatives = df_vs_rep.filter(~pl.col("Query_file").is_in(added))
        candidates = set(candidates_for_representatives['Query_file'])
        print('members candidate for new rep clusters:', len(candidates))

        df_all_vs_all = pl.from_arrow(_read_triangle_output_as_parquet(file_vs_self, sep='\t'))
        df_all_vs_all = df_all_vs_all.filter(pl.col(
            "Query_file1").is_in(candidates) & pl.col(
            "Query_file2").is_in(candidates)).filter(pl.col("ANI") >= self.ani_radius)

        exp_method = ExpMethod(df_all_vs_all)
        visited = exp_method.compute_new_clusters(candidates)

        print('new clusters', len(exp_method.clusters))

        self.catalog_new_representatives(exp_method)

        self.rep_exp_counter += 1

    def catalog_new_representatives(self, exp_method):
        if len(exp_method.clusters) > 0:
            genome_files = [GenomeReference.from_file(x) for x in exp_method.clusters.keys()]
            exp_library_file = f'{self.path_to_exp_library}/iteration_{self.rep_exp_counter}.txt'
            ani_library = self.ani_method.build_library(genome_files, exp_library_file)
            self.rep_exp_libraries[self.rep_exp_counter] = ani_library

            self.rep_clusters.update(exp_method.clusters)
