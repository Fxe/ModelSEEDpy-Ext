import os
from modelseedpy_ext.utils import progress
from modelseedpy_ext.re.core.genome import REAssembly


def _read_skani_out(filename):
    scores = {}
    with open(filename, 'r') as fh:
        header = fh.readline()
        line = fh.readline()
        while line:
            g1, g2, ani, af_ref, af_query, ref_name, query_name =  line.strip().split('\t')
            p = (g1.split('/')[-1], g2.split('/')[-1])
            if p not in scores:
                scores[p] = (float(ani), float(af_ref), float(af_query))
            line = fh.readline()
    return scores


def _get_max_score(scores):
    member_max_score = {}
    for _p in scores:
        g2 = _p[1]
        if g2 not in member_max_score:
            member_max_score[g2] = (0, None)
        score, af1, af2 = scores[_p]
        if score > member_max_score[g2][0]:
            member_max_score[g2] = (score, _p[0])
    return member_max_score


class RepGenomeExpansion:

    def __init__(self, ani_radius=95):
        self.ani_method = None
        self.representatives = None
        self.ani_radius = ani_radius

    def _get_best_quality(self, candidates):
        # candidates
        return list(candidates)[0]

    def _pop_representative_circun(self, clusters, candidates, visited):
        # select the best candidate for new representative
        rep = self._get_best_quality(candidates)

        # setup singleton cluster
        clusters[rep] = {rep}

        # get entities in the candidate pool that has ANI range of ani_radius
        clusters[rep] |= set(node_to_dest_score.get(rep, {})) - visited
        visited |= set(clusters[rep])
        return candidates - visited

    @staticmethod
    def max_score(scores):
        member_max_score = {}
        for _p in scores:
            g2 = _p[1]
            if g2 not in member_max_score:
                member_max_score[g2] = (0, None)
            score, af1, af2 = scores[_p]
            if score > member_max_score[g2][0]:
                member_max_score[g2] = (score, _p[0])

        return member_max_score

    def add_to_repsentative_clusters(self, ani_matrix):
        genome_to_representative = {}
        return genome_to_representative

    def aaaa(self, alexey_triag, candidates):
        edges = {}
        for g1 in alexey_triag.d_ani:
            for g2 in alexey_triag.d_ani[g1]:
                if g1 != g2 and g1 in candidates and g2 in candidates:
                    v = alexey_triag.d_ani[g1][g2]
                    score = float(v[0])
                    if score >= 95:
                        _p = (g1, g2)
                        _p_rev = (g2, g1)
                        if _p not in edges and _p_rev not in edges:
                            edges[_p] = score
        node_to_dest_score = {}
        for (g1, g2), score in edges.items():
            if g1 not in node_to_dest_score:
                node_to_dest_score[g1] = {}
            if g2 not in node_to_dest_score:
                node_to_dest_score[g2] = {}
            node_to_dest_score[g1][g2] = score
            node_to_dest_score[g2][g1] = score

    def expand(self, genome_ids, ani_circun=95):
        # run ani
        ani_matrix = self.ani_method.run(genome_ids, self.representatives)

        # add members to representative
        genome_to_representative = self.add_to_repsentative_clusters(ani_matrix)

        new_rep_candidates = set(genome_ids) - set(genome_to_representative)
        candidates = set(new_rep_candidates)

        clusters = {}
        visited = set()

        while len(candidates) > 0:
            candidates = self._pop_representative_circun(clusters, candidates, visited)
        pass


class NameThisLater:

    def __init__(self, genome_refs, kbase_api):
        self.kbase_api = kbase_api
        self.genome_refs = frozenset(genome_refs)

        self.genomes = None
        self.assemblies = None
        self.contigs = None
        self.genome_to_contig_h = None
        self.contig_filename = None

    @staticmethod
    def from_genome_refs(kbase_api, genome_refs):
        return NameThisLater(genome_refs, kbase_api)

    @staticmethod
    def from_genome_set(kbase_api, genome_set_id_or_ref, ws_id=None):
        genomes_set = kbase_api.get_from_ws(genome_set_id_or_ref, ws_id)
        # check type
        if not genomes_set.info.type == 'KBaseSearch.GenomeSet':
            raise ValueError(f'expected KBaseSearch.GenomeSet got {genomes_set.info.type}')

        genome_refs = [kbase_api.get_object_info(ref) for ref in genomes_set.refs]

        return NameThisLater(genome_refs, kbase_api)

    @staticmethod
    def from_workspaces(kbase_api, ws_ids):
        genome_refs = set()
        for ws_id in ws_ids:
            for o in kbase_api.list_workspace(ws_id):
                if o.type == 'KBaseGenomes.Genome':
                    genome_refs.add(o)
        return NameThisLater(genome_refs, kbase_api)

    def run(self):
        self.genomes = {}
        for ref in progress(self.genome_refs):
            genome = self.kbase_api.get_from_ws(str(ref))
            self.genomes[genome.info.id] = genome

        # make sure genome.id is same as object name
        for genome_id, genome in self.genomes.items():
            genome.id = genome_id

        self.assemblies = {}
        for genome_id, genome in progress(self.genomes.items()):
            assembly = self.kbase_api.get_from_ws(genome.assembly_ref)
            self.assemblies[genome_id] = assembly

        self.contigs = {}
        self.contig_filename = {}
        for genome_id, assembly in progress(self.assemblies.items()):
            fasta_handle_ref = assembly.fasta_handle_ref
            filename = f'/scratch/fliu/data/kbase/cache/handle/{fasta_handle_ref}'
            if not os.path.exists(filename):
                self.kbase_api.download_file_from_kbase2(fasta_handle_ref, filename)
            self.contigs[genome_id] = REAssembly.from_fasta(filename)
            self.contig_filename[genome_id] = filename

        self.genome_to_contig_h = {}
        for genome_id in progress(self.genomes):
            genome_contigs = self.contigs[genome_id]
            self.genome_to_contig_h[genome_id] = genome_contigs.hash_value

    def build_library(self, filename):
        """
        Build library file
        :param filename:
        :return:
        """
        h_to_file = {}
        for genome_id, h in self.genome_to_contig_h.items():
            if h not in h_to_file:
                h_to_file[h] = self.contig_filename[genome_id]

        library_files = set(h_to_file.values())

        with open(filename, 'w') as fh:
            for f in library_files:
                fh.write(f + '\n')

        return filename
