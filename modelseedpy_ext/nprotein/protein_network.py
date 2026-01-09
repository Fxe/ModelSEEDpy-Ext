import networkx as nx
from modelseedpy_ext.re.hash_seq import HashSeq


def convert_location(feature_location):
    contig, p0, strand, sz = feature_location
    start = p0
    end = start + sz - 1
    if strand == '-':
        end = p0
        start = end - sz + 1
    return contig, start, end, strand


class SortedContigFeatures:

    def __init__(self, contig_feature: dict):
        """

        :param contig_feature: dict with features sorted by low to high based no starting base position in contig
        """
        self.ctg_feature_sort = contig_feature

    @staticmethod
    def from_genome(genome, fn_location_conv):
        ctg_feature_sort = {}
        for f in genome.features:
            contig_id, start, end, strand = fn_location_conv(f)
            if contig_id not in ctg_feature_sort:
                ctg_feature_sort[contig_id] = {}
            if start not in ctg_feature_sort[contig_id]:
                ctg_feature_sort[contig_id][start] = f
            else:
                raise ValueError(f'more than 1 feature starts at position: {start}')

        return SortedContigFeatures(ctg_feature_sort)


class GenomeFeaturePairFreq:

    def __init__(self, fn_transform):
        self.genomes = {}
        self.sorted_features = {}
        self.pair_freq = {}
        self.fn_transform = fn_transform

    def add_genome(self, genome_id, genome, fn_location_conv):
        self.genomes[genome_id] = genome
        self.sorted_features[genome_id] = SortedContigFeatures.from_genome(genome, fn_location_conv)

        ctg_feature_sort = self.sorted_features[genome_id].ctg_feature_sort

        for contig, data in ctg_feature_sort.items():
            arr_sorted = [data[k] for k in sorted(data)]
            for f1, f2 in zip(arr_sorted, arr_sorted[1:]):
                p = (self.fn_transform(f1), self.fn_transform(f1))
                p_rev = (p[1], p[0])

                if p not in self.pair_freq:
                    self.pair_freq[p] = set()
                self.pair_freq[p].add(genome_id)


class ProteinNetworkFactory:

    def __init__(self):
        self.G = nx.DiGraph()
        pass

    def add_feature(self, feature):
        contig, start, end, strand = convert_location(feature.location[0])
        attributes = {
            'RAST': feature.ontology_terms.get('RAST'),
            'hash_seq': HashSeq(feature.seq).hash_value,
            'start': start,
            'end': end,
            'strand': strand,
            'contig': contig
        }
        self.add_node(feature.id, attributes)

    def get_k_neighbors(self, protein_id, k=3):
        n_back = []
        n_forward = []
        curr = protein_id
        for i in range(k):
            res = self.get_node_incoming_edges_by_type(curr, 'next')
            if len(res) > 0:
                _left, _ = res[0]
                n_back.append(_left)
                curr = _left
            else:
                n_back.append('*')
                break
        curr = protein_id
        for i in range(k):
            res = self.get_node_outgoing_edges_by_type(curr, 'next')
            if len(res) > 0:
                _, _right = res[0]
                n_forward.append(_right)
                curr = _right
            else:
                n_forward.append('*')
                break
        n_back.reverse()
        return n_back, n_forward

    def add_genome_proteins(self, genome):
        loc_d = {}
        for feature in genome.features:
            if len(feature.location) > 1:
                raise ValueError('feature with > 1 location')
            else:
                self.add_feature(feature)
                contig, start, end, strand = convert_location(feature.location[0])
                if start not in loc_d:
                    loc_d[start] = []
                loc_d[start].append((feature, end))
        sorted_pos = sorted(loc_d.keys())
        for i in range(len(sorted_pos)):
            current_pos = sorted_pos[i]
            for feature, end in loc_d[current_pos]:
                for j in range(i+1, len(sorted_pos)):
                    next_pos = sorted_pos[j]
                    stop = False
                    if next_pos > end:
                        stop = True
                        for other_feature, other_end in loc_d[next_pos]:
                            if other_feature.id != feature.id:
                                self.add_edge(feature.id, other_feature.id, 'next')
                    if stop:
                        break
                

    def add_node(self, protein_id, attributes=None):
        if attributes is None:
            attributes = {}
        self.G.add_node(protein_id, **attributes)

    def add_edge(self, protein_id1, protein_id2, edge_type):
        self.G.add_edge(protein_id1, protein_id2, edge_type=edge_type)

    def get_node_by_id(self, protein_id):
        return self.G.nodes[protein_id]
    
    def get_protein_neighbors(self, protein_id):
        return list(self.G.neighbors(protein_id))
    
    def get_node_incoming_edges(self, protein_id):
        return list(self.G.in_edges(protein_id))
    
    def get_node_outgoing_edges(self, protein_id):
        return list(self.G.out_edges(protein_id))
    
    def get_node_incoming_edges_by_type(self, protein_id, edge_type):
        return [e for e in self.G.in_edges(protein_id) if self.G.edges[e]['edge_type'] == edge_type]
    
    def get_node_outgoing_edges_by_type(self, protein_id, edge_type):
        return [e for e in self.G.out_edges(protein_id) if self.G.edges[e]['edge_type'] == edge_type]
    
    def get_protein_neighbors_by_type(self, protein_id, edge_type):
        return [n for n in self.G.neighbors(protein_id) if self.G.edges[protein_id, n]['edge_type'] == edge_type]
    
    def get_protein_neighbors_by_type_and_direction(self, protein_id, edge_type, direction):
        return [n for n in self.G.neighbors(protein_id) if self.G.edges[protein_id, n]['edge_type'] == edge_type and self.G.edges[protein_id, n]['direction'] == direction]
