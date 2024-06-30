import logging
from Bio.Seq import Seq
from modelseedpy_ext.re.etl.transform_graph import TransformGraph, Node
from modelseedpy_ext.re.etl.etl_transform_contig_set import ETLTransformContigSet
from modelseedpy_ext.re.etl.transform_kbase_object import TransformKBaseObjectInfo
from modelseedpy_ext.utils import progress
from modelseedpy_ext.re.utils import count_symbols
from modelseedpy_ext.re.hash_seq import HashSeq, HashSeqList

logger = logging.getLogger(__name__)


class SeqSet:

    def __init__(self, label, label_set, relationship_seq_set, graph):
        self.hash_list = HashSeqList()
        self.seq_nodes = {}
        self.seq_counts = {}
        self.label = label
        self.label_set = label_set
        self.graph = graph
        self.relationship_seq_set = relationship_seq_set

    def build_seq_node(self, seq):
        hash_seq = HashSeq(seq)
        node = Node(hash_seq.hash_value, self.label, data={"size": len(hash_seq), "symbols": count_symbols(hash_seq)})

        return hash_seq, node

    def add(self, seq):
        hash_seq, node = self.build_seq_node(seq)

        if self.graph is not None:
            node = self.graph.add_transform_node2(node)

        self.hash_list.append(hash_seq)
        if node.id not in self.seq_nodes:
            self.seq_nodes[node.id] = node
            self.seq_counts[node.id] = 0
        self.seq_counts[node.id] += 1

        return node

    def build_set_node(self):
        node_seq_set = Node(self.hash_list.hash_value, self.label_set, data={'size': len(self.hash_list)})
        if self.graph is not None:
            node_seq_set = self.graph.add_transform_node2(node_seq_set)

            for i, n in self.seq_nodes.items():
                count = self.seq_counts[i]
                data = {}
                if count > 1:
                    data['copies'] = count
                self.graph.add_transform_edge2(node_seq_set, n, self.relationship_seq_set, data)

        return node_seq_set


class DriverIndexKBaseGenomes:

    def __init__(self, seq_loader_contigs, seq_loader_dna, seq_loader_protein):
        self.seq_loader_contigs = seq_loader_contigs
        self.seq_loader_protein = seq_loader_protein
        self.seq_loader_dna = seq_loader_dna
        self.etl_contig_set = ETLTransformContigSet(seq_loader_contigs)
        self.etl_kbase_info = TransformKBaseObjectInfo()
        self.graph = TransformGraph()
        self.fn_get_assembly_from_genome = None
        self.fn_get_contigs_from_assembly = None
        self.fn_get_contig_by_id = None
        self.ncbi_protein_mapping = {}

    def index_info(self, info):
        node_kb_type = self.graph.add_transform_node2(TransformKBaseObjectInfo.transform_node_type(info))
        node_kb_object = self.graph.add_transform_node2(TransformKBaseObjectInfo.transform_node_object(info))
        self.graph.add_transform_edge2(node_kb_object,
                                       node_kb_type,
                                       "kbase_object_has_type",
                                       data={'type_version': info.type_version})

        return node_kb_object, node_kb_type

    def get_genome_assembly(self, genome):
        return self.fn_get_assembly_from_genome(genome)

    def get_contigs_from_assembly(self, assembly):
        return self.fn_get_contigs_from_assembly(assembly)

    def index_kbase_genomes(self, genomes):
        for kb_genome in progress(genomes):
            node_kb_genome_object, node_kb_type = self.index_info(kb_genome.info)

            kb_assembly = self.get_genome_assembly(kb_genome)

            node_kb_assembly_object, node_kb_type = self.index_info(kb_assembly.info)

            self.graph.add_transform_edge2(node_kb_genome_object,
                                           node_kb_assembly_object,
                                           "kbase_object_has_object_reference",
                                           data={'value': 'assembly_ref'})

            seq_contigs = self.get_contigs_from_assembly(kb_assembly)

            node_contig_set = self.etl_contig_set.transform_to_graph(seq_contigs, self.graph)

            contig_hash_to_ids = {}

            node_seq_set_protein, node_seq_set_dna, ncbi_protein_mapping = \
                self.transform_features(kb_genome.features, seq_contigs, self.graph)

            self.graph.add_transform_edge2(node_kb_assembly_object,
                                           node_contig_set,
                                           "kbase_object_has_re_contig_set")
            self.graph.add_transform_edge2(node_kb_genome_object,
                                           node_seq_set_protein,
                                           "kbase_object_has_re_seq_protein_set")
            self.graph.add_transform_edge2(node_kb_genome_object,
                                           node_seq_set_dna,
                                           "kbase_object_has_re_seq_dna_set")

            self.ncbi_protein_mapping.update(ncbi_protein_mapping)

            self.graph.add_transform_edge2(node_contig_set,
                                           node_seq_set_dna,
                                           "re_contig_set_has_re_dna_set",
                                           data={'value': 'assembly_ref'})

    def transform_features(self, features, seq_contigs, graph=None):

        if graph is None:
            graph = TransformGraph()

        contig_id_to_h = {}
        for contig_feature in seq_contigs.features:
            contig_id_to_h[contig_feature.id] = HashSeq(contig_feature.seq).hash_value

        ncbi_protein_mapping = {}
        seq_set_dna = SeqSet("re_seq_dna", "re_seq_dna_set", "re_seq_dna_set_has_re_seq_dna", graph)
        seq_set_protein = SeqSet("re_seq_protein", "re_seq_protein_set", "re_seq_protein_set_has_re_seq_protein", graph)
        for feature in features:
            node_re_seq_dna = None
            node_re_seq_protein = None
            if len(feature.dna_sequence) > 0:
                self.seq_loader_dna.load(feature.dna_sequence)
                node_re_seq_dna = seq_set_dna.add(feature.dna_sequence)

            if len(feature.protein_translation) > 0:
                self.seq_loader_protein.load(feature.protein_translation)
                node_re_seq_protein = seq_set_protein.add(feature.protein_translation)

            if node_re_seq_dna is not None and node_re_seq_protein is not None:
                seq_translation = str(Seq(feature.dna_sequence).translate())[:-1]
                is_translation = seq_translation == feature.protein_translation
                if is_translation:
                    graph.add_transform_edge2(node_re_seq_dna, node_re_seq_protein,
                                              "re_seq_dna_has_translation")
                else:
                    graph.add_transform_edge2(node_re_seq_dna, node_re_seq_protein,
                                              "re_seq_dna_has_partial_translation")
            if feature.aliases:
                for alias in feature.aliases:
                    if type(alias) == list and len(alias) == 2:
                        db, value = alias
                        db = db.strip()
                        value = value.strip()
                        if value:
                            if db == 'old_locus_tag':
                                node_alias = graph.add_transform_node2(Node(value, 're_gene_old_locus_tag'))
                                graph.add_transform_edge2(node_re_seq_dna, node_alias,
                                                          "re_seq_dna_has_gene_old_locus_tag")
                            elif db == 'locus_tag':
                                node_alias = graph.add_transform_node2(Node(value, 're_gene_locus_tag'))
                                graph.add_transform_edge2(node_re_seq_dna, node_alias,
                                                          "re_seq_dna_has_gene_locus_tag")
                            elif db == 'EC_number':
                                node_alias = graph.add_transform_node2(Node(value, 're_reaction_ec'))
                                graph.add_transform_edge2(node_re_seq_protein, node_alias,
                                                          "re_seq_protein_has_re_reaction_ec")
                            elif db == 'gene':
                                node_alias = graph.add_transform_node2(Node(value, 're_gene_id'))
                                graph.add_transform_edge2(node_re_seq_dna, node_alias,
                                                          "re_seq_dna_has_gene_id")
                            elif db == 'gene_synonym':
                                node_alias = graph.add_transform_node2(Node(value, 're_gene_id'))
                                graph.add_transform_edge2(node_re_seq_dna, node_alias,
                                                          "re_seq_dna_has_gene_id", data={'comment': 'gene_synonym'})
                            elif db == 'alias':
                                if value.startswith('WP_'):
                                    ncbi_protein_mapping[(node_re_seq_protein.id, feature.id)] = value
                            elif db == 'protein_id':
                                ncbi_protein_mapping[(node_re_seq_protein.id, feature.id)] = value
                            else:
                                logger.warning(f'alias: what is this [{db}, {value}]?')
                    else:
                        logger.warning(f'alias: what is this [{alias}]?')

            if node_re_seq_dna is not None:
                match_blocks = []
                dna_placement = []
                node_contig = None
                for loc in feature.location:
                    contig_id, a, s, b = loc
                    if contig_id in seq_contigs.features and 're_contig' in graph.t_nodes and \
                            f're_contig/{contig_id_to_h[contig_id]}' in graph.t_nodes['re_contig']:
                        node_contig = graph.t_nodes['re_contig'][f're_contig/{contig_id_to_h[contig_id]}']

                        start = a
                        end = b
                        if s == '+':
                            end = a + b - 1
                        elif s == '-':
                            start = a - b + 1
                            end = start + b - 1

                        contig_feature = seq_contigs.features.get_by_id(contig_id)
                        if s == '+':
                            match_blocks.append(contig_feature.seq[(start - 1):end])
                        elif s == '-':
                            match_blocks.append(str(Seq(contig_feature.seq[start - 1:end][::-1]).complement()))

                        dna_placement.append([start, end])
                    else:
                        logger.warning(f"Feature {feature.id} could not find contig {loc} ")

                match = ''.join(match_blocks) == feature.dna_sequence
                data = {"location": dna_placement, "match": match}
                graph.add_transform_edge2(node_re_seq_dna, node_contig, "re_contig_has_re_dna_seq", data=data)

        node_seq_set_protein = seq_set_protein.build_set_node()
        node_seq_set_dna = seq_set_dna.build_set_node()

        return node_seq_set_protein, node_seq_set_dna, ncbi_protein_mapping
