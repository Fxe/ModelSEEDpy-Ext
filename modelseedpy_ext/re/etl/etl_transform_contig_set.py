import hashlib
import logging
from modelseedpy_ext.re.hash_seq import HashSeqList, HashSeq
from modelseedpy_ext.re.etl.etl_transform_graph import ETLTransformGraph
from modelseedpy_ext.re.etl.transform_graph import TransformGraph, Node
from modelseedpy_ext.re.utils import count_symbols
from modelseedpy import MSGenome

logger = logging.getLogger(__name__)


class ETLTransformContigSet:

    def __init__(self, store_contigs=None):
        self.store_contigs = store_contigs

    @staticmethod
    def yield_contigs(contigs: list | MSGenome):
        if type(contigs) == MSGenome:
            for f in contigs.features:
                yield f.seq
        else:
            for s in contigs:
                if not type(s) == str:
                    raise TypeError(f'found {type(s)} contigs must be str')
                yield s

    def transform_to_graph(self, contig_set: list | MSGenome, graph):
        hash_list = HashSeqList()
        contig_nodes = {}
        contig_counts = {}
        for dna_seq in self.yield_contigs(contig_set):
            seq = HashSeq(dna_seq)
            if self.store_contigs:
                self.store_contigs.load(seq)
            hash_list.append(seq)
            node_contig = Node(seq.hash_value,
                               "re_contig",
                               data={
                                   "size": len(dna_seq),
                                   "symbols": count_symbols(dna_seq)}
                               )
            node_contig = graph.add_transform_node2(node_contig)

            if node_contig.id not in contig_nodes:
                contig_nodes[node_contig.id] = node_contig
                contig_counts[node_contig.id] = 0
            contig_counts[node_contig.id] += 1

        node_contig_set = Node(hash_list.hash_value, "re_contig_set")
        node_contig_set = graph.add_transform_node2(node_contig_set)

        for i, n in contig_nodes.items():
            count = contig_counts[i]
            data = {}
            if count > 1:
                data['copies'] = count
            graph.add_transform_edge2(node_contig_set, n, "re_contig_set_has_contig", data)

        return contig_set

    def transform(self, contig_set: list):

        graph = TransformGraph()

        self.transform_to_graph(contig_set, graph)

        return graph

    def transform_old(self, contig_set):
        graph = TransformGraph()

        hash_list = []
        for dna_seq in contig_set:
            symbols = {}
            for o in dna_seq:
                if o not in symbols:
                    symbols[o] = 0
                symbols[o] += 1
            h = self.dna_store.store_sequence(dna_seq)
            hash_list.append(h)
            graph.add_transform_node(h, "re_contig",
                                     data={"size": len(dna_seq), "symbols": symbols})

        hash_list = sorted(hash_list)
        hash_seq = "_".join(hash_list)
        hash_contig_set = hashlib.sha256(hash_seq.encode("utf-8")).hexdigest()

        node_contig_set = graph.add_transform_node(hash_contig_set, "re_contig_set")
        for node in graph.t_nodes["re_contig"].values():
            graph.add_transform_edge(node_contig_set.id, node.id,
                                 "re_contig_set_has_contig")

        return graph


class ETLTransformContigSetOld(ETLTransformGraph):
    def __init__(self, dna_store):
        super().__init__()
        self.dna_store = dna_store

    def transform(self, contig_set):
        nodes = {}
        edges = {}

        def add_node(node_id, label, data=None):
            if label not in nodes:
                nodes[label] = {}
            if label in nodes:
                _node = self.build_node(node_id, label, data)
                if _node.id not in nodes[label]:
                    nodes[label][_node.id] = _node
                else:
                    # print('dup', _node.id)
                    pass

                return nodes[label][_node.id]

            logger.error("add_node error")
            return None

        def add_edge(node_from, node_to, label, data=None):
            if label not in edges:
                edges[label] = []
            if label in edges:
                _edge = self.transform_edge(node_from, node_to, data)
                edges[label].append(_edge)
                return _edge

            logger.error("add_edge error")
            return None

        hash_list = []
        for dna_seq in contig_set:
            symbols = {}
            for o in dna_seq:
                if o not in symbols:
                    symbols[o] = 0
                symbols[o] += 1
            h = self.dna_store.store_sequence(dna_seq)
            hash_list.append(h)
            add_node(h, "re_contig", {"size": len(dna_seq), "symbols": symbols})

        hash_list = sorted(hash_list)
        hash_seq = "_".join(hash_list)
        hash_contig_set = hashlib.sha256(hash_seq.encode("utf-8")).hexdigest()

        node_contig_set = add_node(hash_contig_set, "re_contig_set", {})
        for node in nodes["re_contig"].values():
            add_edge(node_contig_set, node, "re_contig_set_has_contig", {})

        return nodes, edges
