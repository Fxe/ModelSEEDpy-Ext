import hashlib
import logging
from modelseedpy_ext.re.etl.etl_transform_graph import ETLTransformGraph
from modelseedpy_ext.re.etl.transform_graph import TransformGraph

logger = logging.getLogger(__name__)


class ETLTransformContigSet:
    def __init__(self, dna_store):
        self.dna_store = dna_store
        print('!')

    def transform(self, contig_set):
        G = TransformGraph()

        hash_list = []
        for dna_seq in contig_set:
            symbols = {}
            for o in dna_seq:
                if o not in symbols:
                    symbols[o] = 0
                symbols[o] += 1
            h = self.dna_store.store_sequence(dna_seq)
            hash_list.append(h)
            G.add_transform_node(h, "re_contig", 
                                 data={"size": len(dna_seq), "symbols": symbols})

        hash_list = sorted(hash_list)
        hash_seq = "_".join(hash_list)
        hash_contig_set = hashlib.sha256(hash_seq.encode("utf-8")).hexdigest()

        node_contig_set = G.add_transform_node(hash_contig_set, "re_contig_set")
        for node in G.t_nodes["re_contig"].values():
            G.add_transform_edge(node_contig_set.id, node.id,
                                 "re_contig_set_has_contig")

        return G


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
