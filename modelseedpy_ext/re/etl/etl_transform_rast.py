import logging
from modelseedpy_ext.re.etl.etl_transform_graph import ETLTransformGraph
from modelseedpy.core.msgenome import normalize_role

logger = logging.getLogger(__name__)


class ETLRast(ETLTransformGraph):

    def __init__(self, rast, dna_store, protein_store):
        super().__init__()
        self.dna_store = dna_store
        self.protein_store = protein_store
        self.rast = rast

    def transform(self, rast_annotations):
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

            logger.error('add_node error')
            return None

        def add_edge(node_from, node_to, label, data=None):
            if label not in edges:
                edges[label] = []
            if label in edges:
                _edge = self.transform_edge(node_from, node_to, data)
                edges[label].append(_edge)
                return _edge

            logger.error('add_edge error')
            return None

        for h in rast_annotations:
            node_protein_seq = add_node(h, 're_seq_protein', data=None)
            _value = rast_annotations[h]
            node_annotation = add_node(normalize_role(_value), 'rast', data={'value': _value})
            add_edge(node_protein_seq, node_annotation, 're_seq_protein_has_rast_annotation')
        return nodes, edges
