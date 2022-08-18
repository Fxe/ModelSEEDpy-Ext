import logging
from modelseedpy_ext.re.etl.etl_transform_graph import ETLTransformGraph

logger = logging.getLogger(__name__)


class ETLTransformUniparc(ETLTransformGraph):

    def __init__(self, protein_store):
        super().__init__()
        self.protein_store = protein_store

    def transform(self, entry):
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

        node_uniparc = add_node(entry['uniParcId'], 'uniparc', {
            'uniParcCrossReferences': entry['uniParcCrossReferences']
        })

        node_seq = None
        if 'sequence' in entry and 'value' in entry['sequence']:
            seq_hash = self.protein_store.store_sequence(entry['sequence']['value'])
            node_seq = add_node(seq_hash, 're_seq_protein')
        else:
            logger.error('unable to find sequence')

        if node_seq:
            add_edge(node_uniparc, node_seq, 'uniparc_has_protein_sequence')

        return nodes, edges
