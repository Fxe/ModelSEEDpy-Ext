import logging
from modelseedpy_ext.re.etl.etl_transform_graph import ETLTransformGraph
from modelseedpy.core.msgenome import normalize_role
import hashlib


logger = logging.getLogger(__name__)


class ETLRast(ETLTransformGraph):
    def __init__(self, rast, protein_store):
        super().__init__()
        self.protein_store = protein_store
        self.rast = rast

    @staticmethod
    def sha_hex(s: str):
        return hashlib.sha256(s.encode("utf-8")).hexdigest()

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

        for h in rast_annotations:
            node_protein_seq = add_node(h, "re_seq_protein", data=None)

            _value = rast_annotations[h]
            node_annotation = add_node(
                self.sha_hex(_value), "rast_function", data={"value": _value}
            )

            _nmz_value = normalize_role(_value)
            node_search_function = add_node(
                self.sha_hex(_nmz_value),
                "rast_search_function",
                data={"value": _nmz_value},
            )

            add_edge(
                node_protein_seq, node_annotation, "re_seq_protein_has_rast_function"
            )
            add_edge(
                node_annotation,
                node_search_function,
                "rast_function_has_search_function",
            )

        return nodes, edges
