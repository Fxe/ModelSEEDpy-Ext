import logging
from modelseedpy_ext.re.etl.etl_transform_graph import ETLTransformGraph
from modelseedpy_ext.re.etl.transform_graph import TransformGraph

logger = logging.getLogger(__name__)


class TransformKBaseObjectInfo:
    def __init__(self):
        pass

    def transform(self, info):
        G = TransformGraph()

        node_type = G.add_transform_node(info.type, "kbase_type", data=None)
        node_info = G.add_transform_node(str(info).replace("/", "_"),
                                         "kbase_object",
                                         {"name": info.id,
                                          "user": info.user,
                                          "created": info.created_at,
                                          "metadata": info.metadata,
                                          })

        G.add_transform_edge(node_info.id, node_type.id,
                             "kbase_object_has_type")

        return G


class TransformKBaseObjectInfoOld(ETLTransformGraph):
    def __init__(self):
        super().__init__()

    def transform(self, info):
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

        node_type = add_node(info.type, "kbase_type", data=None)
        node_object = add_node(
            str(info).replace("/", "_"),
            "kbase_object",
            {
                "name": info.id,
                "user": info.user,
                "created": info.created_at,
                "metadata": info.metadata,
            },
        )

        add_edge(node_object, node_type, "kbase_object_has_type")

        return nodes, edges
