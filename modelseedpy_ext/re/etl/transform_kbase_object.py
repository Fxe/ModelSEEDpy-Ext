import logging
from modelseedpy_ext.re.etl.transform_graph import TransformGraph, Node

logger = logging.getLogger(__name__)


class TransformKBaseObjectInfo:

    def __init__(self):
        pass

    @staticmethod
    def transform_node_type(info):
        return Node(info.type, "kbase_type", data=None)

    @staticmethod
    def transform_node_object(info):
        return Node(str(info).replace("/", "_"),
                    "kbase_object",
                    data={
                        "name": info.id,
                        "user": info.user,
                        "created": info.created_at,
                        "metadata": info.metadata,
                    })

    def transform(self, info):
        graph = TransformGraph()

        node_type = graph.add_transform_node2(self.transform_node_type(info))
        node_info = graph.add_transform_node2(self.transform_node_object(info))

        graph.add_transform_edge2(node_info, node_type, "kbase_object_has_type")

        return graph
