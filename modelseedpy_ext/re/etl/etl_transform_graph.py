class Node:
    def __init__(self, key, label, data=None):
        self._key = key
        self.label = label
        self.data = data if data else {}

    @property
    def key(self):
        return self._key.replace(" ", "_")

    @property
    def id(self):
        return f"{self.label}/{self.key}"


class ETLTransformGraph:
    def __init__(self):
        pass

    @staticmethod
    def build_node(node_id, label, data=None):
        return Node(node_id, label, data)

    @staticmethod
    def transform_edge(src: Node, dst: Node, data=None):
        _edge = {
            "_key": "{}:{}".format(src.key, dst.key),
            "_from": src,
            "_to": dst,
        }
        if data:
            _edge.update(data)
        return _edge
