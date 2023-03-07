from networkx import DiGraph, compose


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

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return isinstance(other, Node) and self.id == other.id


class TransformGraph(DiGraph):
    def __init__(self, incoming_graph_data=None):
        super().__init__(incoming_graph_data)

    def concat(self, graph):
        return compose(self, graph)

    def add_transform_edge(self, src, dst, data=None):
        l1 = list(filter(lambda x: x.id == src, self.nodes))
        l2 = list(filter(lambda x: x.id == dst, self.nodes))
        if len(l1) == 1 and len(l2) == 1:
            self.add_edge(l1[0], l2[0], data=data if data else {})

    def add_transform_node(self, node_id, label, data=None):
        self.add_node(Node(node_id, label, data))
