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
        self.t_nodes = {}
        self.t_edges = {}

    def concat(self, graph):
        res = compose(self, graph)
        res.t_nodes.update(self.t_nodes)
        for klass in graph.t_nodes:
            if klass not in res.t_nodes:
                res.t_nodes[klass] = {}
                res.t_nodes[klass].update(graph.t_nodes[klass]) 
            else:
                for k, v in graph.t_nodes[klass].items():
                    res.t_nodes[klass][k] = v
        res.t_edges.update(self.t_edges)
        for klass in graph.t_edges:
            if klass not in res.t_edges:
                res.t_edges[klass] = {}
                res.t_edges[klass].update(graph.t_edges[klass])
            else:
                for k, v in graph.t_edges[klass].items():
                    res.t_edges[klass][k] = v
        #res.t_nodes.update(graph.t_nodes)
        return res 

    def add_transform_edge(self, src, dst, label, data=None):
        l1 = list(filter(lambda x: x.id == src, self.nodes))
        l2 = list(filter(lambda x: x.id == dst, self.nodes))
        if len(l1) == 1 and len(l2) == 1:
            if label not in self.t_edges:
                self.t_edges[label] = {}
            self.t_edges[label][(l1[0], l2[0])] = data
            self.add_edge(l1[0], l2[0], data=data if data else {})

    def add_transform_node(self, node_id, label, data=None):
        if label not in self.t_nodes:
            self.t_nodes[label] = {}
        node = Node(node_id, label, data)
        if node.id not in self.t_nodes[label]:
            self.t_nodes[label][node.id] = node
            self.add_node(node)
        else:
            raise Exception('dup')
        return node
