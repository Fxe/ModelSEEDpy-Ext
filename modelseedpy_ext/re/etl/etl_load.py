import logging
import copy
import hashlib
from modelseedpy_ext.re.etl.transform_graph import TransformGraph

logger = logging.getLogger(__name__)

DEBUG = None


class ETLLoadArangoDb:

    def __init__(self, arango_database):
        self.aql_upsert = """
        FOR d IN @docs
            UPSERT { _key: d._key}
            INSERT d
            UPDATE d IN @@col
        """
        self.db = arango_database

    @staticmethod
    def process_key(node):
        k = node['_key']
        if len(k) > 120:
            h = hashlib.sha256(k.encode("utf-8")).hexdigest()
            node["_key"] = h
            node["_key_long"] = k

    @staticmethod
    def resolve_edge_key(node_from, node_to, node_ids):
        k_from = None
        k_to = None
        if node_from.label in node_ids:
            k_from = node_ids[node_from.label][node_from.key]
        if node_to.label in node_ids:
            k_to = node_ids[node_to.label][node_to.key]
        if k_from and k_to:
            return f"{k_from}:{k_to}"
        else:
            raise ValueError("unable to resolve ids")

    def load_edges(self, edges, node_ids):
        for collection_id in edges:
            if not self.db.hasCollection(collection_id):
                logger.warning(f"create {collection_id}")
                self.db.createCollection(name=collection_id, className="Edges")

            keys = set()
            to_load = []
            for edge in edges[collection_id]:
                node_from = edge["_from"]
                node_to = edge["_to"]
                edge_copy = copy.deepcopy(edge)
                edge_copy["_key"] = self.resolve_edge_key(node_from, node_to, node_ids)
                edge_copy[
                    "_from"
                ] = f"{node_from.label}/{node_ids[node_from.label][node_from.key]}"
                edge_copy[
                    "_to"
                ] = f"{node_to.label}/{node_ids[node_to.label][node_to.key]}"
                if edge_copy["_key"] not in keys:
                    keys.add(edge_copy["_key"])
                    to_load.append(edge_copy)
            if len(to_load) > 0:
                # print(to_load)
                bin_vars = {"docs": to_load, "@col": collection_id}
                self.db.AQLQuery(self.aql_upsert, bindVars=bin_vars)

    @staticmethod
    def load_old(data):
        nodes = data["nodes"]
        edges = data["edges"]
        # node_ids = self.load_nodes(nodes)
        # self.load_edges(edges, node_ids)
        raise Exception("removed")

    def load_collection_nodes(self, G, collection_id):
        inserted_nodes = {}
        if not self.db.hasCollection(collection_id):
            logger.warning(f"create {collection_id}")
            self.db.createCollection(name=collection_id)
        keys = set()
        to_load = []
        for node_id, node in G.t_nodes[collection_id].items():
            node_copy = copy.deepcopy(node.data)
            c, key = node.id.split('/', 1)
            node_copy['_key'] = key
            self.process_key(node_copy)
            # print(node_copy, c, key)
            # print(node_id)
            inserted_nodes[key] = node_copy["_key"]
            if node.key not in keys:
                keys.add(node.key)
                to_load.append(node_copy)
            if len(to_load) > 0:
                bin_vars = {"docs": to_load, "@col": collection_id}
                try:
                    self.db.AQLQuery(self.aql_upsert, bindVars=bin_vars)
                    pass
                except Exception as e:
                    global DEBUG
                    DEBUG = bin_vars
                    raise e
        return inserted_nodes

    def load_collection_edges(self, G: TransformGraph, collection_id, node_ids):
        if not self.db.hasCollection(collection_id):
            logger.warning(f"create {collection_id}")
            self.db.createCollection(name=collection_id, className="Edges")

        keys = set()
        to_load = []
        for src, dst in G.t_edges[collection_id]:
            data = G.get_edge_data(src, dst)['data']
            data = data if data else {}
            node_from = src
            node_to = dst
            edge_copy = copy.deepcopy(data)
            edge_copy["_key"] = self.resolve_edge_key(node_from, node_to, node_ids)
            edge_copy[
                "_from"
            ] = f"{node_from.label}/{node_ids[node_from.label][node_from.key]}"
            edge_copy[
                "_to"
            ] = f"{node_to.label}/{node_ids[node_to.label][node_to.key]}"
            if edge_copy["_key"] not in keys:
                keys.add(edge_copy["_key"])
                to_load.append(edge_copy)
        if len(to_load) > 0:
            # print(to_load)
            bin_vars = {"docs": to_load, "@col": collection_id}
            self.db.AQLQuery(self.aql_upsert, bindVars=bin_vars)

    def load(self, graph: TransformGraph):
        inserted_nodes = {}
        for col in graph.t_nodes:
            inserted_nodes[col] = self.load_collection_nodes(graph, col)
        for col in graph.t_edges:
            self.load_collection_edges(graph, col, inserted_nodes)
