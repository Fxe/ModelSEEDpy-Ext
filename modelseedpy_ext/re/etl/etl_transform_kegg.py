from modelseedpy_ext.re.etl.transform_graph import TransformGraph, Node, NodeHash


class ETLTransformKEGG():

    def __init__(self, kegg_database):
        self.kegg_database = kegg_database

    def wut(self, version, kegg_metadata):
        graph = TransformGraph()
        for k in kegg_metadata:
            kegg_id = k[:-4]
            if len(kegg_id) == 6 and kegg_id[0] == 'D':
                v = kegg_metadata[k]
                node_identifier = Node(f'{kegg_id}_{version}', self.kegg_database, data={'metabolite_id': kegg_id})
                if v[3]:
                    node_identifier.data['formula'] = v[3]
                if v[4]:
                    node_identifier.data['mol_weight'] = v[4]
                    node_identifier.data['exact_mass'] = v[5]
                graph.add_transform_node2(node_identifier)
                if v[2]:
                    node_inchi_key = Node(v[2], 're_biochem_inchi_key', data={'inchi': v[1]})
                    graph.add_transform_node2(node_inchi_key)
                    graph.add_transform_edge(node_inchi_key.id, node_identifier.id,
                                             're_biochem_inchi_key_has_metabolite')
                if v[1]:
                    node_smiles = NodeHash(v[0], 're_biochem_smiles')
                    graph.add_transform_node2(node_smiles)
                    graph.add_transform_edge(node_smiles.id, node_identifier.id, 're_biochem_smiles_has_metabolite')

        return graph

    def load(self, graph, re):
        import logging
        logger = logging.getLogger(__name__)

        for collection_id in graph.t_nodes:
            if not re.db.hasCollection(collection_id):
                logger.warning(f"create {collection_id}")
                re.db.createCollection(name=collection_id)
            col = re.db[collection_id]
            payload = []
            for node in graph.t_nodes[collection_id].values():
                payload.append(node.to_json())
            res = col.importBulk(payload)
            print(collection_id, res)

        for collection_id in graph.t_edges:
            if not re.db.hasCollection(collection_id):
                logger.warning(f"create {collection_id}")
                re.db.createCollection(name=collection_id, className="Edges")
            col = re.db[collection_id]
            payload = []
            for src, dst in graph.t_edges[collection_id]:
                data = graph.get_edge_data(src, dst)['data']
                data = data if data else {}
                node_from = src
                node_to = dst
                data["_key"] = f"{node_from.key}:{node_to.key}"
                data["_from"] = node_from.id
                data["_to"] = node_to.id
                payload.append(data)
            res = col.importBulk(payload)
            print(collection_id, res)
