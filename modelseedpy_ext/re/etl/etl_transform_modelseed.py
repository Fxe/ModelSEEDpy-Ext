import logging
from modelseedpy_ext.re.etl.etl_transform_graph import ETLTransformGraph
from modelseedpy_ext.re.utils import get_value_set_hash
import hashlib

logger = logging.getLogger(__name__)


def get_hash(seq):
    h = hashlib.sha256(seq.encode("utf-8")).hexdigest()
    return h


class ETLTransformModelSEEDCompound(ETLTransformGraph):
    def __init__(self, biochem_service):
        self.biochem_service = biochem_service
        super().__init__()

    def transform(self, compound_data):
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

        node_compound_key = compound_data["id"]
        node_compound_data = dict(compound_data)
        node_compound = add_node(
            node_compound_key, "modelseed_compound", data=node_compound_data
        )

        if "smiles" in compound_data and len(compound_data["smiles"]) > 0:
            try:
                biochem_compound = self.biochem_service.get_molecule(
                    compound_data["smiles"]
                )
                # print(biochem_compound)
                biochem_smiles = biochem_compound["smiles"]
                biochem_inchi = biochem_compound["inchi"]
                biochem_inchi_key = biochem_compound["inchiKey"]
                compound_smiles = compound_data["smiles"]

                node_cdk_smiles = add_node(
                    get_hash(biochem_smiles),
                    "re_biochem_smiles",
                    data={"cdk_inchi_smiles": True, "value": biochem_smiles},
                )
                if compound_smiles != biochem_smiles:
                    node_smiles = add_node(
                        get_hash(compound_smiles),
                        "re_biochem_smiles",
                        data={"value": compound_smiles},
                    )
                    add_edge(
                        node_smiles,
                        node_cdk_smiles,
                        "re_biochem_smiles_equivalent_smiles",
                    )
                    add_edge(
                        node_compound, node_smiles, "modelseed_compound_has_smiles"
                    )
                else:
                    add_edge(
                        node_compound, node_cdk_smiles, "modelseed_compound_has_smiles"
                    )

                node_inchi = add_node(
                    get_hash(biochem_inchi),
                    "re_biochem_inchi",
                    data={"value": biochem_inchi},
                )
                node_inchikey = add_node(
                    biochem_inchi_key, "re_biochem_inchikey", data={}
                )
                add_edge(node_cdk_smiles, node_inchi, "re_biochem_smiles_has_inchi")
                add_edge(node_inchi, node_inchikey, "re_biochem_inchi_has_inchi_key")
                add_edge(
                    node_compound, node_inchikey, "modelseed_compound_has_inchi_key"
                )
                add_edge(node_compound, node_inchi, "modelseed_compound_has_inchi")
            except ValueError as e:
                e
        if (
            "linked_compound" in compound_data
            and compound_data["linked_compound"] != None
        ):
            print(compound_data["linked_compound"])
            pass

        return nodes, edges


class ETLTransformModelSEEDReaction(ETLTransformGraph):
    def __init__(self, dna_store, protein_store):
        super().__init__()
