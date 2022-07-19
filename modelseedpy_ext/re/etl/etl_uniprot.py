import logging
import networkx as nx
from modelseedpy_ext.re.etl.etl_transform_graph import ETLTransformGraph

logger = logging.getLogger(__name__)


class ETLTransformUniprot(ETLTransformGraph):
    
    def __init__(self, seq_store_protein, uniprot_type='sprot'):
        super().__init__()
        self.seq_store_protein = seq_store_protein
        self.uniprot_collection = 'uniprotkb_sprot'
        self.uniprot_accession_collection = 'uniprotkb_accession'
        self.uniprot_subcell = 'uniprotkb_subcell'
        self.eco_term = 'eco_term'
        self.rhea_collection = 'RHEA_reaction'
        self.chebi_collection = 'ChEBI_term'
        self.ec_collection = 'EC_number'
        self.re_seq_protein_collection = 're_seq_protein'
        self.kegg_gene_collection = 'kegg_gene'
        self.alphafolddb_collection = 'alphafolddb'
        self.uniprot_collection_has_ec = f'uniprotkb_{uniprot_type}_has_ec'
        self.uniprot_collection_has_subcell = f'uniprotkb_{uniprot_type}_has_subcell'
        self.uniprot_collection_has_accession = f'uniprotkb_{uniprot_type}_has_accession'
        self.uniprot_collection_has_protein_sequence = f'uniprotkb_{uniprot_type}_has_protein_sequence'
        self.uniprot_collection_has_reference_to_kegg_gene = f'uniprotkb_{uniprot_type}_has_reference_to_kegg_gene'
        self.uniprot_collection_has_reference_to_alphafolddb = f'uniprotkb_{uniprot_type}_has_reference_to_alphafolddb'
        self.uniprot_collection_has_cofactor_chebi = f'uniprotkb_{uniprot_type}_has_cofactor_chebi_term'
        self.uniprot_collection_has_reaction_rhea = f'uniprotkb_{uniprot_type}_has_catalytic_activity_rhea_reaction'
        self.uniprot_collection_has_reaction_ec = f'uniprotkb_{uniprot_type}_has_catalytic_activity_ec_number'
        
    @staticmethod
    def get_cofactor(o):
        res = []
        for comment in o['comment']:
            if comment['type'] == 'cofactor' and 'cofactor' in comment:
                for cofactor in comment['cofactor']:
                    cofactor_compound = [None, cofactor.get('evidence')]
                    for db_ref in cofactor.get('dbReference', []):
                        if db_ref['type'] == 'ChEBI':
                            cofactor_compound[0] = db_ref['id']
                        else:
                            print('get_cofactor ignore', db_ref['type'])
                    if cofactor_compound[0]:
                        res.append(cofactor_compound)
        return res
        
    @staticmethod
    def get_catalytic_activity(o):
        catalytic_activity = []
        for comment in o['comment']:
            if comment['type'] == 'catalytic activity':
                if len(comment['reaction']) == 1:

                    reaction = comment['reaction'][0]
                    catalytic_activity_reaction = [None, reaction.get('evidence'), None]
                    for db_ref in reaction.get('dbReference', []):
                        if db_ref['type'] == 'Rhea':
                            catalytic_activity_reaction[0] = db_ref['id']
                        elif db_ref['type'] == 'EC':
                            catalytic_activity_reaction[2] = db_ref['id']
                        elif db_ref['type'] == 'ChEBI':
                            pass
                        else:
                            print('get_catalytic_activity ignore', db_ref['type'])
                    if catalytic_activity_reaction[0]:
                        catalytic_activity.append(catalytic_activity_reaction)
                else:
                    print('error', o['accession'])
        return catalytic_activity
        
    @staticmethod
    def get_subcellular_location(o):
        comment_location = []
        for comment in o['comment']:
            if comment['type'] == 'subcellular location':
                sl = comment.get('subcellularLocation')
                if sl:
                    for locations in sl:
                        for location in locations['location']:
                            comment_location.append([location['value'], location.get('evidence'), None])
                else:
                    logger.error(f"get_subcellular_location [{o['xml_sourceline']}]")
                        
        return comment_location
                        
        for ref in o['dbReference']:
            if ref['type'] == 'GO':
                location = [None, None, ref['id']]
                location_property = False
                for prop in ref.get('property', []):
                    if prop.get('type', '') == 'project' and prop.get('value', '') == 'UniProtKB-SubCell':
                        location_property = True
                    elif prop.get('type', '') == 'evidence':
                        location[1] = prop.get('value', '')
                    elif prop.get('type', '') == 'term':
                        location[0] = prop.get('value', '')
                if location_property:
                    comment_location.append(location)
                    
        return comment_location

    def transform(self, o):
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
                    #print('dup', _node.id)
                    pass
                    
                return nodes[label][_node.id]
            print('error')
            return None

        def add_edge(node_from, node_to, label, data=None):
            if label not in edges:
                edges[label] = []
            if label in edges:
                _edge = self.transform_edge(node_from, node_to, data)
                edges[label].append(_edge)
                return _edge
            print('error')
            return None

        node_uniprotkb = add_node('_'.join(sorted(o['accession'])), self.uniprot_collection)
        
        #nodes[self.uniprot_collection].append(node_uniprotkb)
        for i in o['accession']:
            node_accession = add_node(i, self.uniprot_accession_collection)
            add_edge(node_uniprotkb, node_accession, self.uniprot_collection_has_accession)
         
        # TODO: load reference data (publications, etc)
        eco_key = {}
        for evidence in o['evidence']:
            if evidence['key'] not in eco_key:
                eco_key[evidence['key']] = {}
            eco_key[evidence['key']][evidence['type']] = {}
            #print(evidence)
            add_node(evidence['type'], self.eco_term)
        #print(eco_key)
        
        subcellular_location = self.get_subcellular_location(o)
        if len(subcellular_location) > 0:
            for sl in subcellular_location:
                node_subcell = add_node(sl[0], self.uniprot_subcell)
                edge = add_edge(node_uniprotkb, node_subcell, self.uniprot_collection_has_subcell)
                if sl[1] in eco_key:
                    edge['eco'] = eco_key[sl[1]]
                
        catalytic_activity = self.get_catalytic_activity(o)
        if len(catalytic_activity) > 0:
            for p in catalytic_activity:
                node_rhea = add_node(p[0], self.rhea_collection)
                edge = add_edge(node_uniprotkb, node_rhea, self.uniprot_collection_has_reaction_rhea)
                if p[1] in eco_key:
                    edge['eco'] = eco_key[p[1]]
                if p[2]:
                    node_ec = add_node(p[2], self.ec_collection)
                    add_edge(node_uniprotkb, node_ec, self.uniprot_collection_has_reaction_ec)
                
        cofactor = self.get_cofactor(o)
        if len(cofactor) > 0:
            for p in cofactor:
                node_chebi = add_node(p[0], self.chebi_collection)
                edge = add_edge(node_uniprotkb, node_chebi, self.uniprot_collection_has_cofactor_chebi)
                if p[1] in eco_key:
                    edge['eco'] = eco_key[p[1]]
            
        for ref in o['dbReference']:
            if ref['type'] == 'KEGG':
                node_ref = add_node(ref['id'], self.kegg_gene_collection)
                add_edge(node_uniprotkb, node_ref, self.uniprot_collection_has_reference_to_kegg_gene)
            elif ref['type'] == 'AlphaFoldDB':
                node_ref = add_node(ref['id'], self.alphafolddb_collection)
                add_edge(node_uniprotkb, node_ref, self.uniprot_collection_has_reference_to_alphafolddb)
            elif ref['type'] == 'EC':
                node_ref = add_node(ref['id'], self.ec_collection)
                add_edge(node_uniprotkb, node_ref, self.uniprot_collection_has_ec)
                
        sequence = o.get('protein_sequence', {}).get('value', '').strip()
        if len(sequence) > 0:
            h = self.seq_store_protein.store_sequence(sequence)
            node_sequence = add_node(h, self.re_seq_protein_collection, {'size': len(sequence)})
            add_edge(node_uniprotkb, node_sequence, self.uniprot_collection_has_protein_sequence)
            
        for copy_key in {'name', 'gene', 'proteinExistence', 'protein', 'comment'}:
            if copy_key in o:
                node_uniprotkb.data[copy_key] = o[copy_key]
        #print(sprot_node)

        return nodes, edges
