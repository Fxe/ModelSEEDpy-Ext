import networkx as nx

class ETLTransformUniprot:
    
    def __init__(self, seq_store_protein):
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
        self.uniprot_collection_has_ec = 'uniprotkb_sprot_has_ec'
        self.uniprot_collection_has_subcell = 'uniprotkb_sprot_has_subcell'
        self.uniprot_collection_has_accession = 'uniprotkb_sprot_has_accession'
        self.uniprot_collection_has_protein_sequence = 'uniprotkb_sprot_has_protein_sequence'
        self.uniprot_collection_has_reference_to_kegg_gene = 'uniprotkb_sprot_has_reference_to_kegg_gene'
        self.uniprot_collection_has_reference_to_alphafolddb = 'uniprotkb_sprot_has_reference_to_alphafolddb'
        self.uniprot_collection_has_cofactor_chebi = 'uniprotkb_sprot_has_cofactor_chebi_term'
        self.uniprot_collection_has_reaction_rhea = 'uniprotkb_sprot_has_catalytic_activity_rhea_reaction'
        self.uniprot_collection_has_reaction_ec = 'uniprotkb_has_catalytic_activity_ec_number'
        
    @staticmethod
    def get_cofactor(o):
        res = []
        for comment in o['comment']:
            if comment['type'] == 'cofactor':
                if len(comment['cofactor']) == 1:

                    cofactor = comment['cofactor'][0]
                    cofactor_compound = [None, cofactor.get('evidence')]
                    for db_ref in cofactor.get('dbReference', []):
                        if db_ref['type'] == 'ChEBI':
                            cofactor_compound[0] = db_ref['id']
                        else:
                            print('get_cofactor ignore', db_ref['type'])
                    if cofactor_compound[0]:
                        res.append(cofactor_compound)
                else:
                    print('error', o['accession'])
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
                for locations in sl:
                    for location in locations['location']:
                        comment_location.append([location['value'], location.get('evidence'), None])
                        
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
    
    @staticmethod
    def transform_edge(src, dst, data=None):
        return {
            '_key': '{}:{}'.format(src['_key'], dst['_key']),
            '_from': src,
            '_to': dst,
        }
    
    def transform(self, o):
        nodes = {
            self.uniprot_collection: [],
            self.re_seq_protein_collection: [],
            self.uniprot_accession_collection: [],
            self.uniprot_subcell: [],
            self.eco_term: [],
            self.chebi_collection: [],
            self.rhea_collection: [],
            self.ec_collection: [],
            self.kegg_gene_collection: [],
            self.alphafolddb_collection: []
        }
        edges = {
            self.uniprot_collection_has_subcell: [],
            self.uniprot_collection_has_accession: [],
            self.uniprot_collection_has_protein_sequence: [],
            self.uniprot_collection_has_reference_to_kegg_gene: [],
            self.uniprot_collection_has_reference_to_alphafolddb: [],
            self.uniprot_collection_has_cofactor_chebi: [],
            self.uniprot_collection_has_reaction_rhea: [],
            self.uniprot_collection_has_ec: [],
            self.uniprot_collection_has_reaction_ec: [],
        }
        sprot_node = {
            '_key': '_'.join(sorted(o['accession']))
        }
        nodes[self.uniprot_collection].append(sprot_node)
        for i in o['accession']:
            node_accession = {'_key': i}
            nodes[self.uniprot_accession_collection].append(node_accession)
            edges[self.uniprot_collection_has_accession].append(
                self.transform_edge(sprot_node, node_accession))
         
        # TODO: load reference data (publications, etc)
        eco_key = {}
        for evidence in o['evidence']:
            if evidence['key'] not in eco_key:
                eco_key[evidence['key']] = {}
            eco_key[evidence['key']][evidence['type']] = {}
            #print(evidence)
            nodes[self.eco_term].append({'_key': evidence['type']})
        #print(eco_key)
        
        subcellular_location = self.get_subcellular_location(o)
        if len(subcellular_location) > 0:
            for sl in subcellular_location:
                node_subcell = {'_key': sl[0]}
                nodes[self.uniprot_subcell].append(node_subcell)
                edge = self.transform_edge(sprot_node, node_subcell)
                if sl[1] in eco_key:
                    edge['eco'] = eco_key[sl[1]]
                edges[self.uniprot_collection_has_subcell].append(edge)
                
        catalytic_activity = self.get_catalytic_activity(o)
        if len(catalytic_activity) > 0:
            for p in catalytic_activity:
                node_rhea = {'_key': p[0]}
                nodes[self.rhea_collection].append(node_rhea)
                edge = self.transform_edge(sprot_node, node_rhea)
                if p[1] in eco_key:
                    edge['eco'] = eco_key[p[1]]
                edges[self.uniprot_collection_has_reaction_rhea].append(edge)
                if p[2]:
                    node_ec = {'_key': p[2]}
                    edges[self.uniprot_collection_has_reaction_ec].append(self.transform_edge(sprot_node, node_ec))
                
        cofactor = self.get_cofactor(o)
        if len(cofactor) > 0:
            for p in cofactor:
                node_chebi = {'_key': p[0]}
                nodes[self.chebi_collection].append(node_chebi)
                edge = self.transform_edge(sprot_node, node_chebi)
                if p[1] in eco_key:
                    edge['eco'] = eco_key[p[1]]
                edges[self.uniprot_collection_has_cofactor_chebi].append(edge)
                
        
        for comment in o['comment']:
            pass
            
        for ref in o['dbReference']:
            if ref['type'] == 'KEGG':
                node_ref = {'_key': ref['id']}
                nodes[self.kegg_gene_collection].append(node_ref)
                edges[self.uniprot_collection_has_reference_to_kegg_gene].append(
                    self.transform_edge(sprot_node, node_ref))
            elif ref['type'] == 'AlphaFoldDB':
                node_ref = {'_key': ref['id']}
                nodes[self.alphafolddb_collection].append(node_ref)
                edges[self.uniprot_collection_has_reference_to_alphafolddb].append(
                    self.transform_edge(sprot_node, node_ref))
            elif ref['type'] == 'EC':
                node_ref = {'_key': ref['id']}
                nodes[self.ec_collection].append(node_ref)
                edges[self.uniprot_collection_has_ec].append(
                    self.transform_edge(sprot_node, node_ref))
                
        sequence = o.get('protein_sequence', {}).get('value', '').strip()
        if len(sequence) > 0:
            h = self.seq_store_protein.store_sequence(sequence)
            node_seq = {'_key': h, 'size': len(sequence)}
            nodes[self.re_seq_protein_collection].append(node_seq)
            edges[self.uniprot_collection_has_protein_sequence].append(
                self.transform_edge(sprot_node, node_seq))
            
        for copy_key in {'name', 'gene', 'proteinExistence', 'protein'}:
            if copy_key in o:
                sprot_node[copy_key] = o[copy_key]
        #print(sprot_node)

        return nodes, edges
