import hashlib
import modelseedpy


class RE:
    
    def __init__(self, db, kbase_api, eutils_api, dna_store, protein_store):
        self.db = db
        self.kbase_api = kbase_api
        self.eutils_api = eutils_api
        self.dna_store = dna_store
        self.protein_store = protein_store
        
    def init_collections(self):
        collections = [
            'ncbi_assembly', 'ncbi_assembly_gcf_acc', 'ncbi_assembly_gca_acc', 
            'ncbi_nuccore', 'ncbi_taxon',
            'uniref_90', 'uniref_100', 'uniref_50',
            'uniprotkb_accession',
            'uniprotkb_sprot', 'uniprotkb_trembl',
            'uniprotkb_subcell',
            'ChEBI_term',
            'RHEA_reaction',
            'GO_terms',
            'eco_term',
            
            'EC_number',
            
            're_contig', 're_contig_set', 're_seq_dna', 're_seq_protein',

            'kegg_gene',
            'alphafolddb',

            'rast_function',
            'rast_search_function',
        ]
        for c in collections:
            if not self.db.hasCollection(c):
                print('create node collection', c)
                self.db.createCollection(name=c)
        
        edges = [
            'uniprotkb_has_uniprotkb_accession',
            'uniprotkb_has_subcell_go_term',
            'uniprotkb_has_cofactor_chebi_term',
            'uniprotkb_has_catalytic_activity_rhea_reaction',
            'uniprotkb_has_catalytic_activity_ec_number',
            
            'ncbi_assembly_gca_to_gcf',
            'ncbi_assembly_has_nuccore',
            'ncbi_assembly_has_gca',
            'ncbi_assembly_has_gcf',
            'uniref_90_has_go_molecular_function',
            'uniref_90_has_go_cellular_component',
            'eco_term_sub_class_of_eco_term',
            
            're_contig_set_has_contig',

            'rast_function_has_sub_function',
            'rast_function_has_search_function',
            're_seq_protein_has_rast_function',
        ]
        for c in edges:
            if not self.db.hasCollection(c):
                print('create edge collection', c)
                self.db.createCollection(name=c, className='Edges')
    
    def create_node(self, key, data, col):
        if type(col) == str:
            col = self.db[col]
        if key not in col:
            doc = col.createDocument()
            doc['_key'] = key
            for k in data:
                doc[k] = data[k]
            doc.save()
        else:
            doc = col[key]
        return doc

    def create_edge(self, node_from, node_to, attrib, col):
        edge_col = col
        if type(col) == str:
            edge_col = self.db[col]
        edge_key = f'{node_from._key}:{node_to._key}'
        if edge_key not in edge_col:
            edge = edge_col.createEdge()
            edge['_from'] = node_from._id
            edge['_to'] = node_to._id
            edge['_key'] = edge_key
            for k in attrib:
                edge[k] = attrib[k]
            edge.save()
        else:
            edge = edge_col[edge_key]
        return edge
    
    def get_contig_from_genome(self, genome, token):
        assembly = kbase.get_from_ws(genome.assembly_ref)
        return get_contig_from_assembly(assembly, token)

    def get_contig_from_assembly(self, assembly, token):
        file_id = assembly.fasta_handle_ref
        file_path = "/home/fliu/KE/data/"
        res = kbase.download_file_from_kbase(token, file_id, file_path)
        return res

    def load_genome(self, genome):
        for f in genome.features:
            seq_protein_store.store_sequence(f.protein_translation)
            self.dna_store.store_sequence(f.dna_sequence)

        handle_ref = get_contig_from_genome(genome, "TQTZSDDJFYE5MX4RCD6BHP3LTMBE4YUY")
        hash_genome = load_genome_from_fna_file(handle_ref, ' ')

        return hash_genome

    def load_mash_scores(self, docs):
        aql_upsert = """
        FOR d IN @docs
            UPSERT { _key: d._key}
            INSERT { _key: d._key, _to: d._to, _from: d._from, score: d.score, _created_at: DATE_NOW(), _updated_at: DATE_NOW()}
            UPDATE { score: d.score, _updated_at: DATE_NOW()} IN @@col
        """
        bin_vars = {
            'docs': docs,
            '@col': 're_contig_set_mash_score'
        }
        self.db.AQLQuery(aql_upsert, bindVars=bin_vars)

    def load_genome_from_fna_file(self, path, split=' '):
        genome = modelseedpy.core.MSGenome.from_fasta(path, split=split)
        hash_list = []
        node_list = []
        for f in genome.features:
            symbols = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
            for o in f.seq:
                if o not in symbols:
                    symbols[o] = 0
                symbols[o] += 1
            h = self.dna_store.store_sequence(f.seq)
            hash_list.append(h)
            node_contig = self.create_node(h, {'size': len(f.seq), 'symbols': symbols}, 're_contig')
            node_list.append(node_contig)

        hash_list = sorted(hash_list)
        hash_seq = '_'.join(hash_list)
        hash_genome = hashlib.sha256(hash_seq.encode('utf-8')).hexdigest()

        node_contig_set = self.create_node(hash_genome, {}, 're_contig_set')
        for node_contig in node_list:
            self.create_edge(node_contig_set, node_contig, {}, 're_contig_set_has_contig')
        return hash_genome
    
    def get_contig_set(self, h):
        doc = None
        if h in self.db['re_contig_set']:
            doc = self.db['re_contig_set'][h]
        if doc is not None:
            aql = """
            LET edges = (
            FOR e IN re_contig_set_has_contig
                FILTER e._from == @contig_set_hash
                RETURN e._to
            )

            FOR doc IN re_contig
                FILTER doc._id IN edges
                RETURN doc
            """
            res = self.db.AQLQuery(aql, bindVars={'contig_set_hash': doc['_id']}, rawResults=True)
            return {
                'set': doc,
                'contigs': [o for o in res],
            }
        return None
