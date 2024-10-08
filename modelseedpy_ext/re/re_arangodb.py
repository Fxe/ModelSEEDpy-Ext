import logging
import hashlib
import modelseedpy
from modelseedpy_ext.re.core.biochem_modelseed import ReBiochemModelSEEDCompound
from modelseedpy_ext.utils import progress

logger = logging.getLogger(__name__)


class RE:

    def __init__(self, db, kbase_api, eutils_api, contigs_store, dna_store, protein_store,
                 contig_val=None, dna_val=None, protein_val=None):
        self.db = db
        self.kbase_api = kbase_api
        self.eutils_api = eutils_api
        self.dna_store = dna_store
        self.contigs_store = contigs_store
        self.protein_store = protein_store

        self.fn_validate_contigs = contig_val
        self.fn_validate_dna = dna_val
        self.fn_validate_protein = protein_val

    def init_collections(self):
        collections = [
            "ncbi_assembly",
            "ncbi_assembly_gcf_acc",
            "ncbi_assembly_gca_acc",
            "ncbi_nuccore",
            "ncbi_taxon",
            "uniref_90",
            "uniref_100",
            "uniref_50",
            "uniprotkb_accession",
            "uniprotkb_sprot",
            "uniprotkb_trembl",
            "uniprotkb_subcell",
            "ChEBI_term",
            "RHEA_reaction",
            "GO_terms",
            "eco_term",
            "EC_number",
            "re_contig",
            "re_contig_set",
            "re_seq_dna",
            "re_seq_protein",
            "kegg_gene",
            "alphafolddb",
            "rast_function",
            "rast_search_function",
        ]
        for c in collections:
            if not self.db.hasCollection(c):
                print("create node collection", c)
                self.db.createCollection(name=c)

        edges = [
            "uniprotkb_has_uniprotkb_accession",
            "uniprotkb_has_subcell_go_term",
            "uniprotkb_has_cofactor_chebi_term",
            "uniprotkb_has_catalytic_activity_rhea_reaction",
            "uniprotkb_has_catalytic_activity_ec_number",
            "ncbi_assembly_gca_to_gcf",
            "ncbi_assembly_has_nuccore",
            "ncbi_assembly_has_gca",
            "ncbi_assembly_has_gcf",
            "uniref_90_has_go_molecular_function",
            "uniref_90_has_go_cellular_component",
            "eco_term_sub_class_of_eco_term",
            "re_contig_set_has_contig",
            "rast_function_has_sub_function",
            "rast_function_has_search_function",
            "re_seq_protein_has_rast_function",
        ]
        for c in edges:
            if not self.db.hasCollection(c):
                print("create edge collection", c)
                self.db.createCollection(name=c, className="Edges")

    def get_all_ids(self, collection_id):
        return self.get_all(collection_id, '_key')

    def get_all(self, collection_id, field):
        aql = f"""
        FOR d IN {collection_id}
            RETURN d.{field}
        """
        return set(self.db.AQLQuery(aql, batchSize=50000, rawResults=True))

    def aql_filter(self, collection_id, filter_field, filter_values):
        aql = f"""
        FOR d IN {collection_id}
            FILTER d.{filter_field} in @list
            RETURN d
        """
        return list(self.db.AQLQuery(aql, batchSize=50000, rawResults=True, bindVars={'list': list(filter_values)}))

    @staticmethod
    def _collect_one_to_one(cursor, d: dict):
        for a, b in cursor:
            if a not in b:
                d[a] = b
            else:
                logger.warning(f'multiple [{a}] -> {d[a]} {b}')

    @staticmethod
    def _collect_one_to_many(cursor, d: dict):
        for a, b in cursor:
            if a not in b:
                d[a] = set()
            d[a].add(b)

    def generate_ec_hierarchy(self):
        from modelseedpy_ext.re.utils import generate_ec_hierarchy

        col_ec_number = self.db['re_reaction_ec']
        ec_numbers = {doc['_key'] for doc in col_ec_number.fetchAll(rawResults=True)}
        ec_lv1 = set()
        ec_lv2 = set()
        ec_lv3 = set()
        ec_lv4 = set()
        for ec in ec_numbers:
            if ec.count('.') == 3:
                a, b, c, d = ec.split('.')
                if b == '-' and c == '-' and d == '-':
                    ec_lv1.add(ec)
                elif c == '-' and d == '-':
                    ec_lv2.add(ec)
                elif d == '-':
                    ec_lv3.add(ec)
                else:
                    ec_lv4.add(ec)
            else:
                print(ec)
        ec_hier = generate_ec_hierarchy(ec_lv1, ec_lv2, ec_lv3, ec_lv4)
        new_ecs = set()
        for l1 in ec_hier:
            if l1 not in ec_numbers:
                new_ecs.add(l1)
            for l2 in ec_hier[l1]:
                if l2 not in ec_numbers:
                    new_ecs.add(l2)
                for l3 in ec_hier[l1][l2]:
                    if l3 not in ec_numbers:
                        new_ecs.add(l3)
                    for l4 in ec_hier[l1][l2][l3]:
                        if l4 not in ec_numbers:
                            new_ecs.add(l4)
        print('new_ecs', len(new_ecs))

        raise NotImplementedError()

    def get_link(self, edges_id, ids, payload_size=10000, rel_type='one_to_one', rev=False):
        aql = f"""
        FOR d IN {edges_id}
            FILTER d._from IN @list
            RETURN [d._from, d._to]
        """
        if rev:
            aql = f"""
            FOR d IN {edges_id}
                FILTER d._to IN @list
                RETURN [d._to, d._from]
            """

        res = {}

        payload = []
        for i in progress(ids):
            payload.append(i)
            if len(payload) > payload_size:
                cursor = self.db.AQLQuery(aql, rawResults=True, bindVars={'list': payload})
                if rel_type == 'one_to_one':
                    self._collect_one_to_one(cursor, res)
                elif rel_type == 'one_to_many':
                    self._collect_one_to_many(cursor, res)
                else:
                    raise ValueError(f'bad rel_type {rel_type} expected one_to_one, one_to_many')

        if len(payload) > 0:
            cursor = self.db.AQLQuery(aql, rawResults=True, bindVars={'list': payload})
            if rel_type == 'one_to_one':
                self._collect_one_to_one(cursor, res)
            elif rel_type == 'one_to_many':
                self._collect_one_to_many(cursor, res)
            else:
                raise ValueError(f'bad rel_type {rel_type} expected one_to_one, one_to_many')

        return res

    def load_graph(self, graph, proxies=None):
        if proxies is None:
            proxies = set()
        from modelseedpy_ext.re.re_seq_loader import ArangoLoaderNodes, ArangoLoaderEdges
        for collection_id in graph.t_nodes:
            if collection_id not in proxies:
                with ArangoLoaderNodes(self.db, collection_id, block_size=50000) as loader:
                    loader.load(graph)
                    print(loader.load_stats)
        for collection_id in proxies:
            with ArangoLoaderNodes(self.db, collection_id, block_size=50000) as loader:
                if self.db.hasCollection(collection_id):
                    loader.loaded_data.update(self.get_all_ids(collection_id))
                loader.load(graph)
                print(loader.load_stats)
        for collection_id in graph.t_edges:
            with ArangoLoaderEdges(self.db, collection_id, block_size=50000) as loader:
                loader.load(graph)
                print(loader.load_stats)

    def create_node(self, key, data, col):
        if type(col) == str:
            col = self.db[col]
        if key not in col:
            doc = col.createDocument()
            doc["_key"] = key
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
        edge_key = f"{node_from._key}:{node_to._key}"
        if edge_key not in edge_col:
            edge = edge_col.createEdge()
            edge["_from"] = node_from._id
            edge["_to"] = node_to._id
            edge["_key"] = edge_key
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

    def get_modelseed_compound(self, cpd_id, shallow=True):
        if self.db.hasCollection('re_biochem_modelseed_compound'):
            col = self.db['re_biochem_modelseed_compound']
            if cpd_id in col:
                return ReBiochemModelSEEDCompound(col['cpd00001'], self, shallow=shallow)

    def load_genome(self, genome):
        for f in genome.features:
            seq_protein_store.store_sequence(f.protein_translation)
            self.dna_store.store_sequence(f.dna_sequence)

        handle_ref = get_contig_from_genome(genome, "TQTZSDDJFYE5MX4RCD6BHP3LTMBE4YUY")
        hash_genome = load_genome_from_fna_file(handle_ref, " ")

        return hash_genome

    def load_mash_scores(self, docs):
        aql_upsert = """
        FOR d IN @docs
            UPSERT { _key: d._key}
            INSERT { _key: d._key, _to: d._to, _from: d._from, score: d.score, _created_at: DATE_NOW(), _updated_at: DATE_NOW()}
            UPDATE { score: d.score, _updated_at: DATE_NOW()} IN @@col
        """
        bin_vars = {"docs": docs, "@col": "re_contig_set_mash_score"}
        self.db.AQLQuery(aql_upsert, bindVars=bin_vars)

    def load_genome_from_fna_file(self, path, split=" "):
        genome = modelseedpy.core.MSGenome.from_fasta(path, split=split)
        hash_list = []
        node_list = []
        for f in genome.features:
            symbols = {"A": 0, "T": 0, "C": 0, "G": 0}
            for o in f.seq:
                if o not in symbols:
                    symbols[o] = 0
                symbols[o] += 1
            h = self.dna_store.store_sequence(f.seq)
            hash_list.append(h)
            node_contig = self.create_node(
                h, {"size": len(f.seq), "symbols": symbols}, "re_contig"
            )
            node_list.append(node_contig)

        hash_list = sorted(hash_list)
        hash_seq = "_".join(hash_list)
        hash_genome = hashlib.sha256(hash_seq.encode("utf-8")).hexdigest()

        node_contig_set = self.create_node(hash_genome, {}, "re_contig_set")
        for node_contig in node_list:
            self.create_edge(
                node_contig_set, node_contig, {}, "re_contig_set_has_contig"
            )
        return hash_genome

    def get_contig_set(self, h):
        doc = None
        if h in self.db["re_contig_set"]:
            doc = self.db["re_contig_set"][h]
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
            res = self.db.AQLQuery(
                aql, bindVars={"contig_set_hash": doc["_id"]}, rawResults=True
            )
            return {
                "set": doc,
                "contigs": [o for o in res],
            }
        return None

    def query_ec_childs(self, ec: str, min_dist=1, max_dist=4):
        _aql_param_ec = f"EC_number/{ec}"
        aql_ec_childs = """
        FOR v, e, p IN @min_dist..@max_dist OUTBOUND @ec GRAPH 'ec_graph'
                RETURN v._id
        """
        params = {"min_dist": min_dist, "max_dist": max_dist, "ec": _aql_param_ec}
        # strip the collection ID amd return results
        return {
            x[10:]
            for x in self.db.AQLQuery(aql_ec_childs, rawResults=True, bindVars=params)
        }

    @staticmethod
    def _get_uniprot_database_pairs(db: str):
        if db == "sprot":
            return "uniprotkb_sprot_has_ec", "uniprotkb_sprot_has_accession"
        elif db == "trembl":
            return "uniprotkb_trembl_has_ec", "uniprotkb_trembl_has_accession"
        else:
            raise ValueError(f"found {db}, db must be either sprot or trembl")

    def query_uniprot_by_ec(
        self, ec: str, db="sprot", min_dist=1, max_dist=4, limit=None
    ):
        _param_uprot_ec, _param_uprot_acc = self._get_uniprot_database_pairs(db)
        _param_limit = "" if limit is None else f"LIMIT {limit}"

        aql_ec_to_uniprot = f"""
        LET ec_nodes = (
            FOR v, e, p IN @min_dist..@max_dist OUTBOUND @ec GRAPH 'ec_graph'
                RETURN v._id
            )
    
        LET uniprot_nodes = (
        FOR e IN {_param_uprot_ec}
            FILTER e._to IN ec_nodes
            {_param_limit}
            RETURN e._from
        )
    
        FOR e IN {_param_uprot_acc}
            FILTER e._from IN uniprot_nodes
            RETURN e._to
        """
        params = {"min_dist": min_dist, "max_dist": max_dist, "ec": f"EC_number/{ec}"}
        return {
            x
            for x in self.db.AQLQuery(
                aql_ec_to_uniprot, rawResults=True, bindVars=params
            )
        }

    def query_uniprot_by_ec_v2(
        self, ec: str, db="sprot", min_dist=1, max_dist=4, limit=None
    ):
        _param_uprot_ec, _param_uprot_acc = self._get_uniprot_database_pairs(db)
        res_ec_childs = self.query_ec_childs(ec, min_dist, max_dist)

        _param_limit = "" if limit is None else f"LIMIT {limit}"

        aql_ec_to_uniprot = f"""
        LET uniprot_nodes = (
        FOR e IN {_param_uprot_ec}
            FILTER e._to == @ec
            {_param_limit}
            RETURN e._from
        )
        
        FOR e IN {_param_uprot_acc}
            FILTER e._from IN uniprot_nodes
            RETURN e._to
        """

        res = {}
        for ec_child in res_ec_childs:
            params = {"ec": f"EC_number/{ec_child}"}
            res[ec_child] = {
                x
                for x in self.db.AQLQuery(
                    aql_ec_to_uniprot, rawResults=True, bindVars=params
                )
            }
        return res
