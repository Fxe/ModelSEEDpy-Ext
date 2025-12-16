from modelseedpy_ext.re.hash_seq import HashSeq


def get_ko2(d):
    kos = set()
    for k, v in d:
        if k == 'ko':
            kos.add(v[5:])
    return kos


def get_ko(d):
    kos = set()
    for k, v in d['ontology']:
        if k == 'ko':
            kos.add(v[5:])
    return kos


def collect_ontology(_doc):
    ontology = []
    ontology.append(['bakta_product', _doc['product']])
    for db_xref in _doc.get('db_xrefs', []):
        if db_xref.startswith('UniRef:UniRef50'):
            ontology.append(['uniref_50', db_xref])
        elif db_xref.startswith('UniRef:UniRef90'):
            ontology.append(['uniref_90', db_xref])
        elif db_xref.startswith('SO:'):
            ontology.append(['so', db_xref])
        elif db_xref.startswith('EC:'):
            ontology.append(['ec', db_xref])
        elif db_xref.startswith('KEGG:'):
            ontology.append(['ko', db_xref])
        elif db_xref.startswith('GO:'):
            ontology.append(['go', db_xref])
        elif db_xref.startswith('COG:'):
            ontology.append(['cog', db_xref])
        elif db_xref.startswith('PFAM:'):
            ontology.append(['pfam', db_xref])
        elif db_xref.startswith('UniRef:UniRef100:'):
            ontology.append(['uniref_100', db_xref])
        else:
            #print(db_xref)
            pass
    return ontology


class ReportFactory:

    def __init__(self, genome):
        self.genome = genome
        self.genome_feature_h = ReportFactory.feature_to_h(self.genome)
        self.genome_feature_bakta = {}

    @staticmethod
    def feature_to_h(g):
        return {f.id: HashSeq(f.seq).hash_value for f in g.features if f.seq}

    def get_annotation_bakta(self, col_bakta):
        cursor = col_bakta.find({'_id': {
            '$in': list(set(self.genome_feature_h.values()))}})
        h_to_bakta = {doc['_id']: doc for doc in cursor}

        for f in self.genome.features:
            annotation = h_to_bakta.get(self.genome_feature_h[f.id], None)
            if annotation:
                self.genome_feature_bakta[f.id] = annotation
            else:
                print('!')

    def build_gene_table(self):
        """
        Table gene
        genome_id, feature_id, len, start, end, strand, uniref_50, uniref_90, ke_pangenome_cluster_id,
        ke_pangenome_cluster_core
        :return:
        """
        features = {}
        for f in self.genome.features:
            ontology = collect_ontology(self.genome_feature_bakta[f.id])
            feature_data = {
                'seq_aa': f.seq,
                'seq_nc': f.dna_sequence,
                'location': f.location,
                'ontology': ontology
            }
            features[f.id] = feature_data

    def build_pangenome_table(self):
        """
        Table Pangenome:
        cluster_id n_members core

        Table Pangenome Member Features
        feature_id
        :return:
        """
        # ani ref genomes
        # ani member genomes
        # pangenome clusters
        # pangenome cluster annotation
        pass

    def build_phenotype_kegg(self, modules, ko_to_feature, ko_to_clusters, ko_to_hybrid_clusters):
        filter_pathways = set(modules.pathway_names.keys())
        kegg_data = {
            'pathway': {}
        }
        for pwy_id in filter_pathways:
            kegg_data['pathway'][pwy_id] = {
                'label': modules.pathway_names[pwy_id],
                'modules': list(modules.pathway_modules[pwy_id])
            }
        module_cov = {}
        for md_id in modules.pathway_module_data:
            md = modules.pathway_module_data[md_id]
            steps = []
            for k in md['orthologs']:
                ko_list = {x: {
                    'genome': list(set(ko_to_feature.get(x, []))),
                    'pangenome_cluster': list(set(ko_to_clusters.get(x, []))),
                    'pangenome_hybrid_cluster': list(set(ko_to_hybrid_clusters.get(x, [])))
                } for x in [x.strip() for x in k[1:-1].split(', ')]}
                steps.append(ko_list)
            module_cov[md_id] = steps
        kegg_data['modules'] = module_cov
        kegg_data['module_label'] = {k: v['name'] for k, v in modules.pathway_module_data.items()}
        pass

    def build_uniref(self):
        uniref_90_to_gene = {}
        uniref_90_to_cluster = {}
        for gene_id, d in report['features'].items():
            for k in [x[1] for x in d['ontology'] if x[0] == 'uniref_90']:
                if k not in uniref_90_to_gene:
                    uniref_90_to_gene[k] = set()
                uniref_90_to_gene[k].add(gene_id)
        for cluster_id, d in report['pangenome']['clusters'].items():
            member_u = {}
            for member_id, ontology in d['members'].items():
                member_u[member_id] = {x[1] for x in ontology if x[0] == 'uniref_90'}
            u_sets = {frozenset(x) for x in member_u.values()}
            if len(u_sets) == 1:
                for v in list(u_sets)[0]:
                    if v not in uniref_90_to_cluster:
                        uniref_90_to_cluster[v] = set()
                    uniref_90_to_cluster[v].add(cluster_id)
        merge_u = {}
        for u in uniref_90_to_gene:
            merge_u[u] = [uniref_90_to_gene[u], uniref_90_to_cluster.get(u, [])]

        member_cluster = {}
        for u, [gene_ids, cluster_ids] in merge_u.items():
            for gene_id in gene_ids:
                member_cluster[gene_id] = {}
                for cluster_id in cluster_ids:
                    member_cluster[gene_id][cluster_id] = [u, report['pangenome']['clusters'][cluster_id]['is_core']]

        uniref_50_to_gene = {}
        uniref_50_to_cluster = {}
        for gene_id, d in report['features'].items():
            for k in [x[1] for x in d['ontology'] if x[0] == 'uniref_50']:
                if k not in uniref_50_to_gene:
                    uniref_50_to_gene[k] = set()
                uniref_50_to_gene[k].add(gene_id)
        for cluster_id, d in report['pangenome']['clusters'].items():
            member_u = {}
            for member_id, ontology in d['members'].items():
                member_u[member_id] = {x[1] for x in ontology if x[0] == 'uniref_50'}
            u_sets = {frozenset(x) for x in member_u.values()}
            if len(u_sets) == 1:
                for v in list(u_sets)[0]:
                    if v not in uniref_50_to_cluster:
                        uniref_50_to_cluster[v] = set()
                    uniref_50_to_cluster[v].add(cluster_id)
        merge_u = {}
        for u in uniref_50_to_gene:
            merge_u[u] = [uniref_50_to_gene[u], uniref_50_to_cluster.get(u, [])]

        for u, [gene_ids, cluster_ids] in merge_u.items():
            for gene_id in gene_ids:
                if gene_id not in member_cluster:
                    member_cluster[gene_id] = {}
                for cluster_id in cluster_ids:
                    if cluster_id not in member_cluster[gene_id]:
                        member_cluster[gene_id][cluster_id] = [u,
                                                               report['pangenome']['clusters'][cluster_id]['is_core']]

        for gene_id in report['features']:
            report['features'][gene_id]['ke_pangenome'] = member_cluster.get(gene_id, {})

    def build(self):
        report = {
            'genome': {
                'kbase_ref': str(self.genome.info),
                'kbase_object': str(self.genome.info.id),
                'features': {}
            },
            'ani': None,  # grow_to_gtdb_ani['KBH_5935729'],
            'pangenome': {
                'info': {},
                'core': {},
                'aux': {},
            },
            'features': self.build_gene_table(),
            'phenotype': {
                'kegg': {

                },
                'grow_prediction': {
                    'gem': {},
                    'gapmind': {},
                    'ML': {},
                }
            }
        }