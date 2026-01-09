import os
from modelseedpy import MSGenome
from modelseedpy_ext.re.hash_seq import HashSeq
from modelseedpy_ext.re.core.genome import ProteinSequence
from modelseedpy_ext.ani.skani import read_search_output_as_parquet
from modelseedpy_ext.re.cdm.query_pangenome import QueryPangenomes
from modelseedpy_ext.re.cdm.core import ClusterSet
from collections import defaultdict


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
        elif db_xref.startswith('UniRef:UniRef100'):
            ontology.append(['uniref_100', db_xref])
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
            ontology.append(['others', db_xref])
    return ontology


class ReportFactory:

    def __init__(self, genome: MSGenome, assembly: MSGenome, pg: QueryPangenomes):
        self.genome = genome
        self.assembly = assembly
        self.genome_feature_h = ReportFactory.feature_to_h(self.genome)
        self.genome_feature_bakta = {}
        self.pg = pg
        if not os.path.exists('./query_genome.faa'):
            pass
        if not os.path.exists('./query_assembly.fna'):
            pass
        self.data_pan_genome = None
        self.c_to_m = None  # Cluster To Members
        self.m_to_c = None  # Member To Cluster

    def write_genome_files(self):
        if not os.path.exists('./query_genome.faa'):
            self.genome.to_fasta('./query_genome.faa')
        if not os.path.exists('./query_assembly.fna'):
            self.assembly.to_fasta('./query_assembly.fna')

    @staticmethod
    def run_ani(query, library, output_file, threads=20):
        import subprocess
        threads = threads
        cmd = [
            '/opt/skani/0.3.1/skani', 'search', "-t", str(threads),
            "-q", query,
            "-d", library,
            "-o", output_file,
        ]
        output = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )
        return output

    @staticmethod
    def run_mmseqs2(threads=40):
        import subprocess
        cmd = [
            '/opt/mmseqs2/15-6f452/bin/mmseqs', 'easy-cluster', "--threads", str(threads),
            '--min-seq-id', '0.0',
            '--cov-mode', '0', '-c', '0.80',
            "./all_proteins.faa", 'mmseqs2', 'mmseqs2'
        ]
        output = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )
        return output

    def get_query_genome_to_pan_genome_clusters(self):
        if not os.path.exists('./mmseqs2_cluster.tsv'):
            program_out = ReportFactory.run_mmseqs2()

        self.c_to_m = defaultdict(set)
        self.m_to_c = {}

        with open("./mmseqs2_cluster.tsv", 'r') as fh:
            for line in fh:
                rep, mem = line.strip().split("\t")
                self.c_to_m[rep].add(mem)
                self.m_to_c[mem] = rep

        pan_c_to_m = {}
        pan_m_to_c = {}
        for cluster_id, members in self.data_pan_genome['cluster_to_genes'].items():
            pan_c_to_m[cluster_id] = [self.data_pan_genome['feature_id_to_protein_h'][m] for m in members]
            for m in pan_c_to_m[cluster_id]:
                pan_m_to_c[m] = cluster_id

        pan_c_to_m_mapped_to_representative = {}
        for cluster_id, members in pan_c_to_m.items():
            pan_c_to_m_mapped_to_representative[cluster_id] = [self.m_to_c[m] for m in members]

        representative_to_pan_cluster_id = {}
        for cluster_id, reps in pan_c_to_m_mapped_to_representative.items():
            unique_reps = set(reps)
            for rep in unique_reps:
                if rep not in representative_to_pan_cluster_id:
                    representative_to_pan_cluster_id[rep] = set()
                representative_to_pan_cluster_id[rep].add(cluster_id)

        return pan_c_to_m, pan_m_to_c, pan_c_to_m_mapped_to_representative, representative_to_pan_cluster_id

    def get_ani(self):
        self.write_genome_files()
        from modelseedpy_ext.ani.skani import _read_search_output_as_parquet

        if not os.path.exists('./ani_kepangenomes_fast.out'):
            program_out = ReportFactory.run_ani('./query_assembly.fna',
                                  '/home/fliu/scratch/data/ani/skani/ke-pangenomes_fast',
                                  './ani_kepangenomes_fast.out')
        df_ani_clade = _read_search_output_as_parquet('./ani_kepangenomes_fast.out').to_pandas()

        def t_ncbi_to_gtdb_id(s):
            a = s.split('_')
            if s.startswith('GCA'):
                return f'GB_{a[0]}_{a[1]}'
            elif s.startswith('GCF'):
                return f'RS_{a[0]}_{a[1]}'

        def q_transform(s):
            return 'user_genome'

        def r_transform(s):
            return t_ncbi_to_gtdb_id(s.split('/')[-1].rsplit('_', 1)[0])

        def ani_transform(df, fn_q_transform, fn_r_transform):
            query_to_ref = {}
            for row_id, d in df.iterrows():
                q = fn_q_transform(d['Query_file'])
                r = fn_r_transform(d['Ref_file'])
                if q not in query_to_ref:
                    query_to_ref[q] = {}
                query_to_ref[q][r] = [d['ANI'], d['Align_fraction_ref'], d['Align_fraction_query']]
            return query_to_ref

        if not os.path.exists('./ani_phenotypes_fast.out'):
            program_out = ReportFactory.run_ani('./query_assembly.fna',
                                  '/home/fliu/scratch/data/ani/skani/cdm_fast_phenotypes',
                                  './ani_phenotypes_fast.out')

        _pq = _read_search_output_as_parquet('./ani_phenotypes_fast.out')

        ani_phenotypes = None

        if _pq:
            df_ani_phenotypes = _pq.to_pandas()

            import pandas as pd
            df_anl_ecoli = pd.read_csv('/home/fliu/CDM/ecoli/data/genomes/metadata.tsv', sep='\t')
            df_pmi = pd.read_csv('/home/fliu/CDM/pmi/data/genomes/metadata_pmi.tsv', sep='\t')
            df_leaf = pd.read_csv('/home/fliu/CDM/pmi/data/genomes/metadata_atleaf.tsv', sep='\t')
            contig_h_to_genome_id = {}
            for row_id, row in df_anl_ecoli.iterrows():
                if row['h'] not in contig_h_to_genome_id:
                    contig_h_to_genome_id[row['h']] = row['genome_id']
                else:
                    raise ValueError('dup')
            for row_id, row in df_pmi.iterrows():
                h = row['hash_contigset']
                if h not in contig_h_to_genome_id:
                    contig_h_to_genome_id[h] = row['genome_id']
                else:
                    raise ValueError('dup')
            for row_id, row in df_leaf.iterrows():
                h = row['hash_contigset']
                if h not in contig_h_to_genome_id:
                    contig_h_to_genome_id[h] = row['genome_id']
                else:
                    raise ValueError('dup')

            def r_transform_phenotypes(s):
                return contig_h_to_genome_id[s.split('/')[-1][:-4]]

            ani_phenotypes = ani_transform(df_ani_phenotypes, q_transform, r_transform_phenotypes)

        return ani_transform(
            df_ani_clade,
            q_transform,
            r_transform
        ), ani_phenotypes

    @staticmethod
    def feature_to_h(g):
        return {f.id: ProteinSequence(f.seq).hash_value for f in g.features if f.seq}

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

    def cluster_query_to_pan_genome(self, cluster: ClusterSet):
        pass

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

    @staticmethod
    def get_query_to_ref(f: str):
        df_ani = read_search_output_as_parquet(f).to_pandas()
        query_to_ref = {}
        for row_id, d in df_ani.iterrows():
            grow_handle = d['Query_file'].split('/')[-1]
            gtdb_ncbi = d['Ref_file'].split('/')[-1].rsplit('_', 1)[0]
            if grow_handle not in query_to_ref:
                query_to_ref[grow_handle] = {}
            query_to_ref[grow_handle][gtdb_ncbi] = [d['ANI'], d['Align_fraction_ref'], d['Align_fraction_query']]
        pass

    def build_pan_table_from_clade_id(self, clade_id, col_bakta):
        clade_members = self.pg.get_clade_members(clade_id)

        pan_member_features = {}
        pan_member_assembly = {}
        for row in clade_members.iter_rows(named=True):
            pan_member_features[row['genome_id']] = self.pg.get_member_genome(row['genome_id'])
            pan_member_assembly[row['genome_id']] = self.pg.get_member_assembly(row['genome_id'])

        clade_gene_clusters = self.pg.get_clade_gene_clusters(clade_id)

        clade_cluster_ids = set(clade_gene_clusters['gene_cluster_id'])

        df_gene_genecluster = self.pg.get_clusters_members(clade_cluster_ids)

        d_gene_to_cluster = {o[0]: o[1] for o in df_gene_genecluster.iter_rows()}
        d_cluster_to_genes = {}
        for k, v in d_gene_to_cluster.items():
            if v not in d_cluster_to_genes:
                d_cluster_to_genes[v] = set()
            d_cluster_to_genes[v].add(k)

        feature_id_to_protein_h = {}
        u_proteins = {}
        for g in pan_member_features.values():
            for f in g.features:
                if f.seq:
                    h_seq = HashSeq(f.seq)
                    h = h_seq.hash_value
                    feature_id_to_protein_h[f.id] = h
                    u_proteins[h] = str(h_seq)

        cursor = col_bakta.find({'_id': {'$in': list(u_proteins)}})
        u_proteins_bakta = {_doc['_id']: _doc for _doc in cursor}

        return {
            'ani_clade_members': None,
            'pan_member_metadata': clade_members,
            'pan_member_features': pan_member_features,
            'pan_member_assembly': pan_member_assembly,
            'pan_gene_clusters': clade_gene_clusters,
            'cluster_to_genes': d_cluster_to_genes,  # cluster to member ids (feature id)
            'u_proteins': u_proteins,  # NR Protein collection of pan-genome
            'feature_id_to_protein_h': feature_id_to_protein_h,  # Feature collection to Protein Hash
            'u_proteins_bakta': u_proteins_bakta,  # NR Protein Bakta annotation
        }

    def build_pan_table(self, fn_ani_ref, fn_ani_clade, col_bakta):
        """
        Table Pangenome:
        cluster_id n_members core

        Table Pangenome Member Features
        feature_id
        :return:
        """

        ref_ani = fn_ani_ref(self.genome)

        member_id_top_hit, ani_values = max(ref_ani.items(), key=lambda x: x[1][0])
        clade_id = self.pg.get_member_representative(member_id_top_hit)
        res = self.build_pan_table_from_clade_id(clade_id, col_bakta)
        res['ani_ref_top_hit'] = [member_id_top_hit, ani_values]
        return res

    def wut(self, pan_data):
        pan_c_to_m = {}
        pan_m_to_c = {}
        for cluster_id, members in pan_data['cluster_to_genes'].items():
            pan_c_to_m[cluster_id] = [pan_data['feature_id_to_protein_h'][m] for m in members]
            for m in pan_c_to_m[cluster_id]:
                pan_m_to_c[m] = cluster_id

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