from modelseedpy import MSGenome
import json


class AnnotationDatabase:

    def __init__(self):
        self.annotation_database = {}
        self.gene_id_to_rast = None
        self.genome_paths = None

    @staticmethod
    def get_genome_object(filename):
        if 'handle' in filename:
            genome = MSGenome.from_fasta(filename)
            return genome
        else:
            genome = MSGenome.from_fasta(filename)
            return genome

    @staticmethod
    def load_cliff_genomes(kbase):
        genomes = {}
        for ws_id in [155805, 163264]:
            for o in kbase.list_objects(ws_id):
                if o[2].startswith('KBaseGenomes.Genome'):
                    genome = kbase.get_from_ws(o[1], ws_id)
                    genomes[genome.info.id] = genome

        return genomes

    def get_gene_id_to_rast(self, path='/scratch/fliu/data/cliff/annotation.json'):
        print("size should be # 97805067")
        with open(path, 'r') as fh:
            annotation_gtdb = json.load(fh)

            self.gene_id_to_rast = {}
            for rast in annotation_gtdb:
                for g_id in annotation_gtdb[rast]:
                    if g_id not in self.gene_id_to_rast:
                        self.gene_id_to_rast[g_id] = set()
                    self.gene_id_to_rast[g_id].add(rast)
            return self.gene_id_to_rast

    def write_master_faa(self):
        feature_to_genome = {}
        with open('/scratch/fliu/data/cliff/ani_mmseqs/master.faa', 'w') as fh:
            for h in self.genome_paths:
                filename = self.genome_paths.get(h)
                if filename:
                    if 'handle' in filename:
                        genome = MSGenome.from_fasta(filename)
                        for f in genome.features:
                            f_id = f.id
                            if f.seq and len(f.seq) > 0:
                                if f.id not in feature_to_genome:
                                    fh.write(f'>{f_id}\n')
                                    fh.write(f'{f.seq}\n')
                                    feature_to_genome[f.id] = h
                                else:
                                    print('conflict', h, feature_to_genome[f.id])
                    else:
                        genome = MSGenome.from_fasta(filename)
                        for f in genome.features:
                            s = f.id.split()
                            f_id = s[0]
                            if f.seq and len(f.seq) > 0:
                                if f_id not in feature_to_genome:
                                    fh.write(f'>{f_id}\n')
                                    fh.write(f'{f.seq}\n')
                                    feature_to_genome[f_id] = h
                                else:
                                    print('conflict', h, feature_to_genome[f_id])
        return feature_to_genome

    def get_genomeOO_v2(self, filename, feature_to_genome):
        if 'handle' in filename:
            genome = MSGenome.from_fasta(filename)
        else:
            genome = MSGenome.from_fasta(filename)
            for f in genome.features:
                s = f.id.split()
                f_id = s[0]
                if f.seq and len(f.seq) > 0:
                    if f_id not in feature_to_genome:
                        fh.write(f'>{f_id}\n')
                        fh.write(f'{f.seq}\n')
                        feature_to_genome[f_id] = h
                    else:
                        print('conflict', h, feature_to_genome[f_id])

    def catalog_genome_mmseqs(self, genome_id, genome, clusters):
        for f in genome.features:
            cluster_id = clusters[f.id]
            if cluster_id not in self.annotation_database:
                self.annotation_database[cluster_id] = {}
            if genome_id not in self.annotation_database[cluster_id]:
                self.annotation_database[cluster_id][genome_id] = {}
            self.annotation_database[cluster_id][genome_id][f'self:{f.id}'] = float(100)


    def catalog_genome(self, c, genomes, G):
        if self.gene_id_to_rast is None:
            raise Exception('run get_gene_id_to_rast')

        genome = genomes[c]
        for f in genome.features:
            annotation_rast = f.ontology_terms.get('RAST', [])
            for annotation_str in annotation_rast:
                if annotation_str not in self.annotation_database:
                    self.annotation_database[annotation_str] = {}
                if c not in self.annotation_database[annotation_str]:
                    self.annotation_database[annotation_str][c] = {}
                self.annotation_database[annotation_str][c][f'self:{f.id}'] = float(100)
            # print(rast)
        for k, v in G[c].items():
            score = v['weight']
            genome_path = self.genome_paths.get(k)
            if genome_path:
                genome = AnnotationDatabase.get_genome_object(genome_path)
                _feature_to_annotation = {}
                if genome:
                    if 'handle' in genome_path:
                        for f in genome.features:
                            f_id = f.id
                            _feature_to_annotation[f_id] = self.gene_id_to_rast.get(f_id, [])
                    else:
                        for f in genome.features:
                            s = f.id.split()
                            f_id = s[0]
                            _feature_to_annotation[f_id] = self.gene_id_to_rast.get(f_id, [])
                for f_id in _feature_to_annotation:
                    for annotation_str in _feature_to_annotation[f_id]:
                        if annotation_str not in self.annotation_database:
                            self.annotation_database[annotation_str] = {}
                        if c not in self.annotation_database[annotation_str]:
                            self.annotation_database[annotation_str][c] = {}
                        self.annotation_database[annotation_str][c][f'{k}:{f_id}'] = float(score)
        return None

    @staticmethod
    def _catalog_features(genome_id, genome, c_annotation_database, feature_to_rep_seq_h, score: float):
        for f in genome.features:
            if len(f.seq) > 0:
                hseq = HashSeq(f.seq)
                f_hash = hseq.hash_value
                cluster_id = feature_to_rep_seq_h[f_hash]
                if cluster_id not in c_annotation_database:
                    c_annotation_database[cluster_id] = {}
                if genome_id not in c_annotation_database[cluster_id]:
                    c_annotation_database[cluster_id][genome_id] = {}
                c_annotation_database[cluster_id][genome_id][f'self:{f.id}'] = float(score)

    def catalog_all(self, centroids, genomes, G):
        from tqdm import tqdm
        for c in tqdm(centroids):
            self.catalog_genome(c, genomes, G)

    def get_path(self, h, metadata_gtdb_gca, metadata_gtdb_gcf, metadata_gene_cals_path):
        fetch = set()
        if h in metadata_gtdb_gca['h_to_genome']:
            for s in metadata_gtdb_gca['h_to_genome'][h]:
                p = s.split('_')
                ncbi_id = p[0] + '_' + p[1]
                if ncbi_id in metadata_gene_cals_path:
                    fetch.add(metadata_gene_cals_path[ncbi_id])
        elif h in metadata_gtdb_gcf['h_to_genome']:
            for s in metadata_gtdb_gcf['h_to_genome'][h]:
                p = s.split('_')
                ncbi_id = p[0] + '_' + p[1]
                if ncbi_id in metadata_gene_cals_path:
                    fetch.add(metadata_gene_cals_path[ncbi_id])
        else:
            return None

        if len(fetch) == 1:
            return fetch.pop()
        else:
            # print(fetch)
            return None

    def get_genome_paths(self, data, g_metadata):
        everything = {}
        for k in data:
            _data = data[k]
            for v in _data['visited']:
                if v not in everything:
                    everything[v] = None

        base_dir = '/scratch/fliu/data/gtdb/genomes/global/cfs/cdirs/kbase/jungbluth/Projects/Project_Pangenome_GTDB/GTDB_r214_by_spcluster/'
        base_cache_dir = '/scratch/fliu/data/kbase/cache/handle/'
        for k in everything:
            if k in g_metadata:
                everything[k] = base_cache_dir + '/' + g_metadata[k][1]
            else:
                p = self.get_path(k)
                if p:
                    everything[k] = base_dir + '/' + p

    def split(self, score_split=85):
        annotation_database_geq_85 = {}
        annotation_database_lo_85 = {}
        for annotation_str in self.annotation_database:
            annotation_set_geq_85 = {}
            annotation_set_lo_85 = {}
            for g_id in self.annotation_database[annotation_str]:
                for g_other, score in self.annotation_database[annotation_str][g_id].items():
                    if score >= score_split:
                        if g_id not in annotation_set_geq_85:
                            annotation_set_geq_85[g_id] = {}
                        annotation_set_geq_85[g_id][g_other] = score
                    else:
                        if g_id not in annotation_set_lo_85:
                            annotation_set_lo_85[g_id] = {}
                        annotation_set_lo_85[g_id][g_other] = score
            if len(annotation_set_geq_85) > 0:
                annotation_database_geq_85[annotation_str] = annotation_set_geq_85
            if len(annotation_set_lo_85) > 0:
                annotation_database_lo_85[annotation_str] = annotation_set_lo_85

        return annotation_database_geq_85, annotation_database_lo_85
