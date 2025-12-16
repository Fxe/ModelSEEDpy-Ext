from modelseedpy_ext.re.hash_seq import HashSeq
import logging

logger = logging.getLogger(__name__)

class GenomeMetadataFeature:

    def __init__(self, feature_id, feature_type, start, end, contig, strand, seq):
        self.feature_id = feature_id
        self.feature_type = feature_type
        self.start = start
        self.end = end
        self.contig = contig
        self.strand = strand
        self.seq = seq

    def __str__(self):
        return f"{self.feature_id} {self.contig} {self.feature_type} {self.start} {self.end} {self.strand} len: {len(self.seq)}"

def _convert_location(feature_location):
    contig, p0, strand, sz = feature_location
    start = p0
    end = start + sz - 1
    if strand == '-':
        end = p0
        start = end - sz + 1
    return contig, start, end, strand
        
class GenomeMetadata:

    def __init__(self, features: dict[str, GenomeMetadataFeature]):
        self.features = features
        
        self.h_to_feature_id = None
        self.feature_id_to_h = None
        self.hash_features()

        self.feature_start = None
        self.feature_end = None
        self.feature_contig = None  
        self.index_features()
        
        self.h_to_annotation = {}
        self.annotation_to_h = {}

        self.phenotype_predictions = {}
        self.phenotype_assertions = {}

    def index_features(self):
        self.feature_start = {}
        self.feature_end = {}
        self.feature_contig = {}
        for f in self.features.values():
            self.feature_start[f.feature_id] = f.start
            self.feature_end[f.feature_id] = f.end
            self.feature_contig[f.feature_id] = f.contig

    def add_annotation(self, h, kind, s):
        if kind not in self.h_to_annotation:
            self.h_to_annotation[kind] = {}
            self.annotation_to_h[kind] = {}
        if s not in self.annotation_to_h[kind]:
            self.annotation_to_h[kind][s] = set()
        self.annotation_to_h[kind][s].add(h)
        self.h_to_annotation[kind][h] = s

    def hash_features(self):
        self.h_to_feature_id = {}
        self.feature_id_to_h = {}
        for f in self.features.values():
            h_seq = HashSeq(f.seq).hash_value
            if h_seq not in self.h_to_feature_id:
                self.h_to_feature_id[h_seq] = set()
            self.h_to_feature_id[h_seq].add(f.feature_id)
            self.feature_id_to_h[f.feature_id] = h_seq

    def get_feature_by_annotation(self, kind, annotation):
        if kind in self.annotation_to_h and annotation in self.annotation_to_h[kind]:
            found_h = self.annotation_to_h[kind][annotation]
            res = []
            for h_seq in found_h:
                for feature_id in self.h_to_feature_id[h_seq]:
                    res.append(self.features[feature_id])
            return res
        return []

    @staticmethod
    def from_genome_genmark(genome):
        features = {}
        for f in genome.features:
            meta = f.description.split()
            contig = meta[0]
            start = int(meta[1])
            end = int(meta[2])
            strand = meta[3]
            features[f.id] = GenomeMetadataFeature(f.id, 'CDS', start, end, contig, strand, f.seq)

        return GenomeMetadata(features)
    
    @staticmethod
    def from_genome_kbase(genome):
        features = {}
        for f in genome.features:
            contig, start, end, strand = _convert_location(f.location[0])
            features[f.id] = GenomeMetadataFeature(f.id, 'CDS', start, end, contig, strand, f.seq)
            if len(f.location) > 1:
                logger.warning(f"Feature {f.id} has multiple locations: {f.location}")

        return GenomeMetadata(features)

def find_features_in_region(gmm_list: dict[str, GenomeMetadata], bp_distance, bp_center, search_contig, annotation_kind = 'RAST'):
    bp_distance = 10000
    bp_center = 249522
    search_contig = '562.61106.con.0001'
    result = []
    for gmm_id, gmm in gmm_list.items():
        for f in gmm.features.values():
            if f.start > bp_center - bp_distance and f.end < bp_center + bp_distance and f.contig == search_contig:
                annotation = gmm.h_to_annotation[annotation_kind][gmm.feature_id_to_h[f.feature_id]]
                start = f.start
                end = f.end
                if f.strand == '-':
                    start = f.end
                    end = f.start
                result.append({
                    'start': start, 
                    'end': end,
                    'name': f.feature_id,
                    'rast': annotation,
                    'genome': gmm_id,
                })
    return result

def find_features_in_region_by_annotation(gmm_list: dict[str, GenomeMetadata], bp_distance, search_contig, annotation_value, annotation_kind = 'RAST'):
    found_features = {}
    for gmm_id, gmm in gmm_list.items():
        found_features[gmm_id] = []
        for feature in gmm.get_feature_by_annotation(annotation_kind, annotation_value):
            found_features[gmm_id].append(feature)

    result = []
    for gmm_id, features in found_features.items():
        for feature in features:
            bp_center = (feature.start + feature.end) // 2
            result.extend(find_features_in_region(gmm_list, bp_distance, bp_center, search_contig, annotation_kind))

    return result
