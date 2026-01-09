from modelseedpy.core.msgenome import MSGenome, read_fasta2
from modelseedpy_ext.re.hash_seq import HashSeqList, HashSeq
from collections import Counter


class ProteinSequence(HashSeq):

    _STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
    _EXTENDED_ONLY = {"U", "O"}  # selenocysteine, pyrrolysine
    _AMBIGUOUS_AA = set("BJZX*")

    def __new__(cls, sequence, strip_ending_star=True):
        obj = super().__new__(cls, sequence, strip_ending_star=strip_ending_star)
        return obj

    @property
    def sequence(self) -> str:
        return str(self)

    def is_standard(self) -> bool:
        """True if the sequence contains only standard amino acids."""
        return bool(self) and set(self) <= self._STANDARD_AA

    def is_extended(self) -> bool:
        """True if all residues are standard or extended (U, O)."""
        allowed = self._STANDARD_AA | self._EXTENDED_ONLY
        return bool(self) and set(self) <= allowed

    def is_ambiguous(self) -> bool:
        """True if the sequence contains any ambiguous residue codes."""
        return any(res in self._AMBIGUOUS_AA for res in self)

    def is_valid(self) -> bool:
        """True if all residues are recognized AA codes (standard, extended, or ambiguous)."""
        allowed = self._STANDARD_AA | self._EXTENDED_ONLY | self._AMBIGUOUS_AA
        return bool(self) and set(self) <= allowed

    def z_compress(self):
        import zlib
        return zlib.compress(str(self).encode("utf-8"))


def convert_kbase_location(feature_location):
    contig, p0, strand, sz = feature_location
    start = p0
    end = start + sz - 1
    if strand == '-':
        end = p0
        start = end - sz + 1
    return contig, start, end, strand


def get_gff_contig_id(f):
    _id = f.id
    contig_id = _id[::-1].split('_', 1)[1][::-1] # reverse split once reverse again
    return contig_id


def get_gff_coord(f):
    _p_start, _p_end, _p_strand, _p_props = f.description[2:].split(' # ')
    start = int(_p_start)
    end = int(_p_end)
    strand = '?'
    if _p_strand == '1':
        strand = '+'
    elif _p_strand == '-1':
        strand = '-'
    else:
        raise ValueError(f'bad strand: {_p_strand}')
    return start, end, strand


def get_gff_location(f):
    contig = get_gff_contig_id(f)
    start, end, strand = get_gff_coord(f)
    return contig, start, end, strand


def get_kbase_location(f):
    locations = f.location
    if len(locations) > 1:
        raise ValueError('!')

    contig, start, end, strand = convert_kbase_location(f.location[0])
    return contig, start, end, strand


def get_coords_genemark(f):
    contig, start, end, strand, props = f.description.split(' ', 4)
    return contig, int(start), int(end), strand


def _from_str(s):
    contig_id, source, feature_type, start, end, score, strand, phase, attr_str = s.strip().split('\t')
    attr_str = attr_str[:-1] if attr_str[-1] == ';' else attr_str
    attr = dict([x.split('=') for x in attr_str.split(';')])
    return GffRecord(contig_id, source, feature_type, int(start), int(end), score, strand, phase, attr)


def _read_gff_features(f):
    if f.endswith('.gz'):
        import gzip

        with gzip.open(f, "rb") as fh:
            features_gff = []
            _data = fh.read().decode("utf-8")
            for line in _data.split('\n'):
                if not line.startswith('#'):
                    if line:
                        features_gff.append(_from_str(line))

            return features_gff
    else:
        with open(f, "r") as fh:
            features_gff = []
            _data = fh.read()
            for line in _data.split('\n'):
                if not line.startswith('#'):
                    if line:
                        features_gff.append(_from_str(line))
            return features_gff


class CDMContigSet:

    def __init__(self, sha256):
        self.sha256 = sha256
        self.contigs = []


class CDMContig:

    def __init__(self, contig_set_id: str, seq: str):
        #self.seq = seq
        self.contig_set_id = contig_set_id
        self.hash = HashSeq(self.seq).hash_value
        self.base_count = dict(Counter(list(self.seq.upper())))
        self.length = len(self.seq)
        self.gc = (self.base_count.get('G', 0) + self.base_count.get('C', 0)) / self.length

        self.names = []

    def __repr__(self):
        return f'len: {self.length}, gc: {self.gc}, base_count: {self.base_count}, names: {self.names}'


class CDMProtein:

    def __init__(self, seq: str):
        _seq = seq
        if _seq[-1] == '*':
            _seq = _seq[:-1]
        self.seq = _seq
        self.hash = HashSeq(self.seq).hash_value
        self.length = len(self.seq)

        self.names = []

    def __repr__(self):
        return f'len: {self.length}, hash: {self.hash}'


class CDMFeature:

    def __init__(self, feature_id: str, contig_hash, start, end, strand, attributes=None):
        self.id = feature_id
        self.contig_hash = contig_hash
        self.start = start
        self.end = end
        self.strand = strand
        self.attributes = {} if attributes is None else attributes

        self.names = []


class GffRecord:

    def __init__(self, contig_id: str, source: str,
                 feature_type,
                 start: int, end: int, score, strand, phase, attr):
        self.contig_id = contig_id
        self.source = source
        self.feature_type = feature_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attr = attr

    def get_attribute_string(self):
        attr_values = []
        for k, v in self.attr.items():
            attr_values.append(f"{k}={v}")
        return ';'.join(attr_values)

    def __str__(self):
        return '\t'.join([str(x) for x in [self.contig_id, self.source, self.feature_type,
                                           self.start, self.end, self.score, self.strand, self.phase,
                                           self.get_attribute_string()]])

    @staticmethod
    def from_str(s):
        contig_id, source, feature_type, start, end, score, strand, phase, attr_str = s.strip().split('\t')
        attr = dict([x.split('=') for x in attr_str.split(';')])
        return GffRecord(contig_id, source, feature_type, int(start), int(end), score, strand, phase, attr)


def write_fna_fai(file_fna_fai, genome_assembly):
    with open(file_fna_fai, 'w') as fh:
        for f in genome_assembly.features:
            _p = f.id.split()
            fh.write('\t'.join([str(x) for x in [_p[0], len(f.seq), 0, 60, 61]]) + '\n')


def write_gff(file_gff, genome, rast=None, builder=None):
    if rast:
        for feature in genome.features:
            feature.ontology_terms['RAST'] = []
        rast.annotate_genome(genome)
    model = None
    if builder:
        model = builder.build('model', annotate_with_rast=False)
    gff_records = []
    for f in genome.features:
        reactions = 'null'
        if model:
            if f.id in model.genes:
                gene = model.genes.get_by_id(f.id)
                reactions = ':'.join({r.id for r in gene.reactions})
        rast_functions = f.ontology_terms.get('RAST')
        if rast_functions is None:
            rast_functions = ['hypothetical protein']
        for rast in rast_functions:
            for l in f.location:
                # print(l)
                start = l[1]
                size = l[3]
                end = start + size - 1
                d_type = 'CDS'
                if reactions != 'null':
                    d_type = 'ModelSEED'
                elif rast == 'hypothetical protein':
                    d_type = 'HP'
                gff_record = GffRecord(l[0], 'test', d_type, start, end, '.', l[2], '.',
                                       {'rast': rast, 'feature_id': f.id, 'gene': rast, 'reactions': reactions})
                gff_records.append(gff_record)

    with open(file_gff, 'w') as fh:
        for r in gff_records:
            fh.write(str(r) + '\n')


class REGenome(MSGenome):

    def __init__(self):
        super().__init__()

    def re(self):
        pass

    def ke(self):
        pass

    def re_hash(self):
        pass


class REAssembly(MSGenome):

    def __init__(self):
        super().__init__()
        self.hash_list = HashSeqList()
        for contig in self.features:
            seq = HashSeq(contig.seq)
            self.hash_list.append(seq)

    def re(self):
        pass

    def ke(self):
        pass

    @staticmethod
    def from_fasta(filename, split=" ", h_func=None):
        genome = REAssembly()
        genome.features += read_fasta2(filename, split, h_func)
        return genome

    @property
    def hash_value(self):
        hl = HashSeqList()
        for contig in self.features:
            seq = HashSeq(contig.seq)
            hl.append(seq)
        return hl.hash_value

    @staticmethod
    def _process_contigs(contigs):
        hash_list = HashSeqList()
        contig_h_d = []
        for contig in contigs.features:
            seq = HashSeq(contig.seq)
            hash_list.append(seq)
            seq_h = seq.hash_value
            contig_h_d.append([seq_h, contig.id, contig.description])
        return {
            'genome_h': hash_list.hash_value,
            'contig_h': contig_h_d
        }
