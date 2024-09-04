from modelseedpy.core.msgenome import MSGenome, read_fasta2
from modelseedpy_ext.re.hash_seq import HashSeqList, HashSeq


class GffRecord:

    def __init__(self, contig_id, source, feature_type, start, end, score, strand, phase, attr):
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
