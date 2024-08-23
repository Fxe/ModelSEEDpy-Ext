from modelseedpy.core.msgenome import MSGenome, read_fasta2
from modelseedpy_ext.re.hash_seq import HashSeqList, HashSeq


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
    def from_fasta2(filename, split=" ", h_func=None):
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
