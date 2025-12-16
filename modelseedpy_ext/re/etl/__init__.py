from modelseedpy_ext.re.hash_seq import HashSeqList, HashSeq
from Bio import SeqUtils
from Bio.Seq import Seq
from Bio import Align


class SomeETL:

    def __init__(self):
        self.feature_auto_id = 0
        self.data_contig_set = {}
        self.data_contig = {}
        self.data_features = {}
        self.data_feature_attributes = {}
        self.data_proteins = {}
        pass

    def aaaaa(self, features_gff, genome_cds, genome_contigs):
        id_contig_set = genome_contigs.hash_value

        if id_contig_set not in self.data_contig_set:
            self.data_contig_set[id_contig_set] = {
                'contigs': len(genome_contigs.features),
                'contig_bp': sum([len(f.seq) for f in genome_contigs.features]),
                'ambiguous_bases': None,
                'gc_content': SeqUtils.gc_fraction(Seq(''.join([f.seq for f in genome_contigs.features]))),
            }
            contig_id_to_hash = {}
            for f in genome_contigs.features:
                id_contig = HashSeq(f.seq).hash_value
                contig_id_to_hash[f.id] = id_contig
                if (id_contig, id_contig_set) not in self.data_contig:
                    self.data_contig[(id_contig, id_contig_set)] = {
                        'length': len(f.seq),
                        'gc_content': SeqUtils.gc_fraction(Seq(f.seq)),
                        'ambiguous_bases': None,
                        'alias': f.id
                    }
                    #if id_contig not in self.data_contig_set[id_contig_set]['contigs']:
                    #    self.data_contig_set[id_contig_set]['contigs'][id_contig] = []
                    #self.data_contig_set[id_contig_set]['contigs'][id_contig].append(f.id)
            feature_auto_id = self._t_genome(features_gff, genome_cds, genome_contigs, id_contig_set, contig_id_to_hash)
        else:
            print('dup', id_contig_set)

    def _t_genome(self, features_gff, genome_cds, genome_contigs, contig_set_hash, contig_id_to_hash):

        proteins_read = set()
        for f_index in range(len(features_gff)):
            f = features_gff[f_index]
            self.data_features[self.feature_auto_id] = {
                'start': f.start,
                'end': f.end,
                'strand': f.strand,
                'feature_type': f.feature_type,
                'phase': f.phase,
                'source': f.source,
                'dna_protein_aligment_score': None,
                'protein_hash': None,
                'contig_hash': contig_id_to_hash[f.contig_id],
                'contig_set_hash': contig_set_hash,
            }
            self.data_feature_attributes[self.feature_auto_id] = dict(f.attr)
            if f.feature_type == 'CDS':
                protein_id = f.attr.get('protein_id')
                if protein_id:
                    proteins_read.add(protein_id)
                    feature_cds = genome_cds.features.get_by_id(protein_id)
                    feature_contig = genome_contigs.features.get_by_id(f.contig_id)
                    seq = Seq(feature_contig.seq[f.start - 1:f.end])
                    if f.strand == '-':
                        seq = seq.reverse_complement()
                    seq_from_dna = str(seq.translate())[:-1]
                    eq = feature_cds.seq == seq_from_dna
                    if not eq and len(seq_from_dna) > 0:
                        try:
                            aligner = Align.PairwiseAligner()
                            res = aligner.align(feature_cds.seq, seq_from_dna)
                            self.data_features[self.feature_auto_id]['dna_protein_aligment_score'] = res.score
                        except ValueError as ex:
                            print('error', f_index)
                            raise ex
                    else:
                        self.data_features[self.feature_auto_id]['dna_protein_aligment_score'] = len(feature_cds.seq)
                        # print(eq, f.contig_id, f.start, f.end, f.strand, len(feature_cds.seq), len(seq_from_dna), res.score)

                    if feature_cds.seq:
                        id_protein_seq = HashSeq(feature_cds.seq).hash_value
                        if id_protein_seq not in self.data_proteins:
                            self.data_proteins[id_protein_seq] = {
                                'length': len(feature_cds.seq),
                                'protein_sequence': feature_cds.seq
                            }
                        #print(self.data_features[self.feature_auto_id])
                        self.data_features[self.feature_auto_id]['protein_hash'] = id_protein_seq
                else:
                    pass
            self.feature_auto_id += 1
        return self.feature_auto_id
