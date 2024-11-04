from modelseedpy_ext.re.hash_seq import HashSeqList, HashSeq


def _t_genome(features_gff, genome_cds, genome_contigs, contig_id_to_hash,
              data_features, data_feature_attributes, data_proteins, feature_auto_id):
    from Bio.Seq import Seq
    from Bio import Align

    proteins_read = set()
    for f in features_gff:
        data_features[feature_auto_id] = {
            'start': f.start,
            'end': f.end,
            'strand': f.strand,
            'feature_type': f.feature_type,
            'phase': f.phase,
            'source': f.source,
            'dna_protein_aligment_score': None,
            'protein_hash': None,
            'contig_hash': contig_id_to_hash[f.contig_id],
        }
        data_feature_attributes[feature_auto_id] = dict(f.attr)
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
                if not eq:
                    aligner = Align.PairwiseAligner()
                    res = aligner.align(feature_cds.seq, seq_from_dna)
                    data_features[feature_auto_id]['dna_protein_aligment_score'] = res.score
                else:
                    data_features[feature_auto_id]['dna_protein_aligment_score'] = len(feature_cds.seq)
                    # print(eq, f.contig_id, f.start, f.end, f.strand, len(feature_cds.seq), len(seq_from_dna), res.score)

                if feature_cds.seq:
                    id_protein_seq = HashSeq(feature_cds.seq).hash_value
                    if id_protein_seq not in data_proteins:
                        data_proteins[id_protein_seq] = {
                            'length': len(feature_cds.seq),
                            'protein_sequence': feature_cds.seq
                        }
                    print(data_features[feature_auto_id])
                    data_features[feature_auto_id]['protein_hash'] = id_protein_seq
            else:
                pass
        feature_auto_id += 1
    return feature_auto_id
