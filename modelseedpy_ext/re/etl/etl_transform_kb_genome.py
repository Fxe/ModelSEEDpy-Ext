import hashlib
import logging
from Bio.Seq import Seq
from modelseedpy_ext.re.etl.etl_transform_graph import ETLTransformGraph

def locate_feature_dna_sequence_in_contig(f, contigs):
    ll = []
    for l in f.data['location']:
        contig_id, start, direction, seq_size = l
        if contig_id in contigs:
            off_set = seq_size
            sub_str_start = start - 1    
            sub_str_end = sub_str_start + off_set
            if direction == '+':
                contig_sub_string = contigs[contig_id][0][sub_str_start:(sub_str_start + seq_size)]
                if contig_sub_string == f.dna_sequence:
                    ll.append([contig_id, f.id, sub_str_start, sub_str_start + seq_size, False, direction])
                #print(f.id, contig_sub_string[:10], contig_sub_string[-10:])
                #print(f.id, f.dna_sequence[:10], f.dna_sequence[-10:])
            elif direction == '-':
                contig_sub_string = contigs[contig_id][1][(sub_str_start - seq_size + 1):(sub_str_start + 1)]
                if contig_sub_string == f.dna_sequence[::-1]:
                    ll.append([contig_id, f.id, sub_str_start - seq_size + 1, sub_str_start + 1, True, direction])
                #print(f.id, contig_sub_string[:10], contig_sub_string[-10:])
                #print(f.id, f.dna_sequence[::-1][:10], f.dna_sequence[::-1][-10:])
            else:
                print('not sure what this means:', direction)
            
            #print(f.id, l, f.functions, contig_sub_string == f.dna_sequence)
        else:
            print('contig not found:', contig_id)
    return ll

class ETLTransformKBaseGenome(ETLTransformGraph):
    
    def __init__(self, dna_store, protein_store):
        self.dna_store = dna_store
        self.protein_store = protein_store
    
    def transform(self, kb_genome, contigs):
        nodes = {}
        edges = {}

        def add_node(node_id, label, data=None):
            if label not in nodes:
                nodes[label] = {}
            if label in nodes:
                _node = self.build_node(node_id, label, data)
                if _node.id not in nodes[label]:
                    nodes[label][_node.id] = _node
                else:
                    # print('dup', _node.id)
                    pass

                return nodes[label][_node.id]

            logger.error('add_node error')
            return None

        def add_edge(node_from, node_to, label, data=None):
            if label not in edges:
                edges[label] = []
            if label in edges:
                _edge = self.transform_edge(node_from, node_to, data)
                edges[label].append(_edge)
                return _edge

            logger.error('add_edge error')
            return None
        
        _contigs = {}
        _contigs_hash = {}
        for contig in contigs.features:
            _contigs[contig.id] = [contig.seq, str(Seq(contig.seq).complement())]
            _contigs_hash[contig.id] = self.dna_store.store_sequence(contig.seq)
        
        
        print(len(_contigs))
        print(_contigs_hash)
        
        for f in kb_genome.features:
            seq_dna = f.dna_sequence
            seq_protein = f.protein_translation
            h_dna = self.dna_store.store_sequence(seq_dna)
            h_protein = self.protein_store.store_sequence(seq_protein)
            node_seq_dna = add_node(h_dna, 're_seq_dna', data=None)
            node_seq_protein = add_node(h_protein, 're_seq_protein', data=None)
            
            seq_translation = str(Seq(seq_dna).translate())[:-1]
            if seq_translation == seq_protein:
                add_edge(node_seq_dna, node_seq_protein, 're_seq_dna_has_translation')
            else:
                add_edge(node_seq_dna, node_seq_protein, 're_seq_dna_has_partial_translation')
        
        return nodes, edges