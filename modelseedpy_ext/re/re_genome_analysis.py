import logging
import itertools
from modelseedpy import MSGenome
from Bio.Seq import Seq

logger = logging.getLogger(__name__)


class KEGenome:

    def __init__(self, genome, assembly=None, contigs=None, rast=None, kb_re=None,
                 load_arango=None, transform_rast=None):
        self.kb_re = kb_re
        self.seq_protein_store = kb_re.protein_store
        self.seq_dna_store = kb_re.dna_store
        self.genome = genome
        self.assembly = assembly
        self.rast = rast
        self.p_annotation = {}
        self.contig_data = contigs
        self.contigs = {}
        self.load_arango = load_arango
        self.transform_rast = transform_rast
        if self.contig_data:
            for contig in self.contig_data.features:
                self.contigs[contig.id] = [contig.seq, str(Seq(contig.seq).complement())]

        self.d_hash = {}
        self.p_hash = {}
        for o in self.genome.features:
            if o.dna_sequence:
                hash_d = self.seq_dna_store.store_sequence(o.dna_sequence)
                self.d_hash[o.id] = hash_d
            if o.protein_translation:
                hash_p = self.seq_protein_store.store_sequence(o.protein_translation)
                self.p_hash[o.id] = hash_p
            elif len(o.dna_sequence) % 3 == 0:
                try:
                    protein_translation = Seq(o.dna_sequence).translate()
                    hash_p = self.seq_protein_store.store_sequence(protein_translation)
                    self.p_hash[o.id] = hash_p
                except ValueError as e:
                    logger.debug(f"unable to translate {e}")

        self.loc_pos = {}
        self.gene_pairs = set()
        self.f_pairs = {}
        self.null_pairs = {}
        self.anno_pairs = {}
        self.null_annotation = {'rast/hypotheticalprotein', 'rast/null', None}

    def run(self):
        self.fetch_annotation()
        self.loc_pos = self.get_loc_pos()
        self.gene_pairs = set()
        for contig_id in self.loc_pos:
            pairs = self.get_pairs(self.loc_pos[contig_id])
            self.gene_pairs |= pairs
        self.f_pairs = self.get_f_pairs(self.gene_pairs)
        self.null_pairs, self.anno_pairs = self.get_pp(self.f_pairs)

    @staticmethod
    def load_from_kbase(kbase, object_id, ws, kb_re, rast, load_arango, transform_rast):
        print('get genome')
        genome = kbase.get_from_ws(object_id, ws)
        print('get assembly')
        assembly = kbase.get_from_ws(genome.assembly_ref)
        print('get assembly handle_ref:', assembly.fasta_handle_ref)
        kbase.download_file_from_kbase2(assembly.fasta_handle_ref,
                                        '/home/fliu/kbase/cache/handle/' + assembly.fasta_handle_ref)
        print('read contigs')
        ms_contigs = MSGenome.from_fasta('/home/fliu/kbase/cache/handle/' + assembly.fasta_handle_ref, split=' ')

        return KEGenome(genome, assembly, ms_contigs, rast, kb_re, load_arango, transform_rast)

    def _fetch_annotation_arango(self):
        aql = """
        FOR doc_protein IN re_seq_protein_has_rast_annotation
            FILTER doc_protein._from IN @h
            RETURN [doc_protein._from, doc_protein._to]
        """
        h_list = list({'re_seq_protein/' + x[1] for x in self.p_hash.items()})
        res = self.kb_re.db.AQLQuery(aql, rawResults=True, bindVars={"h": h_list})
        for h, annotation in res:
            if h not in self.p_annotation:
                self.p_annotation[h] = annotation

        return h_list

    def fetch_annotation(self):
        h_list = self._fetch_annotation_arango()
        self.etl_rast(h_list)
        self._fetch_annotation_arango()

        return self.p_annotation

    def etl_rast(self, h_list, batch_size=20):
        q = {}
        for i in h_list:
            if i not in self.p_annotation:
                q[i] = self.seq_protein_store.get_sequence(i.split('/')[-1])
                pass
            if len(q) >= batch_size:
                p_features = [{'id': x[0], "protein_translation": x[1]} for x in q.items()]
                res = self.rast.f(p_features)
                _rast_annotations = {}
                for o in res[0]['features']:
                    _function_value = 'null' if 'function' not in o else o['function']
                    _rast_annotations[o['id'].split('/')[-1]] = _function_value
                _nodes, _edges = self.transform_rast.transform(_rast_annotations)
                self.load_arango.load({
                    'nodes': _nodes,
                    'edges': _edges
                })
                q = {}
        if len(q) > 0:
            p_features = [{'id': x[0], "protein_translation": x[1]} for x in q.items()]
            res = self.rast.f(p_features)
            _rast_annotations = {}
            for o in res[0]['features']:
                _function_value = 'null' if 'function' not in o else o['function']
                _rast_annotations[o['id'].split('/')[-1]] = _function_value
            _nodes, _edges = self.transform_rast.transform(_rast_annotations)
            self.load_arango.load({
                'nodes': _nodes,
                'edges': _edges
            })

    def get_loc_pos(self, skip_seq_match=False):
        loc_pos = {}
        for f in self.genome.features:
            ll = locate_feature_dna_sequence_in_contig(f, self.contigs, skip_seq_match)
            contig_id = [x[0] for x in ll]
            if len(contig_id) == 1:
                contig_id = contig_id[0]
            else:
                raise ValueError('!!!' + str(contig_id))
            if contig_id not in loc_pos:
                loc_pos[contig_id] = {}
            if len(ll) == 0:
                print(ll, f.id)
                break
            else:
                p_min = None
                p_max = None
                for o in ll:
                    if p_min is None:
                        p_min = o[2]
                    elif p_min > o[2]:
                        p_min = o[2]
                    if p_min > o[3]:
                        p_min = o[3]
                    if p_max is None:
                        p_max = o[2]
                    elif p_max < o[2]:
                        p_max = o[2]
                    if p_max < o[3]:
                        p_max = o[3]
                if p_min >= p_max:
                    raise ValueError(f'assertion fail for p_min < pmax: ({p_min}, {p_max})')
                if p_min not in loc_pos[contig_id]:
                    loc_pos[contig_id][p_min] = []
                loc_pos[contig_id][p_min].append([ll, p_min, p_max, [], []])

        return loc_pos

    def get_pairs(self, loc_pos):
        p_gene_id_pairs = set()
        s_pmin = sorted(loc_pos)

        def get_gene(d1):
            res = {}
            for loc in d1:
                loc_sites, p_min, p_max, member_l, member_r = loc
                for loc_site in loc_sites:
                    gene_id = loc_site[1]
                    if gene_id not in res:
                        res[gene_id] = []
                    res[gene_id].append(loc_site)
            return res

        for i in range(len(s_pmin) - 1):
            d1 = loc_pos[s_pmin[i]]
            d2 = loc_pos[s_pmin[i + 1]]
            g1 = get_gene(d1)
            g2 = get_gene(d2)
            somelists = [set(g1), set(g2)]
            for element in itertools.product(*somelists):
                p_gene_id_pairs.add(element)
            # print(set(g1), set(g2))

        return p_gene_id_pairs

    def get_f_pairs(self, p_gene_id_pairs):
        f_pairs = {}
        for g1, g2 in p_gene_id_pairs:

            if g1 in self.p_hash:
                f1 = self.p_annotation['re_seq_protein/' + self.p_hash[g1]]
            else:
                f1 = None
            if g2 in self.p_hash:
                f2 = self.p_annotation['re_seq_protein/' + self.p_hash[g2]]
            else:
                f2=None

            f_pair = (f1, f2)
            if f_pair not in f_pairs:
                f_pairs[f_pair] = 0
            f_pairs[f_pair] += 1

        return f_pairs

    def get_pp(self, f_pairs1):
        null_pairs1 = {}
        anno_pairs1 = {}
        for p in f_pairs1:
            p1, p2 = p
            if p1 in self.null_annotation or p2 in self.null_annotation:
                null_pairs1[p] = f_pairs1[p]
            else:
                anno_pairs1[p] = f_pairs1[p]
        return null_pairs1, anno_pairs1


def locate_feature_dna_sequence_in_contig(f, contigs, skip_seq_match=False):
    ll = []
    contig_sub_strings = []

    for l in f.location:
        contig_id, start, direction, seq_size = l
        if contig_id in contigs:
            off_set = seq_size
            sub_str_start = start - 1
            sub_str_end = sub_str_start + off_set
            if direction == '+':
                contig_sub_string = contigs[contig_id][0][sub_str_start:(sub_str_start + seq_size)]
                contig_sub_strings.append(contig_sub_string)
                # if contig_sub_string == f.dna_sequence:
                ll.append([contig_id, f.id, sub_str_start, sub_str_start + seq_size, False, direction, l])
                # print(f.id, contig_sub_string[:10], contig_sub_string[-10:])
                # print(f.id, f.dna_sequence[:10], f.dna_sequence[-10:])
            elif direction == '-':
                contig_sub_string = contigs[contig_id][1][(sub_str_start - seq_size + 1):(sub_str_start + 1)]
                contig_sub_string = contig_sub_string[::-1]
                contig_sub_strings.append(contig_sub_string)
                ll.append([contig_id, f.id, sub_str_start - seq_size + 1, sub_str_start + 1, True, direction, l])
                # print(f.id, contig_sub_string[:10], contig_sub_string[-10:])
                # print(f.id, f.dna_sequence[::-1][:10], f.dna_sequence[::-1][-10:])
            else:
                print('not sure what this means:', direction)

            # print(f.id, l, f.functions, contig_sub_string == f.dna_sequence)
        else:
            print('contig not found:', contig_id)
    if skip_seq_match or ''.join(contig_sub_strings) == f.dna_sequence:
        return ll
    return []


# ll = locate_feature_dna_sequence_in_contig(f, contigs)
# ll

