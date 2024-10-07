from modelseedpy_ext.utils import progress
import pyarrow as pa
import pyarrow.parquet as pq


def read_ani():
    def read_ani_out(f):
        ani_score = {}
        with open(f, 'r') as fh:
            l = fh.readline()
            while (l):
                _p = l.split('\t')
                g1, g2, score, k1, k2 = _p
                if g1 not in ani_score:
                    ani_score[g1] = {}
                ani_score[g1][g2] = (score, k1, k2)
                l = fh.readline()
        return ani_score

    ani_score_gca = read_ani_out('./ani_gtdb_gca2.txt')
    ani_score_gcf = read_ani_out('./ani_gtdb_gcf2.txt')
    ani_score_self = read_ani_out('./ani_self2.txt')


def read_mmseqs(filename='./data/mmseqs/output/mmseqs_cluster.tsv'):
    feature_to_rep_seq_h = {}
    ct = 0
    with open(filename, 'r') as fh:
        l = fh.readline()
        a, b = [s.strip() for s in l.split('\t')]
        while l:
            if a:
                feature_to_rep_seq_h[b] = a
            l = fh.readline()
            if l:
                a, b = [s.strip() for s in l.split('\t')]
            ct += 1


def write(feature_to_rep_seq_h, output_file='/scratch/fliu/data/cliff/db/mmseqs_clusters.parquet'):
    col_seq_id = []
    col_cluster_id = []
    for s_id, c_id in progress(feature_to_rep_seq_h.items()):
        col_seq_id.append(bytes(bytearray.fromhex(s_id)))
        col_cluster_id.append(bytes(bytearray.fromhex(c_id)))
    table = pa.Table.from_arrays([pa.array(col_cluster_id), pa.array(col_seq_id)],
                                 names=['cluster_id', 'seq_id'])
    pq.write_table(table, output_file)
