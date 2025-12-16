from modelseedpy_ext.utils import progress
import numpy as np
import os
import json


class Reads:

    def __init__(self, sra_id: str, org: str):
        self.sra_id = sra_id
        self.org = org
        self.bam_file_cache_directory = f'/scratch/fliu/data/ANME/IMG/bam_split_{self.sra_id}_{org}/'
        self.genome_fna = None
        self.alignment_block_file = f'{self.bam_file_cache_directory}/ablock.json'

    def index(self):
        _files = {x for x in os.listdir(self.bam_file_cache_directory) if x.endswith('.bam')}
        if len(_files) == 1:
            _bam_file = list(_files)[0]
            _in = f'{self.bam_file_cache_directory}/{_bam_file}'
            _out = f'{self.bam_file_cache_directory}/{_bam_file}.index'
            cmd = ['/opt/samtools/samtools', 'index', _in, '-o', _out]

            print(' '.join(cmd))

    def count(self):
        import pysam
        # genome_fasta_files = org_fasta[org]
        # o_out = f'/scratch/fliu/data/ANME/IMG/coverm_out/{s}_{org}.out'
        if os.path.exists(self.alignment_block_file) and False:
            print('skip', self.alignment_block_file)
        else:
            _files = {x for x in os.listdir(self.bam_file_cache_directory) if x.endswith('.bam')}
            if len(_files) == 1:
                _bam_file = list(_files)[0]
                _bam = f'{self.bam_file_cache_directory}/{_bam_file}'
                _bam_index = f'{self.bam_file_cache_directory}/{_bam_file}.index'
                # _big_seq = concat_seq(self.genome_fna)
                sam_file = pysam.AlignmentFile(_bam, "rb", index_filename=_bam_index)
                blocks_array = _get_read_cov(sam_file)
                sam_file.close()
                with open(self.alignment_block_file, 'w') as fh:
                    fh.write(json.dumps(blocks_array))

    def aaa(self, _assembly):
        bam_file_cache_directory = f'/scratch/fliu/data/ANME/IMG/bam_split_{s}_{org}/'
        _aligment_block_file = f'{bam_file_cache_directory}/ablock.json'
        _aligment_arr_cov_file = f'{bam_file_cache_directory}/ctg_arr_cov.json'
        print(_aligment_block_file)
        with open(_aligment_block_file, 'r') as fh:
            ablock = json.load(fh)
            ctg_arr_cov = _get_arr_cov(_assembly, ablock)
            with open(_aligment_arr_cov_file, 'w') as fh_out:
                fh_out.write(json.dumps({k: v.tolist() for k, v in ctg_arr_cov.items()}))


def split_reads(filename_fastq, filename_fastq_1, filename_fastq_2):
    with open(filename_fastq, 'r') as fh:
        with open(filename_fastq_1, 'w') as fh_1:
            with open(filename_fastq_2, 'w') as fh_2:
                line_header = fh.readline()
                while line_header:
                    _header_parts = line_header.split(' ')
                    end = int(_header_parts[1][0])
                    line_seq = fh.readline()
                    line_sep = fh.readline()
                    line_qua = fh.readline()
                    if end == 1:
                        fh_1.write('_'.join(_header_parts))
                        fh_1.write(line_seq)
                        fh_1.write(line_sep)
                        fh_1.write(line_qua)
                    elif end == 2:
                        fh_2.write('_'.join(_header_parts))
                        fh_2.write(line_seq)
                        fh_2.write(line_sep)
                        fh_2.write(line_qua)
                    else:
                        raise ValueError(f'unable to determine if read 1 or read 2: {line_header}')

                    line_header = fh.readline()


def concat_seq(file):
    seqs = ""
    with open(file, 'r') as fh:
        line = fh.readline()
        while line:
            if line:
                if line.startswith('>'):
                    pass
                else:
                    line = line.strip()
                    seqs += line
            line = fh.readline()
    return seqs


def _get_read_cov(sam_file, save_seq=False):
    blocks_array = []
    for read in progress(sam_file.fetch()):
        blocks = read.get_blocks()

        if len(blocks) != 0:
            _block = [list(x) for x in blocks]
            _tags = [list(x) for x in read.tags]
            if save_seq:
                read_block = [read.qname, _block, _tags, read.is_paired, read.reference_name, read.seq]
            else:
                read_block = [read.qname, _block, _tags, read.is_paired, read.reference_name]
            blocks_array.append(read_block)
    return blocks_array


def _get_arr_cov(genome, blocks_array):
    contig_cov = {}
    for f in genome.features:
        contig_cov[f.id] = np.zeros(len(f.seq))
    for br in progress(blocks_array):
        _genome_id, _contig_id = br[4].split('~')
        arr_cov = contig_cov[_contig_id]
        for b0, b1 in br[1]:
            _dist = b1 - b0
            for i in range(_dist):
                arr_cov[b0 + i] += 1
    return contig_cov


def _get_arr_cov_with_bigseq(seqs, blocks_array):
    arr_cov = np.zeros(len(seqs))
    for br in progress(blocks_array):
        for b0, b1 in br[1]:
            _dist = b1 - b0 + 1
            for i in range(_dist):
                arr_cov[b0 + i] += 1
    return arr_cov


def xxx(sra, org_names, org_fasta):
    from modelseedpy import MSGenome
    import json
    d_assembly = {}
    for s in sra:
        for org in org_names:
            genome_fasta_files = org_fasta[org]
            bam_file_cache_directory = f'/scratch/fliu/data/ANME/IMG/bam_split_{s}_{org}/'
            _aligment_block_file = f'{bam_file_cache_directory}/ablock.json'
            _aligment_arr_cov_file = f'{bam_file_cache_directory}/ctg_arr_cov.json'
            _assembly = d_assembly.get(genome_fasta_files)
            if _assembly is None:
                print('load', genome_fasta_files)
                _assembly = MSGenome.from_fasta2(genome_fasta_files)
                d_assembly[genome_fasta_files] = _assembly
            print(_aligment_block_file)
            with open(_aligment_block_file, 'r') as fh:
                ablock = json.load(fh)
                ctg_arr_cov = _get_arr_cov(_assembly, ablock)
                with open(_aligment_arr_cov_file) as fh_out:
                    fh_out.write(json.dumps({k: v.tolist() for k, v in ctg_arr_cov.items()}))


def _cal_score_func(arr_cov, file):
    cov_len = len(arr_cov)
    score_func = {}
    gene_func = {}
    scores_unmerged = []
    with open(file, 'r') as fh:
        line = fh.readline()
        while line:
            gene_id, _t, _u, pos_start, pos_end, score, strand, flag, data = line.strip().split('\t')
            _data = {i.split('=')[0]: i.split('=')[1] for i in data.split(';')}
            _product = _data.get('product', None)
            # print(pos_start, pos_end)

            if int(pos_end) > cov_len:
                print('protein out of bounds', pos_end, 'max len', cov_len)
            cov = int(min(arr_cov[int(pos_start): int(pos_end)]))
            scores_unmerged.append(cov)
            if cov not in score_func:
                score_func[cov] = set()
            gene_func[gene_id] = _product
            score_func[cov].add(gene_id)
            # print(gene_id, cov, _product)
            line = fh.readline()
    return score_func, gene_func, scores_unmerged


def _cal_score_func2(ctg_arr_cov, file):
    gene_score = {}
    score_func = {}
    gene_func = {}
    scores_unmerged = []
    with open(file, 'r') as fh:
        line = fh.readline()
        while line:
            location_id, _t, _u, pos_start, pos_end, score, strand, flag, data = line.strip().split('\t')
            arr_cov = ctg_arr_cov[location_id]
            cov_len = len(arr_cov)
            _data = {i[0]: i[1] for i in [x.split('=') for x in data.split(';')]}
            _product = _data.get('product', None)
            # print(pos_start, pos_end)

            if int(pos_end) > cov_len:
                print('protein out of bounds', pos_end, 'max len', cov_len)
            cov = int(min(arr_cov[int(pos_start): int(pos_end)]))
            scores_unmerged.append(cov)
            gene_id = _data['ID']

            if cov not in score_func:
                score_func[cov] = set()
            gene_func[gene_id] = _product
            score_func[cov].add(gene_id)
            gene_score[gene_id] = cov
            # print(gene_id, cov, _product)
            line = fh.readline()
    return score_func, gene_func, scores_unmerged, gene_score

