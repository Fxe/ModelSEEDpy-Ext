from modelseedpy_ext.utils import progress
import numpy as np


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
    for read in tqdm(sam_file.fetch()):
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
