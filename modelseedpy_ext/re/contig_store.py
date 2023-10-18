import os
import zlib
import hashlib
from modelseedpy.core.msgenome import parse_fasta_str, MSGenome


def get_sequence_hash(seq):
    h = hashlib.sha256(seq.encode("utf-8")).hexdigest()
    return h


def hash_contig_set(contig_set: MSGenome):
    contigs = {}
    hash_list = []
    for contig in contig_set.features:
        seq = contig.seq
        if seq is None:
            raise ValueError('invalid contig')

        seq = seq.upper()
        h = get_sequence_hash(seq)
        hash_list.append(h)
        contigs[h] = seq

    contig_set_hash = None

    hash_list = sorted(hash_list)
    hash_seq = "_".join(hash_list)
    res = hashlib.sha256(hash_seq.encode("utf-8")).hexdigest()

    return res


class ContigLibrary:
    BLOCK_PREFIX = 'block_'

    def __init__(self, path, mongo_database, max_block_size=5000):
        if not os.path.exists(path) or not os.path.isdir(path):
            raise ValueError(f'invalid path: {path}')
        self.path = path
        self._current_block = None
        self.max_block_size = max_block_size
        self.mongo_collection = mongo_database['contig_store']

    def get_contig_file(self, h):
        d = self.mongo_collection.find_one({'_id': h})
        if d:
            return d['p']

        return None

    def get_contig_set(self, h):

        # get folder
        filepath = self.get_contig_file(h)
        if filepath:
            with open(filepath, 'rb') as fh:
                contig_set = MSGenome()
                contig_set.features += parse_fasta_str(zlib.decompress(fh.read()).decode('utf-8'))
                return contig_set
        else:
            return None

    def _is_current_block_free(self):
        if self._current_block is None:
            return False
        contents = self.block_content(self._current_block)
        return len(contents) < self.max_block_size

    @staticmethod
    def to_fasta(genome, l=80, fn_header=None):
        for feature in genome.features:
            h = f">{feature.id}\n"
            if fn_header:
                h = fn_header(feature)
            yield h
            lines = [
                feature.seq[i: i + l] + "\n" for i in range(0, len(feature.seq), l)
            ]
            for line in lines:
                yield line

    def _get_block_path(self, block_id):
        return f'{self.path}/{self.BLOCK_PREFIX}{block_id}'

    def put_contig_set(self, genome):
        h = hash_contig_set(genome)
        # print(h)
        filepath = self.get_contig_file(h)
        # print(filepath)
        if filepath:
            return filepath

        if not self._is_current_block_free():
            self._current_block = self._locate_free_block_id()

        _block_path = self._get_block_path(self._current_block)
        # print('current_block', _block_path)

        if not os.path.exists(_block_path):
            os.mkdir(_block_path)

        # write file
        filepath = f'{self.path}/{self.BLOCK_PREFIX}{self._current_block}/{h}.zl'
        # print('write file', filepath)
        data = ''
        for l in self.to_fasta(genome):
            data += l
        with open(filepath, 'wb') as fh:
            fh.write(zlib.compress(data.encode("utf-8"), zlib.Z_BEST_COMPRESSION))
        self.mongo_collection.insert_one({'_id': h, 'p': filepath})
        return filepath

    def block_content(self, block_id):
        content = set()
        p = self._get_block_path(block_id)
        for o in os.listdir(p):
            if o.endswith('.zl'):
                content.add(o)
        return content

    def _locate_free_block_id(self):
        max_id = None
        # print('finding free block')
        for i in os.listdir(self.path):
            if i.startswith(self.BLOCK_PREFIX) and os.path.isdir(self.path + '/' + i):
                block_id = int(i[len(self.BLOCK_PREFIX):])
                if max_id is None or block_id > max_id:
                    max_id = block_id
                contents = self.block_content(block_id)
                # print(len(contents), self.max_block_size)
                if self._is_current_block_free():
                    return block_id

        if max_id is None:
            return 0
        else:
            return max_id + 1
