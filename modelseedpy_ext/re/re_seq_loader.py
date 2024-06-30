import logging
from modelseedpy_ext.re.hash_seq import HashSeq


logger = logging.getLogger(__name__)


class SeqLoader:

    def __init__(self, store, block_size=500, fn_validation=None):
        self.block_size = block_size
        self.loaded_sequences = set()
        self.staging = set()
        self.store = store
        self.fn_validation = fn_validation

    def load(self, seq):
        valid = True
        h = HashSeq(seq)
        if self.fn_validation:
            valid = self.fn_validation(h)

        if valid:
            hash_value = h.hash_value
            if hash_value not in self.loaded_sequences:
                self.staging.add(h)
                if len(self.staging) >= self.block_size:
                    logger.debug(f"Load {len(self.staging)} sequences to storage.")
                    self.load_to_storage()
            else:
                logger.debug(f"skip loaded sequence hash {hash_value}")

    def load_to_storage(self):
        res = self.store.store_sequences(self.staging)
        self.loaded_sequences |= set(res)
        self.staging.clear()
        logger.debug(f"Loaded {len(self.loaded_sequences)} sequences.")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if len(self.staging) > 0:
            self.load_to_storage()
        self.store = None
        self.staging = None
        self.loaded_sequences = None
        logger.debug(f"Exit")
