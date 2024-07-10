import logging
from modelseedpy_ext.re.hash_seq import HashSeq
from modelseedpy_ext.utils import progress


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


class Loader:

    def __init__(self, store, block_size=500, fn_validation=None):
        self.block_size = block_size
        self.loaded_data = set()
        self.fn_load_storage = None
        self.fn_transform_data = None
        self.staging = set()
        self.store = store
        self.fn_validation = fn_validation

    def in_cache(self, data):
        return data in self.loaded_data

    def add_data_to_staging(self, data):
        self.staging.add(data)

    def load(self, data):
        valid = True

        if self.fn_transform_data:
            data = self.fn_transform_data(data)
        if self.fn_validation:
            valid = self.fn_validation(data)

        if valid:
            if not self.in_cache(data):
                self.add_data_to_staging(data)
                if len(self.staging) >= self.block_size:
                    logger.debug(f"Load {len(self.staging)} data to storage.")
                    self.load_to_storage()
            else:
                logger.debug(f"skip data in cache")

    def clear_stating(self):
        self.staging.clear()

    def cache_loaded_data(self, loaded_data):
        self.loaded_data |= set(loaded_data)

    def load_to_storage(self):
        res = self.fn_load_storage(self.staging)
        self.cache_loaded_data(res)
        self.clear_stating()
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


class AbstractArangoLoader(Loader):

    def __init__(self, db, collection_id, block_size=500):
        super().__init__(db, block_size=block_size)
        self._init_collection(collection_id)
        self.collection = self.store[collection_id]
        self.collection_id = collection_id
        self.staging = []
        self.load_stats = {
            'error_true': 0,
            'error_false': 0
        }

    def in_cache(self, data):
        return data['_key'] in self.loaded_data

    def add_data_to_staging(self, data):
        self.staging.append(data)

    def cache_loaded_data(self, loaded_data):
        self.loaded_data |= {x['_key'] for x in loaded_data}

    def _update_load_stats(self, stats: dict):
        for k in stats:
            if k == 'error':
                if stats[k]:
                    self.load_stats['error_true'] += 1
                else:
                    self.load_stats['error_false'] += 1
            else:
                if k not in self.load_stats:
                    self.load_stats[k] = 0
                self.load_stats[k] += stats[k]

    def _init_collection(self, collection_id):
        pass

    def load_to_storage(self):
        res = self.collection.importBulk(self.staging)
        self._update_load_stats(res)

        self.cache_loaded_data(self.staging)
        self.clear_stating()
        logger.debug(f"Loaded {len(self.loaded_data)} objects.")


class ArangoLoaderNodes(AbstractArangoLoader):

    def __init__(self, db, collection_id, block_size=500):
        super().__init__(db, collection_id, block_size=block_size)
        self.fn_transform_data = ArangoLoaderNodes.node_to_json

    def _init_collection(self, collection_id):
        if not self.store.hasCollection(collection_id):
            logger.warning(f"create {collection_id}")
            self.store.createCollection(name=collection_id)

    @staticmethod
    def node_to_json(n):
        return n.to_json()

    def load(self, graph):
        for k, node in progress(graph.t_nodes[self.collection_id].items()):
            super(ArangoLoaderNodes, self).load(node)


class ArangoLoaderEdges(AbstractArangoLoader):

    def __init__(self, db, collection_id, block_size=500):
        super().__init__(db, collection_id, block_size=block_size)
        self.fn_transform_data = ArangoLoaderEdges.edge_to_json

    def _init_collection(self, collection_id):
        if not self.store.hasCollection(collection_id):
            logger.warning(f"create {collection_id}")
            self.store.createCollection(name=collection_id, className="Edges")

    @staticmethod
    def edge_to_json(t: tuple):
        """
        Covert tuple with source (src), destination (dst) and edge data
        :param t: tuple with src, dst, data
        :return:
        """
        src, dst, data = t
        data = data if data else {}
        node_from = src
        node_to = dst
        data["_key"] = f"{node_from.key}:{node_to.key}"
        data["_from"] = node_from.id
        data["_to"] = node_to.id
        return data

    def load(self, graph):
        for src, dst in progress(graph.t_edges[self.collection_id]):
            data = graph.get_edge_data(src, dst)['data']
            super(ArangoLoaderEdges, self).load((src, dst, data))
