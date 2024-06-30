from modelseedpy_ext.seq_store_mongo import load_dna_seq_store_mongo, load_protein_seq_store_mongo
from modelseedpy_ext.re.re_arangodb import RE
from pyArango.connection import Connection


def init_re(mongo_host='mongodb://localhost:27017/',
            mongo_database='database',
            arango_url='http://127.0.0.1:8529/',
            arango_user=None,
            arango_pwd=None):

    seq_protein_store = load_protein_seq_store_mongo(mongo_host, mongo_database)
    seq_dna_store = load_dna_seq_store_mongo(mongo_host, mongo_database)
    conn_sequoia = Connection(username=arango_user, password=arango_pwd, arangoURL=arango_url)
    arango_database = conn_sequoia['test']
    kb_re = RE(arango_database, None, None, seq_dna_store, seq_protein_store)
    kb_re.init_collections()
    return kb_re
