from modelseedpy_ext.re.re_arangodb import RE
from modelseedpy_ext.seq_store_mongo import SeqStoreMongo
import pymongo
from pyArango.connection import Connection


def init_re(mongo_host='mongodb://localhost:27017/',
            mongo_database='database',
            arango_url='http://127.0.0.1:8529/',
            arango_user=None,
            arango_pwd=None,
            timeout=30):

    # seq_protein_store = load_protein_seq_store_mongo(mongo_host, mongo_database)
    # seq_dna_store = load_dna_seq_store_mongo(mongo_host, mongo_database)
    conn_sequoia = Connection(username=arango_user, password=arango_pwd, arangoURL=arango_url, timeout=timeout)
    arango_database = conn_sequoia['test']

    mongo_client = pymongo.MongoClient(mongo_host)[mongo_database]
    store_contigs = SeqStoreMongo(mongo_client, "seq_contigs")
    store_contigs.CHARSET_VALIDATION = {
        "A",
        "C",
        "G",
        "T",
        # "U", # Adenine, Cytosine, Guanine, Thymine, Uracil
        "R",  # G/A Guanine / Adenine (puRine)
        "Y",  # C/T Cytosine / Thymine (pyYimidine)
        "K",  # G/T Keto
        "M",  # A/C aMino
        "S",  # Guanine or cytosine
        "W",  # Adenine or thymine
        "B",  # Guanine / Thymine / Cytosine
        "D",  # Guanine / Adenine / Thymine
        "H",  # A or C or T not-G, H follows G in the alphabet
        "V",  # G or C or A not-T (not-U), V follows U
        "N",  # Adenine / Guanine / Cytosine / Thymine
    }
    store_dna = SeqStoreMongo(mongo_client, "seq_dna")
    store_dna.CHARSET_VALIDATION = {
        "A",
        "C",
        "G",
        "T",
        # "U", # Adenine, Cytosine, Guanine, Thymine, Uracil
        "R",  # G/A Guanine / Adenine (puRine)
        "Y",  # C/T Cytosine / Thymine (pyYimidine)
        "K",  # G/T Keto
        "M",  # A/C aMino
        "S",  # Guanine or cytosine
        "W",  # Adenine or thymine
        "B",  # Guanine / Thymine / Cytosine
        "D",  # Guanine / Adenine / Thymine
        "H",  # A or C or T not-G, H follows G in the alphabet
        "V",  # G or C or A not-T (not-U), V follows U
        "N",  # Adenine / Guanine / Cytosine / Thymine
    }
    store_protein = SeqStoreMongo(mongo_client, "seq_protein")
    store_protein.CHARSET_VALIDATION = {
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
        # "O",  # Pyrrolysine
        "U",  # Selenocysteine
        "B",  # Aspartate or Asparagine
        "Z",  # Glutamate or Glutamine
        "J",  # Leucine or Isoleucine
        "X",  # Any
        # "*"
    }

    def validate_protein(s):
        tokens = set(s) - store_protein.CHARSET_VALIDATION
        if len(tokens) > 0:
            raise ValueError(f'Invalid tokens {tokens}')
        return True

    def validate_dna(s):
        tokens = set(s) - store_dna.CHARSET_VALIDATION
        if len(tokens) > 0:
            raise ValueError(f'Invalid tokens {tokens}')
        return True

    def validate_contigs(s):
        tokens = set(s.upper()) - store_contigs.CHARSET_VALIDATION
        if len(tokens) > 0:
            raise ValueError(f'Invalid tokens {tokens}')
        return True

    kb_re = RE(arango_database, None, None, store_contigs, store_dna, store_protein,
               validate_contigs, validate_dna, validate_protein)

    return kb_re
