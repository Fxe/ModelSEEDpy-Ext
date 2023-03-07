import hashlib
import pymongo
import zlib


class SeqStoreMongo:
    def __init__(self, mongo_database, col_name: str):
        self.CHARSET_VALIDATION = {}
        self.db = mongo_database
        self.col_name = col_name

    @property
    def collection(self):
        return self.db[self.col_name]

    def validate_sequence(self, s):
        v = set(s) - self.CHARSET_VALIDATION
        if len(v) == 0:
            return True
        raise ValueError("Invalid input - unaccepted characters:" + ", ".join(v))

    @staticmethod
    def get_sequence_hash(seq):
        h = hashlib.sha256(seq.encode("utf-8")).hexdigest()
        return h

    def store_sequence(self, seq):
        self.validate_sequence(seq)
        h = self.get_sequence_hash(seq)
        if self.collection.find_one({"_id": h}) is None:
            bin_data = zlib.compress(seq.encode("utf-8"))
            self.collection.insert_one({"_id": h, "z_seq": bin_data})
        return h

    def store_sequences(self, seqs: list, batch_size=20):
        docs = []
        res = []
        for seq in seqs:
            self.validate_sequence(seq)
            h = self.get_sequence_hash(seq)
            bin_data = zlib.compress(seq.encode("utf-8"))
            docs.append({"_id": h, "z_seq": bin_data})
            res.append(h)

        skip = set()
        for doc in self.collection.find({"_id": {"$in": [doc["_id"] for doc in docs]}}):
            skip.add(doc["_id"])

        docs_size = len(docs)
        chunks = int(len(docs) / batch_size) + 1
        for k in range(chunks):
            batch = []
            for i in range(k * batch_size, (k + 1) * batch_size):
                if i < docs_size:
                    if docs[i]["_id"] not in skip:
                        skip.add(docs[i]["_id"])
                        batch.append(docs[i])
            if len(batch) > 0:
                self.collection.insert_many(batch)

        return res

    def get_sequence(self, hash_val):
        doc = self.collection.find_one({"_id": hash_val})
        if doc:
            bin_data = doc["z_seq"]
            return zlib.decompress(bin_data).decode("utf-8")

        return None


def load_dna_seq_store_mongo(host, database):
    mongo_client = pymongo.MongoClient(host)
    s = SeqStoreMongo(mongo_client[database], "seq_dna")
    # K	Guanine / Thymine
    # M	Adenine / Cytosine
    # S	Guanine / Cytosine
    # W	Adenine / Thymine
    # B	Guanine / Thymine / Cytosine
    # D	Guanine / Adenine / Thymine
    # H	Adenine / Cytosine / Thymine
    # V	Guanine / Cytosine / Adenine
    s.CHARSET_VALIDATION = {
        "A",
        "C",
        "G",
        "T",
        "U",  # Adenine, Cytosine, Guanine, Thymine, Uracil
        "R",
        "Y",  # Guanine / Adenine (purine), Cytosine / Thymine (pyrimidine)
        "K",
        "M",
        "S",
        "W",
        "B",
        "D",
        "H",
        "V",
        "N",  # Adenine / Guanine / Cytosine / Thymine
    }
    return s


def load_protein_seq_store_mongo(host, database):
    mongo_client = pymongo.MongoClient(host)
    s = SeqStoreMongo(mongo_client[database], "seq_protein")
    s.CHARSET_VALIDATION = {
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
        "O",  # Pyrrolysine
        "U",  # Selenocysteine
        "B",  # Aspartate or Asparagine
        "Z",  # Glutamate or Glutamine
        "J",  # Leucine or Isoleucine
        "X",  # Any
    }
    return s
