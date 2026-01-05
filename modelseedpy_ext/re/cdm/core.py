import uuid
from modelseedpy_ext.re.hash_seq import HashSeq


def generate_uuid():
    return str(uuid.uuid4())


def generate_cdm_uuid(prefix):
    _uuid = generate_uuid()
    return f'{prefix}-{_uuid}'


class ClusterSet:

    def __init__(self, clusters: dict):
        self.members = {}
        self.clusters = clusters
        self.m_to_c = {}
        for c, members in self.clusters.items():
            for m in members:
                self.m_to_c[m] = c

    def get_cluster_by_member(self, member_id):
        return self.m_to_c[member_id]

    def get_cluster_members(self, cluster_id):
        return self.clusters[cluster_id]

    @staticmethod
    def from_mmseqs2(filename_cluster_tsv):
        from collections import defaultdict

        c_to_m = defaultdict(set)

        with open(filename_cluster_tsv) as fh:
            for line in fh:
                rep, mem = line.strip().split("\t")
                c_to_m[rep].add(mem)

        return ClusterSet(c_to_m)


class ProteinSequence(HashSeq):

    _STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
    _EXTENDED_ONLY = {"U", "O"}  # selenocysteine, pyrrolysine
    _AMBIGUOUS_AA = set("BJZX")

    def __new__(cls, sequence, strip_ending_star=True):
        obj = super().__new__(cls, sequence, strip_ending_star=strip_ending_star)
        return obj

    @property
    def sequence(self) -> str:
        return str(self)

    def is_standard(self) -> bool:
        """True if the sequence contains only standard amino acids."""
        return bool(self) and set(self) <= self._STANDARD_AA

    def is_extended(self) -> bool:
        """True if all residues are standard or extended (U, O)."""
        allowed = self._STANDARD_AA | self._EXTENDED_ONLY
        return bool(self) and set(self) <= allowed

    def is_ambiguous(self) -> bool:
        """True if the sequence contains any ambiguous residue codes."""
        return any(res in self._AMBIGUOUS_AA for res in self)

    def is_valid(self) -> bool:
        """True if all residues are recognized AA codes (standard, extended, or ambiguous)."""
        allowed = self._STANDARD_AA | self._EXTENDED_ONLY | self._AMBIGUOUS_AA
        return bool(self) and set(self) <= allowed

    def z_compress(self):
        import zlib
        return zlib.compress(str(self).encode("utf-8"))
