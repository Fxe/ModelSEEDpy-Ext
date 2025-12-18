import uuid
from modelseedpy_ext.re.hash_seq import HashSeq


def convert_kbase_location(feature_location):
    contig, p0, strand, sz = feature_location
    start = p0
    end = start + sz - 1
    if strand == '-':
        end = p0
        start = end - sz + 1
    return contig, start, end, strand


def generate_uuid():
    return str(uuid.uuid4())


def generate_cdm_uuid(prefix):
    _uuid = generate_uuid()
    return f'{prefix}-{_uuid}'


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
