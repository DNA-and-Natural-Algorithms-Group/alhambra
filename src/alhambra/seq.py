# Utility functions for sequences
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple, Union

_L_TO_N = {
    "A": frozenset((0,)),
    "B": frozenset((1, 2, 3)),
    "C": frozenset((1,)),
    "D": frozenset((0, 2, 3)),
    "G": frozenset((2,)),
    "H": frozenset((0, 1, 3)),
    "K": frozenset((2, 3)),
    "M": frozenset((0, 1)),
    "N": frozenset((0, 1, 2, 3)),
    "R": frozenset((0, 2)),
    "S": frozenset((1, 2)),
    "T": frozenset((3,)),
    "V": frozenset((0, 1, 2)),
    "W": frozenset((0, 3)),
    "Y": frozenset((1, 3)),
}

_N_TO_L = {v: i for i, v in _L_TO_N.items()}


# X as synonym for N is allowed but not reversible
_L_TO_N["X"] = frozenset((0, 1, 2, 3))

# Additional convenience: uppercase letters
_L_TO_N |= {k.lower(): v for k, v in _L_TO_N.items()}

_WC = {k: _N_TO_L[frozenset(3 - v for v in s)] for k, s in _L_TO_N.items() if k}
_WC |= {
    k: _N_TO_L[frozenset(3 - v for v in s)].lower()
    for k, s in _L_TO_N.items()
    if k.islower()
}
_WC["X"] = "X"
_WC["x"] = "x"

_PUNC = "-+ \n\t"

_WC_WITH_PUNC = _WC | {x: x for x in _PUNC}

_AMBBASES = frozenset(_L_TO_N.keys() - {"a", "c", "g", "t", "A", "C", "G", "T"})

_VALID_NTS = _L_TO_N.keys()

_VALID_WITH_PUNC = frozenset(_WC_WITH_PUNC.keys())


def revcomp(seq_str: str) -> str:
    """Return the reverse complement of a string sequence.  Preserves
    whitespace and capitalization.  Does not do any sequence checking."""
    return "".join(_WC_WITH_PUNC[c] for c in reversed(seq_str))


def is_null(seq_str: str | None, _check_seq: bool = True) -> bool:
    """Return True if a sequence consists only of Ns (or Xs), or is empty.
    Return False otherwise."""
    if seq_str is None:
        return True
    if len(seq_str) == 0:
        return True
    if _check_seq:
        check_seq(seq_str)  # fixme: do we need this?
    return set(seq_str.lower()).issubset("nx" + _PUNC)


def is_definite(seq_str: str | None, _check_seq: bool = True) -> bool:
    """Return True if a sequence consists only of defined bases.  Return
    False otherwise.  If blank, return False.
    """
    if seq_str is None:
        return False
    if len(seq_str) == 0:
        return False
    if _check_seq:
        check_seq(seq_str)
    if set(seq_str).issubset(_PUNC):
        return False
    return set(seq_str).issubset("acgtACGT" + _PUNC)


def check_seq(seq_str: str):
    """Raises an error (`InvalidSequence`) if the sequence has invalid elements."""
    if not set(seq_str).issubset(_VALID_WITH_PUNC):
        invalids = [(i, v) for i, v in enumerate(seq_str) if v not in _VALID_WITH_PUNC]
        raise InvalidSequence(seq_str, invalids)


def count_ambiguous(seq_str: str, _check_seq: bool = True) -> int:
    """Return the number of ambiguous bases in a sequence."""
    if _check_seq:
        check_seq(seq_str)
    return sum(1 for x in seq_str if x in _AMBBASES)


def dna_length(seq_str: str, _check_seq: bool = True) -> int:
    """Return the length of a sequence, stripping whitespace.  This does not handle
    extended labels, etc."""
    check_seq(seq_str)
    return sum(1 for c in seq_str if c not in _PUNC)


def merge(seq1: str, seq2: str, _check_seq: bool = True) -> str:
    """Merge two sequences together, returning a single sequence that
    represents the constraint of both sequences.  If the sequences
    can't be merged, raise a MergeConflictError."""

    if _check_seq:
        check_seq(seq1)
        check_seq(seq2)

    if dna_length(seq1) != dna_length(seq2):
        raise MergeConflictError(seq1, seq2, "length", len(seq1), len(seq2))

    # Whitespace and capitalization of the *first* sequence is preserved.
    out = []
    i1 = 0
    i2 = 0
    while i1 < len(seq1):
        c1 = seq1[i1]
        c2 = seq2[i2]
        if c1 in _PUNC:
            out.append(seq1[i1])
            i1 += 1
            continue
        if c2 in _PUNC:
            i2 += 1
            continue
        try:
            cn = _N_TO_L[_L_TO_N[c1] & _L_TO_N[c2]]
        except KeyError as e:
            if e.args[0] == frozenset():
                raise MergeConflictError(seq1, seq2, (i1, i2), c1, c2) from None
            else:
                raise e
        if c1.isupper():
            out.append(cn.upper())
        else:
            out.append(cn.lower())
        i1 += 1
        i2 += 1
    return "".join(out)


@dataclass
class InvalidSequence(ValueError):
    sequence: str
    invalids: List[Tuple[int, str]]


@dataclass
class MergeConflictError(ValueError):
    """
    Merge of items failed because of conflicting information.
    Arguments are (item1, item2, location or property, value1, value2)
    """

    item1: str
    item2: str
    location: Union[str, Tuple[int, int], int]
    value1: Union[str, int]
    value2: Union[str, int]


class MergeConflictsError(ValueError):
    """
    Merge of multiple items failed because individual merges
    raised MergeConflictErrors.
    Arguments are ([list of MergeConflictErrors])
    """


@dataclass(frozen=True, init=False)
class Seq:
    seq_str: str

    def __init__(self, seq_str: Union[Seq, str]) -> None:
        if isinstance(seq_str, Seq):
            object.__setattr__(self, "seq_str", seq_str.seq_str)
            # fixme: could we pass through instead?
        else:
            check_seq(seq_str)
            object.__setattr__(self, "seq_str", seq_str)

    @property
    def revcomp(self) -> Seq:
        return Seq(revcomp(self.seq_str))

    @property
    def is_null(self) -> bool:
        return is_null(self.seq_str, _check_seq=False)

    @property
    def is_definite(self) -> bool:
        return is_definite(self.seq_str, _check_seq=False)

    @property
    def dna_length(self) -> int:
        return dna_length(self.seq_str, _check_seq=False)

    @property
    def n_ambiguous(self) -> int:
        return count_ambiguous(self.seq_str, _check_seq=False)

    def __add__(self, other: Seq | str) -> Seq:
        return Seq(self.seq_str + Seq(other).seq_str)

    def __radd__(self, other: Seq | str) -> Seq:
        return Seq(Seq(other).seq_str + self.seq_str)

    def merge(self, other: Union[Seq, str, None]) -> Seq:
        if other is None:
            return self
        if not isinstance(other, Seq):
            other = Seq(other)
        return Seq(merge(self.seq_str, other.seq_str, _check_seq=False))

    def __or__(self, other: Union[Seq, str, None]) -> Seq:
        return self.merge(other)

    def __ror__(self, other: Union[Seq, str, None]) -> Seq:
        return self.merge(other)

    def __str__(self) -> str:
        return self.seq_str

    def __repr__(self) -> str:
        return repr(self.seq_str)

    def __len__(self) -> int:
        return self.dna_length

    @property
    def base_str(self) -> str:
        return "".join(c for c in self.seq_str if c not in _PUNC)

    def __getitem__(self, ix: Union[int, slice]) -> Seq:
        raise NotImplementedError("Seq getitem needs to handle whitespace.")
