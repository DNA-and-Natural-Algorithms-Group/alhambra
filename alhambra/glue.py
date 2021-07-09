from __future__ import annotations
from dataclasses import dataclass

from .seq import Seq
from typing import (
    Any,
    ClassVar,
    Generic,
    Iterable,
    List,
    Literal,
    Mapping,
    MutableMapping,
    Optional,
    Protocol,
    Type,
    TypeVar,
    Union,
)
from warnings import warn
import copy
from enum import Enum
import collections.abc
from .classes import UpdateListD


T = TypeVar("T")


class Use(Enum):
    UNUSED = 0
    NULL = 1
    INPUT = 2
    OUTPUT = 3
    BOTH = 4
    PERMANENT = 5

    def invert(self) -> Use:
        return Use([0, 1, 3, 2, 4, 5][self.value])

    def merge(self, other: Use | None) -> Use:
        if other is None:
            return self
        raise NotImplementedError

    def __or__(self, other: Use | None) -> Use:
        raise NotImplementedError

    def __ror__(self, other: Use | None) -> Use:
        raise NotImplementedError


def merge_items(a: Optional[T], b: Optional[T]) -> Optional[T]:
    if a is None:
        return b
    elif b is None:
        return a
    else:
        assert a == b
        return a


@dataclass(init=False)
class Glue:
    name: Optional[str] = None
    note: Optional[str] = None
    use: Optional[Use] = None
    abstractstrength: Optional[int] = None
    _fields: ClassVar[tuple[tuple[str, str], ...]] = (
        ("name", "name"),
        ("note", "note"),
        ("use", "use"),
        ("abstractstrength", "abstractstrength"),
    )

    def __init__(
        self,
        name: Optional[str] = None,
        note: Optional[str] = None,
        use: Optional[Use] = None,
        abstractstrength: Optional[int] = None,
    ):
        if name and name[-1] == "/":
            name = name[:-1] + "*"
        self.name = name
        self.note = note
        self.use = use
        self.abstractstrength = abstractstrength

    def _into_complement(self):
        if self.name is not None:
            if self.name[-1] == "*":
                self.name = self.name[:-1]
            else:
                self.name = self.name + "*"

    def __repr__(self) -> str:
        p = []
        for f, d in type(self)._fields:
            if (v := getattr(self, f)) != getattr(type(self), f, None):
                p.append(f"{d}={repr(v)}")
        return type(self).__name__ + "(" + ", ".join(p) + ")"

    def copy(self: T) -> T:
        return copy.copy(self)

    def ident(self) -> str:
        if self.name:
            return self.name
        else:
            raise ValueError

    def basename(self) -> str:
        if self.name is None:
            raise ValueError
        if self.name[-1] == "*":
            return self.name[:-1]
        else:
            return self.name

    @property
    def complement(self) -> Glue:
        c = self.copy()
        c._into_complement()
        return c

    def update(self, other: Glue):
        if type(other) == Glue:
            self.name = merge_items(self.name, other.name)
        else:
            raise NotImplementedError

    def merge(self, other: Glue) -> Glue:
        if type(other) == Glue:
            return Glue(merge_items(self.name, other.name))
        else:
            return other.merge(self)

    def __or__(self, other: Glue) -> Glue:
        return self.merge(other)

    def __ior__(self, other: Glue):
        self.update(other)

    def __ror__(self, other: Glue) -> Glue:
        return self.merge(other)

    def as_dict(self) -> dict[str, Any]:
        return {
            k: v for k in ["name", "use", "note"] if (v := getattr(self, k)) is not None
        }

    @staticmethod
    def from_dict(d: dict[str, Any]) -> Glue:
        return glue_factory.from_dict(d)


class GlueFactory:
    types: dict[str, Type[Glue]]

    def __init__(self):
        self.types = {}

    def register(self, c: Type[Glue]):
        self.types[c.__name__] = c

    def from_dict(self, d: str | dict[str, Any]) -> Glue:
        if isinstance(d, str):
            return Glue(d)
        if "type" in d:
            c = self.types[d["type"]]
            del d["type"]
            return c(**d)
        else:
            return Glue(**d)


glue_factory = GlueFactory()


@dataclass(init=False)
class SSGlue(Glue):
    _sequence: Seq
    _fields: ClassVar[tuple[tuple[str, str], ...]] = Glue._fields + (
        ("_sequence", "sequence"),
    )

    def __init__(
        self,
        name: Optional[str] = None,
        length: Union[None, int, str, Seq] = None,
        sequence: Union[None, str, Seq] = None,
        note: Optional[str] = None,
        use: Optional[Use] = None,
    ):
        super().__init__(name, note, use)

        if isinstance(length, int):
            lseq: Seq | None = Seq("N" * length)
        elif isinstance(length, str):
            lseq = Seq(length)
        elif isinstance(length, Seq):
            lseq = length
        else:
            lseq = None

        if sequence is not None:
            if not isinstance(sequence, Seq):
                sequence = Seq(sequence)
            self._sequence = sequence | lseq  # fixme: better error
        elif lseq is not None:
            self._sequence = lseq
        else:
            raise ValueError("Must have at least length or sequence.")

    @property
    def dna_length(self) -> int:
        return self.sequence.dna_length

    @property
    def sequence(self) -> Seq:
        return self._sequence

    @sequence.setter
    def sequence(self, seq: Seq | str | None):  # type: ignore
        if seq is None:
            self._sequence = Seq("N" * self.dna_length)
            return
        elif not isinstance(seq, Seq):
            seq = Seq(seq)
        if self.dna_length is not None:
            assert seq.dna_length == self.dna_length
        self._sequence = seq

    def ident(self) -> str:
        if self.name:
            return super().ident()
        if self.sequence:
            return f"SSGlue_{self.sequence.base_str}"
        else:
            raise ValueError

    def _into_complement(self):
        if self.sequence is not None:
            self.sequence = self.sequence.revcomp
        super()._into_complement()

    def update(self, other: Glue | str):
        if isinstance(other, str):
            other = Glue(other)
        if type(other) is Glue:  # a base glue: we can merge
            self.name = merge_items(self.name, other.name)
        elif type(other) is SSGlue:
            newname = merge_items(self.name, other.name)
            newseq = self.sequence | other.sequence
            self.name = newname
            self.sequence = newseq
        else:
            raise NotImplementedError

    def merge(self, other: Glue | str) -> SSGlue:
        new = self.copy()
        new.update(other)
        return new

    def __repr__(self):
        s = f"SSGlue({repr(self.name)}, {self.dna_length}"
        if not self.sequence.is_null:
            s += f", {repr(self.sequence.seq_str)}"
        return s + ")"

    def as_dict(self) -> dict[str, Any]:
        d = super().as_dict()
        d["type"] = self.__class__.__name__
        d["sequence"] = self.sequence.seq_str
        return d

    @property
    def _shortdesc(self) -> str:
        r = []
        if self.name is not None:
            r.append(self.name)
            if (self.sequence is not None) and not self.sequence.is_null:
                r.append(f"({self.sequence.seq_str})")
        else:
            r.append(self.sequence.seq_str)
        return "".join(r)


glue_factory.register(SSGlue)


class DXGlue(Glue):
    etype: Literal["TD", "DT"]
    _fullseq: Seq
    _fields: ClassVar[tuple[tuple[str, str], ...]] = Glue._fields + (
        ("etype", "etype"),
        ("_fullseq", "fullseq"),
    )

    def _into_complement(self):
        if self.fullseq is not None:
            self.fullseq = self.fullseq.revcomp
        super()._into_complement()

    def __init__(
        self,
        etype,
        name=None,
        length=None,
        sequence: Seq | str | None = None,
        fullseq=None,
        use=None,
        abstractstrength=None,
        note=None,
    ):
        super().__init__(name, note, use, abstractstrength)
        self.etype = etype
        trial_fseq: Optional[Seq] = None
        if length:
            trial_fseq = Seq("N" * (length + 1))
        if sequence:
            if self.etype == "TD":
                if trial_fseq:
                    trial_fseq |= "N" + sequence
                else:
                    trial_fseq = Seq("N" + sequence)
            elif self.etype == "DT":
                if trial_fseq:
                    trial_fseq |= sequence + "N"
                else:
                    trial_fseq = Seq(sequence + "N")
            else:
                raise ValueError(etype)
        if fullseq:
            if trial_fseq:
                trial_fseq |= fullseq
            else:
                trial_fseq = Seq(fullseq)
        if trial_fseq is None:
            raise ValueError("Must have a least some information.")
        self._fullseq = trial_fseq

    @property
    def fseq(self) -> Optional[str]:
        if self._fullseq is not None:
            return self._fullseq.seq_str
        else:
            return None

    @fseq.setter
    def fseq(self, value: str):
        if self._fullseq is not None and len(value) != len(self._fullseq):
            warn("Changing end length")
        self._fullseq = Seq(value)

    @property
    def seq(self) -> Optional[str]:
        if not self.fseq:
            return None
        if self.etype == "TD":
            return self.fseq[1:]
        elif self.etype == "DT":
            return self.fseq[:-1]

    @property
    def comp(self):
        """The complement end sequences of the End, as a string."""
        if not self._fullseq:
            return None
        if self.etype == "TD":
            return self._fullseq.revcomp.base_str[1:]
        elif self.etype == "DT":
            return self._fullseq.revcomp.base_str[:-1]

    def merge(self, other: Glue) -> DXGlue:
        out = self.copy()
        if type(other) not in [Glue, DXGlue]:
            raise ValueError
        for k in ["note", "name", "etype", "abstractstrength"]:
            if (v := getattr(out, k, None)) is not None:
                if (nv := getattr(other, k, None)) is not None:
                    if nv != v:
                        raise ValueError
            else:
                if (nv := getattr(other, k, None)) is not None:
                    setattr(out, k, nv)
        if isinstance(other, DXGlue):
            if out._fullseq:
                out._fullseq = out._fullseq.merge(other._fullseq)
            if out.use and other.use:
                out.use = out.use | other.use
        return out

    def __str__(self):
        if self.fseq and self.seq and self.comp:
            if self.etype == "DT":
                s = self.seq[0] + "-" + self.seq[1:]
                c = self.comp[0] + "-" + self.comp[1:]
            elif self.etype == "TD":
                s = self.seq[:-1] + "-" + self.seq[-1]
                c = self.comp[:-1] + "-" + self.comp[-1]
            else:
                raise ValueError
            return "<dxend {} ({}{}): {} | {}>".format(
                self.name, self.etype, len(self.seq), s, c
            )
        else:
            return "<dxend {} ({})>".format(self.name, getattr(self, "etype", "?"))


SomeGlue = TypeVar("SomeGlue", bound=Glue)


class GlueList(Generic[SomeGlue], UpdateListD[SomeGlue]):
    def merge_complements(self):
        newitems: dict[str, SomeGlue] = {}
        for v in self:
            c = v.complement
            kc = c.ident()
            if kc in self.data:
                self.data[kc] = self.data[kc].merge(c)
            else:
                newitems[kc] = c
        self.data.update(newitems)

    def merge_glue(self, g: SomeGlue) -> SomeGlue:
        if g.ident() in self.data:
            g = self.data[g.ident()].merge(g)
        c = g.complement
        if c.ident() in self.data:
            g = self.data[c.ident()].complement.merge(g)
        return g

    def merge_glue_and_update_list(self, g: SomeGlue) -> SomeGlue:
        kg = g.ident()
        if kg in self.data:
            g = self.data[kg].merge(g)
            self.data[kg] = g
        c = g.complement
        kc = c.ident()
        if kc in self.data:
            c = self.data[kc].merge(c)
            self.data[kc] = c
            g = c.complement
            if kg in self.data:
                self.data[kg] = g
        return g
