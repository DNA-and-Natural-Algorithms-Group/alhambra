from __future__ import annotations

import collections.abc
import copy
from dataclasses import dataclass
from enum import Enum
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
    cast,
)
from warnings import warn

from .classes import UpdateListD
from .seq import Seq

import attrs

T = TypeVar("T")

GlueA = TypeVar("GlueA", bound="Glue")
GlueB = TypeVar("GlueB", bound="Glue")


class Use(int, Enum):
    UNUSED = 0
    NULL = 1
    INPUT = 2
    OUTPUT = 3
    BOTH = 4
    PERMANENT = 5
    UNSET = 6

    @classmethod
    def from_any(cls, other: str | Use | int) -> Use:
        if isinstance(other, str):
            return Use[other.upper()]
        elif isinstance(other, Use):
            return other
        else:
            return Use(other)

    def invert(self) -> Use:
        if self == Use.UNSET:
            raise ValueError
        return Use([0, 1, 3, 2, 4, 5][self.value])

    def merge(self: Use, other: Use) -> Use:
        if self == other:
            return self
        if self == Use.UNSET:
            return other
        elif other == Use.UNSET:
            return self
        elif self == Use.UNUSED:
            return other
        elif other == Use.UNUSED:
            return self
        elif (self == Use.INPUT and other == Use.OUTPUT) or (
            self == Use.OUTPUT and other == Use.INPUT
        ):
            return Use.BOTH
        elif (self == Use.BOTH and 2 <= other.value <= 4) or (
            2 <= self.value <= 4 and self == Use.BOTH
        ):
            return Use.BOTH
        raise ValueError

    def __or__(self, other: Any) -> Use:
        raise NotImplementedError

    def __ror__(self, other: Any) -> Use:
        raise NotImplementedError


def merge_items(a: Optional[T], b: Optional[T]) -> Optional[T]:
    if a is None:
        return b
    elif b is None:
        return a
    else:
        assert a == b
        return a


def _ensure_glue_name(n: str | None) -> str | None:
    "Ensure the glue name uses *, not /, and is reasonable."
    if n is not None:
        if (len(n) >= 1) and n[-1] == "/":
            n = n[:-1] + "*"
    return n


@attrs.define()
class Glue:
    name: Optional[str] = attrs.field(
        converter=_ensure_glue_name, on_setattr=attrs.setters.convert, default=None
    )
    note: Optional[str] = attrs.field(default=None)
    use: Use = attrs.field(default=Use.UNSET)
    abstractstrength: Optional[int] = attrs.field(default=None)
    """The stickydesign type of the glue."""

    @property
    def etype(self) -> str:
        raise NotImplementedError

    def _into_complement(self):
        if self.name is not None:
            if self.is_complement:
                self.name = self.name[:-1]
            else:
                self.name = self.name + "*"

    def copy(self: T) -> T:
        return copy.copy(self)

    def ident(self) -> str:
        if self.name:
            return self.name
        else:
            raise ValueError

    @property
    def type(self) -> str:
        return self.__class__.__name__

    def basename(self) -> str:
        if self.name is None:
            raise ValueError
        if self.is_complement:
            return self.name[:-1]
        else:
            return self.name

    @property
    def is_complement(self: Glue) -> bool:
        if self.name and (self.name[-1] == "*"):
            return True
        return False

    @property
    def complement(self: GlueA) -> GlueA:
        c = self.copy()
        c._into_complement()
        return c

    def update(self, other: Glue):
        if type(other) == Glue:
            self.name = merge_items(self.name, other.name)
            self.abstractstrength = merge_items(
                self.abstractstrength, other.abstractstrength
            )
            self.use = Use.merge(self.use, other.use)
        else:
            raise NotImplementedError

    def merge(self: Glue, other: Glue) -> Glue:
        if type(other) == Glue:
            return Glue(
                merge_items(self.name, other.name),
                abstractstrength=merge_items(
                    self.abstractstrength, other.abstractstrength
                ),
                use=Use.merge(self.use, other.use),
            )
        else:
            return other.merge(self)

    def __or__(self, other: Glue) -> Glue:
        return self.merge(other)

    def __ior__(self, other: Glue):
        self.update(other)

    def __ror__(self, other: Glue) -> Glue:
        return self.merge(other)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            other = Glue(other)
        if not isinstance(other, Glue):
            return False
        return self.ident() == other.ident()

    def to_dict(self) -> dict[str, Any]:
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


@attrs.define(init=False)
class SSGlue(Glue):
    etype: str = "S"
    _sequence: Seq = attrs.field(default=Seq(""))

    def __init__(
        self,
        name: Optional[str] = None,
        length: Union[None, int, str, Seq] = None,
        sequence: Union[None, str, Seq] = None,
        note: Optional[str] = None,
        use: Optional[Use] = None,
    ):
        if use is None:
            use = Use.UNSET
        super(SSGlue, self).__init__(name, note, use)

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

        self.etype = "S"  # FIXME: there should be a better way to do this

    @property
    def dna_length(self) -> int:
        return self.sequence.dna_length

    @property
    def sequence(self) -> Seq:
        return self._sequence

    @sequence.setter
    def sequence(self, seq: Seq | str | None) -> None:  # type: ignore
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

    def to_dict(self) -> dict[str, Any]:
        d = super().to_dict()
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


@attrs.define(init=False)
class DXGlue(Glue):
    etype: Literal["TD", "DT"] = attrs.field(default="TD")
    fullseq: Seq = attrs.field(default=Seq(""))

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
        super(DXGlue, self).__init__(name, note, use, abstractstrength)
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
        self.fullseq = trial_fseq

    @property
    def fseq(self) -> Optional[str]:
        if self.fullseq is not None:
            return self.fullseq.seq_str
        else:
            return None

    @fseq.setter
    def fseq(self, value: str):
        if self.fullseq is not None and len(value) != len(self.fullseq):
            warn("Changing end length")
        self.fullseq = Seq(value)

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
        if not self.fullseq:
            return None
        if self.etype == "TD":
            return self.fullseq.revcomp.base_str[1:]
        elif self.etype == "DT":
            return self.fullseq.revcomp.base_str[:-1]

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
            if out.fullseq:
                out.fullseq = out.fullseq.merge(other.fullseq)
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


class GlueList(Generic[GlueA], UpdateListD[GlueA]):
    def merge_complements(self) -> None:
        newitems: dict[str, GlueA] = {}
        for v in self:
            c = v.complement
            kc = c.ident()
            if kc in self.data:
                self.data[kc] = cast(GlueA, self.data[kc].merge(c))
            else:
                newitems[kc] = c
        self.data.update(newitems)

    def merge_glue(self, g: GlueA) -> GlueA:
        if g.ident() in self.data:
            g = cast(GlueA, self.data[g.ident()].merge(g))
        c = g.complement
        if c.ident() in self.data:
            g = cast(GlueA, self.data[c.ident()].complement.merge(g))
        return g

    def merge_glue_and_update_list(self, g: GlueA) -> GlueA:
        kg = g.ident()
        if kg in self.data:
            g = cast(
                GlueA, self.data[kg].merge(g)
            )  # Glue merges can only constrain type
            self.data[kg] = g
        c = g.complement
        kc = c.ident()
        if kc in self.data:
            c = cast(GlueA, self.data[kc].merge(c))
            self.data[kc] = c
            g = c.complement
            if kg in self.data:
                self.data[kg] = g
        return g

    def to_endarrays(self) -> Any:
        raise NotImplementedError
