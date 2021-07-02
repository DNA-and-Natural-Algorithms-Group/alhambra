from __future__ import annotations
from dataclasses import dataclass
from typing import Any, List, Literal, Optional, ClassVar, Sequence, Tuple, Type, TypeVar, Union
from enum import Enum
import scadnano

from .seq import Seq
from .glue import Glue, SSGlue

Color = str


class D(Enum):
    N = 0
    E = 1
    S = 2
    W = 3


Offset = Tuple[int, int]


@dataclass(frozen=True)
class EdgeLoc:
    direction: D
    position: Offset = (0, 0)


EL = EdgeLoc


class Tile:
    name: Optional[str]
    _edges: Optional[List[Glue]]
    color: Optional[Color]
    stoic: Optional[float]
    note: Optional[str]
    __slots__ = ('name', '_edges', 'color', 'stoic', 'note')

    def __init__(self, edges: Optional[List[Glue]],
                 name: Optional[str] = None,
                 color: Optional[Color] = None,
                 stoic: Optional[float] = None,
                 note: Optional[str] = None) -> None:
        self._edges = edges
        self.name = name
        self.color = color
        self.stoic = stoic
        self.note = note

    @property
    def edges(self):
        return self._edges

    @property
    def rotations(self) -> List[Tile]:
        raise NotImplementedError

    @property
    def xgname(self) -> str:
        if self.name:
            return self.name
        if self._edges:
            return f"{self.__class__.__name__}_"+"_".join(f"{g.xgname}"
                                                          for g in self._edges)
        else:
            raise ValueError

    def as_dict(self, glues_as_refs=False) -> dict[str, Any]:
        b = {k: v for k in ['name', 'edges', 'color', 'stoic', 'note']
             if (v := getattr(self, k)) is not None}
        if self.edges is not None:
            if not glues_as_refs:
                b['edges'] = [x.as_dict() for x in self.edges]
            else:
                # fixme: deal with None
                b['edges'] = [x.name for x in self.edges]
        return b

    @staticmethod
    def from_dict(d: dict[str, Any]) -> Tile:
        return tile_factory.from_dict(d)

    def to_xgrow_dict(self) -> dict[str, Any]:
        d = {}
        if self.color:
            d['color'] = self.color
        if self.stoic:
            d['stoic'] = self.stoic
        d['name'] = self.xgname
        if self._edges:
            d['edges'] = [g.xgname for g in self._edges]
        else:
            raise ValueError
        return d


class SingleTile(Tile):
    """A tile with N, E, S and W edges."""
    ...


class VDupleTile(Tile):
    def to_xgrow_dict(self) -> dict[str, Any]:
        d = super().to_xgrow_dict()
        d['double'] = 'V'
        return d


class HDupleTile(Tile):
    def to_xgrow_dict(self) -> dict[str, Any]:
        d = super().to_xgrow_dict()
        d['double'] = 'H'
        return d


class BaseSSTile(Tile):
    _edges: List[SSGlue]

    def as_dict(self, glues_as_refs=False) -> dict[str, Any]:
        d = super().as_dict(glues_as_refs=glues_as_refs)
        d['type'] = self.__class__.__name__
        d['sequence'] = self.sequence.seq_str
        return d

    @property
    def _base_domains(self) -> ClassVar[List[SSGlue]]:
        ...

    @property
    def _base_edges(self) -> ClassVar[List[SSGlue]]:
        ...

    @property
    def domains(self) -> List[SSGlue]:
        ...

    @property
    def _sequence_length(self) -> int:
        return sum(x.dna_length for x in self._base_domains)

    def __init__(self, edges: Optional[List[Union[str, Glue]]] = None,
                 name: Optional[str] = None,
                 color: Optional[Color] = None,
                 stoic: Optional[float] = None,
                 sequence: Optional[Seq] = None,
                 note: Optional[str] = None):
        super().__init__(None, name, color, stoic, note)

        if edges is None and sequence is None:
            raise ValueError
        if edges is not None:
            self._edges = [bd.merge(g)
                           for bd, g in zip(self._base_edges, edges)]
        else:
            self._edges = [bd.copy() for bd in self._base_edges]
        if sequence is not None:
            self.sequence | sequence
            self.sequence = sequence

    @property
    def sequence(self) -> Seq:
        return Seq("-".join(str(glue.sequence) for glue in self.domains))

    @sequence.setter
    def sequence(self, seq: Seq):
        seq = Seq(seq)
        if seq.dna_length != self._sequence_length:
            raise ValueError

        pos = 0
        base_str = seq.base_str

        # fixme: should we check whitespace?
        for base_domain, domain in zip(self._base_domains, self.domains):
            base_domain.sequence | base_str[pos: pos+base_domain.dna_length]
            domain.sequence = seq.base_str[pos: pos+base_domain.dna_length]
            pos += base_domain.dna_length

    @property
    def edges(self):
        return self._edges  # fixme: should not be mutable here

    @edges.setter
    def edges(self, edges: List[Glue]):
        self._edges = [bd.merge(g) for bd, g in zip(self._base_edges, edges)]

    def set_edge(self, i: int, glue: Glue):
        self._edges[i] = self._base_edges[i].merge(glue)

    def __repr__(self) -> str:
        s = []
        s.append(self.__class__.__name__)
        s.append("(")
        if self.name is not None:
            s.append(f"name={repr(self.name)}, ")
        s.append("edges=[" + ", ".join(repr(x) for x in self.edges) + "]")
        s.append(")")
        return "".join(s)

    def __str__(self) -> str:
        s = [f"<{self.__class__.__name__}: "]
        if self.name is not None:
            s.append(self.name + " | ")
        s.append(", ".join(g._shortdesc for g in self.edges))
        s.append(">")
        return "".join(s)


class BaseSSTSingle(SingleTile, BaseSSTile):
    @property
    def domains(self) -> List[SSGlue]:
        return [self.edges[i] for i in [1, 0, 3, 2]]

    @property
    def _base_edges(self) -> List[SSGlue]:
        return [self._base_domains[i] for i in [1, 0, 3, 2]]

    def to_scadnano(self, design: scadnano.Design,
                    helix: int, offset: int) -> scadnano.Strand:
        s = design.strand(helix, offset + 21)

        for e in self.domains[0:2]:
            s.move(-e.dna_length)
            if e.name is not None:
                s.with_domain_name(e.name)
        s.cross(s.current_helix+1)
        for e in self.domains[2:]:
            s.move(e.dna_length)
            if e.name is not None:
                s.with_domain_name(e.name)

        if self.name is not None:
            s.with_name(self.name)

        s.with_sequence(self.sequence.base_str)
        return s.strand


class FlatishSingleTile9(BaseSSTSingle):
    _base_domains: ClassVar[List[SSGlue]] = [
        SSGlue(length=x) for x in [12, 9, 11, 10]]
    _scadnano_offsets = ((-1, -12), (-1, 9), (1, 11), (1, -10))


class FlatishSingleTile10(BaseSSTSingle):
    _base_domains: ClassVar[List[SSGlue]] = [
        SSGlue(length=x) for x in [11, 10, 12, 9]]
    _scadnano_offsets = ((-1, -11), (-1, 10), (1, 12), (1, -9))

# @dataclass(init=False)
# class FlatishHDupleTile12_E2(HDupleTile, BaseFlatishTile):
#     _base_edges: ClassVar[List[FlatishGlue]] = [
#         FlatishGlue(length=x) for x in [12, 9, 11, 10]]
#     _fixed: ClassVar[List[Optional[str]]] = [
#         None, 8*'T', None, None, None, 8*'T', None, None]


# @dataclass(init=False)
# class FlatishHDupleTile11_E2(HDupleTile, BaseFlatishTile):
#     _base_edges: ClassVar[List[FlatishGlue]] = [
#         FlatishGlue(length=x) for x in [12, 9, 11, 10]]
#     _fixed: ClassVar[List[Optional[str]]] = [
#         None, 8*'T', None, None, None, 8*'T', None, None]

_STANDARD_LOOP = SSGlue("loop", 8*'T')


def _add_domain_from_glue(s: scadnano.StrandBuilder, g: SSGlue, d: Literal[1, -1]):
    s.move(g.dna_length*d)
    if g.name is not None:
        s.with_domain_name(g.name)
    return s


def _add_loopout_from_glue(s: scadnano.StrandBuilder, g: SSGlue, d: Literal[1, -1]):
    s.loopout(s.current_helix+d, g.dna_length)
    if g.name is not None:
        s.with_domain_name(g.name)
    return s


T = TypeVar('T')


def _reorder(seq: Sequence[T], ord: Sequence[int]) -> Sequence[T]:
    return [seq[i] for i in ord]


class FlatishVDupleTile10_E2(VDupleTile, BaseSSTile):
    _base_domains: ClassVar[List[SSGlue]] = [
        SSGlue(length=12),
        _STANDARD_LOOP,
        SSGlue(length=11),
        SSGlue(length=10),
        SSGlue(length=12),
        _STANDARD_LOOP,
        SSGlue(length=11),
        SSGlue(length=10)]
    _base_edges = _reorder(_base_domains, [3, 2, 0, 7, 6, 4])
    _scadnano_offsets = ((-1, -11), (-1, 10), (0, 21),
                         (2, 23), (2, 1), (1, -9))

    @property
    def domains(self):
        e = self.edges
        return [
            e[2], _STANDARD_LOOP, e[1], e[0], e[5], _STANDARD_LOOP, e[4], e[3]
        ]

    def to_scadnano(self, design: scadnano.Design,
                    helix: int, offset: int) -> scadnano.Strand:
        s = design.strand(helix, offset + 33)
        domiter = iter(self.domains)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_loopout_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), -1)
        s.cross(s.current_helix+1)
        _add_domain_from_glue(s, next(domiter), +1)
        _add_loopout_from_glue(s, next(domiter), +1)
        _add_domain_from_glue(s, next(domiter), +1)
        _add_domain_from_glue(s, next(domiter), +1)

        if self.name is not None:
            s.with_name(self.name)

        s.with_sequence(self.sequence.base_str)
        return s.strand


class FlatishVDupleTile9_E2(VDupleTile, BaseSSTile):
    _base_domains: ClassVar[List[SSGlue]] = [
        SSGlue(length=11),
        _STANDARD_LOOP,
        SSGlue(length=12),
        SSGlue(length=9),
        SSGlue(length=11),
        _STANDARD_LOOP,
        SSGlue(length=12),
        SSGlue(length=9)]
    _base_edges = _reorder(_base_domains, [3, 2, 0, 7, 6, 4])
    _scadnano_offsets = ((-1, -12), (-1, 9), (0, 21),
                         (2, 23), (2, 2), (1, -10))

    @property
    def domains(self):
        e = self.edges
        return [
            e[2], _STANDARD_LOOP, e[1], e[0], e[5], _STANDARD_LOOP, e[4], e[3]
        ]

    def to_scadnano(self, design: scadnano.Design,
                    helix: int, offset: int) -> scadnano.Strand:
        s = design.strand(helix, offset + 32)
        domiter = iter(self.domains)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_loopout_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), -1)
        s.cross(s.current_helix+1)
        _add_domain_from_glue(s, next(domiter), +1)
        _add_loopout_from_glue(s, next(domiter), +1)
        _add_domain_from_glue(s, next(domiter), +1)
        _add_domain_from_glue(s, next(domiter), +1)

        if self.name is not None:
            s.with_name(self.name)

        s.with_sequence(self.sequence.base_str)
        return s.strand


class FlatishHDupleTile9_E(HDupleTile, BaseSSTile):
    _base_domains: ClassVar[List[SSGlue]] = [
        SSGlue(length=11),
        SSGlue(length=10),
        _STANDARD_LOOP,
        SSGlue(length=9),
        SSGlue(length=11),
        SSGlue(length=10),
        _STANDARD_LOOP,
        SSGlue(length=9)]
    _base_edges = _reorder(_base_domains, [3, 1, 0, 7, 5, 4])

    @property
    def domains(self):
        e = self.edges
        return [
            e[2], e[1], _STANDARD_LOOP, e[0], e[5], e[4], _STANDARD_LOOP, e[3]
        ]

    def to_scadnano(self, design: scadnano.Design,
                    helix: int, offset: int) -> scadnano.Strand:
        s = design.strand(helix, offset + 30)
        domiter = iter(self.domains)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_loopout_from_glue(s, next(domiter), +1)
        _add_domain_from_glue(s, next(domiter), -1)
        s.cross(s.current_helix+1)
        _add_domain_from_glue(s, next(domiter), +1)
        _add_domain_from_glue(s, next(domiter), +1)
        _add_loopout_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), +1)

        if self.name is not None:
            s.with_name(self.name)

        s.with_sequence(self.sequence.base_str)
        return s.strand


class FlatishHDupleTile10_E(HDupleTile, BaseSSTile):
    _base_domains: ClassVar[List[SSGlue]] = [
        SSGlue(length=12),
        SSGlue(length=9),
        _STANDARD_LOOP,
        SSGlue(length=10),
        SSGlue(length=12),
        SSGlue(length=9),
        _STANDARD_LOOP,
        SSGlue(length=10)]
    _base_edges = _reorder(_base_domains, [3, 1, 0, 7, 5, 4])

    @property
    def domains(self):
        e = self.edges
        return [
            e[2], e[1], _STANDARD_LOOP, e[0], e[5], e[4], _STANDARD_LOOP, e[3]
        ]

    def to_scadnano(self, design: scadnano.Design,
                    helix: int, offset: int) -> scadnano.Strand:
        s = design.strand(helix, offset + 31)
        domiter = iter(self.domains)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_loopout_from_glue(s, next(domiter), +1)
        _add_domain_from_glue(s, next(domiter), -1)
        s.cross(s.current_helix+1)
        _add_domain_from_glue(s, next(domiter), +1)
        _add_domain_from_glue(s, next(domiter), +1)
        _add_loopout_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), +1)

        if self.name is not None:
            s.with_name(self.name)

        s.with_sequence(self.sequence.base_str)
        return s.strand


class TileFactory:
    types: dict[str, Type[Tile]]

    def __init__(self):
        self.types = {}

    def register(self, c: Type[Tile], n: str = None):
        self.types[n if n is not None else c.__name__] = c

    def from_dict(self, d: dict[str, Any]) -> Tile:
        if 'type' in d:
            c = self.types[d['type']]
            del(d['type'])
            return c(**d)
        elif 'structure' in d:
            c = self.types[d['structure']]
            del(d['structure'])
            return c(**d)
        else:
            return Tile(**d)


tile_factory = TileFactory()
for ttype in [Tile, FlatishSingleTile10, FlatishSingleTile9, FlatishHDupleTile10_E,
              FlatishHDupleTile9_E, FlatishVDupleTile10_E2, FlatishVDupleTile9_E2]:
    tile_factory.register(ttype)
