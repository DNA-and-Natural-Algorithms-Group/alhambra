from __future__ import annotations
from .glue import Glue
from .newseed import Seed, seed_factory
from typing import Any, ClassVar, Literal, Sequence, Type, TypeVar, Union, cast
import scadnano
from .newtile import (
    BaseSSTSingle,
    BaseSSTile,
    HDupleTile,
    SSGlue,
    Tile,
    VDupleTile,
    tile_factory,
)
import xgrow.tileset as xgt

__all__ = [
    "FlatishHSeed9",
    "FlatishHDupleTile9_E",
    "FlatishHDupleTile10_E",
    "FlatishVDupleTile9_E2",
    "FlatishVDupleTile10_E2",
    "FlatishSingleTile9",
    "FlatishSingleTile10",
]


def _add_domain_from_glue(
    s: scadnano.StrandBuilder[Any, Any], g: SSGlue, d: Literal[1, -1]
):
    s.move(g.dna_length * d)
    if g.name is not None:
        s.with_domain_name(g.name)
    return s


def _add_loopout_from_glue(
    s: scadnano.StrandBuilder[Any, Any], g: SSGlue, d: Literal[1, -1]
):
    s.loopout(s.current_helix + d, g.dna_length)
    if g.name is not None:
        s.with_domain_name(g.name)
    return s


_STANDARD_LOOP = SSGlue("loop", 8 * "T")


T = TypeVar("T")


def _reorder(seq: Sequence[T], ord: Sequence[int]) -> list[T]:
    return [seq[i] for i in ord]


class FlatishSingleTile9(BaseSSTSingle):
    _base_domains: ClassVar[list[SSGlue]] = [SSGlue(length=x) for x in [12, 9, 11, 10]]
    _scadnano_offsets = ((-1, -12), (-1, 9), (1, 11), (1, -10))


class FlatishSingleTile10(BaseSSTSingle):
    _base_domains: ClassVar[list[SSGlue]] = [SSGlue(length=x) for x in [11, 10, 12, 9]]
    _scadnano_offsets = ((-1, -11), (-1, 10), (1, 12), (1, -9))


class FlatishVDupleTile10_E2(VDupleTile, BaseSSTile):
    _base_domains: ClassVar[list[SSGlue]] = [
        SSGlue(length=12),
        _STANDARD_LOOP,
        SSGlue(length=11),
        SSGlue(length=10),
        SSGlue(length=12),
        _STANDARD_LOOP,
        SSGlue(length=11),
        SSGlue(length=10),
    ]
    _base_edges = _reorder(_base_domains, [3, 2, 0, 7, 6, 4])
    _scadnano_offsets = ((-1, -11), (-1, 10), (0, 21), (2, 23), (2, 1), (1, -9))

    @property
    def domains(self) -> list[SSGlue]:
        e = self.edges
        return [e[2], _STANDARD_LOOP, e[1], e[0], e[5], _STANDARD_LOOP, e[4], e[3]]  # type: ignore

    def to_scadnano(
        self, design: scadnano.Design, helix: int, offset: int
    ) -> scadnano.Strand:
        s = design.strand(helix, offset + 33)
        domiter = iter(self.domains)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_loopout_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), -1)
        s.cross(s.current_helix + 1)
        _add_domain_from_glue(s, next(domiter), 1)
        _add_loopout_from_glue(s, next(domiter), 1)
        _add_domain_from_glue(s, next(domiter), 1)
        _add_domain_from_glue(s, next(domiter), 1)

        if self.name is not None:
            s.with_name(self.name)

        s.with_sequence(self.sequence.base_str)
        return s.strand


class FlatishVDupleTile9_E2(VDupleTile, BaseSSTile):
    _base_domains: ClassVar[list[SSGlue]] = [
        SSGlue(length=11),
        _STANDARD_LOOP,
        SSGlue(length=12),
        SSGlue(length=9),
        SSGlue(length=11),
        _STANDARD_LOOP,
        SSGlue(length=12),
        SSGlue(length=9),
    ]
    _base_edges = _reorder(_base_domains, [3, 2, 0, 7, 6, 4])
    _scadnano_offsets = ((-1, -12), (-1, 9), (0, 21), (2, 23), (2, 2), (1, -10))

    @property
    def domains(self):
        e = self.edges
        return [e[2], _STANDARD_LOOP, e[1], e[0], e[5], _STANDARD_LOOP, e[4], e[3]]

    def to_scadnano(
        self, design: scadnano.Design, helix: int, offset: int
    ) -> scadnano.Strand:
        s = design.strand(helix, offset + 32)
        domiter = iter(self.domains)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_loopout_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), -1)
        s.cross(s.current_helix + 1)
        _add_domain_from_glue(s, next(domiter), 1)
        _add_loopout_from_glue(s, next(domiter), 1)
        _add_domain_from_glue(s, next(domiter), 1)
        _add_domain_from_glue(s, next(domiter), 1)

        if self.name is not None:
            s.with_name(self.name)

        s.with_sequence(self.sequence.base_str)
        return s.strand


class FlatishHDupleTile9_E(HDupleTile, BaseSSTile):
    _base_domains: list[SSGlue] = [
        SSGlue(length=11),
        SSGlue(length=10),
        _STANDARD_LOOP,
        SSGlue(length=9),
        SSGlue(length=11),
        SSGlue(length=10),
        _STANDARD_LOOP,
        SSGlue(length=9),
    ]
    _base_edges = _reorder(_base_domains, [3, 1, 0, 7, 5, 4])

    @property
    def domains(self):
        e = self.edges
        return [e[2], e[1], _STANDARD_LOOP, e[0], e[5], e[4], _STANDARD_LOOP, e[3]]

    def to_scadnano(
        self, design: scadnano.Design, helix: int, offset: int
    ) -> scadnano.Strand:
        s = design.strand(helix, offset + 30)
        domiter = iter(self.domains)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_loopout_from_glue(s, next(domiter), 1)
        _add_domain_from_glue(s, next(domiter), -1)
        s.cross(s.current_helix + 1)
        _add_domain_from_glue(s, next(domiter), 1)
        _add_domain_from_glue(s, next(domiter), 1)
        _add_loopout_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), 1)

        if self.name is not None:
            s.with_name(self.name)

        s.with_sequence(self.sequence.base_str)
        return s.strand


class FlatishHDupleTile10_E(HDupleTile, BaseSSTile):
    _base_domains: ClassVar[list[SSGlue]] = [
        SSGlue(length=12),
        SSGlue(length=9),
        _STANDARD_LOOP,
        SSGlue(length=10),
        SSGlue(length=12),
        SSGlue(length=9),
        _STANDARD_LOOP,
        SSGlue(length=10),
    ]
    _base_edges = _reorder(_base_domains, [3, 1, 0, 7, 5, 4])

    @property
    def domains(self):
        e = self.edges
        return [e[2], e[1], _STANDARD_LOOP, e[0], e[5], e[4], _STANDARD_LOOP, e[3]]

    def to_scadnano(
        self, design: scadnano.Design, helix: int, offset: int
    ) -> scadnano.Strand:
        s = design.strand(helix, offset + 31)
        domiter = iter(self.domains)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), -1)
        _add_loopout_from_glue(s, next(domiter), 1)
        _add_domain_from_glue(s, next(domiter), -1)
        s.cross(s.current_helix + 1)
        _add_domain_from_glue(s, next(domiter), 1)
        _add_domain_from_glue(s, next(domiter), 1)
        _add_loopout_from_glue(s, next(domiter), -1)
        _add_domain_from_glue(s, next(domiter), 1)

        if self.name is not None:
            s.with_name(self.name)

        s.with_sequence(self.sequence.base_str)
        return s.strand


for ttype in [
    FlatishHDupleTile10_E,
    FlatishHDupleTile9_E,
    FlatishVDupleTile10_E2,
    FlatishVDupleTile9_E2,
    FlatishSingleTile10,
    FlatishSingleTile9,
]:
    tile_factory.register(ttype)


T_FHS9 = TypeVar("T_FHS9", bound="FlatishHSeed9")


class FlatishHSeed9(Seed):
    adapter_tiles: list[tuple[Glue | str, FlatishSingleTile9]]

    def __init__(self, adapter_tiles=[]):
        self.adapter_tiles = list(adapter_tiles)

    def to_dict(self, glues_as_refs=False) -> dict:
        d: dict[str, Any] = {}
        d["adapter_tiles"] = [
            [str(g), t.to_dict()] for g, t in self.adapter_tiles  # FIXME
        ]
        d["type"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls: Type[T_FHS9], d: dict) -> T_FHS9:
        return cls([(g, Tile.from_dict(t)) for g, t in d["adapter_tiles"]])

    def to_xgrow(
        self, self_complementary_glues=False
    ) -> tuple[list[xgt.Tile], list[xgt.Bond], xgt.InitState]:

        xgtiles = []
        locs: list[tuple[int, int, str]] = []
        bonds = [xgt.Bond("seed", 10)]

        xgtiles.append(
            xgt.Tile([0, "seed", "seed", "seed"], "seed", stoic=0, color="white")
        )
        x = 1
        ybase = 1
        for y_offset in range(0, 2 * len(self.adapter_tiles)):
            if y_offset % 2:
                locs.append((x, ybase + y_offset, "seed"))
            else:
                adapt = y_offset // 2
                aname = f"adapterNW_{adapt}"
                aglue = self.adapter_tiles[adapt][0]
                if isinstance(aglue, Glue):
                    if self_complementary_glues:
                        aglue = aglue.basename()
                    else:
                        aglue = aglue.ident()
                atile = xgt.Tile(
                    [0, "seed", aglue, "seed"], aname, stoic=0, color="green"
                )
                xgtiles.append(atile)
                locs.append((x, ybase + y_offset, aname))

        x = 2
        ybase = 2
        for adapt, (south, tile) in enumerate(self.adapter_tiles):
            if tile.name:
                aname = "adapterSE_" + tile.name
            else:
                aname = f"adapterSE_{adapt}"
            if self_complementary_glues:
                edges = ["seed"] + [e.basename() for e in tile._edges[1:]]
            else:
                edges = ["seed"] + [e.ident() for e in tile._edges[1:]]
            xgtiles.append(
                xgt.Tile(
                    cast(list[Union[str, int]], edges),
                    name=aname,
                    stoic=0,
                    color="green",
                )
            )
            locs.append((x, ybase + 2 * adapt, aname))

        return xgtiles, bonds, xgt.InitState(locs)


seed_factory.register(FlatishHSeed9)


def flatgrid_hofromxy(
    x: int, y: int, start_helix: int, start_o: int, p: Literal[9, 10] = 9
) -> tuple[int, int]:
    if p == 9:
        pn = 0
    elif p == 10:
        pn = 1
    else:
        raise ValueError
    sx = (pn + y) % 2
    sy = (pn) % 2
    return (
        start_helix - y + x,
        start_o
        + 23 * (x // 2)
        + 19 * (y // 2)
        + (11 + sx) * (x % 2)
        + (9 + sy) * (y % 2),
    )
