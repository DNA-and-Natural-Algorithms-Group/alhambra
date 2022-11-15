"""
Tiles, seeds, glues, and lattices for the 'flatish' tile motif.
"""
from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    ClassVar,
    Literal,
    Optional,
    Sequence,
    Type,
    TypeVar,
    Union,
    cast,
)
import attrs

import numpy as np
import scadnano

from alhambra.grid import (
    AbstractLattice,
    LatticeSupportingScadnano,
    AbstractLatticeSupportingScadnano,
    ScadnanoLattice,
    lattice_factory,
)
from alhambra.seq import Seq
from alhambra.tilesets import XgrowGlueOpts

from .glues import Glue, GlueList
from .seeds import Seed, _convert_adapts, seed_factory, DiagonalSESeed
from .tiles import (
    BaseSSTSingleWithExtensions,
    BaseSSTile,
    BaseSSTSingle,
    Color,
    HDupleTile,
    SSGlue,
    Tile,
    VDupleTile,
    tile_factory,
    TileList,
    _scadnano_color,
)

if TYPE_CHECKING:
    import xgrow.tileset as xgt

__all__ = [
    "FlatishHSeed9",
    "FlatishVSeed9",
    "FlatishHDupleTile9_E",
    "FlatishHDupleTile10_E",
    "FlatishVDupleTile9_E2",
    "FlatishVDupleTile10_E2",
    "FlatishSingleTile9",
    "FlatishSingleTile10",
]


def _add_domain_from_glue(
    s: scadnano.StrandBuilder[Any, Any], g: SSGlue, d: Literal[1, -1]
) -> scadnano.StrandBuilder:
    s.move(g.dna_length * d)
    s.with_domain_name(g.ident())
    if g.sequence.is_definite:
        s.with_domain_sequence(g.sequence.base_str)
    return s


def _add_loopout_from_glue(
    s: scadnano.StrandBuilder[Any, Any], g: SSGlue, d: Literal[1, -1]
) -> scadnano.StrandBuilder:
    s.loopout(s.current_helix + d, g.dna_length)
    s.with_domain_name(g.ident())
    return s


_STANDARD_LOOP = SSGlue(name="flatish_loop8", sequence=8 * "T")


T = TypeVar("T")


def _reorder(seq: Sequence[T], ord: Sequence[int]) -> list[T]:
    return [seq[i] for i in ord]


class FlatishSingleTile9(BaseSSTSingle):
    "Flatish single tile, with domains (5'→3') of 12, 9, 11, and 10 nt.  North edge is 9nt."
    _base_domains: ClassVar[list[SSGlue]] = [SSGlue(length=x) for x in [12, 9, 11, 10]]
    _scadnano_offsets = ((-1, -12), (-1, 9), (1, 11), (1, -10))
    _scadnano_5p_offset = (0, 21)


class FlatishSingleTile10(BaseSSTSingle):
    "Flatish single tile, with domains (5'→3') of 11, 10, 12, and 9 nt. North edge is 10nt."
    _base_domains: ClassVar[list[SSGlue]] = [SSGlue(length=x) for x in [11, 10, 12, 9]]
    _scadnano_offsets = ((-1, -11), (-1, 10), (1, 12), (1, -9))
    _scadnano_5p_offset = (0, 21)


class FlatishSingleTile9WithExtensions(FlatishSingleTile9, BaseSSTSingleWithExtensions):
    ...


class FlatishSingleTile10WithExtensions(
    FlatishSingleTile10, BaseSSTSingleWithExtensions
):
    ...


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
    _scadnano_5p_offset = (1, 33)

    @property
    def domains(self) -> list[SSGlue]:
        e = self.edges
        return [e[2], _STANDARD_LOOP.copy(), e[1], e[0], e[5], _STANDARD_LOOP.copy(), e[4], e[3]]  # type: ignore

    def to_scadnano(
        self, design: scadnano.Design, helix: int, offset: int
    ) -> scadnano.Strand:
        s = design.draw_strand(
            helix + self._scadnano_5p_offset[0], offset + self._scadnano_5p_offset[1]
        )
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

        if self.color is not None:
            s.with_color(_scadnano_color(self.color))

        s.with_sequence(self.sequence.base_str)
        return s.strand

    def _input_neighborhood_domains(
        self,
    ) -> Sequence[tuple[Sequence[str], Sequence[str]]]:
        return []


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
    _scadnano_5p_offset = (1, 32)

    @property
    def domains(self):
        e = self.edges
        return [
            e[2],
            _STANDARD_LOOP.copy(),
            e[1],
            e[0],
            e[5],
            _STANDARD_LOOP.copy(),
            e[4],
            e[3],
        ]

    def to_scadnano(
        self, design: scadnano.Design, helix: int, offset: int
    ) -> scadnano.Strand:
        s = design.draw_strand(
            helix + self._scadnano_5p_offset[0], offset + self._scadnano_5p_offset[1]
        )
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

        if self.color is not None:
            s.with_color(_scadnano_color(self.color))

        s.with_sequence(self.sequence.base_str)
        return s.strand

    def _input_neighborhood_domains(
        self,
    ) -> Sequence[tuple[Sequence[str], Sequence[str]]]:
        return []


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
    _scadnano_5p_offset = (-1, 30)

    @property
    def domains(self):
        e = self.edges
        return [
            e[2],
            e[1],
            _STANDARD_LOOP.copy(),
            e[0],
            e[5],
            e[4],
            _STANDARD_LOOP.copy(),
            e[3],
        ]

    def to_scadnano(
        self, design: scadnano.Design, helix: int, offset: int
    ) -> scadnano.Strand:
        s = design.draw_strand(
            helix + self._scadnano_5p_offset[0], offset + self._scadnano_5p_offset[1]
        )
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

        if self.color is not None:
            s.with_color(_scadnano_color(self.color))

        s.with_sequence(self.sequence.base_str)
        return s.strand

    def _input_neighborhood_domains(
        self,
    ) -> Sequence[tuple[Sequence[str], Sequence[str]]]:
        return []


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
    _scadnano_5p_offset = (-1, 31)

    @property
    def domains(self):
        e = self.edges
        return [
            e[2],
            e[1],
            _STANDARD_LOOP.copy(),
            e[0],
            e[5],
            e[4],
            _STANDARD_LOOP.copy(),
            e[3],
        ]

    def to_scadnano(
        self, design: scadnano.Design, helix: int, offset: int
    ) -> scadnano.Strand:
        s = design.draw_strand(
            helix + self._scadnano_5p_offset[0], offset + self._scadnano_5p_offset[1]
        )
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

        if self.color is not None:
            s.with_color(_scadnano_color(self.color))

        s.with_sequence(self.sequence.base_str)
        return s.strand

    def _input_neighborhood_domains(
        self,
    ) -> Sequence[tuple[Sequence[str], Sequence[str]]]:
        return []


for ttype in [
    FlatishHDupleTile10_E,
    FlatishHDupleTile9_E,
    FlatishVDupleTile10_E2,
    FlatishVDupleTile9_E2,
    FlatishSingleTile10,
    FlatishSingleTile9,
    FlatishSingleTile10WithExtensions,
    FlatishSingleTile9WithExtensions,
]:
    tile_factory.register(cast(Type[Tile], ttype))


T_FHS9 = TypeVar("T_FHS9", bound="FlatishHSeed9")


class FlatishHSeed9(Seed):
    """Flatish origami seed."""

    adapter_tiles: list[tuple[Glue | str, FlatishSingleTile9]]

    def __init__(
        self, adapter_tiles: Sequence[tuple[Glue | str, FlatishSingleTile9]] = tuple()
    ):
        self.adapter_tiles = list(adapter_tiles)

    def to_dict(self, glues_as_refs: bool = False) -> dict:
        d: dict[str, Any] = {}
        d["adapter_tiles"] = [
            [str(g.name if isinstance(g, Glue) else g), t.to_dict()]
            for g, t in self.adapter_tiles  # FIXME
        ]
        d["type"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls: Type[T_FHS9], d: dict) -> T_FHS9:
        dat: tuple[str, dict] = d["adapter_tiles"]
        return cls(
            [(Glue(g), cast(FlatishSingleTile9, Tile.from_dict(t))) for g, t in dat]
        )  # FIXME: should check tiles

    def to_xgrow(
        self,
        gluenamemap: Callable[[str], str] = lambda x: x,
        offset: tuple[int, int] | None = None,
        xgtiles: Optional[Sequence[xgt.Tile]] = None,
    ) -> tuple[list[xgt.Tile], list[xgt.Bond], xgt.InitState]:
        import xgrow.tileset as xgt

        if offset is None:
            offset = (0, 0)

        xgtiles = []
        locs: list[tuple[int, int, str]] = []
        bonds = [xgt.Bond("seed", 100)]

        xgtiles.append(
            xgt.Tile([0, "seed", "seed", "seed"], "seed", stoic=0, color="white")
        )
        x = 2 + offset[0]
        ybase = 1 + offset[1]
        for y_offset in range(0, 2 * len(self.adapter_tiles)):
            if y_offset % 2:
                locs.append((x, ybase + y_offset, "seed"))
            else:
                adapt = y_offset // 2
                aname = f"adapterNW_{adapt}"
                aglue = self.adapter_tiles[adapt][0]
                if isinstance(aglue, Glue):
                    aglue = gluenamemap(aglue.ident())
                atile = xgt.Tile(
                    [0, "seed", aglue, "seed"], aname, stoic=0, color="green"
                )
                xgtiles.append(atile)
                locs.append((x, ybase + y_offset, aname))

        x = 3 + offset[0]
        ybase = 2 + offset[1]
        for adapt, (_, tile) in enumerate(self.adapter_tiles):
            if tile.name:
                aname = "adapterSE_" + tile.name
            else:
                aname = f"adapterSE_{adapt}"
            edges = ["seed"] + [gluenamemap(e.ident()) for e in tile._edges[1:]]
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


class FlatishVSeed9(Seed):
    """Flatish origami seed (vertical)."""

    adapter_tiles: list[tuple[Glue | str, FlatishSingleTile9]]

    def __init__(
        self, adapter_tiles: Sequence[tuple[Glue | str, FlatishSingleTile9]] = tuple()
    ):
        self.adapter_tiles = list(adapter_tiles)

    def to_dict(self, glues_as_refs: bool = False) -> dict:
        d: dict[str, Any] = {}
        d["adapter_tiles"] = [
            [str(g), t.to_dict()] for g, t in self.adapter_tiles  # FIXME
        ]
        d["type"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d: dict) -> FlatishVSeed9:
        return cls(
            [
                (Glue(g), cast(FlatishSingleTile9, Tile.from_dict(t)))
                for g, t in d["adapter_tiles"]
            ]
        )

    def to_xgrow(
        self,
        gluenamemap: Callable[[str], str] = lambda x: x,
        offset: tuple[int, int] | None = None,
        xgtiles: Optional[Sequence[xgt.Tile]] = None,
    ) -> tuple[list[xgt.Tile], list[xgt.Bond], xgt.InitState]:
        import xgrow.tileset as xgt

        if offset is None:
            offset = (0, 0)

        xgtiles = []
        locs: list[tuple[int, int, str]] = []
        bonds = [xgt.Bond("seed", 100)]

        xgtiles.append(
            xgt.Tile(["seed", "seed", "seed", 0], "seed", stoic=0, color="white")
        )
        ybase = 2 + offset[1]
        xbase = 1 + offset[0]
        for x_offset in range(0, 2 * len(self.adapter_tiles)):
            if x_offset % 2:
                locs.append((xbase + x_offset, ybase, "seed"))
            else:
                adapt = x_offset // 2
                aname = f"adapterNW_{adapt}"
                aglue = self.adapter_tiles[adapt][0]
                if isinstance(aglue, Glue):
                    aglue = gluenamemap(aglue.ident())
                atile = xgt.Tile(
                    ["seed", aglue, "seed", 0], aname, stoic=0, color="green"
                )
                xgtiles.append(atile)
                locs.append((xbase + x_offset, ybase, aname))

        ybase = 3 + offset[1]
        xbase = 2 + offset[0]
        for adapt, (_, tile) in enumerate(self.adapter_tiles):
            if tile.name:
                aname = "adapterSE_" + tile.name
            else:
                aname = f"adapterSE_{adapt}"
            edges = [gluenamemap(e.ident()) for e in tile._edges[:-1]] + ["seed"]
            xgtiles.append(
                xgt.Tile(
                    cast(list[Union[str, int]], edges),
                    name=aname,
                    stoic=0,
                    color="green",
                )
            )
            locs.append((xbase + 2 * adapt, ybase, aname))

        return xgtiles, bonds, xgt.InitState(locs)


seed_factory.register(FlatishVSeed9)


@attrs.define()
class FlatishNECornerSeed(Seed):
    hadapts: Sequence[tuple[Glue | str, FlatishSingleTile9]]
    corner: Sequence[FlatishSingleTile9]
    vadapts: Sequence[tuple[FlatishSingleTile9, Glue | str]]
    sim_offset: tuple[int, int] = (0, 0)

    def to_dict(self, glues_as_refs: bool = False) -> dict:
        d: dict[str, Any] = {}
        d["hadapts"] = [[str(g), t.to_dict()] for g, t in self.hadapts]  # FIXME
        d["corner"] = [t.to_dict() for t in self.corner]
        d["vadapts"] = [[t.to_dict(), str(g)] for t, g in self.vadapts]  # FIXME
        if self.sim_offset != (0, 0):
            d["sim_offset"] = self.sim_offset
        d["type"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d: dict) -> FlatishNECornerSeed:
        args = {
            "hadapts": [
                (Glue(g), cast(FlatishSingleTile9, Tile.from_dict(t)))
                for g, t in d["hadapts"]
            ],
            "corner": [
                cast(FlatishSingleTile9, Tile.from_dict(t)) for t in d["corner"]
            ],
            "vadapts": [
                (cast(FlatishSingleTile9, Tile.from_dict(t)), Glue(g))
                for t, g in d["vadapts"]
            ],
        }
        if "sim_offset" in d:
            args["sim_offset"] = d["sim_offset"]
        return cls(**args)  # type: ignore

    def to_xgrow(
        self,
        gluenamemap: Callable[[str], str] = lambda x: x,
        offset: tuple[int, int] | None = None,
        xgtiles: Optional[Sequence[xgt.Tile]] = None,
    ) -> tuple[list[xgt.Tile], list[xgt.Bond], xgt.InitState]:
        import xgrow.tileset as xgt

        if offset is None:
            offset = self.sim_offset

        xgtiles = []
        locs: list[tuple[int, int, str]] = []
        bonds = [xgt.Bond("seed", 100)]

        xgtiles.append(
            xgt.Tile(["seed", "seed", "seed", "seed"], "seed", stoic=0, color="white")
        )

        # Horizontal adapters
        x = 2 + offset[0]
        ybase = 1 + offset[1]
        for y_offset in range(0, 2 * len(self.hadapts)):
            if y_offset % 2:
                locs.append((x, ybase + y_offset, "seed"))
            else:
                adapt = y_offset // 2
                aname = f"adapterNW_{adapt}"
                aglue = self.hadapts[adapt][0]
                if isinstance(aglue, Glue):
                    aglue = gluenamemap(aglue.ident())
                atile = xgt.Tile(
                    [0, "seed", aglue, "seed"], aname, stoic=0, color="green"
                )
                xgtiles.append(atile)
                locs.append((x, ybase + y_offset, aname))

        x = 3 + offset[0]
        ybase = 2 + offset[1]
        for adapt, (_, tile) in enumerate(self.hadapts):
            if tile.name:
                aname = "adapterSE_" + tile.name
            else:
                aname = f"adapterSE_{adapt}"
            edges = ["seed"] + [gluenamemap(e.ident()) for e in tile._edges[1:]]

            # On the last tile, the E/1 glue is also seed
            if adapt == len(self.hadapts) - 1:
                edges[1] = "seed"

            xgtiles.append(
                xgt.Tile(
                    cast(list[Union[str, int]], edges),
                    name=aname,
                    stoic=0,
                    color="green",
                )
            )
            locs.append((x, ybase + 2 * adapt, aname))

        # Corner.  Our start x,y is +1,+1 from the last adapter tile
        xbase = locs[-1][0] + 1
        ybase = locs[-1][1] + 1

        for cn, tile in enumerate(self.corner):
            if tile.name:
                aname = "adaptC_SW_" + tile.name
            else:
                aname = f"adaptC_SW_{cn}"
            edges = [gluenamemap(e.ident()) for e in tile._edges]
            edges[0] = "seed"
            edges[1] = "seed"
            xgtiles.append(
                xgt.Tile(
                    cast(list[Union[str, int]], edges),
                    name=aname,
                    stoic=0,
                    color="green",
                )
            )
            locs.append((xbase + cn - 1, ybase + cn, "seed"))
            locs.append((xbase + cn, ybase + cn, aname))

        # Vertical adapters.  Start tile is at +1, +1 from last corner adapter
        xbase = locs[-1][0] + 1
        ybase = locs[-1][1] + 1

        for adapt_n, (tile, aglue) in enumerate(self.vadapts):
            # NW mimics tile
            if tile.name:
                aname = "adapterNW_" + tile.name
            else:
                aname = f"adapterNW_V_{adapt_n}"

            edges = [gluenamemap(e.ident()) for e in tile._edges]

            edges[1] = "seed"

            # On the first tile, the N/0 glue is also seed
            if adapt_n == 0:
                edges[0] = "seed"

                # We also add a seed tile at x-1
                locs.append((xbase - 1, ybase, "seed"))

            xgtiles.append(
                xgt.Tile(
                    cast(list[Union[str, int]], edges),
                    name=aname,
                    stoic=0,
                    color="green",
                )
            )
            locs.append((xbase + 2 * adapt_n, ybase, aname))

            locs.append((xbase + 2 * adapt_n, ybase + 1, "seed"))

            if isinstance(aglue, Glue):
                aglue = gluenamemap(aglue.ident())

            xgtiles.append(
                xgt.Tile(
                    ["seed", 0, "seed", aglue],
                    name=f"adapterSE_V_{adapt_n}",
                    stoic=0,
                    color="green",
                )
            )
            locs.append((xbase + 2 * adapt_n + 1, ybase + 1, f"adapterSE_V_{adapt_n}"))

        return xgtiles, bonds, xgt.InitState(locs)


seed_factory.register(FlatishNECornerSeed)


class FlatishLattice(AbstractLatticeSupportingScadnano):
    # FIXME: seed should be flatish
    """A lattice of flatish tiles.  Position 0,0 is a 9nt-north tile."""

    def to_scadnano_lattice(self) -> ScadnanoLattice:
        sclat = ScadnanoLattice()
        for ix, t in np.ndenumerate(self.grid):
            if not t:
                continue
            scpos = flatgrid_hofromxy(ix[0], ix[1], self.grid.shape[1], 0)
            sclat[scpos] = t

        if (
            self.seed is not None
            and hasattr(self.seed, "to_scadnano")
            and hasattr(self.seed, "ho_from_seed_offset")
        ):
            sclat.seed = self.seed
            sclat.seed_position = self.seed.ho_from_seed_offset(
                self.seed_offset, self.grid.shape[1]
            )

        return sclat


lattice_factory.register(FlatishLattice)


class FlatishDiagonalSESeed10(DiagonalSESeed):
    adapters: Sequence[tuple[SSGlue, SSGlue]] = attrs.field(
        converter=_convert_adapts, on_setattr=attrs.setters.convert
    )
    _lattice: ClassVar[Type[AbstractLatticeSupportingScadnano]] = FlatishLattice

    def _calculate_valid_offset(self, target: tuple[int, int]) -> tuple[int, int]:
        """
        Given a target seed offset, return an offset that will work with a FlatishLattice.
        """

        if ((target[0] + target[1] + len(self.adapters)) % 2) == 0:
            return target
        else:
            return (target[0], target[1] + 1)

    def ho_from_seed_offset(
        self, seed_offset: tuple[int, int], gridysize: int
    ) -> tuple[int, int]:
        seed_offset = self._calculate_valid_offset(seed_offset)

        return flatgrid_hofromxy(
            seed_offset[0], seed_offset[1] + len(self.adapters) - 1, gridysize, 0
        )

    def to_xgrow(
        self,
        gluenamemap: Callable[[str], str] = lambda x: x,
        offset: tuple[int, int] | None = None,
        xgtiles: Optional[Sequence[xgt.Tile]] = None,
    ) -> tuple[list[xgt.Tile], list[xgt.Bond], xgt.InitState]:

        if offset is None:
            offset = (0, 0)

        offset = self._calculate_valid_offset(offset)

        return super().to_xgrow(gluenamemap, offset)

    def update_details(
        self, glues: GlueList[Glue], tiles: TileList[Tile] | None = None
    ) -> None:
        self.adapters = [
            (glues.merge_glue(g1), glues.merge_glue(g2)) for g1, g2 in self.adapters  # type: ignore
        ]

    def to_scadnano(
        self, design: scadnano.Design, helix: int, offset: int
    ) -> list[scadnano.Strand]:
        """
        For this seed, the `helix` and `offset` refer to the position of the NE-most (ie, most Northerly) "fake" tile that the seed represents.  This corresponds to the tile position N of the NE-most tile that attaches to the seed.
        """

        # apse strands.  The first one here starts 22nt east of the starting offset.
        apse = []
        bpse = []

        for i, gs in enumerate(self.adapters):
            s = design.draw_strand(helix + 2 * i, offset + 21 + 2 * i)
            _add_domain_from_glue(s, gs[0], -1)
            s.move(-31)
            s.cross(s.current_helix + 1)
            s.move(33)
            _add_domain_from_glue(s, gs[1], 1)
            apse.append(s)

        alen = len(self.adapters)

        for i in range(0, alen):
            s = design.draw_strand(helix + 1 + 2 * i, offset - 37 + 2 * i)
            s.move(16)
            s.cross(s.current_helix - 1)
            s.move(-16)
            bpse.append(s)

        scaffold = design.draw_strand(
            helix + 1 + 2 * alen - 2, offset + 12 + 2 * alen - 2
        )
        for _ in range(alen):
            scaffold.move(-49)
            scaffold.cross(scaffold.current_helix - 1)
            scaffold.move(47)
            scaffold.cross(scaffold.current_helix - 1)

        return apse + bpse + [scaffold]


class FlatishDiagonalSESeed9(DiagonalSESeed):
    adapters: Sequence[tuple[SSGlue, SSGlue]] = attrs.field(
        converter=_convert_adapts, on_setattr=attrs.setters.convert
    )
    _lattice = FlatishLattice

    def _calculate_valid_offset(self, target: tuple[int, int]) -> tuple[int, int]:
        """
        Given a target seed offset, return an offset that will work with a FlatishLattice.
        """

        if ((target[0] + target[1] + len(self.adapters)) % 2) == 1:
            return target
        else:
            return (target[0], target[1] + 1)

    def ho_from_seed_offset(self, seed_offset: tuple[int, int], gridysize):
        seed_offset = self._calculate_valid_offset(seed_offset)

        return flatgrid_hofromxy(
            seed_offset[0], seed_offset[1] + len(self.adapters) - 1, gridysize, 0
        )

    def to_xgrow(
        self,
        gluenamemap: Callable[[str], str] = lambda x: x,
        offset: tuple[int, int] | None = None,
        xgtiles: Optional[Sequence[xgt.Tile]] = None,
    ) -> tuple[list[xgt.Tile], list[xgt.Bond], xgt.InitState]:

        if offset is None:
            offset = (0, 0)

        offset = self._calculate_valid_offset(offset)

        return super().to_xgrow(gluenamemap, offset)

    def update_details(
        self, glues: GlueList[Glue], tiles: TileList[Tile] | None = None
    ) -> None:
        self.adapters = [
            (glues.merge_glue(g1), glues.merge_glue(g2)) for g1, g2 in self.adapters  # type: ignore
        ]

    def to_scadnano(
        self, design: scadnano.Design, helix: int, offset: int
    ) -> list[scadnano.Strand]:
        """
        For this seed, the `helix` and `offset` refer to the position of the NE-most (ie, most Northerly) "fake" tile that the seed represents.  This corresponds to the tile position N of the NE-most tile that attaches to the seed.
        """

        # apse strands.  The first one here starts 22nt east of the starting offset.
        apse = []
        bpse = []

        for i, gs in enumerate(self.adapters):
            s = design.draw_strand(helix + 2 * i, offset + 21 + 2 * i)
            _add_domain_from_glue(s, gs[0], -1)
            s.move(-31).with_name(f"a_{i}")
            s.cross(s.current_helix + 1)
            s.move(33)
            _add_domain_from_glue(s, gs[1], 1)
            apse.append(s)

        alen = len(self.adapters)

        for i in range(0, alen):
            s = design.draw_strand(helix + 1 + 2 * i, offset - 38 + 2 * i)
            s.move(16).with_name(f"b_{i}")
            s.cross(s.current_helix - 1)
            s.move(-16)
            bpse.append(s)

        scaffold = design.draw_strand(
            helix + 1 + 2 * alen - 2, offset + 12 + 2 * alen - 3
        )
        for _ in range(alen):
            scaffold.move(-49).with_name("scaffold")
            scaffold.cross(scaffold.current_helix - 1)
            scaffold.move(47)
            scaffold.cross(scaffold.current_helix - 1)

        return apse + bpse + [scaffold]


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


seed_factory.register(FlatishDiagonalSESeed10)
seed_factory.register(FlatishDiagonalSESeed9)
