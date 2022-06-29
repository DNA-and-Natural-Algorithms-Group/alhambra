from __future__ import annotations

import copy
import logging
import warnings
from dataclasses import dataclass, field
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Collection,
    Iterable,
    Literal,
    Mapping,
    Optional,
    Type,
    TypeVar,
    Union,
    cast,
)
from typing_extensions import TypeAlias

import drawSvg_svgy as draw
import numpy as np


from numpy import isin

from alhambra.classes import Serializable
from alhambra.grid import (
    AbstractLattice,
    Lattice,
    LatticeSupportingScadnano,
    ScadnanoLattice,
    _skip_polyT_and_inertname,
    lattice_factory,
)

# from . import fastreduceD as fastreduce
from .glues import Glue, GlueList
from .seeds import Seed, seed_factory
from .tiles import (
    D,
    EdgeLoc,
    SupportsGuards,
    Tile,
    TileList,
    TileSupportingScadnano,
    tile_factory,
)

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    import scadnano
    import xgrow.parseoutput
    import xgrow.tileset as xgt

_gl = {
    EdgeLoc(D.N, (0, 0)): (0, 0, 10, 0),
    EdgeLoc(D.E, (0, 0)): (10, 0, 10, 10),
    EdgeLoc(D.S, (0, 0)): (0, 10, 10, 10),
    EdgeLoc(D.W, (0, 0)): (0, 0, 0, 10),
    EdgeLoc(D.N, (0, 1)): (10, 0, 20, 0),
    EdgeLoc(D.E, (1, 0)): (10, 10, 10, 20),
    EdgeLoc(D.E, (0, 1)): (20, 0, 20, 10),
    EdgeLoc(D.S, (1, 0)): (0, 20, 10, 20),
    EdgeLoc(D.S, (0, 1)): (10, 10, 20, 10),
    EdgeLoc(D.W, (1, 0)): (0, 10, 0, 20),
}

T = TypeVar("T")

XgrowGlueOpts: TypeAlias = Literal["self-complementary", "perfect"]


@dataclass(init=False)
class TileSet(Serializable):
    "Class representing a tileset, whether abstract or sequence-level."
    tiles: TileList[Tile]
    glues: GlueList[Glue]
    seeds: dict[str | int, Seed]
    lattices: dict[str | int, Lattice]
    guards: dict[str | int, list[str]]
    params: dict

    def __init__(
        self,
        tiles: Iterable[Tile] = tuple(),
        glues: Iterable[Glue] = tuple(),
        seeds: Mapping[str | int, Seed] | None = None,
        *,
        lattices: Mapping[str | int, Lattice] | None = None,
        guards: Mapping[str | int, list[str]] = dict(),
        params: dict | None = None,
    ) -> None:
        self.tiles = TileList(tiles)
        self.glues = GlueList(glues)
        self.seeds = dict(seeds) if seeds else dict()
        self.lattices = dict(lattices) if lattices else dict()
        self.guards = dict(guards)
        if params is not None:
            self.params = params
        else:
            self.params = dict()

    ### XGROW METHODS

    def run_xgrow(
        self: TileSet,
        to_lattice: bool = True,
        _include_out: bool = False,
        glues: XgrowGlueOpts = "perfect",
        seed: str | int | Seed | None | Literal[False] = None,
        seed_offset: tuple[int, int] = (0, 0),
        **kwargs: Any,
    ) -> Any:  # FIXME
        """Run the tilesystem in Xgrow."""
        import xgrow
        import xgrow.parseoutput
        from xgrow.parseoutput import XgrowOutput

        xgrow_tileset = self.to_xgrow(seed=seed, seed_offset=seed_offset, glues=glues)

        if not to_lattice:
            return (xgrow.run(xgrow_tileset, **kwargs),)

        out = cast(
            XgrowOutput,
            xgrow.run(
                xgrow_tileset,
                outputopts="array",
                **kwargs,
            ),
        )

        assert out.tiles is not None

        newarray = np.full_like(out.tiles[1:-1, 1:-1], "", dtype=object)

        for ix, xti in np.ndenumerate(out.tiles[1:-1, 1:-1]):
            if xti == 0:
                continue
            if xti > len(xgrow_tileset.tiles):
                continue
            tile_name = xgrow_tileset.tiles[int(xti) - 1].name
            if tile_name in self.tiles.data.keys():
                newarray[ix] = tile_name

        if seed is None:
            if len(self.seeds) > 0:
                seed = next(iter(self.seeds.values()))
        elif seed is False:
            seed = None
        elif isinstance(seed, (int, str)):
            seed = self.seeds[cast(Union[int, str], x)]

        if seed is not None and hasattr(seed, "_lattice"):
            lattice_type: Type[AbstractLattice] = seed._lattice
        else:
            lattice_type = AbstractLattice

        a = lattice_type(newarray, cast(Optional[Seed], seed), seed_offset)

        if not _include_out:
            return a
        else:
            return a, out

    def to_xgrow(
        self,
        glues: XgrowGlueOpts = "perfect",
        seed: str | int | Seed | None | Literal[False] = None,
        seed_offset: tuple[int, int] = (0, 0),
    ) -> xgt.TileSet:
        "Convert Alhambra TileSet to an XGrow TileSet"
        import xgrow.tileset as xgt

        self.tiles.refreshnames()
        self.glues.refreshnames()
        tiles = [t.to_xgrow(glues) for t in self.tiles]

        allglues = self.allglues

        # FIXME
        # bonds = [g.to_xgrow(self_complementary_glues) for g in self.glues]
        if glues == "self-complementary":
            bonds = [
                xgt.Bond(g.name, 0)
                for g in allglues
                if g.name
                and ("null" in g.name or "inert" in g.name or "hairpin" in g.name)
            ]
            xg_glues = []
        elif glues == "perfect":
            bonds = [xgt.Bond(g.name, 0) for g in allglues]
            bonds.extend(
                xgt.Bond(g.complement.name, 0)
                for g in allglues
                if g.complement.name not in allglues
            )
            xg_glues = [
                xgt.Glue(
                    g.name,
                    g.complement.name,
                    g.abstractstrength if g.abstractstrength is not None else 1,
                )
                for g in allglues
            ]
        else:
            raise ValueError("Unknown xgrow glue option", glues)

        if seed is None:
            if self.seeds:
                seed = next(iter(self.seeds.values()))
            else:
                seed = False
        if seed is not False and (isinstance(seed, str) or isinstance(seed, int)):
            seed = self.seeds[seed]
        if seed is False:
            seed_tiles: list[xgt.Tile] = []
            seed_bonds: list[xgt.Bond] = []
            initstate = None
        else:
            seed_tiles, seed_bonds, initstate = cast(Seed, seed).to_xgrow(
                glues, offset=seed_offset
            )

        xgrow_tileset = xgt.TileSet(
            seed_tiles + tiles, seed_bonds + bonds, initstate=initstate, glues=xg_glues
        )

        return xgrow_tileset

    def _to_xgrow_dict(self) -> dict:
        """DEPRECATED: to xgrow dict"""
        return self.to_xgrow().to_dict()

    def summary(self):
        """Returns a short summary line about the TileSet"""
        self.tiles.refreshnames()
        self.glues.refreshnames()
        # self.check_consistent()
        info = {
            "ntiles": len(self.tiles),
            "nrt": len([x for x in self.tiles if not x.is_fake]),
            "nft": len([x for x in self.tiles if x.is_fake]),
            "nends": len(self.glues),
            "ntends": len(self.tiles.glues_from_tiles()),
            "tns": " ".join(x.name for x in self.tiles if x.name),
            "ens": " ".join(x.name for x in self.glues if x.name)
            # if ("info" in self.keys() and "name" in self["info"].keys())
            # else "",
        }
        tun = sum(1 for x in self.tiles if x.name is None)
        if tun > 0:
            info["tns"] += " ({} unnamed)".format(tun)
        eun = sum(1 for x in self.glues if x.name is None)
        if eun > 0:
            info["ens"] += " ({} unnamed)".format(eun)
        if info["nft"] > 0:
            info["nft"] = " (+ {} fake)".format(info["nft"])
        else:
            info["nft"] = ""
        return "TileSet: {nrt} tiles{nft}, {nends} ends, {ntends} ends in tiles.\nTiles: {tns}\nEnds:  {ens}".format(
            **info
        )

    def __repr__(self) -> str:
        self.tiles.refreshnames()
        self.glues.refreshnames()
        return f"TileSet({len(self.tiles)} tiles, {len(self.glues)} glues)"

    def __str__(self) -> str:
        return self.summary()

    @classmethod
    def from_scadnano(
        cls: Type[TileSet], des: scadnano.Design, ret_fails: bool = False
    ) -> TileSet:
        """Create TileSet from Scadnano Design."""
        import scadnano

        ts = cls()
        tiles: TileList[TileSupportingScadnano] = TileList()
        ts.glues = GlueList()
        positions: dict[tuple[int, int], str] = {}

        for strand in des.strands:
            try:
                t, o = tile_factory.from_scadnano(strand, return_position=True)
            except ValueError:
                warnings.warn(f"Failed to import strand {strand.name}.")
            except NotImplementedError:
                warnings.warn(f"Failed to import strand {strand.name}.")
            else:
                positions[o] = t.ident()
                if t.name in tiles.data.keys():
                    if t != tiles[t.ident()]:
                        warnings.warn(f"Skipping unequal duplicate strand {t.name}.")
                else:
                    tiles.add(t)

        ts.tiles = cast(TileList[Tile], tiles)
        ts.lattices = {0: cast(Lattice, ScadnanoLattice(positions))}
        return ts

    def to_scadnano(self, lattice: LatticeSupportingScadnano = None) -> scadnano.Design:
        """Export TileSet (with lattice) as Scadnano Design."""
        import scadnano

        self.tiles.refreshnames()
        self.glues.refreshnames()
        if hasattr(lattice, "seed") and hasattr(lattice.seed, "update_details"):
            lattice.seed.update_details(self.allglues, self.tiles)
        if lattice:
            return lattice.to_scadnano(self)
        for tlattice in self.lattices:
            if isinstance(tlattice, LatticeSupportingScadnano):
                return tlattice.to_scadnano(self)
        raise ValueError

    def to_dict(self) -> dict:
        d: dict[str, Any] = {}
        self.tiles.refreshnames()
        self.glues.refreshnames()
        allglues = self.glues | self.tiles.glues_from_tiles()
        refglues = set(allglues.data.keys())  # FIXME

        if self.tiles:
            d["tiles"] = [t.to_dict(refglues=refglues) for t in self.tiles.aslist()]
        if allglues:
            d["glues"] = [g.to_dict() for g in allglues.aslist()]
        if self.seeds:
            d["seeds"] = {k: v.to_dict() for k, v in self.seeds.items()}
        if self.lattices:
            d["lattices"] = {k: v.asdict() for k, v in self.lattices.items()}
        if self.guards:
            d["guards"] = self.guards
        if self.params:
            d["params"] = self.params.copy()
        return d

    @classmethod
    def from_dict(cls: Type[TileSet], d: dict) -> TileSet:
        ts = cls()
        ts.tiles = TileList(Tile.from_dict(x) for x in d.get("tiles", []))
        ts.glues = GlueList(Glue.from_dict(x) for x in d.get("glues", []))
        ts.seeds = {k: seed_factory.from_dict(v) for k, v in d.get("seeds", {}).items()}
        ts.guards = {k: v for k, v in d.get("guards", {}).items()}
        if not ts.seeds and "seed" in d:
            ts.seeds = {0: seed_factory.from_dict(d["seed"])}
        ts.lattices = {
            k: lattice_factory.from_dict(v) for k, v in d.get("lattices", {}).items()
        }
        if "params" in d:
            ts.params = copy.deepcopy(d["params"])
        return ts

    def _serialize(self) -> Any:
        return self.to_dict()

    @classmethod
    def _deserialize(cls, input: Any) -> TileSet:
        return cls.from_dict(input)

    @property
    def allglues(self) -> GlueList:
        return self.tiles.glues_from_tiles() | self.glues

    @property
    def alldomains(self) -> GlueList:
        return self.tiles.domains_from_tiles() | self.glues

    def lattice_tiles(
        self,
        lattice: AbstractLattice | int | str | np.ndarray,
        *,
        x: int | slice | None = None,
        y: int | slice | None = None,
        copy: bool = False,
    ) -> list[Tile]:
        """Return a list of (unique) tiles in a lattice, potentially taking a slice of the lattice.

        Parameters
        ----------
        lattice : AbstractLattice | int | str
            Lattice or reference to a lattice in the tileset
        x : int | slice | None, optional
            index in the lattice, by default None
        y : int | slice | None, optional
            index in the lattice, by default None
        copy : bool, optional
            return copies if True (useful for creating a new set) or tiles in the set if False (useful for modifying the set), by default False
        """

        if isinstance(lattice, (int, str)):
            lattice = cast(AbstractLattice, self.lattices[lattice])
        elif not isinstance(lattice, AbstractLattice):
            lattice = AbstractLattice(lattice)

        tilenames = np.unique(lattice.grid[x, y])

        if copy:
            return [self.tiles[t].copy() for t in tilenames]
        else:
            return [self.tiles[t] for t in tilenames]

    def create_guards_square(
        self,
        lattice: AbstractLattice,
        square_size: int,
        init_x: int = 0,
        init_y: int = 0,
        skip: Callable[[Glue], bool] = _skip_polyT_and_inertname,
    ) -> list[str]:
        glues: set[str] = set()
        for xi in range(init_x, lattice.grid.shape[0], square_size):
            glues.update(
                g.ident()
                for tile in lattice.grid[xi, :]
                if isinstance(self.tiles[tile], SupportsGuards)
                for g in self.tiles[tile].create_guards("S")
                if not skip(g)
            )
        for yi in range(init_y, lattice.grid.shape[1], square_size):
            glues.update(
                g.ident()
                for tile in lattice.grid[:, yi]
                if isinstance(self.tiles[tile], SupportsGuards)
                for g in self.tiles[tile].create_guards("E")
                if not skip(g)
            )
        return list(glues)

    def create_abstract_diagram(
        self,
        lattice: AbstractLattice | str | int | np.ndarray | None,
        filename=None,
        scale=1,
        guards: Collection[str] | str | int = tuple(),
        seed: str | bool | Seed = True,
        seed_offset: tuple[int, int] = (0, 0),
        **options,
    ):
        """Create an SVG layout diagram from a lattice.

        This currently uses the abstract diagram bases to create the
        layout diagrams.

        Parameters
        ----------

        xgrowarray : ndarray or dict
            Xgrow output.  This may be a numpy array of
            an xgrow state, obtained in some way, or may be the 'array'
            output of xgrow.run.

        filename : string
            File name / path of the output file.

        """

        if isinstance(lattice, str) or isinstance(lattice, int):
            lt = self.lattices[lattice]
            assert isinstance(lt, AbstractLattice)
            lattice = lt
        elif lattice is None:
            lt = next(iter(self.lattices.values()))
            assert isinstance(lt, AbstractLattice)
            lattice = lt
        elif not isinstance(lattice, AbstractLattice):
            lattice = AbstractLattice(lattice)

        if isinstance(guards, str) or isinstance(guards, int):
            guards = self.guards[guards]

        d = draw.Drawing(200, 200)

        svgtiles = {}

        for tile in self.tiles:
            svgtiles[tile.name] = tile.abstract_diagram(**options)

        minxi = 10000
        minyi = 10000
        maxxi = 0
        maxyi = 0

        for (yi, xi), tn in np.ndenumerate(lattice.grid):
            if not tn in svgtiles.keys():
                continue
            minxi = min(minxi, xi)
            minyi = min(minyi, yi)
            maxxi = max(maxxi, xi)
            maxyi = max(maxyi, yi)
            d.append(draw.Use(svgtiles[tn], xi * 10, yi * 10))

        if seed is True:
            try:
                seed = next(iter(self.seeds.values()))
            except StopIteration:
                seed = False
        elif isinstance(seed, str):
            seed = self.seeds[seed]

        if len(guards) > 0:
            for (yi, xi), tn in np.ndenumerate(lattice.grid):
                if tn == "":
                    continue
                t = self.tiles[tn]
                for g, pos in zip(t.edges, t.edge_locations):  # FIXME: deal with duples
                    if g.complement.ident() in guards:
                        d.append(
                            draw.Line(
                                xi * 10 + _gl[pos][0],
                                yi * 10 + _gl[pos][1],
                                xi * 10 + _gl[pos][2],
                                yi * 10 + _gl[pos][3],
                                stroke="black",
                                stroke_width=2.5,
                            )
                        )
                        d.append(
                            draw.Line(
                                xi * 10 + _gl[pos][0],
                                yi * 10 + _gl[pos][1],
                                xi * 10 + _gl[pos][2],
                                yi * 10 + _gl[pos][3],
                                stroke="red",
                                stroke_width=1.0,
                            )
                        )

        d.viewBox = (
            minxi * 10,
            minyi * 10,
            (2 + maxxi - minxi) * 10,
            (2 + maxyi - minyi) * 10,
        )

        d.pixelScale = 3

        if filename:
            d.saveSvg(filename)
        else:
            return d

    def reduce_tiles(
        self,
        preserve=("s22", "ld"),
        tries=10,
        threads=1,
        returntype="equiv",
        best=1,
        key=None,
        initequiv=None,
    ):
        """
        Apply tile reduction algorithm, preserving some set of properties, and using a multiprocessing pool.

        Parameters
        ----------
        tileset: TileSet
            The system to reduce.

        preserve: a tuple or list of strings, optional
            The properties to preserve.  Currently supported are 's1' for first order
            sensitivity, 's2' for second order sensitivity, 's22' for two-by-two sensitivity,
            'ld' for small lattice defects, and 'gs' for glue sense (to avoid spurious
            hierarchical attachment).  Default is currently ('s22', 'ld').

        tries: int, optional
            The number of times to run the algorithm.

        threads: int, optional
            The number of threads to use (using multiprocessing).

        returntype: 'TileSet' or 'equiv' (default 'equiv')
            The type of object to return.  If 'equiv', returns an array of glue equivalences
            (or list, if best != 1) that can be applied to the tileset with apply_equiv, or used
            for further reduction.  If 'TileSet', return a TileSet with the equiv already applied
            (or a list, if best != 1).

        best: int or None, optional
            The number of systems to return.  If 1, the result will be returned
            directly; if k > 1, a list will be returned of the best k results (per cmp);
            if k = None, a list of *all* results will be returned, sorted by cmp. (default 1)

        key: function (ts, equiv1, equiv2) -> some number/comparable
            A comparison function for equivs, to sort the results. FIXME: documentation needed.
            Default (if None) here is to sort by number of glues in the system, regardless of number
            of tiles.

        initequiv: equiv
            If provided, the equivalence array to start from.  If None, start from the tileset without
            any merged glues.

        Returns
        -------
        reduced: single TileSet or equiv, or list
            The reduced system/systems
        """
        raise NotImplementedError
        # from fastreduceD import fastreduce
        # return fastreduce.reduce_tiles(
        #     self, preserve, tries, threads, returntype, best, key, initequiv
        # )

    def reduce_ends(
        self,
        preserve=["s22", "ld"],
        tries=10,
        threads=1,
        returntype="equiv",
        best=1,
        key=None,
        initequiv=None,
    ):
        """
        Apply end reduction algorithm, preserving some set of properties, and using a multiprocessing pool.

        Parameters
        ----------
        tileset: TileSet
            The system to reduce.

        preserve: a tuple or list of strings, optional
            The properties to preserve.  Currently supported are 's1' for first order
            sensitivity, 's2' for second order sensitivity, 's22' for two-by-two sensitivity,
            'ld' for small lattice defects, and 'gs' for glue sense (to avoid spurious
            hierarchical attachment).  Default is currently ('s22', 'ld').

        tries: int, optional
            The number of times to run the algorithm.

        threads: int, optional
            The number of threads to use (using multiprocessing).

        returntype: 'TileSet' or 'equiv' (default 'equiv')
            The type of object to return.  If 'equiv', returns an array of glue equivalences
            (or list, if best != 1) that can be applied to the tileset with apply_equiv, or used
            for further reduction.  If 'TileSet', return a TileSet with the equiv already applied
            (or a list, if best != 1).

        best: int or None, optional
            The number of systems to return.  If 1, the result will be returned
            directly; if k > 1, a list will be returned of the best k results (per cmp);
            if k = None, a list of *all* results will be returned, sorted by cmp. (default 1)

        key: function (ts, equiv1, equiv2) -> some number/comparable
            A comparison function for equivs, to sort the results. FIXME: documentation needed.
            Default (if None) here is to sort by number of glues in the system, regardless of number
            of tiles.

        initequiv: equiv
            If provided, the equivalence array to start from.  If None, start from the tileset without
            any merged glues.

        Returns
        -------
        reduced: single TileSet or equiv, or list
            The reduced system/systems
        """
        raise NotImplementedError
        # from fastreduceD import fastreduce
        # return fastreduce.reduce_ends(
        #     self, preserve, tries, threads, returntype, best, key, initequiv
        # )

    def latticedefects(self, direction="e", depth=2, pp=True, rotate=False):
        """
        Calculate and show possible small lattice defect configurations.
        """
        raise NotImplementedError
        # from . import latticedefect

        # return latticedefect.latticedefects(
        #     self, direction=direction, depth=depth, pp=pp, rotate=rotate
        # )

    # FIXME: disabled temporarily for mypy main branch
    from ._tilesets_dx import (  # type: ignore
        dx_plot_adjacent_regions,
        dx_plot_se_hists,
        dx_plot_se_lv,
        dx_plot_side_strands,
    )  # type: ignore

    from .nuad import tileset_to_nuad_design as to_nuad_design  # type: ignore

    from .nuad import update_nuad_design as update_nuad_design # type: ignore

    from .nuad import load_nuad_design as load_nuad_design  # type: ignore

    def apply_equiv(self, equiv):
        """
        Apply an equivalence array (from, eg, `TileSet.reduce_ends` or `TileSet.reduce_tiles`).

        Parameters
        ----------
        equiv : ndarray
            An equivalence array, *for this tileset*, generated by reduction functions.

        Returns
        -------
        TileSet
            A tileset with the equivalence array, and thus the reduction, applied.
        """
        return fastreduce._FastTileSet(self).applyequiv(self, equiv)

    def check_consistent(self):
        """Check the TileSet consistency.

        Check a number of properties of the TileSet for consistency.
        In particular:

           * Each tile must pass Tile.check_consistent()
           * TileSet.ends and TileSet.tiles.endlist() must not contain conflicting
             ends or end sequences.
           * If there is a seed:
               * It must be of an understood type (it must be in seeds.seedtypes)
               * All adapter locations must be valid.
               * The seed must pass its check_consistent and check_sequence.
        """
        # * END LIST The end list itself must be consistent.
        # ** Each end must be of understood type
        # ** Each end must have a valid sequence or no sequence
        # ** There must be no more than one instance of each name
        # ** WARN if there are ends with no namecounts
        # * TILE LIST
        # ** each tile must be of understood type (must parse)
        # ** ends in the tile list must be consistent (must merge)
        # ** there must be no more than one tile with each name
        # self.tiles.check_consistent()
        endsfromtiles = self.tiles.glues_from_tiles()

        # ** WARN if any end that appears does not have a complement used or vice versa
        # ** WARN if there are tiles with no name
        # * TILE + END
        # ** The tile and end lists must merge validly
        # (checks sequences, adjacents, types, complements)
        self.glues | endsfromtiles

        # ** WARN if tilelist has end references not in ends
        # ** WARN if merge is not equal to the endlist
        # ** WARN if endlist has ends not used in tilelist
        # * ADAPTERS / SEEDS
        # SEED stuff was here

    def copy(self):
        """Return a full (deep) copy of the TileSet"""
        return copy.deepcopy(self)
