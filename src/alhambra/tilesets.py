from __future__ import annotations
from abc import ABC, abstractmethod

import copy
from io import TextIOWrapper
import logging
from os import PathLike
from pathlib import Path
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
from typing_extensions import Self, TypeAlias

import numpy as np

from alhambra.seq import MergeConflictError

from . import drawing

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
from .glues import Glue, GlueList, SSGlue
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

import logging

log = logging.getLogger(__name__)

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    import scadnano
    import xgrow.parseoutput
    import xgrow.tileset as xgt
    import stickydesign as sd

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


class XgrowGlueOpts(ABC):
    def get_xgrow_gse(self, tileset: "TileSet") -> float | None:
        return None

    @abstractmethod
    def calculate_gses(
        self, tileset: "TileSet"
    ) -> tuple[list[xgt.Bond], list[xgt.Glue]]:
        ...

    def glue_name_map(self, tileset: "TileSet") -> Callable[[str], str]:
        return lambda x: x

    @classmethod
    def from_str(cls, glueopts: str) -> XgrowGlueOpts:
        if glueopts == "self-complementary":
            return GrowSelfComplementaryGlues()
        elif glueopts == "perfect":
            return GrowPerfectGlues()
        elif glueopts == "orthogonal":
            return GrowOrthogonalGlues()
        elif glueopts == "full":
            return GrowFullGlues()
        else:
            raise ValueError(f"Unknown glueopts: {glueopts}")


class GrowSelfComplementaryGlues(XgrowGlueOpts):
    def get_xgrow_gse(self, tileset: "TileSet") -> float | None:
        return None

    def glue_name_map(self, tileset: "TileSet") -> Callable[[str], str]:
        return lambda x: x[:-1] if x.endswith("*") else x

    def calculate_gses(
        self, tileset: "TileSet"
    ) -> tuple[list[xgt.Bond], list[xgt.Glue]]:
        allglues = tileset.allglues
        bonds = [
            xgt.Bond(g.name, 0)
            for g in allglues
            if g.name and ("null" in g.name or "inert" in g.name or "hairpin" in g.name)
        ]
        return bonds, []


SD_ENERGETICS_CLASSES: dict[str, "sd.EnergeticsBasic"] = {}
HAS_SD_ACCEL: bool = False


def _generate_stickydesign_energetic_classes() -> None:
    global SD_ENERGETICS_CLASSES
    global HAS_SD_ACCEL

    if SD_ENERGETICS_CLASSES:
        return
    import stickydesign as sd

    SD_ENERGETICS_CLASSES = {"SSGlue": sd.EnergeticsBasic, "DXGlue": sd.EnergeticsDAOE}
    HAS_SD_ACCEL = sd.energetics_basic.ACCEL


class GrowPerfectGlues(XgrowGlueOpts):
    def get_xgrow_gse(self, tileset: "TileSet") -> float | None:
        return None

    def calculate_gses(
        self, tileset: "TileSet"
    ) -> tuple[list[xgt.Bond], list[xgt.Glue]]:
        import xgrow.tileset as xgt

        allglues = tileset.allglues
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
        return bonds, xg_glues


@dataclass
class GrowOrthogonalGlues(XgrowGlueOpts):
    temperature: float | None = None
    alpha: float | None = None

    def get_xgrow_gse(self, tileset: "TileSet") -> float | None:
        return 1.0

    def calculate_gses(
        self, tileset: "TileSet"
    ) -> tuple[list[xgt.Bond], list[xgt.Glue]]:
        import xgrow.tileset as xgt
        import stickydesign as sd

        alpha = tileset.params["alpha"] if self.alpha is None else self.alpha
        temperature = (
            tileset.params["temperature"]
            if self.temperature is None
            else self.temperature
        )

        allglues = tileset.allglues

        _generate_stickydesign_energetic_classes()

        R = 1.9872041e-3  # kcal/mol/K
        T_in_K = temperature + 273.15
        RT = R * T_in_K

        xgbonds = [xgt.Bond(g.name, 0) for g in allglues]
        xgbonds.extend(
            xgt.Bond(g.complement.name, 0)
            for g in allglues
            if g.complement.name not in allglues
        )

        sd_energetics = {
            k: v(temperature=temperature) for k, v in SD_ENERGETICS_CLASSES.items()
        }

        sg: dict[
            tuple[str, str, int], list[SSGlue]
        ] = {}  # FIXME: should support DXGlue
        for g in allglues:
            k = (g.__class__.__name__, g.etype, g.dna_length)
            sg[k] = sg.get(k, [])
            sg[k].append(g)

        sge = {
            k: (
                [x.ident() for x in v],
                sd.endarray([x.sequence.base_str.lower() for x in v], k[1]),
            )
            for k, v in sg.items()
        }

        bonds: dict[str, float] = {}
        for k, v in sge.items():
            m = sd_energetics[k[0]].matching_uniform(v[1]) / RT + alpha
            bonds |= {n: x for n, x in zip(v[0], m) if not n.endswith("*")}

        return xgbonds, [xgt.Glue(n, n + "*", v) for n, v in bonds.items()]


@dataclass
class GrowFullGlues(XgrowGlueOpts):
    temperature: float | None = None
    threshold: float = 0.0
    alpha: float | None = None

    def get_xgrow_gse(self, tileset: "TileSet") -> float | None:
        return 1.0

    def calculate_gses(
        self, tileset: "TileSet"
    ) -> tuple[list[xgt.Bond], list[xgt.Glue]]:
        import xgrow.tileset as xgt
        import stickydesign as sd

        alpha = tileset.params["alpha"] if self.alpha is None else self.alpha
        temperature = (
            tileset.params["temperature"]
            if self.temperature is None
            else self.temperature
        )

        allglues = tileset.allglues

        _generate_stickydesign_energetic_classes()
        if not HAS_SD_ACCEL:
            warnings.warn(
                "StickyDesign acceleration not available: full glue calculations may be very slow."
            )

        R = 1.9872041e-3  # kcal/mol/K
        T_in_K = temperature + 273.15
        RT = R * T_in_K

        bonds = [xgt.Bond(g.name, 0) for g in allglues]
        bonds.extend(
            xgt.Bond(g.complement.name, 0)
            for g in allglues
            if g.complement.name not in allglues
        )

        sd_energetics = {
            k: v(temperature=temperature) for k, v in SD_ENERGETICS_CLASSES.items()
        }

        sg: dict[tuple[str, str, int], list[SSGlue]] = {}
        for g in allglues:
            k = (g.__class__.__name__, g.etype, g.dna_length)
            sg[k] = sg.get(k, [])
            sg[k].append(g)

        sge = {
            k: (
                [x.ident() for x in v],
                sd.endarray([x.sequence.base_str.lower() for x in v], k[1]),
            )
            for k, v in sg.items()
        }

        gluelinks: dict[tuple[str, str], float] = {}
        for k, v in sge.items():
            glm = (
                sd_energetics[k[0]]
                .uniform(
                    np.repeat(v[1], v[1].shape[0], 0), np.tile(v[1], (v[1].shape[0], 1))
                )
                .reshape((v[1].shape[0], v[1].shape[0]))
                / RT
                + alpha
            )

            for i in range(0, len(v[0])):
                for j in range(0, len(v[0])):
                    gluelinks[(v[0][i], v[0][j])] = glm[i, j]

        return bonds, [
            xgt.Glue(n[0], n[1], v) for n, v in gluelinks.items() if v > self.threshold
        ]


from alhambra_mixes import Mix, Q_, nM, ureg
from alhambra_mixes.units import _parse_conc_required, _ratio
from alhambra_mixes.logging import log as log_mix


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

    @classmethod
    def from_mix(
        cls: Type[TileSet],
        mix: Mix,
        tilesets_or_lists: TileSet | TileList | Iterable[TileSet | TileList],
        *,
        seed: bool | Seed = False,
        base_conc: ureg.Quantity | str = Q_(100.0, nM),
    ) -> TileSet:
        """
        Given some :any:`TileSet`\ s, or lists of :any:`Tile`\ s from which to
        take tiles, generate an TileSet from the mix.
        """

        from .tiles import BaseSSTile

        base_conc = _parse_conc_required(base_conc)

        newts = cls()

        if isinstance(tilesets_or_lists, (TileList, TileSet)):
            tilesets_or_lists = [tilesets_or_lists]

        for name, row in mix.all_components().iterrows():
            new_tile = None
            for tl_or_ts in tilesets_or_lists:
                try:
                    if isinstance(tl_or_ts, TileSet):
                        tile = tl_or_ts.tiles[name]
                    else:
                        tile = tl_or_ts[name]
                    new_tile = tile.copy()
                    if isinstance(new_tile, BaseSSTile) and (
                        (seq := getattr(row["component"], "sequence", None)) is not None
                    ):
                        try:
                            new_tile.sequence |= seq
                        except MergeConflictError as e:
                            # FIXME: we should keep track of all of these and have
                            # an error at the end.
                            raise ValueError(
                                f"Component {name} has sequence that does not match "
                                f"tileset tile {tile.name}. First has {seq}, second has {tile.sequence}."
                            ) from None

                    new_tile.stoic = float(
                        _ratio(Q_(row["concentration_nM"], nM), base_conc)
                    )
                    newts.tiles.add(new_tile)
                    break
                except KeyError:
                    pass
            if new_tile is None:
                log_mix.warn(f"Component {name} not found in tile lists.")

        firstts = next(iter(tilesets_or_lists))
        assert isinstance(firstts, TileSet)

        if seed is True:
            newts.seeds["default"] = firstts.seeds["default"]
        elif seed is False:
            pass
        elif isinstance(seed, Seed):
            newts.seeds["default"] = seed
        elif isinstance(seed, str):
            newts.seeds["default"] = firstts.seeds[seed]

        if len(newts.tiles) == 0:
            raise ValueError("No mix components match tiles.")

        return newts

    ### XGROW METHODS

    def run_xgrow(
        self: TileSet,
        to_lattice: bool = True,
        _include_out: bool = False,
        glues: XgrowGlueOpts | str | None = None,
        seed: str | int | Seed | None | Literal[False] = None,
        seed_offset: tuple[int, int] | None = None,
        xgrow_seed: tuple[int, int, int | str] | None = None,
        **kwargs: Any,
    ) -> Any:  # FIXME
        """Run the tilesystem in Xgrow."""
        import xgrow
        import xgrow.parseoutput
        from xgrow.parseoutput import XgrowOutput

        xgrow_tileset = self.to_xgrow(
            seed=seed, seed_offset=seed_offset, glue_handling=glues
        )

        if xgrow_seed is not None:
            kwargs["seed"] = xgrow_seed

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
            seed = self.seeds[cast(Union[int, str], seed)]

        if seed is not None and hasattr(seed, "_lattice"):
            lattice_type: Type[AbstractLattice] = seed._lattice  # type: ignore
        else:
            lattice_type = AbstractLattice

        a = lattice_type(newarray, cast(Optional[Seed], seed), seed_offset)

        if not _include_out:
            return a
        else:
            return a, out

    def to_rgrow(
        self,
        glue_handling: XgrowGlueOpts | str | None = None,
        seed: str | int | Seed | None | Literal[False] = None,
        seed_offset: tuple[int, int] | None = None,
        **kwargs,
    ):
        import rgrow as rg
        import xgrow.tileset as xgt

        d = self.to_rgrow_dict(
            glue_handling=glue_handling, seed=seed, seed_offset=seed_offset, **kwargs
        )

        return rg.TileSet.from_dict(d)  # type: ignore
        # (FIXME: rgrow needs a fix here)

    def to_rgrow_dict(
        self,
        glue_handling: XgrowGlueOpts | str | None = None,
        seed: str | int | Seed | None | Literal[False] = None,
        seed_offset: tuple[int, int] | None = None,
        **kwargs,
    ):
        import rgrow as rg
        import xgrow.tileset as xgt

        d = self.to_xgrow(
            glue_handling=glue_handling, seed=seed, seed_offset=seed_offset
        ).to_dict()

        d["options"] = d.pop("xgrowargs")
        if "initstate" in d:
            d["options"]["seed"] = d.pop("initstate")

        if "k" in d["options"]:
            del d["options"]["k"]  # k value for xgrow is modified; FIXME

        for k in [
            "concentration",
            "alpha",
            "temperature",
            "model",
            "chunk_handling",
            "chunk_size",
            "canvas-type",
            "fission",
            "threshold",
            "size",
        ]:
            if k in self.params:
                d["options"][k] = self.params[k]

        for k, v in kwargs.items():
            if v is not None:
                d["options"][k] = v
            else:
                d["options"].pop(k, None)

        return d

    def to_xgrow(
        self,
        glue_handling: XgrowGlueOpts | str | None = None,
        seed: str | int | Seed | None | Literal[False] = None,
        seed_offset: tuple[int, int] | None = None,
    ) -> xgt.TileSet:
        "Convert Alhambra TileSet to an XGrow TileSet"
        import xgrow.tileset as xgt

        if glue_handling is None:
            glue_handling = self.params.get("glue-handling", "perfect")

        glue_handling = (
            XgrowGlueOpts.from_str(glue_handling)
            if isinstance(glue_handling, str)
            else cast(XgrowGlueOpts, glue_handling)
        )

        gluenamemap = glue_handling.glue_name_map(self)

        self.tiles.refreshnames()
        self.glues.refreshnames()
        tiles = [t.to_xgrow(gluenamemap) for t in self.tiles]

        bonds, xg_glues = glue_handling.calculate_gses(self)

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
                gluenamemap, offset=seed_offset, xgtiles=tiles
            )

        xgrow_tileset = xgt.TileSet(
            seed_tiles + tiles, seed_bonds + bonds, initstate=initstate, glues=xg_glues
        )

        ggse = glue_handling.get_xgrow_gse(self)
        if ggse is not None:
            xgrow_tileset.xgrowargs.Gse = ggse

        gconc = self.params.get("concentration", None)
        galpha = self.params.get("alpha", None)
        if gconc is not None and galpha is not None:
            xgrow_tileset.xgrowargs.Gmc = galpha - np.log(gconc * 1e-9)
            xgrow_tileset.xgrowargs.k = 1e6 * np.exp(galpha)

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

    def to_scadnano(
        self, lattice: LatticeSupportingScadnano | None = None
    ) -> scadnano.Design:
        """Export TileSet (with lattice) as Scadnano Design."""
        import scadnano

        self.tiles.refreshnames()
        self.glues.refreshnames()
        if lattice is not None:
            if hasattr(lattice, "seed") and hasattr(lattice.seed, "update_details"):
                lattice.seed.update_details(self.allglues, self.tiles)  # type: ignore
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
        try:
            ts.guards = {k: v for k, v in d.get("guards", {}).items()}
        except AttributeError:
            log.warning("Failed to load guards.")
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

        d = drawing.Drawing(600, 600)

        svgtiles = {}

        for tile in self.tiles:
            svgtiles[tile.name] = tile.abstract_diagram(**options)
            d.defs.append(svgtiles[tile.name])

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
            d.elements.append(drawing.Use(svgtiles[tn], xi * 10, yi * 10))

        if seed is True:
            try:
                seed = next(iter(self.seeds.values()))
            except StopIteration:
                seed = False
        elif isinstance(seed, str):
            seed = self.seeds[seed]

        # if len(guards) > 0:
        #     for (yi, xi), tn in np.ndenumerate(lattice.grid):
        #         if tn == "":
        #             continue
        #         t = self.tiles[tn]
        #         for g, pos in zip(t.edges, t.edge_locations):  # FIXME: deal with duples
        #             if g.complement.ident() in guards:
        #                 d.elements.append(
        #                     drawing.Line(
        #                         xi * 10 + _gl[pos][0],
        #                         yi * 10 + _gl[pos][1],
        #                         xi * 10 + _gl[pos][2],
        #                         yi * 10 + _gl[pos][3],
        #                         stroke="black",
        #                         stroke_width=2.5,
        #                     )
        #                 )
        #                 d.elements.append(
        #                     drawing.Line(
        #                         xi * 10 + _gl[pos][0],
        #                         yi * 10 + _gl[pos][1],
        #                         xi * 10 + _gl[pos][2],
        #                         yi * 10 + _gl[pos][3],
        #                         stroke="red",
        #                         stroke_width=1.0,
        #                     )
        #                 )

        d.viewBox = (
            minxi * 10,
            minyi * 10,
            (2 + maxxi - minxi) * 10,
            (2 + maxyi - minyi) * 10,
        )

        # d.pixelScale = 3

        if filename:
            d.save_svg(filename)
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

    from .nuad import update_nuad_design as update_nuad_design  # type: ignore

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
        raise NotImplementedError
        # return fastreduce._FastTileSet(self).applyequiv(self, equiv)

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

    @classmethod
    def from_file(
        cls,
        path_or_stream: TextIOWrapper | str | PathLike[str],
        format: Literal["json", "yaml", None] = None,
    ) -> "TileSet":

        if isinstance(path_or_stream, str) or isinstance(path_or_stream, PathLike):
            p = Path(path_or_stream)
            if format is None:
                if p.suffix == ".json":
                    format = "json"
                elif p.suffix in [".yaml", ".yml"]:
                    format = "yaml"

            stream = p.open("r")

        if format is None:
            log.warning("No format specified, trying json, then yaml.")
            try:
                return cls.from_json(stream)
            except:
                stream.seek(0)
                return cls.from_yaml(stream)

        if format == "json":
            return cls.from_json(stream)
        elif format == "yaml":
            return cls.from_yaml(stream)
        else:
            raise ValueError(f"Unknown format {format}")

    def to_file(self, path_or_stream: str | PathLike[str] | TextIOWrapper):
        if isinstance(path_or_stream, (str, PathLike)):
            stream = open(path_or_stream, "w")
        else:
            stream = path_or_stream

        return self.to_json(stream)
