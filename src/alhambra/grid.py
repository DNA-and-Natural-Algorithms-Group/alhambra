from __future__ import annotations

import dataclasses
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Callable, Literal, Protocol, Type, TypeVar, cast

import numpy as np
import scadnano

from alhambra.glues import Glue, SSGlue
from alhambra.seeds import Seed, SeedSupportingScadnano
from alhambra.seq import Seq
from alhambra.tiles import D, SupportsGuards, Tile, TileSupportingScadnano

if TYPE_CHECKING:
    from alhambra.tilesets import TileSet


class Lattice(ABC):
    @abstractmethod
    def __getitem__(self, index) -> str | Any:
        ...

    @abstractmethod
    def __setitem__(self, index, v):
        ...

    @property
    @abstractmethod
    def seed(self) -> Seed | None:
        ...

    def asdict(self) -> dict[str, Any]:
        raise NotImplementedError

    @classmethod
    def fromdict(cls, d: dict[str, Any]) -> Lattice:
        raise NotADirectoryError


@dataclass(init=False)
class AbstractLattice(Lattice):
    grid: np.ndarray
    seed: Seed | None = None
    seed_offset: tuple[int, int] = (0, 0)

    def __getitem__(self, index) -> str | Any:
        return AbstractLattice(self.grid[index])

    def __setitem__(self, index, v):
        self.grid[index] = v

    def __init__(
        self,
        v: AbstractLattice | np.ndarray,
        seed: Seed | None = None,
        seed_offset: tuple[int, int] | None = None,
    ) -> None:
        if isinstance(v, AbstractLattice):
            self.grid = v.grid
            self.seed = v.seed
            self.seed_offset = v.seed_offset
        else:
            self.grid = np.array(v)
        if seed is not None:
            self.seed = seed
        if seed_offset is not None:
            self.seed_offset = seed_offset

    def asdict(self) -> dict[str, Any]:
        d: dict[str, Any] = {}
        d["type"] = self.__class__.__name__
        d["grid"] = self.grid.tolist()
        return d

    @property
    def tilenames(self) -> list[str]:
        return list(np.unique(self.grid))

    @classmethod
    def fromdict(cls: Type[AL], d: dict[str, Any]) -> AL:
        return cls(np.array(d["grid"]))

    @classmethod
    def empty(cls, shape):
        return cls(np.full(shape, "", dtype=object))


class LatticeSupportingScadnano(Lattice):
    @abstractmethod
    def to_scadnano_lattice(self) -> ScadnanoLattice:
        ...

    seed: SeedSupportingScadnano | None = None

    def to_scadnano(self, tileset: "TileSet") -> "scadnano.Design":
        tileset.tiles.refreshnames()
        tileset.glues.refreshnames()
        scl = self.to_scadnano_lattice()
        max_helix = max(helix for helix, offset in scl.positions) + 4
        des = scadnano.Design(helices=[scadnano.Helix() for _ in range(0, max_helix)])

        for (helix, offset), tilename in scl.positions.items():
            cast(TileSupportingScadnano, tileset.tiles[tilename]).to_scadnano(
                des, helix, offset
            )

        if scl.seed is not None:
            scl.seed.to_scadnano(des, scl.seed_position[0], scl.seed_position[1])

        return des


class AbstractLatticeSupportingScadnano(AbstractLattice):
    @abstractmethod
    def to_scadnano_lattice(self) -> ScadnanoLattice:
        ...

    seed: SeedSupportingScadnano | None = None

    def to_scadnano(self, tileset: "TileSet") -> "scadnano.Design":
        tileset.tiles.refreshnames()
        tileset.glues.refreshnames()
        scl = self.to_scadnano_lattice()
        max_helix = max(helix for helix, offset in scl.positions) + 4
        des = scadnano.Design(helices=[scadnano.Helix() for _ in range(0, max_helix)])

        for (helix, offset), tilename in scl.positions.items():
            cast(TileSupportingScadnano, tileset.tiles[tilename]).to_scadnano(
                des, helix, offset
            )

        if scl.seed is not None:
            scl.seed.to_scadnano(des, scl.seed_position[0], scl.seed_position[1])

        return des


@dataclass
class ScadnanoLattice(LatticeSupportingScadnano):
    positions: dict[tuple[int, int], str] = field(default_factory=lambda: {})
    seed: SeedSupportingScadnano | None = None
    seed_position: tuple[int, int] = (0, 0)

    def __getitem__(self, index: tuple[int, int]) -> str | None:
        return self.positions[index]

    def __setitem__(self, index: tuple[int, int], v: str):
        self.positions[cast(tuple[int, int], index)] = cast(str, v)

    def findtile(self, tile: str | Tile) -> list[tuple[int, int]]:
        if isinstance(tile, Tile):
            tile = tile.ident()
        return [k for k, v in self.positions.items() if v == tile]

    def to_scadnano_lattice(self) -> ScadnanoLattice:
        return self

    def asdict(self) -> dict[str, Any]:
        raise NotImplementedError

    @classmethod
    def fromdict(cls, d: dict[str, Any]):
        raise NotADirectoryError


AL = TypeVar("AL", bound="AbstractLattice")


def _skip_polyT_and_inertname(glue: Glue) -> bool:
    if "inert" in glue.ident():
        return True
    elif isinstance(glue, SSGlue):
        if frozenset(glue.sequence.base_str) == frozenset("T"):
            return True
    return False


class LatticeFactory:
    types: dict[str, Type[Lattice]]

    def __init__(self):
        self.types = {}

    def register(self, c: Type[Lattice], n: str = None):
        self.types[n if n is not None else c.__name__] = c

    def from_dict(self, d: dict[str, Any]) -> Lattice:
        if "type" in d:
            c = self.types[d["type"]]
            return c.fromdict({k: v for k, v in d.items() if k != "type"})
        else:
            raise ValueError


lattice_factory = LatticeFactory()

lattice_factory.register(AbstractLattice)
lattice_factory.register(ScadnanoLattice)
