from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import dataclasses
from typing import Any, Callable, Literal, Protocol, Type, TypeVar, cast
import numpy as np
import scadnano

from xgrow.tileset import Tile
from alhambra.glues import Glue, SSGlue
from alhambra.seq import Seq

from alhambra.tiles import D, TileSupportingScadnano, SupportsGuards

if False:
    from alhambra.tilesets import TileSet


class Lattice(ABC):
    @abstractmethod
    def __getitem__(self, index) -> str | Any:
        ...

    @abstractmethod
    def __setitem__(self, index, v):
        ...

    def asdict(self) -> dict[str, Any]:
        raise NotImplementedError

    @classmethod
    def fromdict(cls, d: dict[str, Any]) -> Lattice:
        raise NotADirectoryError


class LatticeSupportingScadnano(ABC):
    @abstractmethod
    def to_scadnano_lattice(self) -> ScadnanoLattice:
        ...

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
        return des


@dataclass
class ScadnanoLattice(LatticeSupportingScadnano, Lattice):
    positions: dict[tuple[int, int], str] = field(default_factory=lambda: {})

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


T_AL = TypeVar("T_AL", bound="Type[AbstractLattice]")


def _skip_polyT_and_inertname(glue: Glue) -> bool:
    if "inert" in glue.ident():
        return True
    elif isinstance(glue, SSGlue):
        if frozenset(glue.sequence.base_str) == frozenset("T"):
            return True
    return False


@dataclass(init=False)
class AbstractLattice(Lattice):
    grid: np.ndarray

    def __getitem__(self, index) -> str | Any:
        return self.grid[index]

    def __setitem__(self, index, v):
        self.grid[index] = v

    def __init__(self, v) -> None:
        if isinstance(v, AbstractLattice):
            self.grid = v.grid
        else:
            self.grid = np.array(v)

    def asdict(self) -> dict[str, Any]:
        d: dict[str, Any] = {}
        d["type"] = self.__class__.__name__
        d["grid"] = self.grid.tolist()
        return d

    @classmethod
    def fromdict(cls: Type[T_AL], d: dict[str, Any]) -> T_AL:
        return cls(np.array(d["grid"]))

    @classmethod
    def empty(cls, shape):
        return cls(np.full(shape, "", dtype=object))


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
