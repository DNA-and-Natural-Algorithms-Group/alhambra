from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import dataclasses
from typing import Any, Literal, Protocol, Type, TypeVar, cast
import numpy as np
import scadnano

from xgrow.tileset import Tile

from alhambra.tiles import TileSupportingScadnano

if False:
    from alhambra.tilesets import TileSet


class Lattice(ABC):
    @abstractmethod
    def __getitem__(self, index) -> str | Any:
        ...

    @abstractmethod
    def __setitem__(self, index, v):
        ...


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


T_AL = TypeVar("T_AL", bound="Type[AbstractLattice]")


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
            self.grid = v

    @classmethod
    def empty(cls, shape):
        return cls(np.full(shape, "", dtype=object))
