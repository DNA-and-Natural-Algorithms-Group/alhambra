from __future__ import annotations
from abc import ABC, ABCMeta, abstractmethod
from typing import Any, Type, TypeVar, Union, cast
import xgrow.tileset as xgt
from .newtile import FlatishSingleTile9, Tile
from .glue import Glue


T = TypeVar("T")


class Seed(ABC):
    @abstractmethod
    def to_xgrow(
        self, self_complementary_glues=False
    ) -> tuple[list[xgt.Tile], list[xgt.Bond], xgt.InitState]:
        ...

    @abstractmethod
    def to_dict(self) -> dict:
        ...

    @classmethod
    @abstractmethod
    def from_dict(cls: Type[T], d: dict) -> T:
        ...


class FlatishSeed(Seed):
    ...


class SeedFactory:
    types: dict[str, Type[Seed]]

    def __init__(self):
        self.types = {}

    def register(self, c: Type[Seed], n: str = None):
        self.types[n if n is not None else c.__name__] = c

    def from_dict(self, d: dict[str, Any]) -> Seed:
        if "type" in d:
            c = self.types[d["type"]]
            return c.from_dict(**d)
        else:
            raise ValueError


class FlatishHSeed9(FlatishSeed):
    adapter_tiles: list[tuple[Glue | str, FlatishSingleTile9]]

    def __init__(self, adapter_tiles=[]):
        self.adapter_tiles = list(adapter_tiles)

    def to_dict(self, glues_as_refs=False) -> dict:
        d = {}
        d["adapter_tiles"] = [
            [str(g), t.to_dict()] for g, t in self.adapter_tiles  # FIXME
        ]
        d["type"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls: Type[T], d: dict) -> T:
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


class DXOrigamiSeed(Seed):
    ...


class DX_TallRect(DXOrigamiSeed):
    ...


seed_factory = SeedFactory()
for ttype in [FlatishHSeed9]:
    seed_factory.register(ttype)

seed_factory.register(DX_TallRect, "longrect")
seed_factory.register(DX_TallRect, "tallrect")
seed_factory.register(DX_TallRect)
