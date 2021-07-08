from __future__ import annotations
from dataclasses import dataclass
from alhambra.classes import Serializable
from typing import Any, Iterable, Optional, Type, TypeVar
from .newtile import Tile, TileList
from .glue import Glue, GlueList
from .newseed import Seed, seed_factory
import xgrow.tileset as xgt
import xgrow
import copy

T = TypeVar("T")


@dataclass(init=False)
class TileSet(Serializable):
    tiles: TileList
    glues: GlueList
    seed: Optional[Seed]
    params: dict

    def run_xgrow(self, *args, **kwargs):
        xgrow_tileset = self.to_xgrow()
        return xgrow.run(xgrow_tileset, *args, **kwargs)

    def to_xgrow(self, self_complementary_glues: bool = True) -> xgt.TileSet:

        tiles = [t.to_xgrow(self_complementary_glues) for t in self.tiles]

        # FIXME
        # bonds = [g.to_xgrow(self_complementary_glues) for g in self.glues]
        bonds = []

        if self.seed:
            seed_tiles, seed_bonds, initstate = self.seed.to_xgrow(
                self_complementary_glues
            )
        else:
            seed_tiles = []
            seed_bonds = []
            initstate = None

        xgrow_tileset = xgt.TileSet(
            seed_tiles + tiles, seed_bonds + bonds, initstate=initstate
        )

        return xgrow_tileset

    def to_xgrow_dict(self) -> dict:
        return self.to_xgrow().to_dict()

    def to_dict(self) -> dict:
        d = {}

        allglues = self.glues | self.tiles.glues_from_tiles()
        refglues = set(allglues.data.keys())  # FIXME

        if self.tiles:
            d["tiles"] = [t.to_dict(refglues=refglues) for t in self.tiles.aslist()]
        if allglues:
            d["glues"] = [g for g in allglues.aslist()]
        if self.seed:
            d["seed"] = self.seed.to_dict()
        if self.params:
            d["params"] = self.params.copy()
        return d

    @classmethod
    def from_dict(cls: Type[T], d: dict) -> T:
        ts = cls()
        ts.tiles = TileList(Tile.from_dict(x) for x in d.get("tiles", []))
        ts.glues = GlueList(Glue.from_dict(x) for x in d.get("glues", []))
        if "seed" in d:
            ts.seed = seed_factory.from_dict(d["seed"])
        if "params" in d:
            ts.params = copy.deepcopy(d["params"])
        return ts

    def _serialize(self) -> Any:
        return self.to_dict()

    @classmethod
    def _deserialize(cls, input: Any) -> TileSet:
        return cls.from_dict(input)

    def __init__(
        self,
        tiles: Iterable[Tile] = tuple(),
        glues: Iterable[Glue] = tuple(),
        seed: Seed | None = None,
        params: dict | None = None,
    ) -> None:
        self.tiles = TileList(tiles)
        self.glues = GlueList(glues)
        self.seed = seed
        if params is not None:
            self.params = params
        else:
            self.params = dict()
