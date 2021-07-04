from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Any, Type, TypeVar
import xgrow.tileset as xgt


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


class DXOrigamiSeed(Seed):
    ...


class DX_TallRect(DXOrigamiSeed):
    ...


seed_factory = SeedFactory()

# seed_factory.register(DX_TallRect, "longrect")
# seed_factory.register(DX_TallRect, "tallrect")
# seed_factory.register(DX_TallRect)
