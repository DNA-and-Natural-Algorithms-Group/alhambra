from __future__ import annotations
from abc import ABC, ABCMeta, abstractmethod
import copy
from dataclasses import dataclass
import uuid
import drawSvg_svgy as draw
from typing import (
    Any,
    Collection,
    Generic,
    Iterable,
    List,
    Literal,
    MutableSequence,
    Optional,
    ClassVar,
    Sequence,
    SupportsIndex,
    Tuple,
    Type,
    TypeVar,
    Union,
    cast,
    overload,
)
from enum import Enum
from xgrow.xcolors import xcolors
from . import mixes

try:
    import scadnano
except ImportError:
    pass
from .seq import Seq
from .glues import DXGlue, Glue, GlueFactory, SSGlue, GlueList
from .classes import UpdateListD
import xgrow.tileset as xgt

Color = str

T = TypeVar("T")

import logging

log = logging.getLogger("alhambra")


class D(Enum):
    N = 0
    E = 1
    S = 2
    W = 3

    @property
    def complement(self) -> D:
        return D((self.value + 2) % 4)


Offset = Tuple[int, int]


@dataclass(frozen=True)
class EdgeLoc:
    direction: D
    position: Offset = (0, 0)


EL = EdgeLoc


class EdgeView(MutableSequence[Glue]):
    _edges: List[Glue]
    _tile: "Tile"
    __slots__ = ("_edges", "_tile")

    def __init__(self, _edges: List[Glue], _tile: Tile):
        self._edges = _edges
        self._tile = _tile

    @overload
    def __getitem__(self, k: int) -> Glue:
        ...

    @overload
    def __getitem__(self, k: slice) -> list[Glue]:
        ...

    def __getitem__(self, k: int | slice) -> Glue | list[Glue]:
        return self._edges.__getitem__(k)

    def __setitem__(self, k: int, v: Glue) -> None:
        self._tile.set_edge(k, v)

    def insert(self, index: int, value: Glue) -> None:
        return self._edges.insert(index, value)

    def __delitem__(self, k: int) -> None:
        self._edges.__delitem__(k)

    def __len__(self) -> int:
        return self._edges.__len__()

    def __repr__(self) -> str:
        return self._edges.__repr__()

    def __str__(self) -> str:
        return self._edges.__str__()


@dataclass(init=False)
class Tile:
    name: Optional[str]
    _edges: List[Glue]
    color: Optional[Color]
    stoic: Optional[float]
    note: Optional[str | dict[str, Any]]
    __slots__ = ("name", "_edges", "color", "stoic", "note")

    def __init__(
        self,
        edges: Optional[Iterable[Glue | str]] = None,
        name: Optional[str] = None,
        color: Optional[Color] = None,
        stoic: Optional[float] = None,
        note: Optional[str | dict[str, Any]] = None,
    ) -> None:
        if edges is None:
            raise ValueError
        self._edges = [(g if isinstance(g, Glue) else Glue(g)) for g in edges]
        self.name = name
        self.color = color
        self.stoic = stoic
        self.note = note

    @property
    def structure(self):
        return self.__class__.__name__

    @property
    def edges(self):
        return EdgeView(self._edges, self)

    @edges.setter
    def edges(self, e):
        self._edges = e

    @property
    def edge_directions(self) -> List[D]:
        raise NotImplementedError

    @property
    def edge_locations(self) -> List[EdgeLoc]:
        raise NotImplementedError

    def set_edge(self, i: int, glue: Glue):
        self._edges[i] = glue

    def copy(self: T) -> T:
        return copy.deepcopy(self)

    @property
    def rotations(self) -> List[Tile]:
        raise NotImplementedError

    def ident(self) -> str:
        if self.name:
            return self.name
        else:
            raise ValueError

    @property
    def is_fake(self) -> bool:
        return False

    def merge(self, other) -> Tile:
        if self == other:
            return self
        raise NotImplementedError

    @property
    def xgname(self) -> str:
        if self.name:
            return self.name
        if self._edges:
            return f"{self.__class__.__name__}_" + "_".join(
                f"{g.ident()}" for g in self._edges
            )
        else:
            raise ValueError

    def to_dict(self, refglues: set[str] = set()) -> dict[str, Any]:
        b = {
            k: v
            for k in ["name", "edges", "color", "stoic", "note"]
            if (v := getattr(self, k)) is not None
        }
        if self.edges is not None:
            b["edges"] = [
                (x.name if x.name in refglues else x.to_dict()) for x in self.edges
            ]
            # fixme: deal with None
        return b

    def update_glues(self, gluedict: GlueList):
        if self.edges is not None:
            self.edges = [gluedict.merge_glue(g) for g in self.edges]

    def update_glues_and_list(self, gluedict: GlueList):
        if self.edges is not None:
            self.edges = [gluedict.merge_glue_and_update_list(g) for g in self.edges]

    @staticmethod
    def from_dict(d: dict[str, Any]) -> Tile:
        # FIXME: legacy conversion, changs input
        if "ends" in d:
            assert "edges" not in d
            d["edges"] = d["ends"]
            del d["ends"]
        if "extra" in d:
            if "type" in d:
                d["type"] += "_" + d["extra"]
            else:
                raise ValueError

        return tile_factory.from_dict(d)

    def to_xgrow(self, self_complementary_glues: bool = False) -> xgt.Tile:
        if self._edges and not self_complementary_glues:
            edges = [g.ident() for g in self._edges]
        elif self._edges:
            edges = [g.basename() for g in self._edges]
        else:
            raise ValueError
        return xgt.Tile(
            cast(list[Union[str, int]], edges),
            name=self.xgname,
            stoic=self.stoic,
            color=self.color,
        )

    def abstract_diagram(
        self, tileset=None, draw_names: bool = True, draw_glues: bool = True
    ) -> draw.Group:
        raise NotImplementedError


class TileSupportingScadnano(ABC, Tile):
    @property
    @abstractmethod
    def _scadnano_5p_offset(self) -> tuple[int, int]:
        ...

    @abstractmethod
    def to_scadnano(
        self, design: scadnano.Design, helix: int, offset: int
    ) -> scadnano.Strand:
        ...

    @abstractmethod
    def __init__(
        self,
        edges: Optional[Iterable[Glue | str]] = None,
        name: Optional[str] = None,
        color: Optional[Color] = None,
        stoic: Optional[float] = None,
        note: Optional[str] = None,
        domains: Optional[Iterable[SSGlue]] = None,
    ) -> None:
        ...


class SingleTile(Tile):
    """A tile with N, E, S and W edges."""

    @property
    def edge_directions(self) -> List[D]:
        return [D.N, D.E, D.S, D.W]

    @property
    def edge_locations(self) -> List[EdgeLoc]:
        return [
            EdgeLoc(D.N, (0, 0)),
            EdgeLoc(D.E, (0, 0)),
            EdgeLoc(D.S, (0, 0)),
            EdgeLoc(D.W, (0, 0)),
        ]

    def abstract_diagram(
        self, tileset=None, draw_names: bool = True, draw_glues: bool = True
    ) -> draw.Group:
        if (self.color is not None) and (self.color in xcolors):
            color = xcolors[self.color]
        else:
            color = "gray"
        box = draw.Rectangle(0, 0, 10, 10, fill=color, stroke="black")

        elems: list[draw.DrawingBasicElement] = [box]

        if draw_glues:
            gluetext_locs = [(5, 1, 0), (9, 5, 90), (5, 9, 0), (1, 5, -90)]
            for loc, glue in zip(gluetext_locs, self.edges):
                elems.append(
                    draw.Text(
                        glue.ident(),
                        0.8,
                        loc[0],
                        loc[1],
                        center=True,
                        transform=f"rotate({loc[2]},{loc[0]},{loc[1]})",
                    )
                )

        if self.name is not None and draw_names:
            nametext = draw.Text(self.name, 1.2, 5, 5, center=True, valign="center")
            elems.append(nametext)

        return draw.Group(elems, id=uuid.uuid4().hex)


class VDupleTile(Tile):
    def to_xgrow(self, self_complementary_glues: bool = False) -> xgt.Tile:
        d = super().to_xgrow(self_complementary_glues)
        d.shape = "V"
        return d

    @property
    def edge_directions(self) -> List[D]:
        return [D.N, D.E, D.E, D.S, D.W, D.W]

    @property
    def edge_locations(self) -> List[EdgeLoc]:
        return [
            EdgeLoc(D.N, (0, 0)),
            EdgeLoc(D.E, (0, 0)),
            EdgeLoc(D.E, (1, 0)),
            EdgeLoc(D.S, (1, 0)),
            EdgeLoc(D.W, (1, 0)),
            EdgeLoc(D.W, (0, 0)),
        ]

    def abstract_diagram(
        self, tileset=None, draw_names: bool = True, draw_glues: bool = True
    ) -> draw.Group:
        if (self.color is not None) and (self.color in xcolors):
            color = xcolors[self.color]
        else:
            color = "gray"
        box = draw.Rectangle(0, 0, 10, 20, fill=color, stroke="black")

        elems: list[draw.DrawingBasicElement] = [box]

        if draw_glues:
            gluetext_locs = [
                (5, 1, 0),
                (9, 5, 90),
                (9, 15, 90),
                (5, 19, 0),
                (1, 15, -90),
                (1, 5, -90),
            ]
            for loc, glue in zip(gluetext_locs, self.edges):
                elems.append(
                    draw.Text(
                        glue.ident(),
                        0.8,
                        loc[0],
                        loc[1],
                        center=True,
                        transform=f"rotate({loc[2]},{loc[0]},{loc[1]})",
                    )
                )

        if self.name is not None and draw_names:
            nametext = draw.Text(
                self.name,
                1.2,
                5,
                10,
                transform="rotate(-90, 5, 10)",
                center=True,
                valign="center",
            )
            elems.append(nametext)

        return draw.Group(elems, id=uuid.uuid4().hex)


class HDupleTile(Tile):
    def to_xgrow(self, self_complementary_glues: bool = False) -> xgt.Tile:
        d = super().to_xgrow(self_complementary_glues)
        d.shape = "H"
        return d

    @property
    def edge_directions(self) -> List[D]:
        return [D.N, D.N, D.E, D.S, D.S, D.W]

    @property
    def edge_locations(self) -> List[EdgeLoc]:
        return [
            EdgeLoc(D.N, (0, 0)),
            EdgeLoc(D.N, (0, 1)),
            EdgeLoc(D.E, (0, 1)),
            EdgeLoc(D.S, (0, 1)),
            EdgeLoc(D.S, (0, 0)),
            EdgeLoc(D.W, (0, 0)),
        ]

    def abstract_diagram(
        self, tileset=None, draw_names: bool = True, draw_glues: bool = True
    ) -> draw.Group:
        if (self.color is not None) and (self.color in xcolors):
            color = xcolors[self.color]
        else:
            color = "rgb(150,150,150)"
        box = draw.Rectangle(0, 0, 20, 10, fill=color, stroke="black")

        elems: list[draw.DrawingBasicElement] = [box]

        if draw_glues:
            gluetext_locs = [
                (5, 1, 0),
                (10, 1, 0),
                (19, 5, 90),
                (15, 9, 0),
                (5, 9, 0),
                (1, 5, -90),
            ]
            for loc, glue in zip(gluetext_locs, self.edges):
                elems.append(
                    draw.Text(
                        glue.ident(),
                        0.8,
                        loc[0],
                        loc[1],
                        center=True,
                        transform=f"rotate({loc[2]},{loc[0]},{loc[1]})",
                    )
                )

        if self.name is not None and draw_names:
            nametext = draw.Text(
                self.name,
                1.2,
                10,
                5,
                center=True,
                valign="center",
            )
            elems.append(nametext)

        return draw.Group(elems, id=uuid.uuid4().hex)


class SupportsGuards:
    @abstractmethod
    def create_guards(self, directions: Collection[str | D] = (D.E, D.S)) -> list[Glue]:
        ...


class BaseSSTile(SupportsGuards, TileSupportingScadnano):
    _edges: List[Glue]  # actually SSGlue

    def to_dict(self, refglues: set[str] = set()) -> dict[str, Any]:
        d = super().to_dict(refglues=refglues)
        d["type"] = self.__class__.__name__
        d["sequence"] = self.sequence.seq_str
        return d

    @property
    @abstractmethod
    def _base_domains(self) -> List[SSGlue]:
        ...

    @property
    @abstractmethod
    def _base_edges(self) -> List[SSGlue]:
        ...

    @property
    @abstractmethod
    def domains(self) -> List[SSGlue]:
        ...

    @property
    @abstractmethod
    def edge_directions(self) -> list[str]:
        ...

    @property
    def _sequence_length(self) -> int:
        return sum(x.dna_length for x in self._base_domains)

    def __init__(
        self,
        edges: Optional[List[Union[str, Glue]]] = None,
        name: Optional[str] = None,
        color: Optional[Color] = None,
        stoic: Optional[float] = None,
        sequence: Optional[Seq] = None,
        domains: Optional[List[SSGlue]] = None,
        note: Optional[str] = None,
    ):
        Tile.__init__(self, edges=[], name=name, color=color, stoic=stoic, note=note)
        if edges is None and sequence is None and domains is None:
            raise ValueError
        if edges is not None:
            self._edges = [bd.merge(g) for bd, g in zip(self._base_edges, edges)]
        else:
            self._edges = [bd.copy() for bd in self._base_edges]
        if sequence is not None:
            self.sequence | sequence
            self.sequence = sequence
        elif domains is not None:
            if len(self.domains) != len(domains):
                raise ValueError
            for td, nd in zip(self.domains, domains):
                td |= nd

    @property
    def sequence(self) -> Seq:
        return Seq("-".join(str(glue.sequence) for glue in self.domains))

    @sequence.setter
    def sequence(self, seq: Seq):
        seq = Seq(seq)
        if seq.dna_length != self._sequence_length:
            raise ValueError

        pos = 0
        base_str = seq.base_str

        # fixme: should we check whitespace?
        for base_domain, domain in zip(self._base_domains, self.domains):
            base_domain.sequence | base_str[pos : pos + base_domain.dna_length]  # noqa
            domain.sequence = Seq(
                seq.base_str[pos : pos + base_domain.dna_length]
            )  # noqa
            pos += base_domain.dna_length

    @property
    def edges(self):
        return super().edges  # type: ignore

    @edges.setter
    def edges(self, edges: Sequence[Glue]):
        self._edges = [bd.merge(g) for bd, g in zip(self._base_edges, edges)]

    def set_edge(self, i: int, glue: Glue):
        self._edges[i] = self._base_edges[i].merge(glue)

    def __repr__(self) -> str:
        s: list[str] = []
        s.append(self.__class__.__name__)
        s.append("(")
        if self.name is not None:
            s.append(f"name={repr(self.name)}, ")
        s.append("edges=[" + ", ".join(repr(x) for x in self.edges) + "]")
        s.append(")")
        return "".join(s)

    def __str__(self) -> str:
        s = [f"<{self.__class__.__name__}: "]
        if self.name is not None:
            s.append(self.name + " | ")
        s.append(", ".join(g.ident() for g in self.edges))
        s.append(">")
        return "".join(s)

    def create_guards(self, directions: Collection[str | D] = (D.E, D.S)) -> list[Glue]:
        guards: list[Glue] = []
        directions = set(x if isinstance(x, D) else D[x] for x in directions)
        for ei, d in enumerate(self.edge_directions):
            if d not in directions:
                continue
            else:
                guards.append(self.edges[ei].complement)
        return guards


class BaseSSTSingle(SingleTile, BaseSSTile):
    _edges: List[Glue]

    @property
    def domains(self) -> List[SSGlue]:
        e = self.edges
        return [e[i] for i in [1, 0, 3, 2]]  # type: ignore

    @property
    def _base_edges(self) -> List[SSGlue]:
        return [self._base_domains[i] for i in [1, 0, 3, 2]]

    def to_scadnano(
        self, design: scadnano.Design, helix: int, offset: int
    ) -> scadnano.Strand:
        s = design.strand(helix, offset + 21)

        for e in self.domains[0:2]:
            s.move(-e.dna_length)
            if e.name is not None:
                s.with_domain_name(e.name)
        s.cross(s.current_helix + 1)
        for e in self.domains[2:]:
            s.move(e.dna_length)
            if e.name is not None:
                s.with_domain_name(e.name)

        if self.name is not None:
            s.with_name(self.name)

        s.with_sequence(self.sequence.base_str)
        return s.strand


class TileFactory:
    types: dict[str, Type[Tile]]

    def __init__(self):
        self.types = {}

    def register(self, c: Type[Tile], n: str = None):
        self.types[n if n is not None else c.__name__] = c

    def from_dict(self, d: dict[str, Any]) -> Tile:
        if "edges" in d:
            for i in range(0, len(d["edges"])):
                glue = d["edges"][i]
                if isinstance(glue, dict):
                    glue = Glue.from_dict(glue)
                d["edges"][i] = glue
        if "type" in d:
            c = self.types[d["type"]]
            del d["type"]

            return c(**d)
        elif "structure" in d:
            c = self.types[d["structure"]]
            del d["structure"]
            return c(**d)
        else:
            return Tile(**d)

    @overload
    def from_scadnano(
        self,
        d: "scadnano.Strand" | Iterable["scadnano.Strand"],
        return_position: Literal[True],
    ) -> tuple[TileSupportingScadnano, tuple[int, int]]:
        ...

    @overload
    def from_scadnano(
        self,
        d: "scadnano.Strand" | Iterable["scadnano.Strand"],
        return_position: Literal[False],
    ) -> TileSupportingScadnano:
        ...

    def from_scadnano(
        self,
        d: "scadnano.Strand" | Iterable["scadnano.Strand"],
        return_position: bool = False,
    ) -> tuple[TileSupportingScadnano, tuple[int, int]] | TileSupportingScadnano:
        if isinstance(d, Iterable):
            raise NotImplementedError

        name = d.name
        domain_names = [x.name for x in d.domains]
        domain_seqs = [x.dna_sequence() for x in d.domains]
        domains = [
            SSGlue(name=n, sequence=s) for n, s in zip(domain_names, domain_seqs)
        ]

        t = None
        for tiletype in self.types.values():
            try:
                if issubclass(tiletype, TileSupportingScadnano):
                    t = tiletype(name=name, domains=domains)  # type: ignore
                    break
            except ValueError:
                continue
        if t:
            if not return_position:
                return t
            else:
                domain = d.first_domain()
                return t, (
                    domain.helix - t._scadnano_5p_offset[0],
                    domain.offset_5p() - t._scadnano_5p_offset[1] + 1,
                )
        else:
            raise ValueError


tile_factory = TileFactory()
for ttype in [Tile]:
    tile_factory.register(ttype)


SomeTile = TypeVar("SomeTile", bound=Tile)


class TileList(Generic[SomeTile], UpdateListD[SomeTile]):
    def glues_from_tiles(self) -> GlueList:
        gl = GlueList()
        for tile in self:
            gl |= tile.edges
        return gl

    def apply_mix(self, mix: mixes.Mix, base_conc):

        newlist = TileList()

        for comp, conc in mix.all_comps().items():
            try:
                tile = self[comp]

                new_tile = tile.copy()
                new_tile.stoic = float(conc / base_conc)

                newlist.add(new_tile)
            except KeyError:
                log.warn(f"Component {comp} not found in tile list.")

        if len(newlist) == 0:
            raise ValueError("No mix components match tiles.")

        return newlist


class DAOETile(Tile):
    _edges: List[Glue]  # actually dxglue

    def to_dict(self, refglues: set[str] = set()) -> dict[str, Any]:
        d = super().to_dict(refglues=refglues)
        d["type"] = self.__class__.__name__
        return d


class DAOESingle(SingleTile, DAOETile, metaclass=ABCMeta):
    # @property
    # @abstractmethod
    # def _endlocs(self) -> list[tuple[int, slice]]:
    #    ...

    # @property
    # @abstractmethod
    # def _baseglues(self) -> List[DXGlue]:
    #    ...
    ...


class DAOESingle5p(DAOESingle):
    _baseglues: ClassVar[List[DXGlue]] = [
        DXGlue("TD", length=5),
        DXGlue("TD", length=5),
        DXGlue("DT", length=5),
        DXGlue("DT", length=5),
    ]
    _gluelocs = [
        (0, slice(0, 5)),
        (3, slice(0, 5)),
        (3, slice(21, None)),
        (0, slice(21, None)),
    ]


tile_factory.register(DAOESingle5p, "tile_daoe_5up")


class DAOESingle3p(DAOESingle):
    _baseglues: ClassVar[List[DXGlue]] = [
        DXGlue("DT", length=5),
        DXGlue("DT", length=5),
        DXGlue("TD", length=5),
        DXGlue("TD", length=5),
    ]
    _gluelocs = [
        (0, slice(21, None)),
        (3, slice(21, None)),
        (3, slice(0, 5)),
        (0, slice(0, 5)),
    ]


tile_factory.register(DAOESingle3p, "tile_daoe_3up")


class DAOEHDouble3p(HDupleTile, DAOETile):
    ...


tile_factory.register(DAOEHDouble3p)
tile_factory.register(DAOEHDouble3p, "tile_daoe_doublehoriz_35up")


class DAOEHDouble5p(HDupleTile, DAOETile):
    ...


class DAOEVDouble3p(VDupleTile, DAOETile):
    ...


tile_factory.register(DAOEVDouble3p)
tile_factory.register(DAOEVDouble3p, "tile_daoe_doublevert_35up")


class DAOEVDouble5p(VDupleTile, DAOETile):
    ...
