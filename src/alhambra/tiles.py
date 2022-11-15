from __future__ import annotations

import copy
import uuid
from abc import ABC, ABCMeta, abstractmethod
from dataclasses import dataclass
from enum import Enum
from typing import (
    Any,
    Callable,
    ClassVar,
    Collection,
    Generator,
    Generic,
    Iterable,
    Iterator,
    List,
    Literal,
    MutableMapping,
    Optional,
    Sequence,
    Tuple,
    Type,
    TypeVar,
    Union,
    cast,
    overload,
    TYPE_CHECKING,
)

from xgrow.xcolors import xcolors

from . import drawing

if TYPE_CHECKING:
    from alhambra.tilesets import XgrowGlueOpts, TileSet
    import xgrow.tileset as xgt

from .glues import Use

try:
    import scadnano
except ImportError:
    pass

from .classes import UpdateListD
from .glues import DXGlue, Glue, GlueFactory, GlueList, SSGlue
from .seq import Seq

Color = str

T = TypeVar("T")

import logging

log = logging.getLogger("alhambra")


class D(Enum):
    """Cardinal directions; also edge directions for single tiles."""

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


class EdgeView:
    """A class to ensure that tile edge glue manipulations are handled through the tile."""

    _edges: List[Glue]
    _tile: Tile
    __slots__ = ("_edges", "_tile")

    def __init__(self, _edges: List[Glue], _tile: Tile):
        self._edges = _edges
        self._tile = _tile

    @overload
    def __getitem__(self, k: int | str) -> Glue:
        ...

    @overload
    def __getitem__(self, k: slice) -> list[Glue]:
        ...

    def __getitem__(self, k: int | slice | str) -> Glue | list[Glue]:
        if isinstance(k, str):
            k = self._tile._get_edge_index(k)
        return self._edges.__getitem__(k)

    def __setitem__(self, k: int | str, v: Glue) -> None:
        if isinstance(k, str):
            k = self._tile._get_edge_index(k)
        self._tile.set_edge(k, v)

    def insert(self, index: int, value: Glue) -> None:
        return self._edges.insert(index, value)

    def __delitem__(self, k: int | str) -> None:
        raise NotImplementedError

    def __iter__(self) -> Iterator[Glue]:
        return iter(self._edges)

    def __next__(self) -> Generator[Glue, None, None]:
        for x in self._edges:
            yield x
        raise StopIteration

    def __len__(self) -> int:
        return self._edges.__len__()

    def __repr__(self) -> str:
        return self._edges.__repr__()

    def __str__(self) -> str:
        return self._edges.__str__()


class UseView(MutableMapping[int, Use]):
    """A class to ensure that tile edge use manipulations are handled through the tile."""

    _tile: Tile
    __slots__ = "_tile"

    def __init__(self, _tile: Tile):
        self._tile = _tile

    @overload
    def __getitem__(self, k: int) -> Use:
        ...

    @overload
    def __getitem__(self, k: slice) -> list[Use]:
        ...

    def __getitem__(self, k: int | slice) -> Use | list[Use]:
        x = self._tile._edges[k]
        if isinstance(x, Iterable):
            return [g.use for g in x]
        else:
            return x.use

    def __setitem__(self, k: int, v: Use) -> None:
        self._tile._edges[k].use = v

    def insert(self, index: int, value: Use) -> None:
        raise NotImplementedError

    def __len__(self) -> int:
        return self._tile._edges.__len__()

    def __repr__(self) -> str:
        return repr([g.use for g in self._tile._edges])

    def __str__(self) -> str:
        return str([g.use for g in self._tile._edges])


@dataclass(init=False)
class Tile:
    """Base class for a tile."""

    name: Optional[str]
    _edges: List[Glue]
    color: Optional[Color]
    stoic: Optional[float]
    note: Optional[str | dict[str, Any]]
    fake: bool
    uses: List[List[Use]]  # FIXME: should use list of tuples
    __slots__ = ("name", "_edges", "color", "stoic", "note", "fake", "uses")

    def get_concentration(self, base_concentration: float) -> float:
        if self.stoic is None:
            return base_concentration
        else:
            return base_concentration * self.stoic

    def get_stoic(self, base_concentration: float) -> float:
        """Returns the stoichiometric ratio of this tile to the base concentration."""
        if self.stoic is not None:
            return self.stoic
        else:
            return 1.0

    def __init__(
        self,
        edges: Optional[Iterable[Glue | str]] = None,
        name: Optional[str] = None,
        color: Optional[Color] = None,
        stoic: Optional[float] = None,
        note: Optional[str | dict[str, Any]] = None,
        use: Sequence[Use | int | str] | None = None,
        fake: bool = False,
        uses: Sequence[Sequence[Use | int | str]] | None = None,
    ) -> None:
        if edges is None:
            raise ValueError
        self._edges = [(g if isinstance(g, Glue) else Glue(g)) for g in edges]
        self.uses = []
        if use is not None:
            self.uses += [[Use.from_any(x) for x in use]]
        if uses is not None:
            self.uses += [[Use.from_any(x) for x in u] for u in uses]
        self.name = name
        self.color = color
        self.stoic = stoic
        self.note = note
        self.fake = fake

    @property
    def structure(self):
        return self.__class__.__name__

    @property
    def use(self) -> Sequence[Use]:
        if self.uses:
            # FIXME: combine uses here
            return self.uses[0]
        else:
            return [g.use for g in self._edges]

    @use.setter
    def use(self, uses: Iterable[Use]) -> None:
        # Fixme: make this actually set tile uses as well
        for glue, use in zip(self._edges, uses):
            glue.use = use

    @property
    def edges(self) -> EdgeView:
        return EdgeView(self._edges, self)

    @edges.setter
    def edges(self, e: Iterable[Glue | str]) -> None:
        self._edges = [g if isinstance(g, Glue) else Glue(g) for g in e]

    @property
    def edge_directions(self) -> List[D]:
        raise NotImplementedError

    @property
    def edge_locations(self) -> List[EdgeLoc]:
        raise NotImplementedError

    def set_edge(self, i: int, glue: Glue) -> None:
        self._edges[i] = glue

    def copy(self: T) -> T:
        return copy.deepcopy(self)

    @property
    def rotations(self) -> List[Tile]:
        raise NotImplementedError

    @classmethod
    def _get_edge_index(cls, v: str) -> int:
        raise NotImplementedError

    def ident(self) -> str:
        if self.name:
            return self.name
        else:
            raise ValueError

    @property
    def is_fake(self) -> bool:
        return self.fake

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
            for k in ["name", "edges", "color", "stoic", "note", "uses"]
            if (v := getattr(self, k)) is not None
        }
        if self.edges is not None:
            b["edges"] = [
                (x.name if x.name in refglues else x.to_dict()) for x in self.edges
            ]
            # fixme: deal with None
        return b

    def update_glues(self, gluedict: GlueList[Glue]) -> None:
        if self.edges is not None:
            self._edges = [gluedict.merge_glue(g) for g in self.edges]

    def update_glues_and_list(self, gluedict: GlueList[Glue]) -> None:
        if self.edges is not None:
            self._edges = [gluedict.merge_glue_and_update_list(g) for g in self.edges]

    @staticmethod
    def from_dict(d: dict[str, Any]) -> Tile:
        # FIXME: legacy conversion, changs input
        if "ends" in d:
            assert "edges" not in d
            d["edges"] = d["ends"]
            del d["ends"]
        if "extra" in d:
            if "type" in d:
                d["type"] += "_" + d.pop("extra")
            else:
                raise ValueError

        if "uses" in d:
            d["uses"] = [
                tuple(Use(x) if not isinstance(x, Use) else x for x in y)
                for y in d["uses"]
            ]

        return tile_factory.from_dict(d)

    def to_xgrow(self, gluenamemap: Callable[[str], str] = lambda x: x) -> xgt.Tile:
        import xgrow.tileset as xgt

        edges = [gluenamemap(g.ident()) for g in self._edges]
        return xgt.Tile(
            cast(list[Union[str, int]], edges),
            name=self.xgname,
            stoic=self.stoic,
            color=self.color,
        )

    def abstract_diagram(
        self,
        tileset: TileSet | None = None,
        draw_names: bool = True,
        draw_glues: bool = True,
    ) -> drawing.Group:
        raise NotImplementedError


class TileSupportingScadnano(ABC, Tile):
    """Abstract base class for a tile that allows export to scadnano."""

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

    @classmethod
    def _get_edge_index(cls, v: str) -> int:
        return "NESW".index(v)

    def abstract_diagram(
        self,
        tileset: TileSet | None = None,
        draw_names: bool = True,
        draw_glues: bool = True,
    ) -> drawing.Group:
        if (self.color is not None) and (self.color in xcolors):
            color = xcolors[self.color]
        elif (self.color is not None) and (self.color[0] == "#"):
            color = self.color
        else:
            color = "gray"

        g = drawing.Group(id=uuid.uuid4().hex)

        g.elements.append(drawing.Rectangle(0, 0, 10, 10, fill=color, stroke="black"))

        if draw_glues:
            gluetext_locs = [(5, 1, 0), (9, 5, 90), (5, 9, 0), (1, 5, -90)]
            for loc, glue in zip(gluetext_locs, self.edges):
                g.elements.append(
                    drawing.Text(
                        glue.ident(),
                        0.8,
                        loc[0],
                        loc[1],
                        text_anchor="middle",
                        transform=f"rotate({loc[2]},{loc[0]},{loc[1]})",
                    )
                )

        if self.name is not None and draw_names:
            nametext = drawing.Text(
                self.name, 1.2, 5, 5, text_anchor="middle", valign="center"
            )
            g.elements.append(nametext)

        return g


class VDupleTile(Tile):
    def to_xgrow(self, gluenamemap: Callable[[str], str] = lambda x: x) -> xgt.Tile:
        d = super().to_xgrow(gluenamemap)
        d.shape = "V"
        return d

    @classmethod
    def _get_edge_index(cls, v: str) -> int:
        return ["N", "NE", "SE", "S", "SW", "NW"].index(v)

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
        self,
        tileset: TileSet | None = None,
        draw_names: bool = True,
        draw_glues: bool = True,
    ) -> drawing.Group:
        if (self.color is not None) and (self.color in xcolors):
            color = xcolors[self.color]
        elif (self.color is not None) and self.color[0] == "#":
            color = self.color
        else:
            color = "gray"

        g = drawing.Group(id=uuid.uuid4().hex)

        g.append(drawing.Rectangle(0, 0, 10, 20, fill=color, stroke="black"))

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
                g.append(
                    drawing.Text(
                        glue.ident(),
                        0.8,
                        loc[0],
                        loc[1],
                        text_anchor="middle",
                        transform=f"rotate({loc[2]},{loc[0]},{loc[1]})",
                    )
                )

        if self.name is not None and draw_names:
            nametext = drawing.Text(
                self.name,
                1.2,
                5,
                10,
                transform="rotate(-90, 5, 10)",
                text_anchor="middle",
                valign="center",
            )
            g.append(nametext)

        return g


class HDupleTile(Tile):
    def to_xgrow(self, gluenamemap: Callable[[str], str] = lambda x: x) -> xgt.Tile:
        d = super().to_xgrow(gluenamemap)
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

    @classmethod
    def _get_edge_index(cls, v: str) -> int:
        return ["NW", "NE", "E", "SE", "SW", "W"].index(v)

    def abstract_diagram(
        self,
        tileset: TileSet | None = None,
        draw_names: bool = True,
        draw_glues: bool = True,
    ) -> drawing.Group:
        if (self.color is not None) and (self.color in xcolors):
            color = xcolors[self.color]
        elif (self.color is not None) and (self.color[0] == "#"):
            color = self.color
        else:
            color = "rgb(150,150,150)"

        g = drawing.Group(id=uuid.uuid4().hex)

        g.append(drawing.Rectangle(0, 0, 20, 10, fill=color, stroke="black"))

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
                g.append(
                    drawing.Text(
                        glue.ident(),
                        0.8,
                        loc[0],
                        loc[1],
                        text_anchor="middle",
                        transform=f"rotate({loc[2]},{loc[0]},{loc[1]})",
                    )
                )

        if self.name is not None and draw_names:
            nametext = drawing.Text(
                self.name,
                1.2,
                10,
                5,
                text_anchor="middle",
                valign="center",
            )
            g.append(nametext)

        return g


class SupportsGuards:
    @abstractmethod
    def create_guards(self, directions: Collection[str | D] = (D.E, D.S)) -> list[Glue]:
        ...


class BaseSSTile(SupportsGuards, TileSupportingScadnano):
    """Base class for single-stranded tiles."""

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

    def _input_neighborhood_domains(
        self,
    ) -> Sequence[Tuple[Sequence[str], Sequence[str]]]:
        raise NotImplementedError

    @property
    @abstractmethod
    def edge_directions(self) -> list[D]:
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
        domains: Optional[Sequence[SSGlue]] = None,
        note: Optional[str] = None,
        use: Optional[Sequence[Use]] = None,
        uses: Optional[Sequence[Sequence[Use]]] = None,
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
        if use is not None:
            for e, u in zip(self.edges, use, strict=True):
                e.use = u
        if uses is not None:
            self.uses = list(list(u) for u in uses)
        else:
            self.uses = []

    @property
    def sequence(self) -> Seq:
        return Seq("-".join(str(glue.sequence) for glue in self.domains))

    @sequence.setter
    def sequence(self, seq: Seq) -> None:
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
    def edges(self, edges: Sequence[Glue]) -> None:
        self._edges = [bd.merge(g) for bd, g in zip(self._base_edges, edges)]

    def set_edge(self, i: int, glue: Glue) -> None:
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


def _scadnano_color(color: Optional[str]) -> Optional[scadnano.Color]:
    if color is None:
        return None
    elif color[0] == "#":
        return scadnano.Color(hex_string=color)
    else:
        xc = xcolors[color]
        ci = [int(x) for x in xc[4:-1].split(",")]
        return scadnano.Color(*ci)


class BaseSSTSingle(SingleTile, BaseSSTile):
    """Base class for a standard-orientation SST single tile."""

    _edges: List[Glue]

    @property
    def domains(self) -> List[SSGlue]:
        e = self.edges
        return [e[i] for i in [1, 0, 3, 2]]  # type: ignore

    @property
    def _base_edges(self) -> List[SSGlue]:
        return [self._base_domains[i] for i in [1, 0, 3, 2]]

    def _input_neighborhood_domains(
        self,
    ) -> Sequence[Tuple[Sequence[str], Sequence[str]]]:
        if not self.uses:
            return []

        strands = []

        for use in self.uses:
            if use == [Use.INPUT, Use.OUTPUT, Use.OUTPUT, Use.INPUT]:  # NW
                strands.append(
                    (
                        ("W", "N"),
                        (
                            self.edges["W"].complement.ident(),
                            "algo_fake_spacer",
                            self.edges["N"].complement.ident(),
                        ),
                    )
                )
            elif use == [Use.INPUT, Use.INPUT, Use.OUTPUT, Use.OUTPUT]:  # NE
                strands.append(
                    (
                        ("N", "E"),
                        (
                            self.edges["N"].complement.ident(),
                            "algo_fake_spacer",
                            self.edges["E"].complement.ident(),
                        ),
                    )
                )
            elif use == [Use.OUTPUT, Use.OUTPUT, Use.INPUT, Use.INPUT]:  # SW
                strands.append(
                    (
                        ("S", "W"),
                        (
                            self.edges["S"].complement.ident(),
                            "algo_fake_spacer",
                            self.edges["W"].complement.ident(),
                        ),
                    )
                )
            elif use == [Use.OUTPUT, Use.INPUT, Use.INPUT, Use.OUTPUT]:  # SE
                strands.append(
                    (  # FIXME: needed to avoid pseudoknot
                        ("S", "E"),
                        (
                            self.edges["S"].complement.ident(),
                            "algo_fake_spacer",
                            self.edges["E"].complement.ident(),
                        ),
                    )
                )

        return strands

    def to_scadnano(
        self, design: scadnano.Design, helix: int, offset: int
    ) -> scadnano.Strand:
        s = design.draw_strand(helix, offset + 21)

        for e in self.domains[0:2]:
            s.move(-e.dna_length)
            s.with_domain_name(e.ident())
        s.cross(s.current_helix + 1)
        for e in self.domains[2:]:
            s.move(e.dna_length)
            s.with_domain_name(e.ident())

        if self.name is not None:
            s.with_name(self.name)

        if self.color is not None:
            s.with_color(_scadnano_color(self.color))

        # We generally don't want to assign complements here because (a) we will be assigning
        # sequences for every tile, and (b) in some cases, we will have intentional mismatches.
        s.with_sequence(self.sequence.base_str, assign_complement=False)
        return s.strand


Domain = SSGlue


@dataclass(init=False)
class BaseSSTSingleWithExtensions(BaseSSTSingle):
    mod5: Domain | None = None
    mod3: Domain | None = None

    def to_dict(self, refglues: set[str] = set()) -> dict[str, Any]:
        d = super().to_dict(refglues=refglues)
        d["type"] = self.__class__.__name__
        d["mod5"] = self.mod5.to_dict() if self.mod5 else None
        d["mod3"] = self.mod3.to_dict() if self.mod3 else None
        return d

    @property
    def domains(self) -> List[SSGlue]:
        d = super().domains
        if self.mod5 is not None:
            d.insert(0, self.mod5)
        if self.mod3 is not None:
            d.append(self.mod3)
        return d

    @property
    def _sequence_length(self) -> int:
        noext = sum(x.dna_length for x in self._base_domains)
        if self.mod5 is not None:
            noext += self.mod5.dna_length
        if self.mod3 is not None:
            noext += self.mod3.dna_length
        return noext

    @property
    def sequence(self) -> Seq:
        return Seq("-".join(str(glue.sequence) for glue in self.domains))

    @sequence.setter
    def sequence(self, seq: Seq) -> None:
        seq = Seq(seq)
        if seq.dna_length != self._sequence_length:
            raise ValueError

        pos = 0
        base_str = seq.base_str

        _base_domains = self._base_domains.copy()
        if self.mod5 is not None:
            _base_domains.insert(0, self.mod5)
        if self.mod3 is not None:
            _base_domains.append(self.mod3)

        # fixme: should we check whitespace?
        for base_domain, domain in zip(_base_domains, self.domains):
            base_domain.sequence | base_str[pos : pos + base_domain.dna_length]  # noqa
            domain.sequence = Seq(
                seq.base_str[pos : pos + base_domain.dna_length]
            )  # noqa
            pos += base_domain.dna_length

    def __init__(
        self,
        edges: Optional[list[Union[str, Glue]]] = None,
        name: Optional[str] = None,
        color: Optional[Color] = None,
        stoic: Optional[float] = None,
        sequence: Optional[Seq] = None,
        domains: Optional[list[Domain]] = None,
        note: Optional[str] = None,
        use: Optional[Sequence[Use]] = None,
        uses: Optional[Sequence[Sequence[Use]]] = None,
        mod5: Seq | str | int | dict | None = None,
        mod3: Seq | str | int | dict | None = None,
    ):
        # Don't deal with the sequence just yet.
        super().__init__(edges, name, color, stoic, None, domains, note, use, uses)

        # Now add the mod5 and mod3 domains, if they exist:
        if isinstance(mod5, (str, Seq)):
            self.mod5 = Domain(sequence=mod5)
        elif isinstance(mod5, int):
            self.mod5 = Domain(length=mod5)
        elif isinstance(mod5, dict):
            g = Glue.from_dict(mod5)
            assert isinstance(g, SSGlue)
            self.mod5 = g
        elif mod5 is None:
            self.mod5 = None
        else:
            raise TypeError(f"mod5 must be str, Seq, int, or None, not {type(mod5)}")

        if isinstance(mod3, (str, Seq)):
            self.mod3 = Domain(sequence=mod3)
        elif isinstance(mod3, int):
            self.mod3 = Domain(length=mod3)
        elif isinstance(mod3, dict):
            g = Glue.from_dict(mod3)
            assert isinstance(g, SSGlue)
            self.mod3 = g
        elif mod3 is None:
            self.mod3 = None
        else:
            raise TypeError(f"mod3 must be str, Seq, int, or None, not {type(mod3)}")

        if sequence is not None:
            self.sequence = sequence


class SST10_5S(BaseSSTSingle):
    "Single SST, with domains (5'→3') of 11, 10, 10, and 11 nt. North edge is 10nt. 5' is S, 3' is E."
    _base_domains: ClassVar[list[SSGlue]] = [SSGlue(length=x) for x in [11, 10, 10, 11]]
    _scadnano_offsets = ((-1, -11), (-1, 10), (1, 10), (1, -11))
    _scadnano_5p_offset = (1, 21)

    @property
    def domains(self) -> List[SSGlue]:
        e = self.edges
        return [e[i] for i in [2, 3, 0, 1]]  # type: ignore

    @property
    def _base_edges(self) -> List[SSGlue]:
        return [self._base_domains[i] for i in [2, 3, 0, 1]]


class SST11_5S(BaseSSTSingle):
    "Single SST, with domains (5'→3') of 10, 11, 11, and 10 nt. North edge is 11nt. 5' is S, 3' is E."
    _base_domains: ClassVar[list[SSGlue]] = [SSGlue(length=x) for x in [10, 11, 11, 10]]
    _scadnano_offsets = ((-1, -10), (-1, 11), (1, 11), (1, -10))
    _scadnano_5p_offset = (1, 21)

    @property
    def domains(self) -> List[SSGlue]:
        e = self.edges
        return [e[i] for i in [2, 3, 0, 1]]  # type: ignore

    @property
    def _base_edges(self) -> List[SSGlue]:
        return [self._base_domains[i] for i in [2, 3, 0, 1]]


class SST10(BaseSSTSingle):
    "Single SST, with domains (5'→3') of 11, 10, 10, and 11 nt. North edge is 10nt. 5' is E, 3' is S."
    _base_domains: ClassVar[list[SSGlue]] = [SSGlue(length=x) for x in [11, 10, 10, 11]]
    _scadnano_offsets = ((-1, -11), (-1, 10), (1, 10), (1, -11))
    _scadnano_5p_offset = (1, 21)


class SST11(BaseSSTSingle):
    "Single SST, with domains (5'→3') of 10, 11, 11, and 10 nt. North edge is 11nt. 5' is E, 3' is S."
    _base_domains: ClassVar[list[SSGlue]] = [SSGlue(length=x) for x in [10, 11, 11, 10]]
    _scadnano_offsets = ((-1, -10), (-1, 11), (1, 11), (1, -10))
    _scadnano_5p_offset = (1, 21)


class TileFactory:
    types: dict[str, Type[Tile]]

    def __init__(self):
        self.types = {}

    def register(self, c: Type[Tile], n: Optional[str] = None) -> None:
        self.types[n if n is not None else c.__name__] = c

    def from_dict(self, d: dict[str, Any]) -> Tile:
        if "ends" in d:  # old ends format
            d["edges"] = d.pop("ends")
        if "fullseqs" in d:  # old fullseqs format
            d["sequences"] = d.pop("fullseqs")
        if "input" in d:  # old input format
            d["uses"] = [[3 - x for x in d.pop("input")]]  # FIXME
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
for ttype in [Tile, SST10, SST10_5S, SST11, SST11_5S]:
    tile_factory.register(ttype)

tile_factory.register(SST10, "SST10_5E")
tile_factory.register(SST11, "SST11_5E")

SomeTile = TypeVar("SomeTile", bound=Tile)


class TileList(Generic[SomeTile], UpdateListD[SomeTile]):
    def glues_from_tiles(self) -> GlueList:
        gl: GlueList[Glue] = GlueList()
        for tile in self:
            gl |= tile.edges
        return gl

    def domains_from_tiles(self) -> GlueList:
        gl: GlueList[Glue] = GlueList()
        for tile in self:
            gl |= tile.domains
        return gl


class DAOETile(Tile):
    _edges: List[Glue]  # actually dxglue
    _strand_sequences: List[str] | None

    def to_dict(self, refglues: set[str] = set()) -> dict[str, Any]:
        d = super().to_dict(refglues=refglues)
        d["type"] = self.__class__.__name__
        return d

    def __init__(
        self,
        name: str,
        sequence: str | None = None,
        sequences: list[str] | None = None,
        edges: List[Glue | str] | None = None,
        label: str | None = None,
        color: str | None = None,
        **kwargs,
    ):
        if sequence:
            assert not sequences
            self._strand_sequences = sequence.split("+")
        elif sequences:
            self._strand_sequences = sequences
        else:
            self._strand_sequences = None

        if label:
            log.warning("label not currently handled")

        super().__init__(name=name, edges=edges, color=color, **kwargs)


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
tile_factory.register(DAOESingle5p, "tile_daoe_5up_2h")


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
tile_factory.register(DAOESingle3p, "tile_daoe_3up_2h")


class DAOEHDouble3p(HDupleTile, DAOETile):
    ...


tile_factory.register(DAOEHDouble3p)
tile_factory.register(DAOEHDouble3p, "tile_daoe_doublehoriz_35up")
tile_factory.register(DAOEHDouble3p, "tile_daoe_doublehoriz_35up_2h3h")  # FIXME
tile_factory.register(DAOEHDouble3p, "tile_daoe_doublehoriz_35up_1h2i")  # FIXME


class DAOEHDouble5p(HDupleTile, DAOETile):
    ...


class DAOEVDouble3p(VDupleTile, DAOETile):
    ...


tile_factory.register(DAOEVDouble3p)
tile_factory.register(DAOEVDouble3p, "tile_daoe_doublevert_35up")
tile_factory.register(DAOEVDouble3p, "tile_daoe_doublevert_35up_4h5h")  # FIXME


class DAOEVDouble5p(VDupleTile, DAOETile):
    ...
