# type: ignore

from __future__ import annotations

from copy import copy
from dataclasses import dataclass

# from turtle import st
from typing import (
    Any,
    Iterable,
    Literal,
    Protocol,
    Sequence,
    Type,
    TypeVar,
    cast,
    overload,
    TYPE_CHECKING,
)

if TYPE_CHECKING:
    from .tilesets import TileSet

import numpy as np

if False:
    from .tilesets import TileSet

import logging
from random import shuffle

import numpy.typing as npt

# from ..alhambra_old.alhambra_old import fastlatticedefect as fld
# from ..alhambra_old.alhambra_old import util
from .glues import Glue, GlueList, Use
from .tiles import HDupleTile, SingleTile, Tile, TileList, VDupleTile

TMap = dict[tuple[int, int], set[tuple[int, int, int]]]


@dataclass
class FTile:
    name: str
    used: bool
    glues: npt.NDArray[np.int64]
    use: npt.NDArray[np.int64]
    color: Any
    structure: Type[Tile]
    dfake: int
    sfake: int


@dataclass
class FTilesArray:
    color: Any
    use: np.ndarray[Any, np.dtype[np.int64]]
    glues: np.ndarray[Any, np.dtype[np.int64]]
    name: np.ndarray[Any, np.dtype[np.str0]]
    used: np.ndarray[Any, np.dtype[np.bool8]]
    structure: np.ndarray[Any, np.dtype[np.str0]]
    dfake: np.ndarray[Any, np.dtype[np.int64]]
    sfake: np.ndarray[Any, np.dtype[np.int64]]

    def __len__(self) -> int:
        return len(self.name)

    @classmethod
    def from_ftiles(cls, ftiles: Iterable[FTile]) -> FTilesArray:

        color = np.array([x.color for x in ftiles])
        name = np.array([x.name for x in ftiles])
        glues = np.array([x.glues for x in ftiles])
        use = np.array([x.use for x in ftiles])
        used = np.array([x.used for x in ftiles], dtype=bool)
        structure = np.array([x.structure for x in ftiles])
        dfake = np.array([x.dfake for x in ftiles])
        sfake = np.array([x.sfake for x in ftiles])

        return FTilesArray(
            color=color,
            name=name,
            glues=glues,
            use=use,
            used=used,
            structure=structure,
            dfake=dfake,
            sfake=sfake,
        )


TAU = 2

other = [2, 3, 0, 1]


class TileWithUse(Protocol):
    @property
    def use(self) -> list[int]:
        ...


# Use is now going to have 4 possible values: Null, Input, Output, Both, Permanent
uU = 0
uN = 1
uI = 2
uO = 3
uB = 4
uP = 5
# null, input, output, both, permanent
invertuse = [0, 1, 3, 2, 4, 5]
usedict = {
    "U": uU,
    "N": uN,
    "I": uI,
    "O": uO,
    "B": uB,
    "u": uU,
    "n": uN,
    "i": uI,
    "o": uO,
    "b": uB,
}

Equiv = npt.NDArray[np.int64]


class FGlueList:
    def __init__(self, glues: GlueList[Glue]):
        name = []
        strength = []
        type = []
        complement = []
        use = []
        self.tonum = {}

        # Ensure every glue has a complement
        newglues = GlueList()
        for g in glues:
            if g.complement.ident() not in glues.data.keys():
                newglues.add(g.complement)
        glues.update(newglues)

        for i, g in enumerate(glues):
            g = cast(Glue, g)
            name.append(g.ident())
            type.append(g.type)
            self.tonum.update({g.ident(): i})
            strength.append(g.abstractstrength if g.abstractstrength is not None else 1)
            if g.use is None:
                use.append(uB)
            else:
                use.append(g.use.value)

        for gn in name:
            complement.append(name.index(gn + "*" if gn[-1] != "*" else gn[:-1]))

        self.name = np.array(name)
        self.strength = np.array(strength)
        self.type = np.array(type)
        self.complement = np.array(complement)
        self.use = np.array(use)

    def blankequiv(self) -> Equiv:
        return np.arange(0, len(self.name))

    def iseq(self, equiv: Equiv, a: int, b: int) -> bool:
        return equiv[a] == equiv[b]

    def domerge(self, equiv: Equiv, a: int, b: int, preserveuse: bool = False) -> Equiv:
        if self.type[a] != self.type[b]:
            raise ValueError("structure")
        elif self.strength[a] != self.strength[b]:
            raise ValueError("strength")
        elif equiv[a] == equiv[self.complement[b]]:
            raise ValueError("self-comp")
        elif preserveuse and (self.use[a] != self.use[b]):
            raise ValueError("use")
        else:
            equiv = copy(equiv)
            newg, oldg = sorted((equiv[a], equiv[b]))
            newc, oldc = sorted((equiv[self.complement[a]], equiv[self.complement[b]]))
            equiv[equiv == oldg] = newg
            equiv[equiv == oldc] = newc
        return equiv


# U_UNUSED = -1
# U_INPUT = 1
# U_OUTPUT = 0


class FTileList:
    stiles: FTilesArray
    htiles: FTilesArray
    vtiles: FTilesArray

    def __init__(self, tiles: TileList[Tile], gluelist: FGlueList):
        self.tiles: list[FTile] = []
        self.totile: dict[str, FTile] = {}
        for t in tiles:
            glues = np.array([gluelist.tonum[x.ident()] for x in t.edges])
            if t.is_fake:
                continue
            if t.use is None:
                # if "input" in t.keys():
                #    use = np.array([[uO, uI][int(x)] for x in t["input"]])
                #    used = True
                # else:
                used = False
                use = np.array([uU for _ in t.edges])
            else:
                used = True
                use = np.array([x.value if x is not None else uB for x in t.use])
            # color = "label" in t.keys()
            self.tiles.append(
                FTile(
                    name=t.ident(),
                    color=False,  # FIXME
                    use=use,
                    glues=glues,
                    used=used,
                    structure=t.__class__,
                    dfake=0,
                    sfake=0,
                )
            )
            assert t.name not in self.totile.keys()
            self.totile[t.name] = self.tiles[-1]
        stiles = []
        htiles = []
        vtiles = []

        for t in self.tiles:
            if issubclass(t.structure, SingleTile):
                stiles.append(t)
            elif issubclass(t.structure, HDupleTile):
                stiles += _ffakesingle(t, gluelist)
                htiles.append(t)
            elif issubclass(t.structure, VDupleTile):
                stiles += _ffakesingle(t, gluelist)
                vtiles.append(t)
            else:
                raise NotImplementedError
        self.stiles = FTilesArray.from_ftiles(stiles)

        for i in range(0, len(self.stiles)):
            x, y = _ffakedouble_n(i, self.stiles, gluelist)
            htiles += x
            vtiles += y

        self.htiles = FTilesArray.from_ftiles(htiles)
        self.vtiles = FTilesArray.from_ftiles(vtiles)


RSEL = (2, 3, 0, 1)
FTI = (1, 0, 1, 0)
FTS = (VDupleTile, HDupleTile, VDupleTile, HDupleTile)


def _fdg(dir, gs1, gs2):
    if dir == 0:
        return [gs2[0], gs2[1], gs1[1], gs1[2], gs1[3], gs2[3]]
    elif dir == 1:
        return [gs1[0], gs2[0], gs2[1], gs2[2], gs1[2], gs1[3]]
    elif dir == 2:
        return [gs1[0], gs1[1], gs2[1], gs2[2], gs2[3], gs1[3]]
    elif dir == 3:
        return [gs2[0], gs1[0], gs1[1], gs1[2], gs2[2], gs2[3]]
    else:
        raise ValueError(dir)


@overload
def _ffakedouble_n(
    tn: int,
    sta,
    gluelist: FGlueList,
    outputonly: bool = True,
    *,
    dir4: Literal[True],
    equiv=None,
) -> tuple[list[FTile], list[FTile], list[FTile], list[FTile]]:
    ...


@overload
def _ffakedouble_n(
    tn: int,
    sta,
    gluelist: FGlueList,
    outputonly: bool = True,
    *,
    dir4: Literal[False] = False,
    equiv=None,
) -> tuple[list[FTile], list[FTile]]:
    ...


def _ffakedouble_n(
    tn: int,
    sta,
    gluelist: FGlueList,
    outputonly: bool = True,
    *,
    dir4: bool = False,
    equiv=None,
) -> tuple[list[FTile], list[FTile]] | tuple[
    list[FTile], list[FTile], list[FTile], list[FTile]
]:
    if equiv is None:
        equiv = gluelist.blankequiv()
    if not dir4:
        faketiles: tuple[list[FTile], list[FTile]] | tuple[
            list[FTile], list[FTile], list[FTile], list[FTile]
        ] = ([], [])
        fti = FTI
    else:
        faketiles = (
            [],
            [],
            [],
            [],
        )
        fti = (0, 1, 2, 3)
    if sta.sfake[tn]:
        ddir = np.flatnonzero(gluelist.tonum["fakedouble"] == sta.glues[tn])[0]
    for dir in range(0, 4):
        if (sta.use[tn, dir] == uI) and outputonly:
            continue
        if outputonly:
            oti = np.nonzero(
                (
                    equiv[gluelist.complement[sta.glues[tn, dir]]]
                    == equiv[sta.glues[:, RSEL[dir]]]
                )
                & (sta.use[:, RSEL[dir]] == uI)
                & (sta.used)
            )
        else:
            oti = np.nonzero(
                (
                    equiv[gluelist.complement[sta.glues[tn, dir]]]
                    == equiv[sta.glues[:, RSEL[dir]]]
                )
            )
        for i in oti[0]:
            faketiles[fti[dir]].append(
                FTile(
                    color=False,
                    used=True,
                    name=sta.name[tn] + "_{}_".format(dir) + sta.name[i],
                    structure=FTS[dir],
                    glues=np.array(_fdg(dir, sta.glues[tn, :], sta.glues[i])),
                    use=np.array(_fdg(dir, sta.use[tn, :], sta.use[i])),
                    dfake=dir + 1,
                    sfake=(sta.sfake[tn] or sta.sfake[i]),
                )
            )
        if sta.sfake[tn] and (dir == ddir):  # type: ignore
            faketiles[fti[dir]].append(
                FTile(
                    color=False,
                    used=True,
                    name=sta.name[tn]
                    + "_{}_".format(dir)
                    + sta.name[tn + sta.sfake[tn]],
                    structure=FTS[dir],
                    glues=np.array(
                        _fdg(
                            dir,
                            sta.glues[tn, :],
                            sta.glues[tn + sta.sfake[tn]],
                        )
                    ),
                    use=np.array(
                        _fdg(
                            dir,
                            sta.use[tn, :],
                            sta.use[tn + sta.sfake[tn]],
                        )
                    ),
                    dfake=dir + 1,
                    sfake=True,
                )
            )
    return faketiles


# THESE SELECTORS ARE 1-INDEXED!!! 0 corresponds to fake double tile bond.
HSEL = ((1, 0, 5, 6), (2, 3, 4, 0))
VSEL = ((1, 2, 0, 6), (0, 3, 4, 5))


def _ffakesingle(ftile: FTile, gluelist: FGlueList) -> Sequence[FTile]:
    # FIXME: should be more generalized.  Currently only tau=2
    if issubclass(ftile.structure, VDupleTile):
        sel = HSEL
    elif issubclass(ftile.structure, HDupleTile):
        sel = VSEL
    else:
        raise NotImplementedError

    # Start by making the tiles, then change around the inputs
    fdb = gluelist.tonum["fakedouble"]
    glues: list[list[int]] = [[([fdb] + list(ftile.glues))[x] for x in y] for y in sel]
    fuse = [[([uP] + list(ftile.use))[x] for x in y] for y in sel]
    use = []
    used = []
    names = [ftile.name + "_fakedouble_a", ftile.name + "_fakedouble_b"]
    for gu in zip(glues, fuse):
        if sum(gluelist.strength[g] for g, u in zip(*gu) if u == uI) >= TAU:
            use.append(gu[1])
            used.append(True)
        else:
            use.append(gu[1])
            used.append(False)
    return [
        FTile(
            color=ftile.color,
            use=u,
            glues=np.array(g),
            name=n,
            used=ud,
            structure=SingleTile,
            dfake=0,
            sfake=sfo,
        )
        for u, g, n, ud, sfo in zip(use, glues, names, used, [1, -1])
    ]


class _FastTileSet:
    gluelist: FGlueList
    tilelist: FTileList

    def __init__(self, tilesystem: "TileSet"):
        self.gluelist = FGlueList(
            tilesystem.allglues
            + [
                # End({"name": "hp", "type": "hairpin", "strength": 0}),
                Glue("fakedouble", "fakedouble", abstractstrength=0, use=Use.PERMANENT)
            ]
        )
        self.tilelist = FTileList(
            tilesystem.tiles
            # + sum([x.named_rotations() for x in tilesystem.tiles], TileList()) # FIXME
            ,
            self.gluelist,
        )

    def applyequiv(self, ts: TileSet, equiv) -> TileSet:
        ts = ts.copy()
        alreadythere: list[tuple[list[Glue], bool]] = []
        for tile in ts.tiles:
            tile.edges = [
                self.gluelist.name[equiv[self.gluelist.tonum[e.ident()]]]
                for e in tile.edges
            ]
            if (tile.edges, False) in alreadythere:  # FIXME: False was label
                tile.fake = True
                continue
            rs = [tile]  # + tile.rotations FIXME
            alreadythere.extend((list(t.edges), False) for t in rs)
        # FIXME
        # if "seed" in ts.keys():
        #    for t in ts.seed["adapters"]:
        #        if "ends" in t.keys():
        #            t["ends"] = [
        #                self.gluelist.name[equiv[self.gluelist.tonum[e]]]
        #                for e in t["ends"]
        #            ]
        # ts["info"] = ts.get("info", dict())
        # ts["info"]["fgluemerge"] = ts["info"].get("fgluemerge", list())
        # ts["info"]["fgluemerge"].append([int(x) for x in equiv])
        return ts

    def togluemergespec(self, ts, equiv):
        gms = util.GlueMergeSpec()
        for i in range(0, len(equiv)):
            if i != equiv[i]:
                gms.add(self.gluelist.name[equiv[i]], self.gluelist.name[i])
        return gms


def ptins(
    fts: _FastTileSet, equiv: Equiv | None = None, tau=2
) -> list[list[np.ndarray[Any, np.dtype[np.int64]]]]:
    """Calculate potential tile attachments to input neighborhoods"""
    ptins = []
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
    for ta in [fts.tilelist.stiles, fts.tilelist.htiles, fts.tilelist.vtiles]:
        ptin = []
        ti: int
        for ti in np.arange(0, len(ta.used))[ta.used]:
            # Only iterate through used tiles
            isel = ta.use[ti] == uI  # used edges
            gsel = equiv[ta.glues[ti, isel]] == equiv[ta.glues[:, isel]]
            # matching glues
            strs = np.sum(
                fts.gluelist.strength[ta.glues[ti, isel]] * gsel, axis=1
            )  # strengths of matching
            matches = np.nonzero((strs >= tau) & (ta.dfake == 0))
            # indices of matching
            # (excludes fake doubles, which can't actually attach
            # (because they are
            # actually two singles) to their own local neighborhoods)
            ptin.append(matches)
        ptins.append(ptin)
    return ptins


def gmatch(fts, g1, g2, equiv=None):
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
    return equiv[g1] == equiv[g2]


def gcomp(fts, g1, g2, equiv=None):
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
    return equiv[fts.gluelist.complement[g1]] == equiv[g2]


def findmovetiles(fts, tilei, direction, allin, allout):
    use = fts.tilelist.stiles.use[tilei, direction]
    if use == uN:
        return {}
    if use == uP:
        otherdoublei = tilei + fts.tilelist.stiles.sfake[tilei]
        return {(otherdoublei, allin, allout, other[direction])}
    if use == uI:
        tilematches = np.nonzero(
            gcomp(
                fts,
                fts.tilelist.stiles.glues[:, other[direction]],
                fts.tilelist.stiles.glues[tilei, direction],
            )
            & (
                (fts.tilelist.stiles.use[:, other[direction]] == uO)
                | (fts.tilelist.stiles.use[:, other[direction]] == uB)
            )
        )[0]
        return {(tmatch, allin, False, other[direction]) for tmatch in tilematches}
    if use == uO:
        tilematches = np.nonzero(
            gcomp(
                fts,
                fts.tilelist.stiles.glues[:, other[direction]],
                fts.tilelist.stiles.glues[tilei, direction],
            )
            & (
                (fts.tilelist.stiles.use[:, other[direction]] == uI)
                | (fts.tilelist.stiles.use[:, other[direction]] == uB)
            )
        )[0]
        return {(tmatch, False, allout, other[direction]) for tmatch in tilematches}
    if use == uB:
        tilematches = np.nonzero(
            gcomp(
                fts,
                fts.tilelist.stiles.glues[:, other[direction]],
                fts.tilelist.stiles.glues[tilei, direction],
            )
            & (
                (fts.tilelist.stiles.use[:, other[direction]] == uI)
                | (fts.tilelist.stiles.use[:, other[direction]] == uO)
                | (fts.tilelist.stiles.use[:, other[direction]] == uB)
            )
        )[0]
        return {(tmatch, allin, allout, other[direction]) for tmatch in tilematches}


def _2go_moveandfill(
    fts: _FastTileSet,
    t: int,
    i: int,
    allin=True,
    allout=True,
    pdir=None,
    x=0,
    y=0,
    old_d=None,
) -> TMap:
    if old_d is None:
        d: dict[tuple[int, int], set[Any]] = dict()
    else:
        d = old_d

    s: set[tuple[int, int, int]] = d.get((x, y), set())
    d[(x, y)] = s
    for edge in (0, 1, 2, 3):
        if allin:
            continue
        elif edge == pdir:
            if (fts.tilelist.stiles.use[t, edge] == uI) & (not allout):
                s.add((edge, fts.tilelist.stiles.glues[t, edge], t))
        elif fts.tilelist.stiles.use[t, edge] == uI:
            s.add((edge, fts.tilelist.stiles.glues[t, edge], t))
    if i == 0:
        return d
    for edge in (0, 1, 2, 3):
        if edge == pdir:
            continue
        dx, dy = directions[edge]
        moves = findmovetiles(fts, t, edge, allin, allout)
        # assert moves is not None  # FIXME
        if moves is None:
            continue  # FIXME
        for newt, newallin, newallout, newpdir in moves:
            d = _2go_moveandfill(
                fts, newt, i - 1, newallin, newallout, newpdir, x + dx, y + dy, d
            )
    return d


def _2go_findtrialmoves(fts: _FastTileSet, tilei, direction, equiv=None):
    use = fts.tilelist.stiles.use[tilei, direction]
    if use == uN:
        return {}
    if use == uP:
        otherdoublei = tilei + fts.tilelist.stiles.sfake[tilei]
        return {
            (otherdoublei, 0, other[direction])
        }  # other tile, (is not a new tile -> 0), prev direction
    else:
        tilematches = np.nonzero(
            gcomp(
                fts,
                fts.tilelist.stiles.glues[:, other[direction]],
                fts.tilelist.stiles.glues[tilei, direction],
                equiv,
            )
        )[0]
        return {(tmatch, 1, other[direction]) for tmatch in tilematches}


directions = [(0, -1), (1, 0), (0, 1), (-1, 0)]


def _2go_checkandmove(fts, ct, x, y, i, tmap: TMap, exclude=tuple(), equiv=None):

    for direction in (0, 1, 2, 3):
        if direction in exclude:
            continue
        if (fts.tilelist.stiles.use[ct, direction] == uP) | (
            fts.tilelist.stiles.use[ct, direction] == uN
        ):
            continue
        for d, mg, ot in tmap.get((x, y), set()):
            if d != direction:
                continue  # FIXME optimize out
            if gmatch(fts, mg, fts.tilelist.stiles.glues[ct, direction], equiv):
                return True, ot, ct, x, y
    # Move to next, if we still have moves
    if i == 0:
        return False, None, None, None, None  # reached the end
    for direction in (0, 1, 2, 3):
        if direction in exclude:
            continue
        tm = _2go_findtrialmoves(fts, ct, direction, equiv)
        dx, dy = directions[direction]
        for newct, inccount, newexclude in tm:
            # print("checkandmove({}, {}, {}, {},tmap, ({}))".format(newct, x+dx, y+dy, i-inccount, newexclude))
            r = _2go_checkandmove(
                fts, newct, x + dx, y + dy, i - inccount, tmap, (newexclude,), equiv
            )
            if r[0] is not False:
                return r
    return False, None, None, None, None


@overload
def is_2go_nn(
    fts: _FastTileSet,
    tn: int,
    un: int,
    equiv: np.ndarray[Any, np.dtype[np.int64]],
    tau: int = 2,
    tmaps: dict[int, TMap] | None = None,
    *,
    retall: Literal[True],
    also22go: bool = False,
) -> tuple[bool, tuple[Any, Any] | None]:
    ...


@overload
def is_2go_nn(
    fts: _FastTileSet,
    tn: int,
    un: int,
    equiv: np.ndarray[Any, np.dtype[np.int64]],
    tau: int = 2,
    tmaps: dict[int, TMap] | None = None,
    *,
    retall: Literal[False] = False,
    also22go: bool = False,
) -> bool:
    ...


def is_2go_nn(
    fts: _FastTileSet,
    tn: int,
    un: int,
    equiv: np.ndarray[Any, np.dtype[np.int64]],
    tau: int = 2,
    tmaps: dict[int, TMap] | None = None,
    *,
    retall: bool = False,
    also22go: bool = False,
) -> bool | tuple[bool, tuple[Any, Any] | None]:

    if also22go:
        raise NotImplementedError()

    # Before starting, check to see if t and u are actually the same:
    if np.all(
        equiv[fts.tilelist.stiles.glues[tn]] == equiv[fts.tilelist.stiles.glues[un]]
    ):
        # if not also22go:
        if not retall:
            return False
        else:
            return False, None
        # else:
        #    if not retall:
        #        return False, False
        #    else:
        #        return False, False, None, None

    if tmaps is not None:
        tmap = tmaps[tn]
    else:
        tmap = _2go_moveandfill(fts, tn, 2)

    res = _2go_checkandmove(fts, un, 0, 0, 1, tmap=tmap, exclude=tuple(), equiv=equiv)
    if res[0] is False:
        if not retall:
            return False
        else:
            return False, None
    else:
        if not retall:
            return True
        else:
            return True, (res[1], res[2])


def gen_2go_maps(fts: _FastTileSet) -> dict[int, TMap]:
    sel: np.ndarray[Any, np.dtype[np.int64]] = np.flatnonzero(fts.tilelist.stiles.used)
    maps: dict[int, TMap] = dict()
    for t in sel:
        maps[t] = _2go_moveandfill(fts, t, 2)
    return maps


def gen_2go_profile(fts: _FastTileSet, equiv=None, tmaps=None, also22go: bool = False):
    if tmaps is None:
        tmaps = gen_2go_maps(fts)
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
    sens1s = ptins(fts, tau=1, equiv=equiv)[0]
    sens2s = []
    sens22s = []
    sel = np.flatnonzero(fts.tilelist.stiles.used)
    for tn, uns in enumerate(sens1s):
        t = sel[tn]
        x2 = []
        x22 = []
        for un in uns[0]:
            if also22go:
                s2, s22 = is_2go_single_nn(  # FIXME: this was still 2go_single; why?
                    fts, t, un, equiv, tau=2, also22go=True
                )
            else:
                s2 = is_2go_nn(fts, t, un, equiv, tmaps=tmaps, also22go=False)
            if s2:
                x2.append(un)
            if also22go and s22:
                x22.append(un)
        sens2s.append(x2)
        if also22go:
            sens22s.append(x22)

    if also22go:
        return sens2s, sens22s
    else:
        return sens2s


def is_2go_equiv(
    fts: _FastTileSet, equiv: Equiv | None = None, tmaps=None, origsens=None
):
    if tmaps is None:
        tmaps = gen_2go_maps(fts)
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
    if origsens is None:
        origsens = gen_2go_profile(fts, tmaps=tmaps)
    sens1s = ptins(fts, tau=1, equiv=equiv)[0]
    sel = np.flatnonzero(fts.tilelist.stiles.used)

    for tn, (uns, os) in enumerate(zip(sens1s, origsens)):
        t = sel[tn]
        for un in uns[0]:
            r, p = is_2go_nn(fts, t, un, equiv, tmaps=tmaps, retall=True)
            if r and (un not in os):
                # Next if checks to see if sensitive tile is actually identical
                # to an already sensitive tile:
                # FIXME: this would not be necessary if we stored whether a
                # tile was a duplicate.  However, that would break the
                # glue-equiv-only method...
                if not np.any(
                    np.all(
                        equiv[fts.tilelist.stiles.glues[os]]
                        == equiv[fts.tilelist.stiles.glues[un]],
                        axis=1,
                    )
                ):
                    return (
                        False,
                        (
                            fts.tilelist.stiles.name[fts.tilelist.stiles.used][tn],
                            fts.tilelist.stiles.name[un],
                        ),
                        [fts.tilelist.stiles.name[x] for x in p],
                    )
    return True, None, None


def is_22go_equiv(fts, equiv=None, ins2go=None, orig22go=None):
    if ins2go is None:
        ins2go = gen_2go_single_ins(fts, tau=2)
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
    if orig22go is None:
        _, orig22go = gen_2go_profile(fts, also22go=True)
    sens1s = ptins(fts, tau=1, equiv=equiv)[0]

    for tn, (uns, in2go, os22) in enumerate(zip(sens1s, ins2go, orig22go)):
        for un in uns[0]:
            r2, r22, p2, p22 = is_2go_single_nn(
                fts, tn, un, equiv, tau=2, in2go=in2go, also22go=True, retall=True
            )
            if r22 and (un not in os22):
                # Next if checks to see if sensitive tile is actually identical
                # to an already sensitive tile:
                # FIXME: this would not be necessary if we stored whether a
                # tile was a duplicate.  However, that would break the
                # glue-equiv-only method...
                if not np.any(
                    np.all(
                        equiv[fts.tilelist.stiles.glues[os22]]
                        == equiv[fts.tilelist.stiles.glues[un]],
                        axis=1,
                    )
                ):
                    return (
                        False,
                        (
                            fts.tilelist.stiles.name[fts.tilelist.stiles.used][tn],
                            fts.tilelist.stiles.name[un],
                        ),
                        p22,
                    )
    return True, None, None


def fta_to_ft(stiles, un, sel=None):
    if sel is None:
        sel = np.ones_like(stiles.used, dtype=bool)
    return FTile(
        color=stiles.color[sel][un],
        use=stiles.use[sel, :][un, :],
        glues=stiles.glues[sel, :][un, :],
        name=stiles.name[sel][un],
        used=stiles.used[sel][un],
        structure=stiles.structure[sel][un],
        dfake=stiles.dfake[sel][un],
        sfake=stiles.sfake[sel][un],
    )


def isatamequiv(fts, equiv, initptins=None):
    if initptins is None:
        initptins = ptins(fts, fts.gluelist.blankequiv())
    npt = ptins(fts, equiv)
    for x, y, tl in zip(
        initptins, npt, [fts.tilelist.stiles, fts.tilelist.htiles, fts.tilelist.vtiles]
    ):
        for xx, yy in zip(x, y):
            if len(xx[0]) == len(yy[0]) and (np.all(xx[0] == yy[0])):  # No change
                continue
            elif len(xx[0]) == 1:  # Deterministic start
                mm = ~np.all((equiv[tl.glues[xx]] == equiv[tl.glues[yy]]), axis=1) | ~(
                    tl.color[xx] == tl.color[yy]
                )
                if np.any(mm):
                    return False, (tl.name[xx[0][0]], tl.name[yy[0][mm][0]])
            elif len(xx[0]) == 0:
                return False, None
            else:
                raise NotImplementedError(xx)
    return True, None


def tilemerge(fts, equiv, t1, t2, preserveuse=False):
    if t1.structure.name != t2.structure.name:
        raise ValueError
    if t1.color != t2.color:
        raise ValueError
    for g1, g2 in zip(t1.glues, t2.glues):
        equiv = fts.gluelist.domerge(equiv, g1, g2, preserveuse=preserveuse)
    return equiv


def _findpotentialtilemerges(fts, equiv):
    ppairs = []
    for ti in range(0, len(fts.tilelist.tiles)):
        t1 = fts.tilelist.tiles[ti]
        if not t1.used:
            continue
        ppairs.extend(
            (t1, t)
            for t in fts.tilelist.tiles[ti:]
            if (t1.color == t.color)
            and (t1.structure.name == t.structure.name)
            and (t1.name != t.name)
        )
    shuffle(ppairs)
    return ppairs


def _findpotentialgluemerges(fts, equiv):
    ppairs = []
    for g1 in np.arange(0, len(fts.gluelist.strength)):
        g2s = (
            g1
            + 1
            + np.nonzero(
                (fts.gluelist.strength[g1 + 1 :] == fts.gluelist.strength[g1])
                & (fts.gluelist.structure[g1 + 1 :] == fts.gluelist.structure[g1])
            )[0]
        )
        ppairs.extend((g1, g2) for g2 in g2s)
    shuffle(ppairs)
    return ppairs


def _recfix(
    fts,
    equiv,
    tp,
    initptins,
    check2go=False,
    tmaps=None,
    orig2go=None,
    orig22go=None,
    check22go=False,
    checkld=False,
    chain=None,
    preserveuse=False,
):
    log = logging.getLogger(__name__)
    if chain is None:
        chain = []
    equiv = tilemerge(
        fts,
        equiv,
        fts.tilelist.totile[tp[0]],
        fts.tilelist.totile[tp[1]],
        preserveuse=preserveuse,
    )
    ae, badpair = isatamequiv(fts, equiv, initptins=initptins)
    if check2go and ae:
        ae, badpair, _ = is_2go_equiv(fts, equiv, tmaps=tmaps, origsens=orig2go)
    if ae and check22go:
        ae, badpair, _ = is_22go_equiv(fts, equiv, tmaps=tmaps, orig22go=orig22go)
    if ae and checkld:
        ld_e = fld.latticedefects(fts, "e", equiv=equiv)
        ld_w = fld.latticedefects(fts, "w", equiv=equiv)
        if (len(ld_e) > 0) or (len(ld_w) > 0):
            raise ValueError
    if ae:
        log.debug("Recfix succeeds with {}, {}, {}".format(tp[0], tp[1], chain))
        return equiv
    elif badpair is None:
        raise ValueError
    else:
        chain.append(tp)
        return _recfix(
            fts,
            equiv,
            badpair,
            initptins,
            check2go,
            tmaps,
            orig2go,
            orig22go=orig22go,
            check22go=check22go,
            checkld=checkld,
            chain=chain,
            preserveuse=preserveuse,
        )


def _tilereduce(
    fts,
    equiv=None,
    check2go=False,
    initptins=None,
    check22go=False,
    checkld=False,
    preserveuse=False,
):
    log = logging.getLogger(__name__)
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
    todo = _findpotentialtilemerges(fts, equiv)
    if initptins is None:
        initptins = ptins(fts)
    if check2go:
        tmaps = gen_2go_maps(fts)
        origsens = gen_2go_profile(fts, tmaps=tmaps)
    else:
        tmaps = None
        origsens = None
    for todoi, (t1, t2) in enumerate(todo):
        try:
            nequiv = tilemerge(fts, equiv, t1, t2, preserveuse=preserveuse)
        except ValueError:
            continue
        ae, badpair = isatamequiv(fts, nequiv, initptins=initptins)
        if ae and check2go:
            ae, badpair, _ = is_2go_equiv(fts, nequiv, tmaps, origsens=origsens)
        if ae and check22go:
            ae, badpair, _ = is_22go_equiv(
                fts, nequiv, ins2go=ins2go, orig22go=orig22go
            )
        if ae and checkld:
            ld_e = fld.latticedefects(fts, "e", equiv=nequiv)
            ld_w = fld.latticedefects(fts, "w", equiv=nequiv)
            if (len(ld_e) > 0) or (len(ld_w) > 0):
                continue
        if ae:
            equiv = nequiv
            log.debug(
                "Reduced {}, {} ({} of {} pairs done)".format(
                    t1.name, t2.name, todoi, len(todo)
                )
            )
        elif badpair is None:
            continue
        else:
            try:
                equiv = _recfix(
                    fts,
                    equiv,
                    badpair,
                    initptins,
                    check2go,
                    tmaps,
                    origsens,
                    check22go=check22go,
                    checkld=checkld,
                    chain=[(t1.name, t2.name)],
                    preserveuse=preserveuse,
                )
            except ValueError:
                continue
            except KeyError:
                continue
    return equiv


def _gluereduce(
    fts: _FastTileSet,
    equiv: Equiv | None = None,
    check2go: bool = False,
    check22go: bool = False,
    checkld: bool = False,
    initptins=None,
    preserveuse: bool = False,
):
    log = logging.getLogger(__name__)
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
    todo = _findpotentialgluemerges(fts, equiv)
    if initptins is None:
        initptins = ptins(fts)
    if check2go:
        tmaps = gen_2go_maps(fts)
        origsens = gen_2go_profile(fts, tmaps=tmaps)
    else:
        ins2go = None
        origsens = None
        orig22go = None
        tmaps = None
    for todoi, (g1, g2) in enumerate(todo):
        try:
            nequiv = fts.gluelist.domerge(equiv, g1, g2, preserveuse=preserveuse)
        except ValueError:
            continue
        ae, badpair = isatamequiv(fts, nequiv, initptins=initptins)
        if ae and check2go:
            ae, badpair, _ = is_2go_equiv(fts, nequiv, tmaps=tmaps, origsens=origsens)
        if ae and check22go:
            ae, badpair, _ = is_22go_equiv(
                fts, nequiv, ins2go=ins2go, orig22go=orig22go
            )
        if ae and checkld:
            ld_e = fld.latticedefects(fts, "e", equiv=nequiv)
            ld_w = fld.latticedefects(fts, "w", equiv=nequiv)
            if (len(ld_e) > 0) or (len(ld_w) > 0):
                continue
        if ae:
            equiv = nequiv
            log.debug(
                "Glue reduction: {}, {} ({}/{} done)".format(
                    fts.gluelist.name[g1], fts.gluelist.name[g2], todoi, len(todo)
                )
            )
        elif badpair is None:
            continue
        else:
            try:
                equiv = _recfix(
                    fts,
                    equiv,
                    badpair,
                    initptins,
                    check2go,
                    tmaps,
                    origsens,
                    check22go=check22go,
                    checkld=checkld,
                    chain=[(fts.gluelist.name[g1], fts.gluelist.name[g2])],
                    preserveuse=preserveuse,
                )
            except ValueError:
                continue
            except KeyError:
                continue
    return equiv


def _single_reduce_tiles(p):
    feq, ts, fts, params = p
    # starttime = time.time()
    initptins = ptins(fts)
    e = _tilereduce(fts, equiv=feq, initptins=initptins, **params)
    # endtime = time.time()
    # te = fts.applyequiv(ts, e)
    # nt = len([y for y in te.tiles if 'fake' not in y.keys()])
    # ng = len(te.allends)

    # a = isatamequiv(fts, e)
    # g2 = is_2go_equiv(fts, e)[0]
    # g22 = is_22go_equiv(fts, e)[0]
    # lde = len(fld.latticedefects(fts, equiv=e, direction='e')) == 0
    # ldw = len(fld.latticedefects(fts, equiv=e, direction='w')) == 0
    # print('TR: Found {}t, {}g: aTAM {} 2GO {} 22GO {} LD {}/{} in {}s'.format(
    #    nt, ng, a, g2, g22, lde, ldw, endtime - starttime))
    # te.to_file(nb + 'step1-{}t{}g.yaml'.format(nt, ng))
    return e


def _single_reduce_ends(p):
    equiv, ts, fts, params = p
    # starttime = time.time()
    initptins = ptins(fts)
    e = _gluereduce(fts, equiv=equiv, initptins=initptins, **params)
    # FIXME: all this is code for checking things and printing logs.
    # endtime = time.time()
    # te = fts.applyequiv(ts, e)
    # nt = len([y for y in te.tiles if 'fake' not in y.keys()])
    # ng = len(te.allends)
    # a = isatamequiv(fts, e)
    # g2 = is_2go_equiv(fts, e)[0]
    # g22 = is_22go_equiv(fts, e)[0]
    # lde = len(fld.latticedefects(fts, equiv=e, direction='e')) == 0
    # ldw = len(fld.latticedefects(fts, equiv=e, direction='w')) == 0
    # print('ER: Found {}t, {}g: aTAM {} 2GO {} 22GO {} LD {}/{} in {}s'.format(
    #    nt, ng, a, g2, g22, lde, ldw, endtime - starttime))
    # te.to_file(nb + '-step2-{}t{}g.yaml'.format(nt, ng))
    return e


def reduce_tiles(
    tileset,
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
    reduced: single TileSet
        The reduced system/systems
    """

    fts = _FastTileSet(tileset)

    if key is None:
        key = lambda x: len(
            np.unique(x)
        )  # number of unique numbers in equiv, equivalent
        # to number of glues in reduced system. FIXME?

    # FIXME: could do a better job here
    params = {
        "check2go": "s2" in preserve,
        "check22go": "s22" in preserve,
        "checkld": "ld" in preserve,
        "preserveuse": "gs" in preserve,
    }

    if initequiv is None:
        initequiv = fts.gluelist.blankequiv()

    if threads > 1:
        from multiprocessing import Pool

        with Pool(threads) as pool:
            equivs = pool.map(
                _single_reduce_tiles, [[initequiv, tileset, fts, params]] * tries
            )
    else:
        equivs = [
            _single_reduce_tiles(x)
            for x in ([[initequiv, tileset, fts, params]] * tries)
        ]

    equivs.sort(key=key)

    if returntype == "TileSet":
        equivs = [tileset.apply_equiv(equiv) for equiv in equivs]

    if best == 1:
        return equivs[0]
    else:
        return equivs[0:best]


def reduce_ends(
    tileset,
    preserve=("s22", "ld"),
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
    reduced: single TileSet
        The reduced system/systems
    """

    fts = _FastTileSet(tileset)

    if key is None:
        key = lambda x: len(
            np.unique(x)
        )  # number of unique numbers in equiv, equivalent
        # to number of glues in reduced system. FIXME?

    # FIXME: could do a better job here
    params = {
        "check2go": "s2" in preserve,
        "check22go": "s22" in preserve,
        "checkld": "ld" in preserve,
        "preserveuse": "gs" in preserve,
    }

    if initequiv is None:
        initequiv = fts.gluelist.blankequiv()

    if threads > 1:
        from multiprocessing import Pool

        with Pool(threads) as pool:
            equivs = pool.map(
                _single_reduce_ends, [[initequiv, tileset, fts, params]] * tries
            )
    else:
        equivs = [
            _single_reduce_ends(x)
            for x in ([[initequiv, tileset, fts, params]] * tries)
        ]

    equivs.sort(key=key)

    if returntype == "TileSet":
        equivs = [tileset.apply_equiv(equiv) for equiv in equivs]

    if best == 1:
        return equivs[0]
    else:
        return equivs[0:best]
