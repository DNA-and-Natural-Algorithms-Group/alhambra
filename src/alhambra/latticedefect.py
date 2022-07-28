import itertools
import re

from .sensitivitynew import _fakesingles  # type: ignore
from .tiles import TileList
from .util import GlueMergeSpec, comp

_OE = [2, 3, 0, 1]
_ENDER = [3, 2, 1, 0]


def _generic_branch(direction, tiles, tile, n, f=False, gms=GlueMergeSpec([])):
    # Generic branch finder.  MUST USE ALL SINGLES (EG, FAKESINGLE)
    branches = []
    if n == 0:
        return [([tile], tile["ends"][_ENDER[direction]])]
    for tn in tiles:
        if gms.eq(tile["ends"][direction], comp(tn["ends"][_OE[direction]])):
            branches.extend(
                ([tile] + x, y) for x, y in _generic_branch(direction, tiles, tn, n - 1)
            )
    return branches


def _latticedefect_tile(tiles, tile, direction="e", n=2, gms=GlueMergeSpec([])):
    """With n tiles in each branch, can tile t form a lattice defect?"""
    d1, d2 = {"e": (1, 2), "w": (0, 3), "n": (0, 1), "s": (2, 3)}[direction]
    b1 = _generic_branch(d1, tiles, tile, n, gms=gms)
    b2 = _generic_branch(d2, tiles, tile, n, gms=gms)
    neighborhoods = itertools.product(b1, b2)
    res = []
    for n in neighborhoods:
        res.extend(
            (n, tile)
            for tile in tiles
            if (
                gms.eq(n[0][1], comp(tile["ends"][_OE[_ENDER[d1]]]))
                and gms.eq(n[1][1], comp(tile["ends"][_OE[_ENDER[d2]]]))
            )
        )
    return res


def _ppld(res):
    """Pretty-print lattice defects, to some extent.  I is the initial
    tile, A/B are the two branches, and F is the tile that
    attaches to those branches.
    """
    return [
        "I:"
        + n[0][0][0]["name"]
        + " "
        + "A:["
        + ",".join(t["name"] for t in n[0][0][1:])
        + "] "
        + "B:["
        + ",".join(t["name"] for t in n[1][0][1:])
        + "] "
        + "F:"
        + tt["name"]
        for n, tt in res
    ]


def latticedefects(
    ts, direction="e", depth=2, pp=True, rotate=False, gms=GlueMergeSpec([])
):
    if depth < 2:
        raise ValueError("Depth cannot be less than 2, received {}.".format(depth))
    tiles = _fakesingles(ts.tiles)
    rtiles = _fakesingles(
        TileList([x for x in ts.tiles if "fake" not in x.keys()])
        + sum([x.rotations for x in ts.tiles if "fake" not in x.keys()], TileList())
    )
    if rotate:
        tll = rtiles
    else:
        tll = tiles
    alldefects = sum(
        (
            _latticedefect_tile(tll, tile, direction=direction, n=depth, gms=gms)
            for tile in tll
        ),
        [],
    )
    if pp:
        return _ppld(alldefects)
    else:
        return alldefects
