import itertools
import re

import numpy as np

from .sensitivitynew import _fakesingles  # type: ignore
from .tiles import TileList
from .util import GlueMergeSpec, comp

_OE = np.array([2, 3, 0, 1])
_ENDER = np.array([3, 2, 1, 0])


def _branchtiles(gl, st, tns, direction, equiv=None, onlyusedandfs=True):
    """For an array of tile number chains (newest-to-oldest along dim 1)
    in a FTilesArray, find the tiles that match the direction (with equiv).
    Return a new chain array with the new tiles added."""
    if tns is None:
        return None
    newtnl = []
    for tn in tns:
        sel = (
            equiv[gl.complement[st.glues[tn[0], direction]]]
            == equiv[st.glues[:, _OE[direction]]]
        )
        if onlyusedandfs:
            sel = sel & (st.used)
        ntns = np.flatnonzero(sel)
        ftns = np.empty((len(ntns), len(tn) + 1), dtype="int")
        ftns[:, 0] = ntns
        ftns[:, 1:] = tn[None, :]
        newtnl.append(ftns)
    if len(newtnl) == 0:
        return None
    return np.concatenate(newtnl, axis=0)


def _latticedefect_tile(gl, st, tn, direction="e", n=2, equiv=None):
    """With n tiles in each branch, can tile number tn form
    a lattice defect?"""
    if equiv is None:
        equiv = gl.blankequiv()
    d1, d2 = {"e": (1, 2), "w": (0, 3), "n": (0, 1), "s": (2, 3)}[direction]
    b1 = np.array([[tn]], dtype="int")
    b2 = b1
    for _ in range(0, n):
        b1 = _branchtiles(gl, st, b1, d1, equiv=equiv)
        b2 = _branchtiles(gl, st, b2, d2, equiv=equiv)

    if (b1 is None) or (b2 is None):
        return None
    tiledirs = _OE[_ENDER[[d1, d2]]]
    # shape here is:
    tilegluecomps = equiv[gl.complement[st.glues[:, tiledirs]]]

    c1g = st.glues[b1[:, 0], _ENDER[d1]]
    c2g = st.glues[b2[:, 0], _ENDER[d2]]

    # dims here become (CLOSER, BRANCH)
    s1 = tilegluecomps[:, 0, None] == c1g[None, :]
    s2 = tilegluecomps[:, 1, None] == c2g[None, :]

    # dims here become (CLOSER, BRANCH1, BRANCH2)
    ts = (s1[:, :, None] & s2[:, None, :]) & (st.used)[:, None, None]

    r = np.nonzero(ts)

    if len(r[0]) == 0:
        return None
    else:
        return r


def latticedefects(fts, direction="e", depth=2, rotate=False, equiv=None):
    if equiv is None:
        equiv = fts.blankequiv()
    st = fts.tilelist.stiles
    gl = fts.gluelist

    itns = np.flatnonzero(st.used)

    ldl = []
    for tn in itns:
        r = _latticedefect_tile(gl, st, tn, direction, n=depth, equiv=equiv)
        if r is None:
            continue
        ldl.append((tn, r))
    return ldl
