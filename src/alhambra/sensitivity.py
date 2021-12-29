import logging as log

import numpy as np

secondclasses = {
    (False, True): "O",
    (False, False): "U",
    (True, True): "U",
    (True, False): "U",
}
firstsortorder = ["F", "M", "C", "D"]
secondsortorder = ["f", "b", "r", "d"]
flipdir = [2, 3, 0, 1]


def senspairs(tileset, **params):

    # Begin by changing the tileset into a form useful for our calculations.

    tiles = []
    for tile in tileset["tiles"]:

        # Ignore non-single tiles.
        if len(tile["ends"]) != 4:
            log.warning("Ignoring non-single tile {0}".format(tile["name"]))
            continue

        # Make the actual tile list.
        ends = np.array(tile["ends"])
        inputs = np.array(tile["input"])
        tiles.append((tile["name"], ends, inputs))

    senslist = []

    # Now consider each tile:
    for An, Ae, Ai in tiles:
        senslist += _find_pairs_from_tile(An, Ae, Ai, tiles, **params)

    # Now do something with the senslist.

    return senslist


def _find_pairs_from_tile(An, Ae, Ai, tiles, **params):
    """
    Given a tile of name `An`, ends `Ae`, and IO `Ai`, and a list of tiles `tiles`,
    find all sensitive end type pairs.
    """

    senslist = []

    # We need to iterate over tiles (again!):
    for Bn, Be, Bi in tiles:

        # Skip if the tiles are identical:
        if (Be == Ae).all():
            if An == Bn:
                continue
            else:
                log.warning("{0} and {1} have *identical* ends.".format(An, Bn))

        # Now we need to iterate over each non-inert common side:
        for side in (np.nonzero(np.ravel((Ae == Be) & (Ai != -1) & (Bi != -1))))[0]:
            senslist += _pairs_commonside(An, Ae, Ai, Bn, Be, Bi, side, tiles, **params)

    return senslist


def _pairs_commonside(An, Ae, Ai, Bn, Be, Bi, common, tiles, order=2, **params):
    """
    Given two tiles and a common side, find all sensitive pairs.
    """

    senslist = []

    # We need to loop over every side except the common one.
    for mismatch in range(0, 4):

        # Skip the common mismatch
        if mismatch == common:
            continue

        # Skip any inert mismatch
        if Ai[mismatch] == -1 or Bi[mismatch] == -1:
            continue

        # Skip, but warn, about any second common mismatch
        if Ae[mismatch] == Be[mismatch]:
            log.debug("Tiles {} and {} have two common sides.".format(An, Bn))
            continue

        # Now we know we have a sensitive pair, and we just need to determine type.
        # First order type at this point is entirely IO dependent.
        # Second order will require more loops...

        # This is the first-order sensitivity type.
        if (Ai[common]) >= 1 and (Ai[mismatch]) >= 1:
            senstype = "O1"
        else:
            senstype = "U1"

        # Now we need to find second-order
        if order == 2:
            st1, st2 = (
                senstype,
                _second_order_type(
                    An, Ae, Ai, Bn, Be, Bi, common, mismatch, tiles, senstype
                ),
            )
            if st1 == "O1" and st2 == "O2":
                sens = ("2GO", "2NGO")
            elif st1 == "O1" and st2 == "U2":
                sens = ("1GO", "2NGO")
            elif st1 == "O1" and st2 == "":
                sens = ("1GO", "1NGO")
            elif st1 == "U1" and st2 == "":
                sens = ("0GO", "1NGO")
            elif st1 == "U1" and st2 == "U2":
                sens = ("0GO", "2NGO")
            else:
                raise "Impossible situation {},{}".format(st1, st2)

        senslist.append((Ae[mismatch], Be[mismatch], sens))

    return senslist


def _second_order_type(An, Ae, Ai, Bn, Be, Bi, common, mismatch, tiles, firsttype):

    types = set()
    # We need to consider every non-mismatch, non-common side of A:
    for internal in range(0, 4):
        if internal == common:
            continue
        if internal == mismatch:
            continue

        # Now we need to consider every tile that fits there:
        for A2n, A2e, A2i in tiles:
            if not (A2e[flipdir[internal]] == comp(Ae[internal])):
                continue
            # And then consider every tile that corresponds...
            for B2n, B2e, B2i in tiles:
                if B2e[flipdir[internal]] != comp(Be[internal]):
                    continue
                # And now we consider every common side of that except the bond it has with the
                # first tile...
                for secondedge in (
                    np.nonzero(np.ravel((A2e == B2e) & (A2i != -1) & (B2i != -1)))
                )[0]:
                    if secondedge == flipdir[internal]:
                        continue
                    types.add(secondclasses[(Ai[internal] >= 1, A2i[secondedge] >= 1)])

    if len(types) == 0:
        return ""
    elif "O" in types and firsttype == "O1":
        return "O2"
    else:
        return "U2"


def consolidate_pairs(senslist, onlytop=False, **params):
    # We start out with a [ (e1, e2, type) ] list of sensitivity types.
    # We need to consolidate it.

    # To begin, we'll make a dict...

    sensdict = {"1GO": set(), "2GO": set(), "1NGO": set(), "2NGO": set()}

    for end1, end2, senstype in senslist:
        # First we need to sort the ends:
        if end1 > end2:
            end2, end1 = end1, end2
        if params["comcomp"] >= 1:
            if end1[-1] == "/" and end2[-1] == "/":
                end1 = end1[:-1]
                end2 = end2[:-1]
        if params["comcomp"] == 2:
            if end1[-1] == "/" and end2[-1] != "/":
                end1 = end1[:-1]
                end2 = end2 + "/"
        for s in senstype:
            if s == "0GO":
                continue
            sensdict[s].add((end1, end2))
            if s == "2GO":
                sensdict["1GO"].add((end1, end2))
            if s == "2NGO":
                sensdict["1NGO"].add((end1, end2))

    if onlytop:
        sensdict["1NGO"].difference_update(sensdict["1GO"])
        sensdict["1NGO"].difference_update(sensdict["2NGO"])
        if not sensdict["1GO"].issubset(sensdict["2NGO"]):
            log.warning(
                "1GO is not a subset of 2NGO. Contains {}.".format(
                    sensdict["1GO"].difference(sensdict["2NGO"])
                )
            )
        sensdict["2NGO"].difference_update(sensdict["1GO"])
        sensdict["1GO"].difference_update(sensdict["2GO"])

    return sensdict


def comp(edge):
    if not (type(edge) is np.ndarray):
        if edge[-1] == "/":
            return edge[:-1]
        else:
            return edge + "/"
    if type(edge) is np.ndarray:
        return np.array([comp(x) for x in edge])
