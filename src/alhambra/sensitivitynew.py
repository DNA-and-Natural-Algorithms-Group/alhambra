# type: ignore
# New sensitivity code.
# from .tilestructures import tile_daoe_single
from collections import Counter

from .tiles import HDupleTile, Tile, TileList, VDupleTile
from .util import comp


class FakeSingle(Tile):
    pass


class CounterSet(Counter):
    def add(self, val, valthrowaway):
        self[val] += 1


class SetSet(set):
    def add(self, val, valthrowaway):
        set.add(self, val)


class PathSet(dict):
    def add(self, key, val):
        if key not in self.keys():
            self[key] = set()
        self[key].add(val)


def _fakesingle(tile: Tile):
    if not isinstance(tile, (HDupleTile, VDupleTile)):
        return TileList([tile])
    ret = TileList()
    c = ["", "*"]
    cn = ["_db_a", "_db_b"]
    for qq, es in enumerate(tile.structure.singleends):
        ft = FakeSingle()
        ft.structure._endlocs = [None, None, None, None]
        ft["ends"] = []
        ft["name"] = tile.name + cn[qq]
        if "input" in tile.keys():
            ft["input"] = []
        for i in es:
            if i is not None:
                ft["ends"].append(tile["ends"][i])
                ft.structure._endtypes.append(tile.structure._endtypes[i])
                if "input" in tile.keys():
                    ft["input"].append(tile["input"][i])
            else:
                ft["ends"].append(tile.name + "_db" + c[qq])
                ft.structure._endtypes.append("fakedouble")
                if "input" in tile.keys():
                    ft["input"].append(0)
        ret.append(ft)
    return ret


def _fakesingles(tiles):
    return sum((_fakesingle(x) for x in tiles), TileList())


_rev = [2, 3, 0, 1]


def sensitivity_classes(tileset, count=False, _maxorder=2):

    singles = _fakesingles(tileset.tiles)
    rotatedsingles = singles + _fakesingles(
        sum([x.rotations for x in tileset.tiles], TileList())
    )

    if count == "paths":
        sclass = PathSet
    elif count:
        sclass = CounterSet
    else:
        sclass = SetSet

    spairs = {"1GO": sclass(), "1NGO": sclass(), "2NGO": sclass(), "2GO": sclass()}

    if _maxorder >= 3:
        spairs["22NGO"] = sclass()
        spairs["22GO"] = sclass()
        maxc = "22GO"
    else:
        maxc = "2GO"

    for t1 in singles:
        for i in range(0, 4):
            if t1.structure._endtypes[i] in {"fakedouble", "hairpin"}:
                continue
            for t2 in rotatedsingles:
                if t2.edges[i] != t1.edges[i]:
                    continue
                for j in range(0, 4):
                    if i == j:
                        continue
                    if frozenset((t2["ends"][j], t1["ends"][j])) in spairs[maxc]:
                        continue
                        # don't go through the rest of the process
                        # if the pair is already in the "worst"
                        # class.
                    if t1["ends"][j] == t2["ends"][j]:
                        continue
                    if t1.structure._endtypes[j] != t2.structure._endtypes[j]:
                        continue
                    if t1.structure._endtypes[j] in {"fakedouble", "hairpin"}:
                        continue
                    spairs["1NGO"].add(
                        frozenset((t2["ends"][j], t1["ends"][j])), (t1.name, t2.name)
                    )
                    spairs["1NGO"].add(
                        frozenset((comp(t2["ends"][j]), comp(t1["ends"][j]))),
                        (t1.name, t2.name),
                    )
                    go1 = False
                    if t1["input"][i] and t1["input"][j]:
                        go1 = True
                        spairs["1GO"].add(
                            frozenset((t2["ends"][j], t1["ends"][j])),
                            (t1.name, t2.name),
                        )
                        spairs["1GO"].add(
                            frozenset((comp(t2["ends"][j]), comp(t1["ends"][j]))),
                            (t1.name, t2.name),
                        )
                    for k in set(range(0, 4)) - {i, j}:
                        if (t1.structure._endtypes[k] == "hairpin") or (
                            t2.structure._endtypes[k] == "hairpin"
                        ):
                            continue
                        ec1 = comp(t1["ends"][k])
                        ec2 = comp(t2["ends"][k])
                        kc = _rev[k]
                        for t12 in singles:
                            if t12["ends"][kc] != ec1:
                                continue
                            for t22 in rotatedsingles:
                                if t22["ends"][kc] != ec2:
                                    continue
                                for m in range(0, 4):
                                    if m == kc:
                                        continue
                                    if t22["ends"][m] != t12["ends"][m]:
                                        continue
                                    if t12.structure._endtypes[m] in {
                                        "fakedouble",
                                        "hairpin",
                                    }:
                                        continue
                                    spairs["2NGO"].add(
                                        frozenset((t2["ends"][j], t1["ends"][j])),
                                        (t1.name, t2.name, t12.name, t22.name),
                                    )
                                    spairs["2NGO"].add(
                                        frozenset(
                                            (comp(t2["ends"][j]), comp(t1["ends"][j]))
                                        ),
                                        (t1.name, t2.name, t12.name, t22.name),
                                    )
                                    if t12["input"][m] and t12["input"][kc] and go1:
                                        spairs["2GO"].add(
                                            frozenset((t2["ends"][j], t1["ends"][j])),
                                            (t1.name, t2.name, t12.name, t22.name),
                                        )
                                        spairs["2GO"].add(
                                            frozenset(
                                                (
                                                    comp(t2["ends"][j]),
                                                    comp(t1["ends"][j]),
                                                )
                                            ),
                                            (t1.name, t2.name, t12.name, t22.name),
                                        )
                                    if _maxorder > 2:
                                        for k2 in set(range(0, 4)) - {
                                            i,
                                            j,
                                            k,
                                        }:  # should just be one number!
                                            # hairpin? then continue.
                                            if (
                                                t1.structure._endtypes[k2] == "hairpin"
                                            ) or (
                                                t2.structure._endtypes[k2] == "hairpin"
                                            ):
                                                continue
                                            # fakedouble on either: then assume we're doomed: (FIXME)
                                            if (
                                                t1.structure._endtypes[k2]
                                                == "fakedouble"
                                            ) or (
                                                t2.structure._endtypes[k2]
                                                == "fakedouble"
                                            ):
                                                spairs["22NGO"].add(
                                                    frozenset(
                                                        (t2["ends"][j], t1["ends"][j])
                                                    ),
                                                    (
                                                        t1.name,
                                                        t2.name,
                                                        t12.name,
                                                        t22.name,
                                                    ),
                                                )
                                                spairs["22NGO"].add(
                                                    frozenset(
                                                        (
                                                            comp(t2["ends"][j]),
                                                            comp(t1["ends"][j]),
                                                        )
                                                    ),
                                                    (
                                                        t1.name,
                                                        t2.name,
                                                        t12.name,
                                                        t22.name,
                                                    ),
                                                )
                                                spairs["22GO"].add(
                                                    frozenset(
                                                        (t2["ends"][j], t1["ends"][j])
                                                    ),
                                                    (
                                                        t1.name,
                                                        t2.name,
                                                        t12.name,
                                                        t22.name,
                                                    ),
                                                )
                                                spairs["22GO"].add(
                                                    frozenset(
                                                        (
                                                            comp(t2["ends"][j]),
                                                            comp(t1["ends"][j]),
                                                        )
                                                    ),
                                                    (
                                                        t1.name,
                                                        t2.name,
                                                        t12.name,
                                                        t22.name,
                                                    ),
                                                )
                                                continue

                                            ec1_3 = comp(t1["ends"][k2])
                                            ec2_3 = comp(t2["ends"][k2])
                                            kc2 = _rev[k2]
                                            for t13 in singles:
                                                if t13["ends"][kc2] != ec1_3:
                                                    continue
                                                for t23 in rotatedsingles:
                                                    if t23["ends"][kc2] != ec2_3:
                                                        continue
                                                    for m2 in range(0, 4):
                                                        if m2 == kc2:
                                                            continue
                                                        if (
                                                            t23["ends"][m2]
                                                            != t13["ends"][m2]
                                                        ):
                                                            continue
                                                        if t13.structure._endtypes[
                                                            m2
                                                        ] in {
                                                            "fakedouble",
                                                            "hairpin",
                                                        }:
                                                            continue
                                                        spairs["22NGO"].add(
                                                            frozenset(
                                                                (
                                                                    t2["ends"][j],
                                                                    t1["ends"][j],
                                                                )
                                                            ),
                                                            (
                                                                t1.name,
                                                                t2.name,
                                                                t12.name,
                                                                t22.name,
                                                                t13.name,
                                                                t23.name,
                                                            ),
                                                        )
                                                        spairs["22NGO"].add(
                                                            frozenset(
                                                                (
                                                                    comp(t2["ends"][j]),
                                                                    comp(t1["ends"][j]),
                                                                )
                                                            ),
                                                            (
                                                                t1.name,
                                                                t2.name,
                                                                t12.name,
                                                                t22.name,
                                                                t13.name,
                                                                t23.name,
                                                            ),
                                                        )
                                                        if (
                                                            t12["input"][m]
                                                            and t12["input"][kc]
                                                            and t13["input"][m2]
                                                            and t13["input"][kc2]
                                                            and go1
                                                        ):
                                                            spairs["22GO"].add(
                                                                frozenset(
                                                                    (
                                                                        t2["ends"][j],
                                                                        t1["ends"][j],
                                                                    )
                                                                ),
                                                                (
                                                                    t1.name,
                                                                    t2.name,
                                                                    t12.name,
                                                                    t22.name,
                                                                    t13.name,
                                                                    t23.name,
                                                                ),
                                                            )
                                                            spairs["22GO"].add(
                                                                frozenset(
                                                                    (
                                                                        comp(
                                                                            t2["ends"][
                                                                                j
                                                                            ]
                                                                        ),
                                                                        comp(
                                                                            t1["ends"][
                                                                                j
                                                                            ]
                                                                        ),
                                                                    )
                                                                ),
                                                                (
                                                                    t1.name,
                                                                    t2.name,
                                                                    t12.name,
                                                                    t22.name,
                                                                    t13.name,
                                                                    t23.name,
                                                                ),
                                                            )

    return spairs
