# New sensitivity code.
from collections import Counter

from .tiles import Tile, TileList
from .tilestructures import tile_daoe_single
from .util import GlueMergeSpec, TileMergeSpec, comp


def _fakesingle(tile):
    if not tile.structure.double:
        return TileList([tile])
    ret = TileList()
    c = ["", "/"]
    cn = ["_db_a", "_db_b"]
    for qq, es in enumerate(tile.structure.singleends):
        ft = Tile()
        ft.structure = tile_daoe_single()
        ft.structure._endlocs = [None, None, None, None]
        ft["ends"] = []
        ft["name"] = tile.name + cn[qq]
        if "input" in tile.keys():
            ft["input"] = []
        ft.structure._endtypes = []
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


def sensitivity_profiles_fakesingles(
    tileset,
    _maxorder=2,
    ms=GlueMergeSpec([]),
    tms=TileMergeSpec([]),
    oldclasses=None,
    checks=[],
    stopfirst=False,
):
    singles = _fakesingles(tileset.tiles)
    rotatedsingles = singles + _fakesingles(
        sum([x.rotations for x in tileset.tiles], TileList())
    )

    spairs = {"1GO": set(), "1NGO": set(), "2NGO": set(), "2GO": set()}
    if _maxorder >= 3:
        spairs["22GO"] = set()
        spairs["22NGO"] = set()

    occ = dict()

    # Convert tile names
    if oldclasses:
        for k, v in oldclasses.items():
            occ[k] = {(tms._cm.get(t1, t1), tms._cm.get(t2, t2)) for t1, t2 in v}

    for t1 in singles:
        # Iterate glues searching for a shared glue
        for i in range(0, 4):
            # If the glue is fake, then continue
            if t1.structure._endtypes[i] in {"fakedouble", "hairpin"}:
                continue
            for t2 in rotatedsingles:
                if not ms.eq(t1["ends"][i], t2["ends"][i]):
                    # not shared glue
                    continue
                if all(ms.eq(x, y) for x, y in zip(t1["ends"], t2["ends"])):
                    # NoIO(t1) = NoIO(t2)
                    continue
                # shared glue and not same tile: search for
                # viable unshared glue
                for j in range(0, 4):
                    if i == j:
                        # same edge
                        continue
                    if ms.eq(t1["ends"][j], t2["ends"][j]):
                        continue
                    if t1.structure._endtypes[j] != t2.structure._endtypes[j]:
                        continue
                    if t1.structure._endtypes[j] in {"fakedouble", "hairpin"}:
                        continue
                    spairs["1NGO"].add((t1.name, t2.name))
                    go1 = False
                    if t1["input"][i] and t1["input"][j]:
                        go1 = True
                        if (
                            stopfirst
                            and "1GO" in checks
                            and (
                                (
                                    tms._cm.get(t1.name, t1.name),
                                    tms._cm.get(t2.name, t2.name),
                                )
                                not in occ["1GO"]
                            )
                        ):
                            return (t1, t2)
                        spairs["1GO"].add((t1.name, t2.name))
                    for k in set(range(0, 4)) - {i, j}:
                        if (t1.structure._endtypes[k] == "hairpin") or (
                            t2.structure._endtypes[k] == "hairpin"
                        ):
                            continue
                        ec1 = comp(t1["ends"][k])
                        ec2 = comp(t2["ends"][k])
                        kc = _rev[k]
                        for t12 in singles:
                            if not ms.eq(t12["ends"][kc], ec1):
                                continue
                            for t22 in rotatedsingles:
                                if not ms.eq(t22["ends"][kc], ec2):
                                    continue
                                for m in range(0, 4):
                                    if m == kc:
                                        continue
                                    if t12.structure._endtypes[m] in {
                                        "fakedouble",
                                        "hairpin",
                                    }:
                                        continue
                                    if not ms.eq(t22["ends"][m], t12["ends"][m]):
                                        continue
                                    spairs["2NGO"].add((t1.name, t2.name))
                                    if t12["input"][m] and t12["input"][kc] and go1:
                                        spairs["2GO"].add((t1.name, t2.name))
                                        if (
                                            stopfirst
                                            and "2GO" in checks
                                            and (
                                                (
                                                    tms._cm.get(t1.name, t1.name),
                                                    tms._cm.get(t2.name, t2.name),
                                                )
                                                not in occ["2GO"]
                                            )
                                        ):
                                            return (t1, t2)
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
                                                spairs["22NGO"].add((t1.name, t2.name))
                                                spairs["22GO"].add((t1.name, t2.name))
                                                if (
                                                    stopfirst
                                                    and "22GO" in checks
                                                    and (
                                                        (
                                                            tms._cm.get(
                                                                t1.name, t1.name
                                                            ),
                                                            tms._cm.get(
                                                                t2.name, t2.name
                                                            ),
                                                        )
                                                        not in occ["22GO"]
                                                    )
                                                ):
                                                    return (t1, t2)
                                                continue

                                            ec1_3 = comp(t1["ends"][k2])
                                            ec2_3 = comp(t2["ends"][k2])
                                            kc2 = _rev[k2]
                                            for t13 in singles:
                                                if not ms.eq(t13["ends"][kc2], ec1_3):
                                                    continue
                                                for t23 in rotatedsingles:
                                                    if not ms.eq(
                                                        t23["ends"][kc2], ec2_3
                                                    ):
                                                        continue
                                                    for m2 in range(0, 4):
                                                        if m2 == kc2:
                                                            continue
                                                        if not ms.eq(
                                                            t23["ends"][m2],
                                                            t13["ends"][m2],
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
                                                            (t1.name, t2.name)
                                                        )
                                                        if (
                                                            t12["input"][m]
                                                            and t12["input"][kc]
                                                            and t13["input"][m2]
                                                            and t13["input"][kc2]
                                                            and go1
                                                        ):
                                                            if (
                                                                stopfirst
                                                                and "22GO" in checks
                                                                and (
                                                                    (
                                                                        tms._cm.get(
                                                                            t1.name,
                                                                            t1.name,
                                                                        ),
                                                                        tms._cm.get(
                                                                            t2.name,
                                                                            t2.name,
                                                                        ),
                                                                    )
                                                                    not in occ["22GO"]
                                                                )
                                                            ):
                                                                return (t1, t2)
                                                            spairs["22GO"].add(
                                                                (t1.name, t2.name)
                                                            )

    if not stopfirst:
        return spairs
    else:
        return False
