import numpy as np
import copy
import time
from typing import Generic, MutableMapping, Protocol, TypeVar

from ruamel.yaml.comments import CommentedSeq
from ruamel.yaml.representer import RoundTripRepresenter

from . import seq

import stickydesign as sd
import stickydesign.stickydesign2 as sd2

DEFAULT_ENERGETICS = sd.EnergeticsDAOE(temperature=33, coaxparams=True)

DEFAULT_MULTIMODEL_ENERGETICS = [
    sd.EnergeticsDAOE(temperature=33, coaxparams="protozanova"),
    sd.EnergeticsDAOE(temperature=33, coaxparams="pyshni"),
    sd.EnergeticsDAOE(temperature=33, coaxparams="peyret"),
    sd.EnergeticsDAOE(temperature=33, coaxparams=False),
]

DEFAULT_SD2_MULTIMODEL_ENERGETICS = [
    sd2.EnergeticsDAOEC(5, temperature=36, coaxparams="protozanova"),
    sd2.EnergeticsDAOEC(5, temperature=36, coaxparams="pyshni"),
    sd2.EnergeticsDAOEC(5, temperature=36, coaxparams="peyret"),
    sd2.EnergeticsDAOEC(5, temperature=36, coaxparams=False),
]

DEFAULT_MM_ENERGETICS_NAMES = ["Prot", "Pysh", "Peyr", "None"]

DEFAULT_REGION_ENERGETICS = sd.EnergeticsBasic(
    temperature=33, coaxparams=False, danglecorr=False
)

MAPPER = map

MPOBJECT = None


def multimap(function, data):
    return MAPPER(function, data)


def setup_multi(method, ncores=None):
    global MPOBJECT
    global MAPPER
    if method == "none":
        MAPPER = map
    elif method == "multiprocessing":
        import multiprocessing

        if not ncores:
            import os

            ncores = os.cpu_count() - 1
        if not isinstance(MPOBJECT, multiprocessing.pool.Pool):
            MPOBJECT = multiprocessing.Pool(ncores)
        MAPPER = MPOBJECT.map
    elif method == "ipyparallel":
        import ipyparallel

        if not isinstance(MPOBJECT, ipyparallel.client.view.LoadBalancedView):
            rc = ipyparallel.Client()
            MPOBJECT = rc.load_balanced_view()
        MAPPER = MPOBJECT.map


class NamedList(CommentedSeq):
    """A class for a list of dicts, where some dicts have a 'name' item
    which should be unique, but others might not.  Indexing works with
    either number or name.  Note that updating dicts may make the list
    inconsistent.

    """

    def __init__(self, x=[]):
        CommentedSeq.__init__(self, x)

    def __getitem__(self, i):
        if isinstance(i, str):
            r = [x for x in self if x.get("name", None) == i]
            if len(r) > 1:
                raise KeyError("There are {} elements named {}.".format(len(r), i))
            elif len(r) == 0:
                raise KeyError("No element named {} found.".format(i))
            else:
                return r[0]
        else:
            return CommentedSeq.__getitem__(self, i)

    def check_consistent(self):
        """Checks that each name appears only once.  On failure, returns a
        ValueError with (message, {failed_name: count}).  Otherwise, return
        with no output.

        """
        names = [v["name"] for v in self if "name" in v.keys()]

        # if we have *no* names, the tests will fail, but we are obviously
        # consistent, so:
        if not names:
            return

        from collections import Counter

        namecounts = Counter(names)

        if max(namecounts.values()) > 1:
            badcounts = {n: v for n, v in namecounts.items() if v > 1}
            raise ValueError("Inconsistent NamedList.", badcounts)

    def __setitem__(self, i, v):
        if isinstance(i, str):
            r = [(ii, x) for ii, x in enumerate(self) if x.get("name", None) == i]
            if len(r) > 1:
                raise KeyError("There are {} elements named {}.".format(len(r), i))
            elif len(r) == 0:
                self.append(v)
            else:
                CommentedSeq.__setitem__(self, r[0][0], v)
        else:
            CommentedSeq.__setitem__(self, i, v)

    def __delitem__(self, i):
        if isinstance(i, str):
            r = [ii for ii, x in enumerate(self) if x.get("name", None) == i]
            if len(r) > 1:
                raise KeyError("There are {} elements named {}.".format(len(r), i))
            elif len(r) == 0:
                raise KeyError("No element named {} found.".format(i))
            else:
                CommentedSeq.__delitem__(self, r[0])
        else:
            return CommentedSeq.__delitem__(self, i)

    def keys(self):
        return [x["name"] for x in self if "name" in x.keys()]


RoundTripRepresenter.add_representer(NamedList, RoundTripRepresenter.represent_list)

RoundTripRepresenter.add_representer(np.str_, RoundTripRepresenter.represent_str)


class ProgressLogger(object):
    def __init__(self, logger, N, seconds_interval=60):
        self.logger = logger
        stime = time.perf_counter()
        self.stime = stime
        self.ltime = stime
        self.li = 0
        self
        self.seconds_interval = seconds_interval
        self.N = N
        self.logger.info("starting {} tasks".format(self.N))

    def update(self, i):
        ctime = time.perf_counter()
        if ctime - self.ltime > self.seconds_interval:
            self.logger.info(
                "finished {}/{}, {} s elapsed, {} s est remaining".format(
                    i,
                    self.N,
                    int(ctime - self.stime),
                    int((self.N - i) * (ctime - self.stime) / i),
                )
            )
            self.ltime = ctime
            self.li = i


def comp(endname):
    "Return the complementary name of a given end (eg, for 'a', return 'a/')"
    if endname[-1] == "/":
        return endname[:-1]
    else:
        return endname + "/"


def base(endname):
    "Return the base name of a given end name (eg, for either 'a' or 'a/', return 'a')"
    if endname[-1] == "/":
        return endname[:-1]
    else:
        return endname


class GlueMergeSpec:
    def __init__(self, ops=[]):
        self._cm = {}
        self._ecs = []
        for op in ops:
            self.add(*op)

    def add(self, a, b):
        ai = self._cm.get(a, -1)
        bi = self._cm.get(b, -1)
        if (ai > -1) and (bi == -1):
            self._ecs[ai].add(b)
            self._ecs[self._cm[comp(a)]].add(comp(b))
            self._rebuild_map()
        if (ai == -1) and (bi > -1):
            self._ecs[bi].add(a)
            self._ecs[self._cm[comp(b)]].add(comp(a))
            self._rebuild_map()
        if (ai == -1) and (bi == -1):
            self._ecs.append({a, b})
            self._ecs.append({comp(a), comp(b)})
            self._rebuild_map()
        if (ai > -1) and (bi > -1):
            if ai == bi:
                return
            if ai == self._cm[comp(b)]:
                raise ValueError
            self._ecs[ai].update(self._ecs[bi])
            self._ecs[self._cm[comp(a)]].update(self._ecs[self._cm[comp(b)]])
            del self._ecs[bi]
            if self._cm[comp(b)] > bi:
                del self._ecs[self._cm[comp(b)] - 1]
            else:
                del self._ecs[self._cm[comp(b)]]
            self._rebuild_map()

    def copyadd(self, a, b):
        v = copy.deepcopy(self)
        v.add(a, b)
        return v

    def _rebuild_map(self):
        for i, s in enumerate(self._ecs):
            for v in s:
                if comp(v) in s:
                    raise ValueError
                self._cm[v] = i

    def eq(self, a, b):
        if a == b:
            return True
        elif (
            (a in self._cm.keys())
            and (b in self._cm.keys())
            and (self._cm[a] == self._cm[b])
        ):
            return True
        else:
            return False


class TileMergeSpec:
    def __init__(self, ops=[]):
        self._cm = {}
        self._ecs = []
        for op in ops:
            self.add(*op)

    def add(self, a, b):
        ai = self._cm.get(a, -1)
        bi = self._cm.get(b, -1)
        if (ai > -1) and (bi == -1):
            self._ecs[ai].add(b)
            self._rebuild_map()
        if (ai == -1) and (bi > -1):
            self._ecs[bi].add(a)
            self._rebuild_map()
        if (ai == -1) and (bi == -1):
            self._ecs.append({a, b})
            self._rebuild_map()
        if (ai > -1) and (bi > -1):
            if ai == bi:
                return
            self._ecs[ai].update(self._ecs[bi])
            del self._ecs[bi]
            self._rebuild_map()

    def copyadd(self, a, b):
        v = copy.deepcopy(self)
        v.add(a, b)
        return v

    def _rebuild_map(self):
        for i, s in enumerate(self._ecs):
            for v in s:
                self._cm[v] = i

    def eq(self, a, b):
        if a == b:
            return True
        elif (
            (a in self._cm.keys())
            and (b in self._cm.keys())
            and (self._cm[a] == self._cm[b])
        ):
            return True
        else:
            return False
