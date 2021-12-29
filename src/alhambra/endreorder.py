# End reordering code to optimize placement.

import random as pyrand
from copy import deepcopy

import numpy as np
import numpy.ma as ma
import numpy.random as random
import stickydesign as sd

from . import anneal
from . import sensitivity as sens


def flatten(seq):
    for item in seq:
        if is_string_like(item) or not iterable(item):
            yield item
        else:
            for subitem in flatten(item):
                yield subitem


def ecomp(x):
    if x[-1] == "/":
        return x[:-1]
    else:
        return x + "/"


class FseqState:
    def __init__(self, seqs=None):
        if not seqs:
            self.seqs = {}
        else:
            self.seqs = seqs

    def copy(self):
        return FseqState({"DT": self.seqs["DT"].copy(), "TD": self.seqs["TD"].copy()})


class cachedarray:
    def __init__(self, func, shape):
        self.arr = ma.masked_all(shape)
        self.func = func

    def __getitem__(self, index):
        if self.arr.mask[index]:
            self.arr[index] = self.func(*index)
        return self.arr[index]


from copy import deepcopy


class FastState:
    def __init__(self, state):
        self.state = state

    def copy(self):
        return deepcopy(self)

    def __getitem__(self, it):
        return self.state[it]


class EndSystemFseq:
    def __init__(
        self,
        tilesys,
        newends=None,
        pairs=None,
        energetics=None,
        inputpairs=False,
        multiscore=False,
    ):
        tilesys = deepcopy(tilesys)
        self.ends = tilesys.ends
        self.tiles = tilesys.tiles
        self.tilesystem = tilesys
        self.ef = energetics
        self.multiscore = multiscore
        if not pairs:
            pairs = sens.consolidate_pairs(
                sens.senspairs(tilesys), comcomp=1, onlytop=True
            )

        if inputpairs:
            inputpairs = [
                tuple(z[:-1] for z, i in zip(x["ends"], x["input"]) if i)
                for x in self.tiles
            ]
            print(inputpairs)
        self.inputpairs = inputpairs

        self.names = {}
        fseqsTD, self.names["TD"] = (
            list(x)
            for x in zip(
                *[
                    [end["fseq"].lower(), end["name"]]
                    for end in self.ends
                    if end["type"] == "TD"
                ]
            )
        )
        fseqsDT, self.names["DT"] = (
            list(x)
            for x in zip(
                *[
                    [end["fseq"].lower(), end["name"]]
                    for end in self.ends
                    if end["type"] == "DT"
                ]
            )
        )
        self.seqs = {}
        self.seqs["TD"] = sd.endarray(fseqsTD, "TD")
        self.seqs["DT"] = sd.endarray(fseqsDT, "DT")
        self.initstate = FastState(
            {
                "DT": np.arange(0, len(self.seqs["DT"])),
                "TD": np.arange(0, len(self.seqs["TD"])),
            }
        )
        self.enlocs = {}
        for i, endn in enumerate(self.names["TD"]):
            self.enlocs[endn] = (i, "TD")
        for i, endn in enumerate(self.names["DT"]):
            self.enlocs[endn] = (i, "DT")

        # ends that can be reordered.
        if newends:
            self.mutableTD = [
                i for i, t in [self.enlocs[x] for x in newends] if t == "TD"
            ]
            self.mutableDT = [
                i for i, t in [self.enlocs[x] for x in newends] if t == "DT"
            ]
        else:
            self.mutableTD = range(0, len(fseqsTD))
            self.mutableDT = range(0, len(fseqsDT))

        self.pairdict = {}
        for pairclass, memberset in pairs.items():
            for x, y in memberset:
                self.pairdict[(x, ecomp(y))] = pairclass
                self.pairdict[(y, ecomp(x))] = pairclass
        tdsh = (len(self.seqs["TD"]), len(self.seqs["TD"]))
        dtsh = (len(self.seqs["DT"]), len(self.seqs["DT"]))

        if not multiscore:
            # Get the mean non-spurious interaction
            self.meangse = 0.5 * (
                np.mean(self.ef.matching_uniform(self.seqs["TD"]))
                + np.mean(self.ef.matching_uniform(self.seqs["DT"]))
            )

            self.mult = {
                "1NGO": np.exp(-2.0 * self.meangse),
                "2NGO": np.exp(-1.65 * self.meangse),
                "1GO": np.exp(-1.5 * self.meangse),
                "2GO": np.exp(-1.1 * self.meangse),
                "I": np.exp(-1.0 * self.meangse),
            }

            self.ecache_cc = {
                "TD": cachedarray(
                    lambda x, y: self.ef.uniform(
                        self.seqs["TD"][x : x + 1].comps,
                        self.seqs["TD"][y : y + 1].comps,
                    ),
                    tdsh,
                ),
                "DT": cachedarray(
                    lambda x, y: self.ef.uniform(
                        self.seqs["DT"][x : x + 1].comps,
                        self.seqs["DT"][y : y + 1].comps,
                    ),
                    dtsh,
                ),
            }
            self.ecache_ce = {
                "TD": cachedarray(
                    lambda x, y: self.ef.uniform(
                        self.seqs["TD"][x : x + 1].comps,
                        self.seqs["TD"][y : y + 1].ends,
                    ),
                    tdsh,
                ),
                "DT": cachedarray(
                    lambda x, y: self.ef.uniform(
                        self.seqs["DT"][x : x + 1].comps,
                        self.seqs["DT"][y : y + 1].ends,
                    ),
                    dtsh,
                ),
            }
            self.ecache_ec = {
                "TD": cachedarray(
                    lambda x, y: self.ef.uniform(
                        self.seqs["TD"][x : x + 1].ends,
                        self.seqs["TD"][y : y + 1].comps,
                    ),
                    tdsh,
                ),
                "DT": cachedarray(
                    lambda x, y: self.ef.uniform(
                        self.seqs["DT"][x : x + 1].ends,
                        self.seqs["DT"][y : y + 1].comps,
                    ),
                    dtsh,
                ),
            }
            self.ecache_ee = {
                "TD": cachedarray(
                    lambda x, y: self.ef.uniform(
                        self.seqs["TD"][x : x + 1].ends, self.seqs["TD"][y : y + 1].ends
                    ),
                    tdsh,
                ),
                "DT": cachedarray(
                    lambda x, y: self.ef.uniform(
                        self.seqs["DT"][x : x + 1].ends, self.seqs["DT"][y : y + 1].ends
                    ),
                    dtsh,
                ),
            }
        else:
            self.ecache_cc = []
            self.ecache_ee = []
            self.ecache_ce = []
            self.ecache_ec = []
            self.mult = []
            self.meangse = []
            for ef in self.ef:
                self.meangse.append(
                    0.5
                    * (
                        np.mean(ef.matching_uniform(self.seqs["TD"]))
                        + np.mean(ef.matching_uniform(self.seqs["DT"]))
                    )
                )

                self.mult.append(
                    {
                        "1NGO": np.exp(-2.0 * self.meangse[-1]),
                        "2NGO": np.exp(-1.65 * self.meangse[-1]),
                        "1GO": np.exp(-1.5 * self.meangse[-1]),
                        "2GO": np.exp(-1.1 * self.meangse[-1]),
                        "I": np.exp(-1.0 * self.meangse[-1]),
                    }
                )

                self.ecache_cc.append(
                    {
                        "TD": cachedarray(
                            lambda x, y: ef.uniform(
                                self.seqs["TD"][x : x + 1].comps,
                                self.seqs["TD"][y : y + 1].comps,
                            ),
                            tdsh,
                        ),
                        "DT": cachedarray(
                            lambda x, y: ef.uniform(
                                self.seqs["DT"][x : x + 1].comps,
                                self.seqs["DT"][y : y + 1].comps,
                            ),
                            dtsh,
                        ),
                    }
                )
                self.ecache_ce.append(
                    {
                        "TD": cachedarray(
                            lambda x, y: ef.uniform(
                                self.seqs["TD"][x : x + 1].comps,
                                self.seqs["TD"][y : y + 1].ends,
                            ),
                            tdsh,
                        ),
                        "DT": cachedarray(
                            lambda x, y: ef.uniform(
                                self.seqs["DT"][x : x + 1].comps,
                                self.seqs["DT"][y : y + 1].ends,
                            ),
                            dtsh,
                        ),
                    }
                )
                self.ecache_ec.append(
                    {
                        "TD": cachedarray(
                            lambda x, y: ef.uniform(
                                self.seqs["TD"][x : x + 1].ends,
                                self.seqs["TD"][y : y + 1].comps,
                            ),
                            tdsh,
                        ),
                        "DT": cachedarray(
                            lambda x, y: ef.uniform(
                                self.seqs["DT"][x : x + 1].ends,
                                self.seqs["DT"][y : y + 1].comps,
                            ),
                            dtsh,
                        ),
                    }
                )
                self.ecache_ee.append(
                    {
                        "TD": cachedarray(
                            lambda x, y: ef.uniform(
                                self.seqs["TD"][x : x + 1].ends,
                                self.seqs["TD"][y : y + 1].ends,
                            ),
                            tdsh,
                        ),
                        "DT": cachedarray(
                            lambda x, y: ef.uniform(
                                self.seqs["DT"][x : x + 1].ends,
                                self.seqs["DT"][y : y + 1].ends,
                            ),
                            dtsh,
                        ),
                    }
                )

    def slowseqs(self, state):
        "Give the state as the slow version would have"
        return {"DT": self.seqs["DT"][state["DT"]], "TD": self.seqs["TD"][state["TD"]]}

    def mutate(self, state):
        # Start by deciding to swap TD or DT ends.
        if random.rand() > 1.0 * len(self.mutableTD) / (
            len(self.mutableDT) + len(self.mutableTD)
        ):
            a = random.choice(self.mutableDT)
            b = random.choice(self.mutableDT)
            state["DT"][[a, b]] = state["DT"][[b, a]]
        else:
            a = random.choice(self.mutableTD)
            b = random.choice(self.mutableTD)
            state["TD"][[a, b]] = state["TD"][[b, a]]

    def score(self, state):

        sc = 0.0
        if self.inputpairs:
            if self.multiscore:
                ip = np.zeros((len(self.ef), len(self.inputpairs)))
            else:
                ip = np.zeros(len(self.inputpairs))
            for i, (xn, yn) in enumerate(self.inputpairs):
                if xn[-1] == "/":
                    xn = xn[:-1]
                if yn[-1] == "/":
                    yn = yn[:-1]
                xi, xt = self.enlocs[xn]
                yi, yt = self.enlocs[yn]
                if self.multiscore:
                    scc = 0
                    for j in range(len(self.ef)):
                        ip[j][i] = np.abs(
                            (
                                self.ecache_ec[j][xt][state[xt][xi], state[xt][xi]]
                                + self.ecache_ec[j][yt][state[yt][yi], state[yt][yi]]
                            )
                        )
                else:
                    ip[i] = np.abs(
                        (
                            self.ecache_ec[xt][state[xt][xi], state[xt][xi]]
                            + self.ecache_ec[yt][state[yt][yi], state[yt][yi]]
                        )
                    )
            if self.multiscore:
                sc += np.mean(
                    np.array([z["I"] for z in self.mult]) * np.exp(np.ptp(ip, axis=1))
                )
            else:
                sc += self.mult["I"] * np.exp(np.ptp(ip))

        for (xn, yn), pairclass in self.pairdict.items():
            # set comp flags
            xc = False
            yc = False
            if xn[-1] == "/":
                xc = True
                xn = xn[:-1]
            if yn[-1] == "/":
                yc = True
                yn = yn[:-1]

            # get end indexes and types
            xi, xt = self.enlocs[xn]
            yi, yt = self.enlocs[yn]
            # print "%s, %s: (%s %s)" % (xn,yn,xt,yt)
            # skip if not same type
            if xt != yt:
                continue
            if not self.multiscore:
                if yc and xc:
                    val = self.ecache_cc[xt][
                        state[xt][xi], state[yt][yi]
                    ]  # self.ef.uniform(state.seqs[xt][xi:xi+1].comps,state.seqs[yt][yi:yi+1].comps)[0]
                elif xc:
                    val = self.ecache_ce[xt][
                        state[xt][xi], state[yt][yi]
                    ]  # self.ef.uniform(state.seqs[xt][xi:xi+1].comps,state.seqs[yt][yi:yi+1].ends)[0]
                elif yc:
                    val = self.ecache_ec[xt][
                        state[xt][xi], state[yt][yi]
                    ]  # state[f.uniform(state.seqs[xt][xi:xi+1].ends,state.seqs[yt][yi:yi+1].comps)[0]
                else:
                    val = self.ecache_ee[xt][
                        state[xt][xi], state[yt][yi]
                    ]  # self.ef.uniform(state.seqs[xt][xi:xi+1].ends,state.seqs[yt][yi:yi+1].ends)[0]

                sc += self.mult[pairclass] * np.exp(val)
            else:
                scc = 0
                for i in range(len(self.ef)):
                    if yc and xc:
                        val = self.ecache_cc[i][xt][
                            state[xt][xi], state[yt][yi]
                        ]  # self.ef.uniform(state.seqs[xt][xi:xi+1].comps,state.seqs[yt][yi:yi+1].comps)[0]
                    elif xc:
                        val = self.ecache_ce[i][xt][
                            state[xt][xi], state[yt][yi]
                        ]  # self.ef.uniform(state.seqs[xt][xi:xi+1].comps,state.seqs[yt][yi:yi+1].ends)[0]
                    elif yc:
                        val = self.ecache_ec[i][xt][
                            state[xt][xi], state[yt][yi]
                        ]  # state[f.uniform(state.seqs[xt][xi:xi+1].ends,state.seqs[yt][yi:yi+1].comps)[0]
                    else:
                        val = self.ecache_ee[i][xt][
                            state[xt][xi], state[yt][yi]
                        ]  # self.ef.uniform(state.seqs[xt][xi:xi+1].ends,state.seqs[yt][yi:yi+1].ends)[0]

                    scc += self.mult[i][pairclass] * np.exp(val)
                sc += scc / len(self.ef)

        return sc


wcd = {
    "a": "t",
    "b": "v",
    "c": "g",
    "d": "h",
    "g": "c",
    "h": "d",
    "k": "m",
    "m": "k",
    "n": "n",
    "s": "s",
    "t": "a",
    "v": "b",
    "w": "w",
}


def wc(seqstr):
    return "".join(wcd[x] for x in reversed(seqstr))
