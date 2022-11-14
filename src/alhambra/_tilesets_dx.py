"""DX Tile TileSet routines."""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import logging
from typing import (
    Optional,
    TYPE_CHECKING,
    Any,
    Callable,
    Iterable,
    Literal,
    Sequence,
    cast,
)
from random import shuffle
import stickydesign.multimodel as multimodel
import numpy as np
import stickydesign

SELOGGER = logging.getLogger(__name__)

from stickydesign.energetics_daoe import EnergeticsDAOE

from .glues import DXGlue, GlueList
from .util import (  # DEFAULT_ENERGETICS,; DEFAULT_SD2_MULTIMODEL_ENERGETICS,
    DEFAULT_MM_ENERGETICS_NAMES,
    DEFAULT_MULTIMODEL_ENERGETICS,
    DEFAULT_REGION_ENERGETICS,
    DEFAULT_ENERGETICS,
)
from . import seq, util

if TYPE_CHECKING:
    from .tilesets import TileSet


def stickydesign_create_dx_glue_sequences(
    self,
    method: Literal["default", "multimodel"] = "default",
    energetics: "stickydesign.Energetics" | None = None,
    trials: int = 100,
    devmethod="dev",
    sdopts: dict[str, Any] | None = None,
    ecpars: dict[str, Any] | None = None,
    listends: bool = False,
) -> tuple[TileSet, Sequence[str]]:
    """Create sticky end sequences for the TileSet, using stickydesign,
    and returning a new TileSet including the ends.


    Parameters
    ----------
    method: {'default', 'multimodel'}
        if 'default', use the default, single-model sequence design.
        If 'multimodel', use multimodel end choice.

    energetics : stickydesign.Energetics
        the energetics instance to use for the design, or list
        of energetics for method='multimodel', in which case the first
        will be the primary.  If None (default), will use
        alhambra.DEFAULT_ENERGETICS, or, if method='multimodel', will use
        alhambra.DEFAULT_MM_ENERGETICS.

    trials : int
        the number of trials to attempt. FIXME

    sdopts : dict
        a dictionary of parameters to pass to stickydesign.easy_ends.

    ecpars : dict
        a dictionary of parameters to pass to the endchooser function
        generator (useful only in limited circumstances).

    listends : bool
        if False, return just the TileSet.  If True, return both the
        TileSet and a list of the names of the ends that were created.

    Returns
    -------
    tileset : TileSet
        TileSet with designed end sequences included.
    new_ends : list
        Names of the ends that were designed.

    """
    if sdopts is None:
        sdopts = {}
    if ecpars is None:
        ecpars = {}
    info: dict[str, Any] = {}
    info["method"] = method
    info["time"] = datetime.now(tz=timezone.utc).isoformat()
    info["sd_version"] = stickydesign.version.__version__

    if not energetics:
        if method == "multimodel":
            all_energetics = DEFAULT_MULTIMODEL_ENERGETICS
        else:
            energetics = DEFAULT_ENERGETICS
    if method == "multimodel" and not isinstance(energetics, Iterable):
        raise ValueError("Energetics must be an iterable for multimodel.")
    elif method == "multimodel":
        all_energetics = cast(list[EnergeticsDAOE], energetics)
        energetics = all_energetics[0]
        info["energetics"] = [str(e) for e in all_energetics]
        info["trails"] = trials
    elif method == "default":
        info["energetics"] = str(energetics)
        energetics = cast(EnergeticsDAOE, energetics)
    else:
        raise ValueError

    # Steps for doing this:

    # Build a list of ends from the endlist in the tileset.  Do this
    # by creating a NamedList, then merging them into it.
    ends: GlueList[DXGlue] = GlueList()

    ends.update(g for g in self.glues if isinstance(g, DXGlue))

    ends.update(g for g in self.tiles.glues_from_tiles() if isinstance(g, DXGlue))

    # Ensure that if there are any resulting completely-undefined ends, they
    # have their sequences removed.
    # for end in ends:
    #    if end.fseq and set(end.fseq) == {'n'}:
    #        del(end.fseq)

    # Build inputs suitable for stickydesign: lists of old sequences for TD/DT,
    # and numbers of new sequences needed.
    oldDTseqs = [
        end.fseq for end in ends if end.etype == "DT" and seq.is_definite(end.fseq)
    ]
    if oldDTseqs:
        oldDTarray = stickydesign.endarray(oldDTseqs, "DT")
    else:
        oldDTarray = None
    oldTDseqs = [
        end.fseq for end in ends if end.etype == "TD" and seq.is_definite(end.fseq)
    ]
    if oldTDseqs:
        oldTDarray = stickydesign.endarray(oldTDseqs, "TD")
    else:
        oldTDarray = None

    newTD = [end for end in ends if end.etype == "TD" and not seq.is_definite(end.fseq)]
    newDT = [end for end in ends if end.etype == "DT" and not seq.is_definite(end.fseq)]

    # Deal with energetics, considering potential old sequences.
    # FIXME: EXPLAIN WHAT THIS ABSTRUSE CODE DOES...
    # TODO: tests needs to test this
    targets = []
    if len(oldDTseqs) == 0 and len(oldTDseqs) == 0:
        targets.append(
            stickydesign.enhist("DT", 5, energetics=energetics)[2]["emedian"]
        )
        targets.append(
            stickydesign.enhist("TD", 5, energetics=energetics)[2]["emedian"]
        )
    if oldDTarray:
        targets.append(energetics.matching_uniform(oldDTarray))
    if oldTDarray:
        targets.append(energetics.matching_uniform(oldTDarray))
    targetint = np.average(targets)

    if any(not seq.is_null(end.fseq) for end in newTD):
        TDtemplates = [end.fseq for end in newTD]
    else:
        TDtemplates = []
    if any(not seq.is_null(end.fseq) for end in newDT):
        DTtemplates = [end.fseq for end in newDT]
    else:
        DTtemplates = []

    if method == "default":
        if TDtemplates or DTtemplates:
            raise NotImplementedError
        # Create new sequences.
        newTDseqs = stickydesign.easyends(
            "TD",
            5,
            number=len(newTD),
            energetics=energetics,
            interaction=targetint,
            **sdopts,
        ).tolist()

        newDTseqs = stickydesign.easyends(
            "DT",
            5,
            number=len(newDT),
            energetics=energetics,
            interaction=targetint,
            **sdopts,
        ).tolist()

    elif method == "multimodel":
        SELOGGER.info(
            "starting multimodel sticky end generation "
            + "of TD ends for {} DT and {} TD ends, {} trials.".format(
                len(newDT), len(newTD), trials
            )
        )

        newTDseqs = []
        pl = util.ProgressLogger(SELOGGER, trials * 2)
        presetavail = None
        for i in range(0, trials):
            endchooserTD = multimodel.endchooser(
                all_energetics, templates=TDtemplates, devmethod=devmethod, **ecpars
            )

            e, presetavail = stickydesign.easyends(
                "TD",
                5,
                number=len(newTD),
                oldends=oldTDseqs,
                energetics=energetics,
                interaction=targetint,
                echoose=endchooserTD,
                _presetavail=presetavail,
                **sdopts,
            )
            newTDseqs.append(e)
            pl.update(i)

        if oldTDarray:
            tvals = [
                [e.matching_uniform(oldTDarray[0:1]) for e in all_energetics]
                * len(newTDseqs)
            ] * len(newTDseqs)
            SELOGGER.debug(tvals[0])
        else:
            tvals = [
                [e.matching_uniform(x[0:1]) for e in all_energetics] for x in newTDseqs
            ]

        endchoosersDT = [
            multimodel.endchooser(
                all_energetics,
                target_vals=tval,
                templates=DTtemplates,
                devmethod=devmethod,
                **ecpars,
            )
            for tval in tvals
        ]

        SELOGGER.info("generating corresponding DT ends")
        newDTseqs = []
        presetavail = None

        for i, echoose in enumerate(endchoosersDT):
            e, presetavail = stickydesign.easyends(
                "DT",
                5,
                number=len(newDT),
                oldends=oldDTseqs,
                energetics=energetics,
                interaction=targetint,
                echoose=echoose,
                _presetavail=presetavail,
                **sdopts,
            )
            newDTseqs.append(e)

            pl.update(i + trials)

        arr = [
            [
                stickydesign.endarray(oldTDseqs + x.tolist(), "TD"),
                stickydesign.endarray(oldDTseqs + y.tolist(), "DT"),
            ]
            for x, y in zip(newTDseqs, newDTseqs)
        ]

        scores = [
            multimodel.deviation_score(list(e), all_energetics, devmethod=devmethod)
            for e in arr
        ]

        sort = np.argsort(scores)

        newTDseqs = newTDseqs[sort[0]].tolist()[len(oldTDseqs) :]
        newDTseqs = newDTseqs[sort[0]].tolist()[len(oldDTseqs) :]
        info["score"] = float(scores[sort[0]])
        info["maxscore"] = float(scores[sort[-1]])
        info["meanscore"] = float(np.mean(scores))

    # FIXME: move to stickydesign
    assert len(newTDseqs) == len(newTD)
    assert len(newDTseqs) == len(newDT)

    # Shuffle the lists of end sequences, to ensure that they're
    # random order, and that ends used earlier in the set are not
    # always better than those used later. But only shuffle if
    # there were no templates:
    if not TDtemplates:
        shuffle(newTDseqs)
    if not DTtemplates:
        shuffle(newDTseqs)

    # Make sure things are consistent if there are templates:
    if TDtemplates:
        for t, s in zip(TDtemplates, newTDseqs):
            seq.merge(t, s)
    if DTtemplates:
        for t, s in zip(DTtemplates, newDTseqs):
            seq.merge(t, s)

    for end, s in zip(newDT, newDTseqs):
        ends[end.ident()].fseq = s
    for end, s in zip(newTD, newTDseqs):
        ends[end.ident()].fseq = s

    # Ensure that the old and new sets have consistent end definitions,
    # and that the tile definitions still fit.
    self.glues.update(ends)
    # self.tiles.glues_from_tiles().update(ends)

    newendnames = [e.name for e in newTD] + [e.name for e in newDT]
    info["newends"] = newendnames

    # Apply new sequences to tile system.
    self.ends = ends
    # if "info" not in self.keys():
    #    self["info"] = {}
    # if "end_design" not in self["info"].keys():
    #    self["info"]["end_design"] = []
    # if isinstance("end_design", dict):  # convert old
    #    self["info"]["end_design"] = [self["info"]["end_design"]]
    # self["info"]["end_design"].append(info)

    return self, newendnames


def dx_plot_se_hists(
    self: TileSet,
    all_energetics: Sequence[stickydesign.Energetics] | None = None,
    energetics_names: Sequence[str] | None = None,
    title: str | None = None,
    **kwargs: Any,
) -> Any:
    """Plot histograms of sticky end energies, using stickydesign.plots.hist_multi.

    Parameters
    ----------

    all_energetics : list of Energetics
        A list of energetics to use.  Defaults to DEFAULT_MULTIMODEL_ENERGETICS.

    energetics_names : list of str
        Names for energetics in all_energetics.  Defaults to DEFAULT_MM_ENERGETICS_NAMES.

    title : str
        Title for the plot.

    **kwargs
        kwargs passed to stickydesign.plots.hist_multi.

    """
    if all_energetics is None:
        all_energetics = DEFAULT_MULTIMODEL_ENERGETICS

    if energetics_names is None:
        energetics_names = DEFAULT_MM_ENERGETICS_NAMES

    ends = self.glues

    if title is None:
        # FIXME
        title = "Title"

    td = stickydesign.endarray(
        [x.fseq for x in ends if isinstance(x, DXGlue) and x.etype == "TD"], "TD"
    )

    dt = stickydesign.endarray(
        [x.fseq for x in ends if isinstance(x, DXGlue) and x.etype == "DT"], "DT"
    )
    import stickydesign.plots as sdplots

    return sdplots.hist_multi(
        [td, dt], all_energetics, energetics_names, title, **kwargs
    )


def dx_plot_se_lv(
    self: TileSet,
    all_energetics: Sequence[stickydesign.Energetics] | None = None,
    energetics_names: Sequence[str] | None = None,
    pltcmd: Optional[Callable[[Any], Any]] = None,
    title: str | None = None,
    **kwargs,
):
    """
    Uses an LV plot to show sticky end energetics.
    """

    if all_energetics is None:
        all_energetics = DEFAULT_MULTIMODEL_ENERGETICS

    if energetics_names is None:
        energetics_names = DEFAULT_MM_ENERGETICS_NAMES
    import stickydesign.plots as sdplots

    m, s = sdplots._multi_data_pandas(
        self.glues.to_endarrays(), all_energetics, energetics_names
    )

    import matplotlib.pyplot as plt
    import seaborn as sns

    if pltcmd is None:
        pltcmd = sns.lvplot

    pltcmd(data=m, **kwargs)  # type: ignore
    pltcmd(data=s, marker="x", **kwargs)  # type: ignore
    if title:
        plt.title(title)
    plt.ylabel("Energy (kcal/mol)")


def dx_plot_adjacent_regions(self: TileSet, energetics=None):
    """
    Plots the strength of double-stranded regions in DX tiles adjacent
    to sticky ends.

    Parameters
    ----------

    energetics : stickydesign.Energetics
        The energetics to use.  Defaults to DEFAULT_REGION_ENERGETICS.
    """

    if energetics is None:
        energetics = DEFAULT_REGION_ENERGETICS

    regions = [t.structure._side_bound_regions(t) for t in self.tiles]
    regions = [[x.lower() for x in y] for y in regions]
    allregions: list[str] = sum(regions, [])
    count: list[list[Counter]] = [[Counter(x) for x in y] for y in regions]
    gc_count = [[x["g"] + x["c"] for x in c] for c in count]
    gc_counts: list[int] = sum(gc_count, [])

    ens = energetics.matching_uniform(stickydesign.endarray(allregions, "DT"))
    from matplotlib import pylab

    pylab.figure(figsize=(10, 4))
    pylab.subplot(121)
    pylab.hist(gc_counts, bins=np.arange(min(gc_counts) - 0.5, max(gc_counts) + 0.5))
    pylab.title("G/C pairs in arms")
    pylab.ylabel("# of 8 nt arms")
    pylab.xlabel("# of G/C pairs")
    pylab.subplot(122)
    pylab.hist(ens)
    pylab.title("ΔG, T=33, no coaxparams/danglecorr")
    pylab.ylabel("# of 8 nt regions")
    pylab.xlabel("stickydesign ΔG")
    pylab.suptitle("8 nt end-adjacent region strengths")


def dx_plot_side_strands(self: TileSet, energetics=None):
    """
    Plots the binding strength of short strands in DX tiles.

    Parameters
    ----------

    energetics : stickydesign.Energetics
        The energetics to use.  Defaults to DEFAULT_REGION_ENERGETICS.
    """

    if energetics is None:
        energetics = DEFAULT_REGION_ENERGETICS

    regions = [t.structure._short_bound_full(t) for t in self.tiles]
    regions = [[x.lower() for x in y] for y in regions]
    allregions: list[str] = sum(regions, [])
    count: list[list[Counter]] = [[Counter(x) for x in y] for y in regions]
    gc_count = [[x["g"] + x["c"] for x in c] for c in count]
    gc_counts: list[int] = sum(gc_count, [])

    ens = energetics.matching_uniform(stickydesign.endarray(allregions, "DT"))
    from matplotlib import pylab

    pylab.figure(figsize=(10, 4))
    pylab.subplot(121)
    pylab.hist(gc_counts, bins=np.arange(min(gc_counts) - 0.5, max(gc_counts) + 0.5))
    pylab.title("G/C pairs in arms")
    pylab.ylabel("# of 8 nt arms")
    pylab.xlabel("# of G/C pairs")
    pylab.subplot(122)
    pylab.hist(ens)
    pylab.title("ΔG, T=33, no coaxparams/danglecorr")
    pylab.ylabel("# of 16 nt regions")
    pylab.xlabel("stickydesign ΔG")
    pylab.suptitle("16 nt arm region strengths")
