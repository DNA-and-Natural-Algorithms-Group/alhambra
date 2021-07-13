from __future__ import annotations
import pkg_resources

from stickydesign.energetics_daoe import EnergeticsDAOE
from alhambra import seq
import collections
from dataclasses import dataclass
import datetime
import warnings
import logging
import numpy as np

import stickydesign as sd

import stickydesign.multimodel as multimodel

from random import shuffle

import fastreduce


from .util import (
    DEFAULT_SD2_MULTIMODEL_ENERGETICS,
    DEFAULT_MM_ENERGETICS_NAMES,
    DEFAULT_REGION_ENERGETICS,
    DEFAULT_MULTIMODEL_ENERGETICS,
    DEFAULT_ENERGETICS,
)
from alhambra import util

SELOGGER = logging.getLogger(__name__)

from numpy import isin
from alhambra.classes import Serializable
from typing import Any, Iterable, Literal, Optional, Type, TypeVar, cast
from .newtile import Tile, tile_factory, TileList, TileSupportingScadnano
from .glue import DXGlue, Glue, GlueList
from .newseed import Seed, seed_factory
import xgrow.tileset as xgt
import xgrow
import copy

try:
    import scadnano
except ImportError:
    pass

T = TypeVar("T")


@dataclass(init=False)
class TileSet(Serializable):
    tiles: TileList[Tile]
    glues: GlueList[Glue]
    seed: Optional[Seed]
    params: dict

    def run_xgrow(self, *args, **kwargs):
        xgrow_tileset = self.to_xgrow()
        return xgrow.run(xgrow_tileset, *args, **kwargs)

    def to_xgrow(self, self_complementary_glues: bool = True) -> xgt.TileSet:

        tiles = [t.to_xgrow(self_complementary_glues) for t in self.tiles]

        # FIXME
        # bonds = [g.to_xgrow(self_complementary_glues) for g in self.glues]
        bonds = []

        if self.seed:
            seed_tiles, seed_bonds, initstate = self.seed.to_xgrow(
                self_complementary_glues
            )
        else:
            seed_tiles = []
            seed_bonds = []
            initstate = None

        xgrow_tileset = xgt.TileSet(
            seed_tiles + tiles, seed_bonds + bonds, initstate=initstate
        )

        return xgrow_tileset

    def summary(self):
        """Returns a short summary line about the TileSet"""
        # self.check_consistent()
        info = {
            "ntiles": len(self.tiles),
            "nrt": len([x for x in self.tiles if not x.is_fake]),
            "nft": len([x for x in self.tiles if x.is_fake]),
            "nends": len(self.glues),
            "ntends": len(self.tiles.glues_from_tiles()),
            "tns": " ".join(x.name for x in self.tiles if x.name),
            "ens": " ".join(x.name for x in self.glues if x.name),
            # "name": " {}".format(self["info"]["name"])
            # if ("info" in self.keys() and "name" in self["info"].keys())
            # else "",
        }
        tun = sum(1 for x in self.tiles if x.name is None)
        if tun > 0:
            info["tns"] += " ({} unnamed)".format(tun)
        eun = sum(1 for x in self.glues if x.name is None)
        if eun > 0:
            info["ens"] += " ({} unnamed)".format(eun)
        if info["nft"] > 0:
            info["nft"] = " (+ {} fake)".format(info["nft"])
        else:
            info["nft"] = ""
        return "TileSet{name}: {nrt} tiles{nft}, {nends} ends, {ntends} ends in tiles.\nTiles: {tns}\nEnds:  {ens}".format(
            **info
        )

    def __str__(self):
        return self.summary()

    def to_xgrow_dict(self) -> dict:
        return self.to_xgrow().to_dict()

    @classmethod
    def from_scadnano(cls: Type[T], des: "scadnano.Design", ret_fails=False) -> T:
        ts = cls()
        tiles: TileList[TileSupportingScadnano] = TileList()
        ts.glues = GlueList()
        positions: dict[tuple[int, int], str] = {}

        for strand in des.strands:
            try:
                t, o = tile_factory.from_scadnano(strand, return_position=True)
            except ValueError:
                warnings.warn(f"Failed to import strand {strand.name}.")
            except NotImplementedError:
                warnings.warn(f"Failed to import strand {strand.name}.")
            else:
                positions[o] = t.ident()
                if t.name in tiles.data.keys():
                    if t != tiles[t.ident()]:
                        warnings.warn(f"Skipping unequal duplicate strand {t.name}.")
                else:
                    tiles.add(t)

        ts.tiles = tiles
        ts.params["scadnano_positions"] = positions
        return ts

    def to_scadnano(self) -> "scadnano.Design":
        max_helix = (
            max(helix for helix, offset in self.params["scadnano_positions"]) + 2
        )
        des = scadnano.Design(helices=[scadnano.Helix() for _ in range(0, max_helix)])

        for (helix, offset), tilename in self.params["scadnano_positions"].items():
            cast(TileSupportingScadnano, self.tiles[tilename]).to_scadnano(
                des, helix, offset
            )

        return des

    def to_dict(self) -> dict:
        d = {}

        allglues = self.glues | self.tiles.glues_from_tiles()
        refglues = set(allglues.data.keys())  # FIXME

        if self.tiles:
            d["tiles"] = [t.to_dict(refglues=refglues) for t in self.tiles.aslist()]
        if allglues:
            d["glues"] = [g for g in allglues.aslist()]
        if self.seed:
            d["seed"] = self.seed.to_dict()
        if self.params:
            d["params"] = self.params.copy()
        return d

    @classmethod
    def from_dict(cls: Type[T], d: dict) -> T:
        ts = cls()
        ts.tiles = TileList(Tile.from_dict(x) for x in d.get("tiles", []))
        ts.glues = GlueList(Glue.from_dict(x) for x in d.get("glues", []))
        if "seed" in d:
            ts.seed = seed_factory.from_dict(d["seed"])
        if "params" in d:
            ts.params = copy.deepcopy(d["params"])
        return ts

    def _serialize(self) -> Any:
        return self.to_dict()

    @classmethod
    def _deserialize(cls, input: Any) -> TileSet:
        return cls.from_dict(input)

    def __init__(
        self,
        tiles: Iterable[Tile] = tuple(),
        glues: Iterable[Glue] = tuple(),
        seed: Seed | None = None,
        params: dict | None = None,
    ) -> None:
        self.tiles = TileList(tiles)
        self.glues = GlueList(glues)
        self.seed = seed
        if params is not None:
            self.params = params
        else:
            self.params = dict()

    def stickydesign_create_dx_glue_sequences(
        self,
        method: Literal["default", "multimodel"] = "default",
        energetics=None,
        trials=100,
        devmethod="dev",
        sdopts={},
        ecpars={},
        listends=False,
    ):
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
        info = {}
        info["method"] = method
        info["time"] = datetime.datetime.now(tz=datetime.timezone.utc).isoformat()
        info["sd_version"] = sd.version.__version__

        if not energetics:
            if method == "multimodel":
                all_energetics = DEFAULT_MULTIMODEL_ENERGETICS
            else:
                energetics = DEFAULT_ENERGETICS
        if method == "multimodel" and not isinstance(energetics, collections.Iterable):
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
            oldDTarray = sd.endarray(oldDTseqs, "DT")
        else:
            oldDTarray = None
        oldTDseqs = [
            end.fseq for end in ends if end.etype == "TD" and seq.is_definite(end.fseq)
        ]
        if oldTDseqs:
            oldTDarray = sd.endarray(oldTDseqs, "TD")
        else:
            oldTDarray = None

        newTD = [
            end for end in ends if end.etype == "TD" and not seq.is_definite(end.fseq)
        ]
        newDT = [
            end for end in ends if end.etype == "DT" and not seq.is_definite(end.fseq)
        ]

        # Deal with energetics, considering potential old sequences.
        # FIXME: EXPLAIN WHAT THIS ABSTRUSE CODE DOES...
        # TODO: tests needs to test this
        targets = []
        if len(oldDTseqs) == 0 and len(oldTDseqs) == 0:
            targets.append(sd.enhist("DT", 5, energetics=energetics)[2]["emedian"])
            targets.append(sd.enhist("TD", 5, energetics=energetics)[2]["emedian"])
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
            newTDseqs = sd.easyends(
                "TD",
                5,
                number=len(newTD),
                energetics=energetics,
                interaction=targetint,
                **sdopts,
            ).tolist()

            newDTseqs = sd.easyends(
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

                e, presetavail = sd.easyends(
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
                    [e.matching_uniform(x[0:1]) for e in all_energetics]
                    for x in newTDseqs
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
                e, presetavail = sd.easyends(
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
                    sd.endarray(oldTDseqs + x.tolist(), "TD"),
                    sd.endarray(oldDTseqs + y.tolist(), "DT"),
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

    @property
    def allglues(self) -> GlueList:
        return self.tiles.glues_from_tiles() | self.glues

    def create_abstract_diagrams(self, filename=None, only_real=True, **options):
        """Write abstract tile diagrams in SVG for the TileSet.

        Parameters
        ----------

        filename : str
            a path or filename for the SVG output.
        only_real: bool
            If True, only create diagrams for tiles that are not fake (ie,
            tiles that actually exist in the system, rather than being provided
            by other tiles as a result of reduction/etc.) (default is True)
        *options
            (currently unused)
        """
        import drawSvg as draw

        d = draw.Drawing(200, 200)

        pos = 0
        for tile in self.tiles:
            if tile.is_fake:
                continue
            group = tile.abstract_diagram(**options)

            group.args["transform"] = "translate({},{})".format(
                (pos % 12) * 22, (pos // 12) * 22
            )
            pos += 2

            d.append(group)

        if filename is not None:
            d.saveSvg(filename)
        else:
            return d

    def create_layout_diagrams(
        self, xgrowarray: np.ndarray, filename=None, scale=1, **options
    ):
        """Create an SVG layout diagram from xgrow output.

        This currently uses the abstract diagram bases to create the
        layout diagrams.

        Parameters
        ----------

        xgrowarray : ndarray or dict
            Xgrow output.  This may be a numpy array of
            an xgrow state, obtained in some way, or may be the 'array'
            output of xgrow.run.

        filename : string
            File name / path of the output file.

        """
        import drawSvg as draw

        d = draw.Drawing(200, 200)

        svgtiles = {}

        if isinstance(xgrowarray, dict):
            if "tiles" in xgrowarray:
                xgrowarray = xgrowarray["tiles"]
            elif "array" in xgrowarray:
                xgrowarray = xgrowarray["array"]["tiles"]

        for tile in self.tiles:
            svgtiles[tile.name] = tile.abstract_diagram(**options)

        tilelist = self.to_xgrow().tiles
        minxi = 10000
        minyi = 10000
        maxxi = 0
        maxyi = 0

        for yi in range(0, xgrowarray.shape[0]):
            for xi in range(0, xgrowarray.shape[1]):
                ti = int(xgrowarray[yi, xi] - 1)
                if ti >= len(tilelist):  # the extra side of a duple
                    continue
                if ti == -1:
                    continue
                tn = tilelist[ti].name
                if not (tn in svgtiles.keys()):
                    continue
                minxi = min(minxi, xi)
                minyi = min(minyi, yi)
                maxxi = max(maxxi, xi)
                maxyi = max(maxyi, yi)
                d.append(draw.Use(svgtiles[tn], xi * 10, yi * 10))

        d.viewBox = (minxi * 10, minyi * 10, (maxxi - minxi) * 10, (maxyi - minyi) * 10)

        if filename:
            d.saveSvg(filename)
        else:
            return d

    def reduce_tiles(
        self,
        preserve=["s22", "ld"],
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
        reduced: single TileSet or equiv, or list
            The reduced system/systems
        """
        return fastreduce.reduce_tiles(
            self, preserve, tries, threads, returntype, best, key, initequiv
        )

    def reduce_ends(
        self,
        preserve=["s22", "ld"],
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
        reduced: single TileSet or equiv, or list
            The reduced system/systems
        """
        return fastreduce.reduce_ends(
            self, preserve, tries, threads, returntype, best, key, initequiv
        )

    def latticedefects(self, direction="e", depth=2, pp=True, rotate=False):
        """
        Calculate and show possible small lattice defect configurations.
        """
        from . import latticedefect

        return latticedefect.latticedefects(
            self, direction=direction, depth=depth, pp=pp, rotate=rotate
        )

    def apply_equiv(self, equiv):
        """
        Apply an equivalence array (from, eg, `TileSet.reduce_ends` or `TileSet.reduce_tiles`).

        Parameters
        ----------
        equiv : ndarray
            An equivalence array, *for this tileset*, generated by reduction functions.

        Returns
        -------
        TileSet
            A tileset with the equivalence array, and thus the reduction, applied.
        """
        return fastreduce._FastTileSet(self).applyequiv(self, equiv)

    def dx_plot_se_hists(
        tileset, all_energetics=None, energetics_names=None, title=None, **kwargs
    ):
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

        if "ends" in tileset.keys():
            ends = tileset.glues
        else:
            ends = tileset

        if title is None:
            # FIXME
            title = "Title"

        td = sd.endarray([x["fseq"] for x in ends if x["type"] == "TD"], "TD")

        dt = sd.endarray([x["fseq"] for x in ends if x["type"] == "DT"], "DT")
        import stickydesign.plots as sdplots

        return sdplots.hist_multi(
            [td, dt], all_energetics, energetics_names, title, **kwargs
        )

    def dx_plot_se_lv(
        self,
        all_energetics=None,
        energetics_names=None,
        pltcmd=None,
        title=None,
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

        import seaborn as sns
        import matplotlib.pyplot as plt

        if pltcmd is None:
            pltcmd = sns.lvplot

        pltcmd(data=m, **kwargs)
        pltcmd(data=s, marker="x", **kwargs)
        if title:
            plt.title(title)
        plt.ylabel("Energy (kcal/mol)")

    def dx_plot_adjacent_regions(tileset, energetics=None):
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

        regions = [t.structure._side_bound_regions(t) for t in tileset.tiles]
        regions = [[x.lower() for x in y] for y in regions]
        allregions = sum(regions, [])
        count = [[Counter(x) for x in y] for y in regions]
        gc_count = [[x["g"] + x["c"] for x in c] for c in count]
        gc_counts = sum(gc_count, [])

        ens = energetics.matching_uniform(sd.endarray(allregions, "DT"))
        from matplotlib import pylab

        pylab.figure(figsize=(10, 4))
        pylab.subplot(121)
        pylab.hist(
            gc_counts, bins=np.arange(min(gc_counts) - 0.5, max(gc_counts) + 0.5)
        )
        pylab.title("G/C pairs in arms")
        pylab.ylabel("# of 8 nt arms")
        pylab.xlabel("# of G/C pairs")
        pylab.subplot(122)
        pylab.hist(ens)
        pylab.title("ΔG, T=33, no coaxparams/danglecorr")
        pylab.ylabel("# of 8 nt regions")
        pylab.xlabel("stickydesign ΔG")
        pylab.suptitle("8 nt end-adjacent region strengths")

    def dx_plot_side_strands(tileset, energetics=None):
        """
        Plots the binding strength of short strands in DX tiles.

        Parameters
        ----------

        energetics : stickydesign.Energetics
            The energetics to use.  Defaults to DEFAULT_REGION_ENERGETICS.
        """

        if energetics is None:
            energetics = DEFAULT_REGION_ENERGETICS

        regions = [t.structure._short_bound_full(t) for t in tileset.tiles]
        regions = [[x.lower() for x in y] for y in regions]
        allregions = sum(regions, [])
        count = [[Counter(x) for x in y] for y in regions]
        gc_count = [[x["g"] + x["c"] for x in c] for c in count]
        gc_counts = sum(gc_count, [])

        ens = energetics.matching_uniform(sd.endarray(allregions, "DT"))
        from matplotlib import pylab

        pylab.figure(figsize=(10, 4))
        pylab.subplot(121)
        pylab.hist(
            gc_counts, bins=np.arange(min(gc_counts) - 0.5, max(gc_counts) + 0.5)
        )
        pylab.title("G/C pairs in arms")
        pylab.ylabel("# of 8 nt arms")
        pylab.xlabel("# of G/C pairs")
        pylab.subplot(122)
        pylab.hist(ens)
        pylab.title("ΔG, T=33, no coaxparams/danglecorr")
        pylab.ylabel("# of 16 nt regions")
        pylab.xlabel("stickydesign ΔG")
        pylab.suptitle("16 nt arm region strengths")

    def check_consistent(self):
        """Check the TileSet consistency.

        Check a number of properties of the TileSet for consistency.
        In particular:

           * Each tile must pass Tile.check_consistent()
           * TileSet.ends and TileSet.tiles.endlist() must not contain conflicting
             ends or end sequences.
           * If there is a seed:
               * It must be of an understood type (it must be in seeds.seedtypes)
               * All adapter locations must be valid.
               * The seed must pass its check_consistent and check_sequence.
        """
        # * END LIST The end list itself must be consistent.
        # ** Each end must be of understood type
        # ** Each end must have a valid sequence or no sequence
        # ** There must be no more than one instance of each name
        # ** WARN if there are ends with no namecounts
        # * TILE LIST
        # ** each tile must be of understood type (must parse)
        # ** ends in the tile list must be consistent (must merge)
        # ** there must be no more than one tile with each name
        # self.tiles.check_consistent()
        endsfromtiles = self.tiles.glues_from_tiles()

        # ** WARN if any end that appears does not have a complement used or vice versa
        # ** WARN if there are tiles with no name
        # * TILE + END
        # ** The tile and end lists must merge validly
        # (checks sequences, adjacents, types, complements)
        self.glues | endsfromtiles

        # ** WARN if tilelist has end references not in ends
        # ** WARN if merge is not equal to the endlist
        # ** WARN if endlist has ends not used in tilelist
        # * ADAPTERS / SEEDS
        # SEED stuff was here

    def copy(self):
        """Return a full (deep) copy of the TileSet"""
        return copy.deepcopy(self)
