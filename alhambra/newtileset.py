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

    def create_abstract_diagrams(self, filename, only_real=True, *options):
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
        from lxml import etree
        import os
        import pkg_resources

        base = etree.parse(
            pkg_resources.resource_stream(
                __name__, os.path.join("seqdiagrambases", "blank.svg")
            )
        )
        baseroot = base.getroot()

        pos = 0
        for tile in self.tiles:
            if tile.is_fake:
                continue
            group, n = tile.abstract_diagram(self)

            group.attrib["transform"] = "translate({},{})".format(
                (pos % 12) * 22, (pos // 12) * 22
            )
            pos += n

            baseroot.append(group)

        base.write(filename)

    def create_layout_diagrams(
        self, xgrowarray: np.ndarray, filename, scale=1, *options
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
        from lxml import etree

        base = etree.parse(
            pkg_resources.resource_stream(
                __name__, os.path.join("seqdiagrambases", "blank.svg")
            )
        )
        baseroot = base.getroot()

        svgtiles = {}

        if isinstance(xgrowarray, dict):
            if "tiles" in xgrowarray:
                xgrowarray = xgrowarray["tiles"]
            elif "array" in xgrowarray:
                xgrowarray = xgrowarray["array"]["tiles"]

        for tile in self.tiles:
            group, n = tile.abstract_diagram(self)
            svgtiles[tile.name] = group

        tilelist = self.to_xgrow_dict()["tiles"]
        tilen = [None] + [x["name"] for x in tilelist]
        firstxi = 10000
        firstyi = 10000
        import copy

        for yi in range(0, xgrowarray.shape[0]):
            for xi in range(0, xgrowarray.shape[1]):
                tn = tilen[xgrowarray[yi, xi]]
                if tn and tn[-5:] == "_left":
                    tn = tn[:-5]
                if tn and tn[-7:] == "_bottom":
                    tn = tn[:-7]
                if not (tn in svgtiles.keys()):
                    continue
                if xi < firstxi:
                    firstxi = xi
                if yi < firstyi:
                    firstyi = yi
                st = copy.deepcopy(svgtiles[tn])
                st.attrib["transform"] = "translate({},{})".format(xi * 10, yi * 10)
                baseroot.append(st)

        base.write(filename)
