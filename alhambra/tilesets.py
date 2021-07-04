from typing import Sequence
import ruamel.yaml as yaml
from ruamel.yaml.comments import CommentedMap
from ruamel.yaml.representer import RoundTripRepresenter
import warnings

import copy
import re

import pkg_resources
import os

from .tiles import TileList
from .ends import EndList

from . import tilestructures
from . import seeds
from . import util
from . import seq
from . import sensitivitynew as sensitivity

# from . import fastreduce # FIXME

# Need to disable this for now
from peppercompiler import compiler as compiler
from peppercompiler.design import spurious_design as spurious_design
from peppercompiler import finish as finish
from peppercompiler.DNA_classes import wc

import numpy as np
import stickydesign as sd

import stickydesign.multimodel as multimodel
from collections import Counter

import collections
from random import shuffle
from datetime import datetime, timezone

import logging

import stickydesign.stickydesign2 as sd2

from .util import (
    DEFAULT_SD2_MULTIMODEL_ENERGETICS,
    DEFAULT_MM_ENERGETICS_NAMES,
    DEFAULT_REGION_ENERGETICS,
    DEFAULT_MULTIMODEL_ENERGETICS,
    DEFAULT_ENERGETICS,
)

SELOGGER = logging.getLogger(__name__)


class TileSet(CommentedMap):
    """A class representing a single DNA tile assembly system.

    This class parses and manipulates Alhambra tilesets.  Note that,
    unless you have already imported or made a dict, you will likely
    want to construct this with :meth:`TileSet.from_file`.

    Parameters
    ----------
    val : dict
        a dict in the Alhambra tileset format.
    """

    def __init__(self, val={}):
        CommentedMap.__init__(self, val)

        self["ends"] = EndList(self.get("ends", []))
        self["tiles"] = TileList(self.get("tiles", []))

    @classmethod
    def from_file(cls, name_or_stream, *args, **kwargs):
        """
        Class method to create a TileSet from a file or stream.

        Loads a YAML stream/file into TileSet.  This uses ruamel.yaml and tries to preserve
        comments, which is supported through `TileSet.to_file`.

        Parameters
        ----------
        name_or_stream : file-like or path-like object
            The path of the file, or a file object, to load.

        Returns
        -------
        TileSet
            The loaded TileSet.
        """
        # Assume a stream:
        if getattr(name_or_stream, "read", None) is None:
            return cls(yaml.round_trip_load(open(name_or_stream, "r"), *args, **kwargs))
        else:
            return cls(yaml.round_trip_load(name_or_stream, *args, **kwargs))

    def to_file(self, name_or_stream):
        """
        Write a TileSet to a file or stream.

        Parameters
        ----------
        name_or_stream : file-like or path-like object
            The path of the file, or a file object, to write to.

        """

        if getattr(name_or_stream, "write", None) is None:
            return yaml.round_trip_dump(self, open(name_or_stream, "w"))
        else:
            return yaml.round_trip_dump(self, name_or_stream)

    @property
    def tiles(self) -> TileList:
        """Property returning the `TileList` of tiles in the TileSet

        Returns
        -------

        TileList
            Tiles listed in the tileset.
        """
        return self["tiles"]

    @tiles.setter
    def tiles(self, val: TileList):
        if isinstance(val, TileList):
            self["tiles"] = val
        else:
            self["tiles"] = TileList(val)

    @property
    def allends(self) -> EndList:
        """All ends in the system, both from ends and tiles.

        Returns
        -------

        EndList

        See Also
        --------

        TileSet.ends
        TileSet.tiles.endlist"""
        return self.ends.merge(self.tiles.glues_from_tiles())

    def ends():
        doc = """The EndList of specified ends in the TileSet (not including
                 ends that are only in tiles.

                 Returns
                 -------

                 EndList

                 See Also
                 --------

                 TileSet.allends

                 TileSet.tiles.endlist"""

        def fget(self):
            return self["ends"]

        def fset(self, value):
            self["ends"] = value

        def fdel(self):
            del self["ends"]

        return locals()

    ends = property(**ends())

    @property
    def seed(self):
        """Information on the tileset's seed"""
        return self.get("seed", None)

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

    def create_sequence_diagrams(tileset, filename, *options):
        """Write sequence tile diagrams in SVG for the TileSet.

        Parameters
        ----------

        filename : str
            a path or filename for the SVG output.
        *options
            (currently unused)
        """

        from lxml import etree
        import pkg_resources
        import os.path

        base = etree.parse(
            pkg_resources.resource_stream(
                __name__, os.path.join("seqdiagrambases", "blank.svg")
            )
        )
        baseroot = base.getroot()
        pos = 150
        for tile in tileset.tiles:
            if tile.is_fake:
                continue
            group = tile.sequence_diagram()

            group.attrib["transform"] = "translate(0,{})".format(pos)
            pos += 150
            baseroot.append(group)

        base.write(filename)

    def create_adapter_sequences(tileset):
        """Generate adapter sequences for the TileSet, returning a new
        TileSet with the generated sequences included.

        Returns
        -------

        TileSet
           TileSet with generated adapter sequences.
        """
        seedclass = seeds.seedtypes[tileset["seed"]["type"]]
        if seedclass.needspepper:
            warnings.warn(
                "This set must have adapter sequences created during\
     regular sequence design. You can ignore this if you just created sequences."
            )
            return tileset
        return seedclass.create_adapter_sequences(tileset)

    def create_layout_diagrams(tileset, xgrowarray, filename, scale=1, *options):
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

        for tile in tileset.tiles:
            group, n = tile.abstract_diagram(tileset)
            svgtiles[tile.name] = group

        tilelist = tileset.generate_xgrow_dict(perfect=True)["tiles"]
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

    @property
    def strand_order_list(self):
        """Return a list of orderable strand sequences for all tile strands.

        These sequences are generated by the Tile.orderableseqs attribute, which
        should include useful strand names and modified bases for labels.  The
        sequences may thus not be valid ACGT-only sequences (eg, they may contain
        /iBiodT/).

        Returns
        -------

        list
            a list of tuple (strand_name, strand_sequence) for the system.
        """
        return [y for x in self.tiles if not x.is_fake for y in x.orderableseqs]

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
        self.tiles.check_consistent()
        endsfromtiles = self.tiles.glues_from_tiles()

        # ** WARN if any end that appears does not have a complement used or vice versa
        # ** WARN if there are tiles with no name
        # * TILE + END
        # ** The tile and end lists must merge validly
        # (checks sequences, adjacents, types, complements)
        self.ends.merge(endsfromtiles)

        # ** WARN if tilelist has end references not in ends
        # ** WARN if merge is not equal to the endlist
        # ** WARN if endlist has ends not used in tilelist
        # * ADAPTERS / SEEDS
        if "seed" in self.keys():
            # ** seeds must be of understood type
            assert self["seed"]["type"] in seeds.seedtypes.keys()
            # ** adapter locations must be valid
            sclass = seeds.seedtypes[self["seed"]["type"]]

            sclass.check_consistent(self)

            sclass.check_sequence(self)
            # ** each adapter must have no sequence or a consistent sequence
            # *** the RH strand must match the associated tile
            # *** the ends in the sequence must match the ends in the endlist
            # *** the LH sequence must be validly binding to both RH and
            #     origami
            # ** each adapter must have valid definition, which means for us:
            # *** if both tile mimic and ends are specified, they must match

    def summary(self):
        """Returns a short summary line about the TileSet"""
        self.check_consistent()
        info = {
            "ntiles": len(self.tiles),
            "nrt": len([x for x in self.tiles if not x.is_fake]),
            "nft": len([x for x in self.tiles if x.is_fake]),
            "nends": len(self.ends),
            "ntends": len(self.tiles.glues_from_tiles()),
            "tns": " ".join(x.name for x in self.tiles if x.name),
            "ens": " ".join(x.name for x in self.ends if x.name),
            "name": " {}".format(self["info"]["name"])
            if ("info" in self.keys() and "name" in self["info"].keys())
            else "",
        }
        tun = sum(1 for x in self.tiles if "name" not in x.keys())
        if tun > 0:
            info["tns"] += " ({} unnamed)".format(tun)
        eun = sum(1 for x in self.ends if "name" not in x.keys())
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

    def copy(self):
        """Return a full (deep) copy of the TileSet"""
        return copy.deepcopy(self)

    def dump(self, stream):
        """Dump the tileset into a stream in YAML format.

        Parameters
        ----------

        stream : stream
            the stream / file object to dump to.
        """

        return yaml.round_trip_dump(self, stream)

    def design_set(
        tileset,
        name="tsd_temp",
        includes=[pkg_resources.resource_filename(__name__, "peppercomps-j1")],
        energetics=None,
        stickyopts={},
        reorderopts={},
        coreopts={},
        keeptemp=False,
    ):
        """Helper function to design sets from scratch, calling the numerous parts of
        tilesetdesigner. You may want to use the tilesetdesigner shell script
        instead.

        As with other functions in tilesetdesigner, this should not clobber inputs.

        Parameters
        ----------

        name : str, optional
            Base name for temporary files (default tsd_temp).

        Returns
        -------
        TileSet
            A copy of the TileSet with added sequences.

        """
        if energetics is None:
            energetics = DEFAULT_ENERGETICS

        if hasattr(tileset, "read"):
            tileset = TileSet.load(tileset)
        else:
            tileset = TileSet(tileset)

        tileset.check_consistent()
        tileset_with_ends_randomorder, new_ends = stickyends.create_sequences(
            tileset, energetics=energetics, **stickyopts
        )
        tileset_with_ends_ordered = stickyends.reorder(
            tileset_with_ends_randomorder,
            newends=new_ends,
            energetics=energetics,
            **reorderopts
        )
        tileset_with_strands = create_strand_sequences(
            tileset_with_ends_ordered, name, includes=includes, **coreopts
        )

        if "guards" in tileset_with_strands.keys():
            tileset_with_strands = create_guard_strand_sequences(tileset_with_strands)

        # FIXME: this is temporary, until we have a better way of deciding.
        if "createseqs" in tileset_with_strands["seed"].keys():
            tileset_with_strands = create_adapter_sequences(tileset_with_strands)

        if not keeptemp:
            os.remove(name + ".fix")
            os.remove(name + ".mfe")
            os.remove(name + ".pil")
            os.remove(name + ".save")
            os.remove(name + ".seqs")
            os.remove(name + ".sys")

        tileset_with_strands.check_consistent()
        return tileset_with_strands

    def create_end_sequences(
        tileset,
        method="default",
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
        info["time"] = datetime.now(tz=timezone.utc).isoformat()
        info["sd_version"] = sd.version.__version__

        if not energetics:
            if method == "multimodel":
                energetics = DEFAULT_MULTIMODEL_ENERGETICS
            else:
                energetics = DEFAULT_ENERGETICS
        if method == "multimodel" and not isinstance(energetics, collections.Iterable):
            raise ValueError("Energetics must be an iterable for multimodel.")
        elif method == "multimodel":
            all_energetics = energetics
            energetics = all_energetics[0]
            info["energetics"] = [str(e) for e in all_energetics]
            info["trails"] = trials
        elif method == "default":
            info["energetics"] = str(energetics)

        # Steps for doing this:

        # Create a copy of the tileset.
        newtileset = tileset.copy()

        # Build a list of ends from the endlist in the tileset.  Do this
        # by creating a NamedList, then merging them into it.
        ends = EndList()

        if newtileset.ends:
            ends.update(newtileset.ends, fail_immediate=False)

        # This is the endlist from the tiles themselves.
        if newtileset.tiles:  # maybe you just want ends?
            # this checks for end/complement usage, and whether any
            # previously-describedends are unused
            # FIXME: implement
            # tilestructures.check_end_usage(newtileset.tiles, ends)

            endlist_from_tiles = newtileset.tiles.glues_from_tiles()

        ends.update(endlist_from_tiles)

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
        oldTDseqs = [
            end.fseq for end in ends if end.etype == "TD" and seq.is_definite(end.fseq)
        ]
        if oldTDseqs:
            oldTDarray = sd.endarray(oldTDseqs, "TD")

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
        if len(oldDTseqs) > 0:
            targets.append(energetics.matching_uniform(oldDTarray))
        if len(oldTDseqs) > 0:
            targets.append(energetics.matching_uniform(oldTDarray))
        targetint = np.average(targets)

        if any(not seq.is_null(end.fseq) for end in newTD):
            TDtemplates = [end.fseq for end in newTD]
        else:
            TDtemplates = None
        if any(not seq.is_null(end.fseq) for end in newDT):
            DTtemplates = [end.fseq for end in newDT]
        else:
            DTtemplates = None

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
                **sdopts
            ).tolist()

            newDTseqs = sd.easyends(
                "DT",
                5,
                number=len(newDT),
                energetics=energetics,
                interaction=targetint,
                **sdopts
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
                    **sdopts
                )
                newTDseqs.append(e)
                pl.update(i)

            if oldTDseqs:
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
                    **ecpars
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
                    **sdopts
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
            ends[end.name].fseq = s
        for end, s in zip(newTD, newTDseqs):
            ends[end.name].fseq = s

        ends.check_consistent()

        # Ensure that the old and new sets have consistent end definitions,
        # and that the tile definitions still fit.
        tileset.ends.merge(ends)
        newtileset.tiles.glues_from_tiles().merge(ends)

        newendnames = [e.name for e in newTD] + [e.name for e in newDT]
        info["newends"] = newendnames

        # Apply new sequences to tile system.
        newtileset.ends = ends
        if "info" not in newtileset.keys():
            newtileset["info"] = {}
        if "end_design" not in newtileset["info"].keys():
            newtileset["info"]["end_design"] = []
        if isinstance("end_design", dict):  # convert old
            newtileset["info"]["end_design"] = [newtileset["info"]["end_design"]]
        newtileset["info"]["end_design"].append(info)

        return newtileset, newendnames

    def _add_info(self, ititle, data):
        if "info" not in self.keys():
            self["info"] = {}
        if ititle not in self["info"].keys():
            self["info"][ititle] = []
        self["info"][ititle].append(data)

    def reorder_ends(
        tileset,
        newends=[],
        hightemp=0.1,
        lowtemp=1e-7,
        steps=45000,
        update=1000,
        energetics=None,
        fsopts={},
    ):
        """Given a tileset dictionary that includes sticky end sequences, reorder these
        to try to optimize error rates.
        """
        from . import endreorder
        from . import anneal

        if energetics is None:
            energetics = DEFAULT_ENERGETICS

        tset = tileset.copy()

        if "info" not in tset.keys():
            tset["info"] = {}

        reordersys = endreorder.EndSystemFseq(
            tset, newends, energetics=energetics, **fsopts
        )

        # FIXME: better parameter control here.
        annealer = anneal.Annealer(reordersys.score, reordersys.mutate)

        newstate = annealer.anneal(
            reordersys.initstate, hightemp, lowtemp, steps, update
        )

        # Now take that new state, and apply it to the new tileset.
        seqs = reordersys.slowseqs(newstate[0])
        for end in tset.ends:
            if end.etype in ["DT", "TD"]:
                eloc = reordersys.enlocs[end["name"]]
                end.fseq = seqs[eloc[1]].tolist()[eloc[0]]

        ri = {}

        ri["score"] = float(reordersys.score(newstate[0]))

        tset["info"]["reorder"] = ri

        # Ensure that only ends in newends moved: that all others remain mergeable:
        if newends:
            old_ends_from_new_set = EndList(
                end for end in tset.ends if end["name"] not in newends
            )
            tileset.ends.merge(old_ends_from_new_set)

        # Ensure system consistency
        tset.check_consistent()
        return tset

    def _create_end_sequences_sd2(
        tileset,
        method="default",
        energetics=None,
        trials=100,
        devmethod="dev",
        sdopts={},
        ecpars={},
        listends=False,
    ):
        """Create sticky end sequences for the TileSet, using stickydesign,
        and returning a new TileSet including the ends. USES STICKYDESIGN2.


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
        info["time"] = datetime.now(tz=timezone.utc).isoformat()
        info["sd_version"] = "stickydesign2-predev"

        if not energetics:
            if method == "multimodel":
                energetics = DEFAULT_SD2_MULTIMODEL_ENERGETICS
            else:
                energetics = DEFAULT_SD2_MULTIMODEL_ENERGETICS[0]
        if method == "multimodel" and not isinstance(energetics, collections.Iterable):
            raise ValueError("Energetics must be an iterable for multimodel.")
        elif method == "multimodel":
            all_energetics = energetics
            energetics = all_energetics[0]
            info["energetics"] = [str(e) for e in all_energetics]
            info["trails"] = trials
        elif method == "default":
            all_energetics = [energetics]
            info["energetics"] = str(energetics)

        # Steps for doing this:
        print(all_energetics)
        # Create a copy of the tileset.
        newtileset = tileset.copy()

        # Build a list of ends from the endlist in the tileset.  Do this
        # by creating a NamedList, then merging them into it.
        ends = EndList()

        if newtileset.ends:
            ends.update(newtileset.ends, fail_immediate=False)

        # This is the endlist from the tiles themselves.
        if newtileset.tiles:  # maybe you just want ends?
            # this checks for end/complement usage, and whether any
            # previously-describedends are unused
            # FIXME: implement
            # tilestructures.check_end_usage(newtileset.tiles, ends)

            endlist_from_tiles = newtileset.tiles.glues_from_tiles()

        ends.update(endlist_from_tiles)

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
            oldDTarray = sd2.EndPairArrayDT(oldDTseqs, "DT")
        else:
            oldDTarray = None
        oldTDseqs = [
            end.fseq for end in ends if end.etype == "TD" and seq.is_definite(end.fseq)
        ]
        if oldTDseqs:
            oldTDarray = sd2.EndPairArrayTD(oldTDseqs, "TD")
        else:
            oldTDarray = None

        newTD = [
            end for end in ends if end.etype == "TD" and not seq.is_definite(end.fseq)
        ]
        newDT = [
            end for end in ends if end.etype == "DT" and not seq.is_definite(end.fseq)
        ]

        # If there are old sequences, use those to get an interaction
        # target for the filter:
        interactions = []
        if len(oldDTseqs) > 0:
            interaction.append(np.mean(energetics.gse(oldDTarray)))
        if len(oldTDseqs) > 0:
            interactions.append(np.mean(energetics.gse(oldTDarray)))
        if len(oldDTseqs) == 0 and len(oldTDseqs) == 0:
            interactions.append(sd2.genpre("DT", 5, energetics)[1])
            interactions.append(sd2.genpre("TD", 5, energetics)[1])
        interaction = np.average(interactions)

        availendsDT, _, sfDT = sd2.genpre(
            "DT", 5, energetics, interaction=interaction, oldends=oldDTarray
        )
        availendsTD, _, sfTD = sd2.genpre(
            "TD", 5, energetics, interaction=interaction, oldends=oldTDarray
        )

        if any(not seq.is_null(end.fseq) for end in newTD):
            raise NotImplementedError
        if any(not seq.is_null(end.fseq) for end in newDT):
            raise NotImplementedError

        newTDseqs = []
        newDTseqs = []

        rejects = 0
        min_score = 1000.0
        SELOGGER.info(
            "starting sticky end generation "
            + "of TD ends for {} DT and {} TD ends, {} trials.".format(
                len(newDT), len(newTD), trials
            )
        )

        for i in range(0, trials):

            if method == "default":
                chooser = sd2.BasicChooser(energetics, interaction)

            elif method == "multimodel":
                chooser = sd2.MultiModelChooser(all_energetics, **ecpars)

            eTD = sd2.easyends(
                "TD",
                5,
                number=len(newTD),
                # oldends=oldTDseqs,
                energetics=energetics,
                interaction=interaction,
                chooser=chooser,
                preavail=availendsTD,
                seqfilter=sfTD,
                **sdopts
            )
            if len(eTD) < len(newTD):
                rejects += 1
                continue

            eDT = sd2.easyends(
                "DT",
                5,
                number=len(newDT),
                # oldends=oldDTseqs,
                energetics=energetics,
                interaction=interaction,
                chooser=chooser,
                preavail=availendsDT,
                seqfilter=sfDT,
                **sdopts
            )
            if len(eDT) < len(newDT):
                rejects += 1
                continue

            # FIXME: old ends
            score = sd2.deviation_score([eTD, eDT], all_energetics, devmethod=devmethod)
            if score < min_score:
                min_score = score
                minTD = eTD
                minDT = eDT
                SELOGGER.info(
                    "Found set: score {} trial {} rejects {}".format(score, i, rejects)
                )

        newTDseqs = minTD.strs
        newDTseqs = minDT.strs
        info["score"] = float(min_score)

        # Shuffle the lists of end sequences, to ensure that they're
        # random order, and that ends used earlier in the set are not
        # always better than those used later. But only shuffle if
        # there were no templates:
        shuffle(newTDseqs)
        shuffle(newDTseqs)

        # Make sure things are consistent if there are templates:
        for end, s in zip(newDT, newDTseqs):
            ends[end.name].fseq = s
        for end, s in zip(newTD, newTDseqs):
            ends[end.name].fseq = s

        ends.check_consistent()

        # Ensure that the old and new sets have consistent end definitions,
        # and that the tile definitions still fit.
        tileset.ends.merge(ends)
        newtileset.tiles.glues_from_tiles().merge(ends)

        newendnames = [e.name for e in newTD] + [e.name for e in newDT]
        info["newends"] = newendnames

        # Apply new sequences to tile system.
        newtileset.ends = ends
        if "info" not in newtileset.keys():
            newtileset["info"] = {}
        if "end_design" not in newtileset["info"].keys():
            newtileset["info"]["end_design"] = []
        if isinstance("end_design", dict):  # convert old
            newtileset["info"]["end_design"] = [newtileset["info"]["end_design"]]
        newtileset["info"]["end_design"].append(info)

        return newtileset, newendnames

    def create_strand_sequences(
        tileset,
        basename="alhambratemp",
        includes=[pkg_resources.resource_filename(__name__, "peppercomps-j1")],
        spurious_pars="verboten_weak=1.5",
        *options
    ):
        """Given a tileset dictionary with sticky ends sequences, create core sequences
        for tiles, using Pepper.

        Parameters
        ----------

        basename : str, optional
            The base name to use in Pepper (FIXME: why is this needed?).

        includes : list of str paths, optional
            The include paths for Pepper.

        spurious_pars : str
            Command line options for SpuriousDesign

        Returns
        -------

        TileSet
            A TileSet with generated strand sequences for all tiles.
        """

        if not basename:
            import uuid

            basename = str(uuid.uuid4()).replace("-", "")
        # FIXME: is this a valid basename?

        newtileset = copy.deepcopy(tileset)

        newtileset._create_pepper_input_files(basename)

        compiler.compiler(
            basename,
            [],
            basename + ".pil",
            basename + ".save",
            fixed_file=basename + ".fix",
            includes=includes,
            synth=True,
        )

        spurious_design.design(
            basename,
            infilename=basename + ".pil",
            outfilename=basename + ".mfe",
            verbose=True,
            struct_orient=True,
            tempname=basename + "-temp",
            extra_pars=spurious_pars,
            findmfe=False,
            cleanup=False,
        )

        if "info" not in newtileset.keys():
            newtileset["info"] = {}
        if "core" not in newtileset["info"].keys():
            newtileset["info"]["core"] = []
        elif isinstance(newtileset["info"]["core"], dict):
            newtileset["info"]["core"] = [newtileset["info"]["core"]]

        with open(basename + "-temp.sp") as f:
            a = f.read()
            cdi = {}
            cdi["basename"] = basename
            cdi["score_verboten"] = float(
                re.findall(r"score_verboten\s+score\s+=\s+([+-]?[\d.,]+)", a)[1]
            )
            cdi["score_spurious"] = float(
                re.findall(r"score_spurious\s+score\s+=\s+([+-]?[\d.,]+)", a)[1]
            )
            cdi["score_bonds"] = float(
                re.findall(r"score_bonds\s+score\s+=\s+([+-]?[\d.,]+)", a)[1]
            )
            cdi["score"] = float(
                re.findall(r"weighted score\s+=\s+([+-]?[\d.,]+)", a)[1]
            )
            cdi["spurious_output"] = re.search(
                r"(?<=FINAL\n\n)[\w\W]+weighted.*", a, re.MULTILINE
            ).group(0)

        cdi["time"] = datetime.now(tz=timezone.utc).isoformat()
        newtileset["info"]["core"].append(cdi)

        finish.finish(
            basename + ".save",
            designname=basename + ".mfe",
            seqsname=basename + ".seqs",
            strandsname=None,
            run_kin=False,
            cleanup=False,
            trials=0,
            time=0,
            temp=27,
            conc=1,
            spurious=False,
            spurious_time=0,
        )  # FIXME: shouldn't need so many options.

        tileset_with_strands = newtileset._load_pepper_output_files(basename)

        # Ensure:
        tileset.ends.merge(
            tileset_with_strands.tiles.glues_from_tiles()
        )  # Ends still fit
        for tile in tileset_with_strands.tiles:
            oldtile = tileset.tiles[tile.name]
            if "fullseqs" in oldtile.keys():
                for old, new in zip(oldtile["fullseqs"], tile["fullseqs"]):
                    seq.merge(old, new)  # old tile sequences remain
            assert oldtile.ends == tile.ends

        # Check that old end sequences remain
        tileset["ends"].merge(tileset_with_strands["ends"])

        return tileset_with_strands

    def _create_pepper_input_files(tileset, basename):
        # Are we creating adapters in Pepper?
        if tileset.seed and seeds.seedtypes[tileset["seed"]["type"]].needspepper:
            seedclass = seeds.seedtypes[tileset["seed"]["type"]]
            createadapts = True
        else:
            createadapts = False

        fixedfile = open(basename + ".fix", "w")
        # We first need to create a fixed sequence list/file for pepper.
        # Add fixed sticky end and adjacent tile sequences.
        for end in tileset["ends"]:
            if "fseq" not in end.keys():
                continue
            seq = end["fseq"][1:-1]
            if end["type"] == "TD":
                adj = end["fseq"][-1]
                cadj = end["fseq"][0]  # FIXME: WAS [1], OFF BY ONE!
            elif end["type"] == "DT":
                adj = end["fseq"][0]
                cadj = end["fseq"][-1]  # FIXME: WAS [1], OFF BY ONE!
            else:
                print("warning! end {} not recognized".format(end["name"]))
            fixedfile.write("signal e_{0} = {1}\n".format(end["name"], seq.upper()))
            fixedfile.write("signal a_{0} = {1}\n".format(end["name"], adj.upper()))
            fixedfile.write("signal c_{0} = {1}\n".format(end["name"], cadj.upper()))
            # If we are creating adapter tiles in Pepper, add origami-determined
            # sequences
        if createadapts:
            for i, core in enumerate(seedclass.cores, 1):
                fixedfile.write("signal origamicore_{0} = {1}\n".format(i, core))

        # Now we'll create the system file in parts.
        importlist = set()
        compstring = ""

        for tile in tileset.tiles:
            if tile.is_fake:
                continue
            e = [[], []]
            for end in tile["ends"]:
                if end == "hp":
                    continue
                    # skip hairpins, etc that aren't designed by stickydesign
                e[0].append("e_" + end.replace("/", "*"))
                if end[-1] == "/":
                    a = "c_" + end[:-1] + "*"
                else:
                    a = "a_" + end
                e[1].append(a)
            s1 = " + ".join(e[0])
            s2 = " + ".join(e[1])
            tiletype = tile["structure"]
            if "extra" in tile.keys():
                tiletype += "_" + tile["extra"]
            compstring += "component {} = {}: {} -> {}\n".format(
                tile["name"], tiletype, s1, s2
            )
            importlist.add(tiletype)
            if "fullseqs" in tile.keys():
                fixedfile.write(
                    "structure {}-tile = ".format(tile["name"])
                    + "+".join([seq.upper() for seq in tile["fullseqs"]])
                    + "\n"
                )

        if createadapts:
            importlist, compstring = seedclass._create_pepper_input_files(
                tileset["seed"], importlist, compstring
            )

        with open(basename + ".sys", "w") as sysfile:
            sysfile.write("declare system {}: ->\n\n".format(basename))
            sysfile.write("import " + ", ".join(importlist) + "\n\n")
            sysfile.write(compstring)

    def _load_pepper_output_files(tileset, basename):
        import re

        # Are we creating adapters in Pepper?
        # if seeds.seedtypes[tileset['seed']['type']].needspepper:
        #     seedclass = seeds.seedtypes[tileset['seed']['type']]
        #     createadapts = True

        tset = copy.deepcopy(tileset)

        seqsstring = open(basename + ".seqs").read()

        for tile in tset.tiles:
            if tile.is_fake:
                continue
            pepperstrands = re.compile(
                "strand " + tile["name"] + "-([^ ]+) = ([^\n]+)"
            ).findall(seqsstring)
            tile["fullseqs"] = tilestructures.order_pepper_strands(pepperstrands)

        # for adapter in tset['seed']['adapters']:
        #    pepperstrands = re.compile('strand ' + adapter['name'] +
        #                               '-([^ ]+) = ([^\n]+)').findall(
        #                                   seqsstring)
        #    adapter['fullseqs'] = tilestructures.order_pepper_strands(
        #        pepperstrands)

        return tset

    def create_guard_strand_sequences(tileset):
        """Given a tileset dictionary with core tile sequences,
        create guard strand sequences.

        Returns
        -------

        TileSet
            A TileSet with generated guard strand sequences.
        """
        tset = tileset.copy()

        for guard in tset["guards"]:
            tile = tset.tiles[guard[0]]
            guard.append(wc(tile["fullseqs"][guard[1] - 1]))

        return tset

    def create_adapter_sequence_diagrams(tileset, filename, *options):
        """Create sequence diagrams of adapters for the seed."""
        from lxml import etree
        import pkg_resources
        import os.path

        base = etree.parse(
            pkg_resources.resource_stream(
                __name__, os.path.join("seqdiagrambases", "blank.svg")
            )
        )
        baseroot = base.getroot()
        pos = 200
        for adapterdef in tileset["seed"]["adapters"]:

            seedclass = seeds.seedtypes[tileset["seed"]["type"]]
            group = seedclass.create_adapter_sequence_diagram(adapterdef)

            group.attrib["transform"] = "translate(0,{})".format(pos)
            pos += 200
            baseroot.append(group)

        base.write(filename)

    def run_xgrow(
        self,
        xgrowparams={},
        perfect=False,
        rotate=False,
        energetics=None,
        ui=False,
        output=None,
        labelsonly=False,
        onlyreal=True,
    ):
        """Run Xgrow for the system.

        Parameters
        ----------

        xgrowparams : dict, optional
            Extra Xgrow parameters.  For example, to change Gse and use 2px block, use
        {'Gse': 9.2, 'block': 2}.

        perfect : bool
            If True, each end binds to its complement by strength 1.  If False,
            stickydesign energetics are used.  Note that TileSets with odd ends
            (eg, with different "ends" that are actually just the same end with
            different adjacent bases, as in the original COUNT system) will likely
            not behave as expected with perfect=True.

        rotate : bool, optional
            If False (Default), tiles are simply converted to Xgrow tiles.  If True,
            rotated/flipped copies of each tile are created in order to simulate valid
            tile orientations other than those intended.

        energetics : stickydesign.Energetics, optional
            The energetics model to use.  If not specified, DEFAULT_ENERGETICS is used.


        output : str or list, optional
            either a string or list, specifying output options.  If a
            string, one of 'final', 'array', or 'trace', corresponding to 'datafile',
            'arrayfile' and 'tracefile', respectively.  If a list of multiple, then
            do those.  These will manage the output, and return the data in usable
            form.

        Returns
        -------
        list or various
            data as requested from output.  See xgrow.run for more information.

        """
        import xgrow

        return xgrow.run(
            self.generate_xgrow_dict(
                perfect=perfect,
                rotate=rotate,
                energetics=energetics,
                labelsonly=labelsonly,
                onlyreal=onlyreal,
            ),
            extraparams=xgrowparams,
            outputopts=output,
            ui=ui,
        )

    # FIXME: NEED TO IMPLEMENT THIS WITH BETTER CODE
    # def sensitivity_classes(ts, count=False, _maxorder=2):
    #    return sensitivity.sensitivity_classes(
    #        ts, count=False, _maxorder=_maxorder)

    def generate_xgrow_dict(
        self,
        perfect: bool = False,
        rotate: bool = False,
        energetics=None,
        labelsonly: bool = False,
        onlyreal: bool = True,
    ):
        """Generate a Xgrow tileset dict.

        Parameters
        ----------

        perfect : bool
            If True, each end binds to its complement by strength 1.  If False,
            stickydesign energetics are used.  Note that TileSets with odd ends
            (eg, with different "ends" that are actually just the same end with
            different adjacent bases, as in the original COUNT system) will likely
            not behave as expected with perfect=True.

        rotate : bool, optional
            If False (Default), tiles are simply converted to Xgrow tiles.  If True,
            rotated/flipped copies of each tile are created in order to simulate valid
            tile orientations other than those intended.

        energetics : stickydesign.Energetics, optional
            The energetics model to use.  If not specified, DEFAULT_ENERGETICS is used.

        Returns
        -------

        dict
            An xgrow-compatible tileset definition, suitable for, eg, xgrow.run.
        """

        # Combine ends and tile-specified adjacents
        newtiles = []
        newends = []
        doubleends = []
        doubles = []
        vdoubleends = []
        vdoubles = []
        ts = copy.deepcopy(self)
        if ts.seed:
            seedtype = seeds.seedtypes[ts.seed["type"]]
            newtiles.append(
                {
                    "name": "origami",
                    "edges": ["origami", "origami", "origami", "origami"],
                    "stoic": 0,
                    "color": "white",
                }
            )

            atiles = [None] * len(seedtype.cores)
            to_use = []
            # If we have use_adapters, use that, otherwise use every adapter:
            if "use_adapters" in ts["seed"]:
                for tilename in ts["seed"]["use_adapters"]:
                    try:
                        tile = [
                            x
                            for x in ts["seed"]["adapters"]
                            if x.get("name") == tilename
                        ][0]
                        to_use.append(tile)
                    except IndexError as e:
                        raise Exception("Can't find {}".format(tilename)) from e
            else:
                to_use = ts["seed"]["adapters"]

            for tile in to_use:
                newtile = {}
                if "ends" in tile.keys():
                    newtile["edges"] = (
                        ["origami"]
                        + [re.sub("/", "_c", x) for x in tile["ends"]]
                        + ["origami"]
                    )
                else:
                    mtile = ts.tiles[tile["tilebase"]]
                    newtile["edges"] = (
                        ["origami"]
                        + [
                            re.sub("/", "_c", x)
                            for x in mtile.ends[
                                seedtype._mimicadapt[mtile.structure.name]["ends"]
                            ]
                        ]
                        + ["origami"]
                    )
                newtile["name"] = tile.get("name", "")
                newtile["stoic"] = 0
                newtile["color"] = "white"
                atiles[tile["loc"] - 1] = newtile
            for tile in atiles:
                if tile:
                    newtiles.append(tile)
                else:
                    newtiles.append(
                        {
                            "name": "emptyadapt",
                            "edges": ["origami", 0, 0, "origami"],
                            "stoic": 0,
                            "color": "white",
                        }
                    )

        if onlyreal:
            ts.tiles = TileList([x for x in ts.tiles if "fake" not in x])

        if rotate:
            rotatedtiles = []
            for tile in ts.tiles:
                # only include rotated tiles that aren't identical (handles symmetric tiles)
                trot = [t for t in tile.rotations if (t["ends"] != tile["ends"])]
                td = []
                for i, tr in enumerate(trot):
                    for j in range(0, i):
                        if tr["ends"] == trot[j]["ends"]:
                            td.append(i)
                            break
                for i in reversed(td):
                    del trot[i]
                rotatedtiles += trot
            ts.tiles += rotatedtiles

        for tile in ts.tiles:
            if "rotation" in tile.keys():
                tile["name"] = tile["name"] + "_rot{}".format(tile["rotation"])
            if re.match("tile_daoe_3up", tile.structure.name) or re.match(
                "tile_daoe_5up", tile.structure.name
            ):
                newtile = {}
                newtile["edges"] = [re.sub("/", "_c", x) for x in tile["ends"]]
                if "name" in tile:
                    newtile["name"] = tile["name"]
                if "conc" in tile:
                    newtile["stoic"] = tile["conc"]
                if ("color" in tile) and (not labelsonly):
                    newtile["color"] = tile["color"]
                elif labelsonly:
                    if "label" in tile.keys():
                        newtile["color"] = "white"
                    else:
                        newtile["color"] = "gray50"
                newtiles.append(newtile)

            if re.match("tile_daoe_doublehoriz", tile.structure.name):
                newtile1 = {}
                newtile2 = {}
                newtile1["edges"] = (
                    [re.sub("/", "_c", x) for x in tile["ends"][0:1]]
                    + [tile["name"] + "_db"]
                    + [re.sub("/", "_c", x) for x in tile["ends"][4:]]
                )
                newtile2["edges"] = [
                    re.sub("/", "_c", x) for x in tile["ends"][1:4]
                ] + [tile["name"] + "_db"]
                newtile1["name"] = tile["name"] + "_left"
                newtile2["name"] = tile["name"] + "_right"

                doubleends.append(tile["name"] + "_db")
                doubles.append((newtile1["name"], newtile2["name"]))

                if "conc" in tile:
                    newtile1["stoic"] = tile["conc"]
                    newtile2["stoic"] = tile["conc"]

                if ("color" in tile) and (not labelsonly):
                    newtile1["color"] = tile["color"]
                    newtile2["color"] = tile["color"]
                elif labelsonly:
                    if "label" in tile.keys():
                        newtile1["color"] = "white"
                        newtile2["color"] = "white"
                    else:
                        newtile1["color"] = "gray50"
                        newtile2["color"] = "gray50"

                newtiles.append(newtile1)
                newtiles.append(newtile2)
            if re.match("tile_daoe_doublevert", tile.structure.name):
                newtile1 = {}
                newtile2 = {}
                newtile1["edges"] = (
                    [re.sub("/", "_c", x) for x in tile["ends"][0:2]]
                    + [tile["name"] + "_db"]
                    + [re.sub("/", "_c", x) for x in tile["ends"][5:]]
                )
                newtile2["edges"] = [tile["name"] + "_db"] + [
                    re.sub("/", "_c", x) for x in tile["ends"][2:5]
                ]
                newtile1["name"] = tile["name"] + "_top"
                newtile2["name"] = tile["name"] + "_bottom"

                vdoubleends.append(tile["name"] + "_db")
                vdoubles.append((newtile1["name"], newtile2["name"]))

                if "conc" in tile:
                    newtile1["stoic"] = tile["conc"]
                    newtile2["stoic"] = tile["conc"]

                if ("color" in tile) and (not labelsonly):
                    newtile1["color"] = tile["color"]
                    newtile2["color"] = tile["color"]
                elif labelsonly:
                    if "label" in tile.keys():
                        newtile1["color"] = "white"
                        newtile2["color"] = "white"
                    else:
                        newtile1["color"] = "gray50"
                        newtile2["color"] = "gray50"

                newtiles.append(newtile1)
                newtiles.append(newtile2)

        newends.append({"name": "origami", "strength": 100})

        for end in doubleends:
            newends.append({"name": end, "strength": 10})
        for end in vdoubleends:
            newends.append({"name": end, "strength": 10})

        # check for whether we need to use perfect:
        if not perfect:
            for end in ts.allends:
                if (end["type"] in {"TD", "DT"}) and "fseq" not in end.keys():
                    perfect = True
                    SELOGGER.warn(
                        "setting perfect=True, because {} has no sequence.".format(
                            end["name"]
                        )
                    )
                    break

        gluelist = []
        if not perfect:
            glueends = {"DT": [], "TD": []}
            for end in ts.allends:
                newends.append({"name": end["name"], "strength": 0})
                newends.append({"name": end["name"] + "_c", "strength": 0})
                if (end["type"] == "TD") or (end["type"] == "DT"):
                    glueends[end["type"]].append((end["name"], end["fseq"]))

            if energetics:
                ef = energetics
            else:
                ef = DEFAULT_ENERGETICS

            eavg = {}
            for t in ["DT", "TD"]:
                names, fseqs = zip(*glueends[t])
                ea = sd.endarray(fseqs, t)
                eavg[t] = np.average(ef.matching_uniform(ea))
            eavg_combined = (eavg["DT"] + eavg["TD"]) / 2.0

            for t in ["DT", "TD"]:
                names, fseqs = zip(*glueends[t])
                allnames = names + tuple(x + "_c" for x in names)
                ea = sd.endarray(fseqs, t)
                ar = sd.energy_array_uniform(ea, ef) / eavg_combined
                for i1, n1 in enumerate(names):
                    for i2, n2 in enumerate(allnames):
                        gluelist.append([n1, n2, max(float(ar[i1, i2]), 0.0)])

        else:
            if "ends" not in ts.keys():
                ts["ends"] = []
            endsinlist = set(e["name"] for e in ts["ends"])
            endsintiles = set()
            for tile in ts.tiles:
                endsintiles.update(
                    re.sub("/", "", e) for e in tile["ends"] if e != "hp"
                )
            for end in ts["ends"] + list({"name": e} for e in endsintiles):
                newends.append({"name": end["name"], "strength": 0})
                newends.append({"name": end["name"] + "_c", "strength": 0})
                gluelist.append([end["name"], end["name"] + "_c", 1.0])

        newends.append({"name": "hp", "strength": 0})

        xga = {}
        xga["doubletiles"] = [list(x) for x in doubles]
        xga["vdoubletiles"] = [list(x) for x in vdoubles]

        xga.update(ts.get("xgrow_options", dict()))

        # if not perfect:
        #    xga['gse_calc_avg'] = eavg_combined

        sts = {"tiles": newtiles, "bonds": newends, "xgrowargs": xga, "glues": gluelist}

        return sts

    def plot_se_hists(
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
            ends = tileset["ends"]
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

    def plot_se_lv(
        self,
        all_energetics=None,
        energetics_names=None,
        pltcmd=None,
        title=None,
        **kwargs
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
            self.ends.to_endarrays(), all_energetics, energetics_names
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

    def plot_adjacent_regions(tileset, energetics=None):
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

    def plot_side_strands(tileset, energetics=None):
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

    def reduce_tiles(
        tileset,
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
            tileset, preserve, tries, threads, returntype, best, key, initequiv
        )

    def reduce_ends(
        tileset,
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
            tileset, preserve, tries, threads, returntype, best, key, initequiv
        )

    def latticedefects(ts, direction="e", depth=2, pp=True, rotate=False):
        """
        Calculate and show possible small lattice defect configurations.
        """
        from . import latticedefect

        return latticedefect.latticedefects(
            ts, direction=direction, depth=depth, pp=pp, rotate=rotate
        )

    def apply_equiv(ts, equiv):
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
        return fastreduce._FastTileSet(ts).applyequiv(ts, equiv)


RoundTripRepresenter.add_representer(TileSet, RoundTripRepresenter.represent_dict)
