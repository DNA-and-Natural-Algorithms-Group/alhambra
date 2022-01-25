from __future__ import annotations

import copy
import warnings

import stickydesign as sd
import stickydesign.stickydesign2 as sd2
from peppercompiler.DNA_classes import wc
from ruamel.yaml.comments import CommentedMap
from ruamel.yaml.representer import RoundTripRepresenter

from . import seq
from .util import (
    DEFAULT_MM_ENERGETICS_NAMES,
    DEFAULT_SD2_MULTIMODEL_ENERGETICS,
    NamedList,
)


class End(CommentedMap):
    """A class representing a single end of some type, with or without a sequence."""

    def __str__(self):
        if self.fseq and self.seq and self.comp:
            if self.etype == "DT":
                s = self.seq[0] + "-" + self.seq[1:]
                c = self.comp[0] + "-" + self.comp[1:]
            elif self.etype == "TD":
                s = self.seq[:-1] + "-" + self.seq[-1]
                c = self.comp[:-1] + "-" + self.comp[-1]
            else:
                raise ValueError(self.etype)
            return "<end {} ({}{}): {} | {}>".format(
                self["name"], self["type"], len(self.seq), s, c
            )
        else:
            return "<end {} ({})>".format(self["name"], self.get("type", "?"))

    @property
    def fseq(self) -> str:
        """The fseq / full sequence of the End, including the adjacent base on the end
        and the complement of the adjacent base on the complement."""
        return self.get("fseq", None)

    @fseq.setter
    def fseq(self, value):
        if self.get("fseq", None) and len(value) != len(self["fseq"]):
            warnings.warn("Changing end length")
        self["fseq"] = value

    @fseq.deleter
    def fseq(self):
        del self["fseq"]

    @property
    def name(self) -> str:
        return self["name"]

    @name.setter
    def name(self, value):
        self["name"] = value

    @property
    def use(self) -> int:
        return self.get("use", 0b0)

    @use.setter
    def use(self, value: int):
        self["use"] = value

    @property
    def strength(self) -> int:
        """The strength of the End (as int)"""
        return self.get("strength", 1)

    @strength.setter
    def strength(self, value: int):
        self["strength"] = value

    @property
    def seq(self):
        """The end sequence (of just the end) of the End, as a string."""
        if not self.fseq:
            return None
        if self.etype == "TD":
            return self.fseq[1:]
        elif self.etype == "DT":
            return self.fseq[:-1]

    @property
    def comp(self):
        """The complement end sequences of the End, as a string."""
        if not self.fseq:
            return None
        if self.etype == "TD":
            return wc(self.fseq[:-1].upper()).lower()
        elif self.etype == "DT":
            return wc(self.fseq[1:].upper()).lower()

    @property
    def etype(self):
        """The end type of the end."""
        return self["type"]

    def merge(self: End, end2: End) -> End:
        """Given ends end1 and end2, assuming they describe the same sticky
        end, merge them into a single end, combining information from each and
        enforcing that the two input ends consistently make a single output
        end.
        """
        # Of things in the ends: fseq might need special care, everything
        # else just needs to match.
        out = copy.deepcopy(self)
        for i, v in end2.items():
            if i in out.keys() and i == "fseq":
                # we merge sequences
                out[i] = seq.merge(out[i], v)
            elif i in out.keys() and i == "use":
                out[i] = out[i] | v
            elif i in out.keys() and out[i] != v:
                # we insisted all others must be equal
                raise ValueError(
                    "end1 has {}={}, end2 has {}={}".format(i, out[i], i, v)
                )
            elif i not in out.keys():
                out[i] = copy.deepcopy(v)
        return out


class EndList(NamedList):
    """A list of End instances, which can be merged together, converted to
    stickydesign.endarray instances, etc."""

    def __init__(self, val=[]):
        """Create an EndList instance.

        Parameters
        ----------

        val : iterable
            an enumerable iterable of objects that can initialize End instances.
        """
        NamedList.__init__(self, val)
        for i, end in enumerate(self):
            self[i] = End(end)

    def to_endarrays(self):
        """Return stickydesign endarrays of each type of (non-hairpin) end.

        Returns
        -------

        list of stickydesign.endarray
            the endarrays of ends in the system.
        """

        endtypes = {x["type"] for x in self}
        endtypes = endtypes - {"hp", "hairpin"}
        return list(
            sd.endarray([x["fseq"] for x in self if x["type"] == y], y)
            for y in endtypes
        )

    @property
    def _epas(self):
        """Return stickydesign2 EndPairArrays"""
        endtypes = {x["type"] for x in self}
        endtypes = endtypes - {"hp", "hairpin"}
        epat = {"TD": sd2.EndPairArrayTD, "DT": sd2.EndPairArrayDT}
        return list(epat[y]([x.fseq for x in self if x["type"] == y]) for y in endtypes)

    @property
    def _epans(self):
        """Return stickydesign2 EndPairArrays"""
        endtypes = {x["type"] for x in self}
        endtypes = endtypes - {"hp", "hairpin"}
        return list(list([x.name for x in self if x["type"] == y]) for y in endtypes)

    def _pandas_data(
        self,
        models=DEFAULT_SD2_MULTIMODEL_ENERGETICS,
        modelnames=DEFAULT_MM_ENERGETICS_NAMES,
    ):
        import stickydesign.stickydesign2.plots as s2pl

        return s2pl._pandas_data(self._epas, self._epans, models, modelnames)

    def update(self: EndList, endlist2: EndList, fail_immediate=False):
        # Check consistency of each NamedList
        self.check_consistent()
        try:
            endlist2.check_consistent()
        except AttributeError:
            pass

        exceptions = []
        errors = []

        for end in endlist2:
            if end.name in self.keys():
                try:
                    self[end.name] = self[end.name].merge(end)
                except ValueError as e:
                    if fail_immediate:
                        raise e
                    else:
                        exceptions.append(e)
                        errors.append((end.name, self[end.name], end))
            else:
                self.append(copy.deepcopy(end))

        if errors:
            errorstring = " ".join([e[0] for e in errors])
            raise ValueError(
                "Errors merging {}".format(errorstring), exceptions, errors, self
            )

        return self

    def merge(self, endlist2: EndList, fail_immediate=False) -> EndList:
        """
        Given end lists `endlist1` and `endlist2`, merge the two lists, using
        `merge_ends` to merge any named ends that are present in both.

        Parameters
        ----------

        endlist1: NamedList of sticky ends.

        endlist2: NamedList OR just a list of sticky ends.  If it just a list,
            which may have multiple copies of the same named sticky end, each
            sticky end in order is merged.

        fail_immediate: (default False) if True, fail immediately on a
        failed merge, passing through the ValueError from the end merge.
        If False, finish merging the two lists, then raise a
        ValueError listing *all* ends that failed to merge.

        in_place: (default False) if True, do merging in place in endlist1.
        Note the merged and added ends from endlist2 will be copies regardless.

        Returns
        -------

        EndList
            a merged list of sticky ends


        Exceptions
        ----------

        ValueError: In the event of a failed merge, when named ends
        cannot be merged.  If fail_immediate is True, then this is passed
        through from merge_ends.  If fail_immediate is False, then the
        ValueError has the following args: ("message", [exceptions],
        [failed_name,failed_end1,failed_end2) ...],out)

        """
        return copy.deepcopy(self).update(endlist2)


RoundTripRepresenter.add_representer(EndList, RoundTripRepresenter.represent_list)
RoundTripRepresenter.add_representer(End, RoundTripRepresenter.represent_dict)
