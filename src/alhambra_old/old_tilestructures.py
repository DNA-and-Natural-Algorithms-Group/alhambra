import copy
import logging
import os.path
import re
import warnings
from abc import ABCMeta, abstractmethod
from typing import Literal, Optional, Union

import cssutils

# Color dictionary for xgrow colors...
import pkg_resources
from lxml import etree

from . import seq as sq
from .edotparen import check_edotparen_consistency, check_edotparen_sequence
from .ends import End, EndList

cssutils.log.setLevel(logging.FATAL)


rgbv = pkg_resources.resource_stream(__name__, os.path.join("data", "rgb.txt"))
xcolors = {
    " ".join(y[3:]): "rgb({},{},{})".format(y[0], y[1], y[2])
    for y in [x.decode().split() for x in rgbv]
}
del rgbv


class TileStructure(object, metaclass=ABCMeta):
    @property
    @abstractmethod
    def edotparen(self) -> Union[str, None]:
        ...

    @property
    @abstractmethod
    def _endlocs(self) -> list[tuple[int, Optional[int], Optional[int]]]:
        ...

    def check_consistent(self):
        if self.edotparen:
            check_edotparen_consistency(self.edotparen)
        else:
            warnings.warn("No edotparen")  # FIXME

    @property
    def numends(self):
        return len(self._endlocs)

    def check_strands(self, tile):
        check_edotparen_sequence(self.edotparen, "+".join(tile.strands))

    @property
    def name(self):  # FIXME: bad!
        return self.__class__.__name__


class tile_daoe(TileStructure):
    @property
    @abstractmethod
    def _endtypes(self) -> list[Literal["DT", "TD", "hairpin"]]:
        ...

    @abstractmethod
    def _seqdiagseqstrings(self, tile) -> list[str]:
        ...

    def sequence_diagram(self, tile):
        ttype = tile.structure.name
        from lxml import etree

        base_svg = etree.parse(
            pkg_resources.resource_stream(
                __name__, os.path.join("seqdiagrambases", ttype + ".svg")
            )
        )

        strings = (
            self._seqdiagseqstrings(tile)
            + [
                e
                for e, t in self.tile_ends(tile)
                if not (t in ("hairpin", "blunt", "inert"))
            ]
            + [tile["name"]]
        )

        texts = base_svg.findall("//{http://www.w3.org/2000/svg}text")

        for text, s in zip(texts, strings):
            if text.getchildren():
                text.getchildren()[0].text = s
            else:
                text.text = s

        return base_svg.xpath(
            "/svg:svg/svg:g", namespaces={"svg": "http://www.w3.org/2000/svg"}
        )[0]

    def tile_ends(self, tile):
        return zip(tile["ends"], self._endtypes)

    def check_sequence(self, tile):
        try:
            check_edotparen_sequence(self.edotparen, "+".join(tile.strands))
        except ValueError as e:
            raise ValueError(
                "{} core is inconsistent.".format(tile["name"]), tile["name"]
            ) from e

    def get_endlist(self, tile):
        es = EndList()
        for (strand, start, end), endtype, endname, endinput in zip(
            self._endlocs,
            self._endtypes,
            tile["ends"],
            tile.get("input", [None] * len(tile["ends"])),
        ):

            if endtype in ("blunt", "inert", "hairpin"):
                continue

            if endtype == "DT":
                if not start:
                    raise ValueError
                sl = slice(start - 1, end)
            elif endtype == "TD":
                if not end:
                    raise ValueError
                sl = slice(start, end + 1)
            else:
                sl = None

            if sl and tile.strands:
                seq = tile.strands[strand][sl]
                if endtype == "DT":
                    seq = (seq + "n").lower()
                elif endtype == "TD":
                    seq = ("n" + seq).lower()
            else:
                seq = None

            isc = False
            if endname[-1] == "/":
                endname = endname[:-1]
                isc = True
                if seq:
                    seq = sq.revcomp(seq)

            e = End({"name": endname, "type": endtype})

            if seq:
                e.fseq = seq

            if endinput is not None:
                if isc:
                    m = 0b11
                else:
                    m = 0b00
                if endinput == -1:
                    use = 0b00 ^ m
                elif endinput == 0:
                    use = 0b10 ^ m
                elif endinput == 1:
                    use = 0b01 ^ m
                e.use = use

            es.update(EndList([e]))  # FIXME
        return es

    def orderableseqs(self, tile):
        seqs = copy.deepcopy(tile.strands)
        assert "label" not in tile.keys()
        return [("{}-{}".format(tile.name, i + 1), seq) for i, seq in enumerate(seqs)]


class tile_daoe_single(tile_daoe):
    """Base class for single DAO-E tiles."""

    _abase = "tile_daoe_single"
    _a_endlocs = ["north", "east", "south", "west"]
    double = False
    singleends = ((0, 1, 2, 3),)
    _dirs = (0, 1, 2, 3)

    def _seqdiagseqstrings(self, tile) -> list[str]:
        s = tile.strands
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16],
            s[1][16:24],
            s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            (s[2][:8] + "--" + s[2][8:16])[::-1],
            s[2][16:24][::-1],
            s[2][24:32],
            (s[2][32:40] + "--" + s[2][40:48]),
            (s[3][:5] + "--" + s[3][5:13])[::-1],
            (s[3][:-6:-1] + "--" + s[3][-6:-14:-1])[::-1],
        ]

    def _short_bound_full(self, tile):
        s = tile.strands
        return [s[0][5:-5], s[3][5:-5]]

    def _side_bound_regions(self, tile):
        s = tile.strands
        return [s[0][5 : 5 + 8], s[0][-5 - 8 : -5], s[3][5 : 5 + 8], s[3][-5 - 8 : -5]]


class tile_daoe_5up(tile_daoe_single):
    _orient = ("3", "5")  # SW, then NE of whole tile
    _endtypes = ["TD", "TD", "DT", "DT"]
    # endlocs is strand, loc, length
    _endlocs = [(0, 0, 5), (3, 0, 5), (3, 21, None), (0, 21, None)]
    # valid edotparen for both 3up and 5up
    edotparen = "5.16(5.+8)16[16{8)+8(16]16}8(+5.16)5."

    @property
    def rotations(self):
        return [
            (tile_daoe_3up, (3, 2, 1, 0)),
            (tile_daoe_5up, (1, 0, 3, 2)),
            (tile_daoe_3up, (2, 3, 0, 1)),
        ]

    def orderableseqs(self, tile):
        seqs = copy.deepcopy(tile.strands)
        if tile.get("label", None) == "both":
            assert seqs[1][16] == "T"
            assert seqs[2][16] == "T"
            seqs[1] = seqs[1][:16] + "/iBiodT/" + seqs[1][17:]
            seqs[2] = seqs[2][:16] + "/iBiodT/" + seqs[2][17:]
            assert re.sub("/iBiodT/", "T", seqs[1]) == tile.strands[1]
            assert re.sub("/iBiodT/", "T", seqs[2]) == tile.strands[2]
        elif "label" in tile.keys():
            raise NotImplementedError
        return [("{}-{}".format(tile.name, i + 1), seq) for i, seq in enumerate(seqs)]


class tile_daoe_3up(tile_daoe_single):
    _orient = ("5", "3")
    _endtypes = ["DT", "DT", "TD", "TD"]
    _endlocs = [(0, 21, None), (3, 21, None), (3, 0, 5), (0, 0, 5)]
    # valid edotparen for both 3up and 5up
    edotparen = "5.16(5.+8)16[16{8)+8(16]16}8(+5.16)5."

    @property
    def rotations(self):
        return [
            (tile_daoe_5up, (3, 2, 1, 0)),
            (tile_daoe_3up, (1, 0, 3, 2)),
            (tile_daoe_5up, (2, 3, 0, 1)),
        ]

    def orderableseqs(self, tile):
        seqs = copy.deepcopy(tile.strands)
        if tile.get("label", None) == "both":
            assert seqs[1][31] == "T"
            assert seqs[2][31] == "T"
            seqs[1] = seqs[1][:31] + "/iBiodT/" + seqs[1][32:]
            seqs[2] = seqs[2][:31] + "/iBiodT/" + seqs[2][32:]
            assert re.sub("/iBiodT/", "T", seqs[1]) == tile.strands[1]
            assert re.sub("/iBiodT/", "T", seqs[2]) == tile.strands[2]
        elif "label" in tile.keys():
            raise NotImplementedError
        return [("{}-{}".format(tile.name, i + 1), seq) for i, seq in enumerate(seqs)]


class tile_daoe_5up_1h(tile_daoe_5up):
    _orient = ("3", "5")
    _endtypes = ["hairpin", "TD", "DT", "DT"]


class tile_daoe_5up_3h(tile_daoe_5up):
    _orient = ("3", "5")
    _endtypes = ["TD", "TD", "hairpin", "DT"]


class tile_daoe_5up_3b(tile_daoe_5up):
    _orient = ("3", "5")
    _endtypes = ["TD", "TD", "hairpin", "DT"]


class tile_daoe_5up_4h(tile_daoe_5up):
    _orient = ("3", "5")
    _endtypes = ["TD", "TD", "DT", "hairpin"]


class tile_daoe_5up_2h(tile_daoe_5up):
    _orient = ("3", "5")
    _endtypes = ["TD", "hairpin", "DT", "DT"]
    _endlocs = [(0, 0, 5), (3, 0, 18), (3, -5, None), (0, -5, None)]

    edotparen = "5.16(5.+8)16[16{8)+8(16]16}8(+5.16)7(4.7)"

    @property
    def rotations(self):
        return [
            (tile_daoe_3up_3h, (3, 2, 1, 0)),
            (tile_daoe_5up_1h, (1, 0, 3, 2)),
            (tile_daoe_3up_4h, (2, 3, 0, 1)),
        ]

    def _short_bound_full(self, tile):
        s = tile.strands
        return [s[0][5:-5], s[3][18:-5]]

    def _side_bound_regions(self, tile):
        s = tile.strands
        return [
            s[0][5 : 5 + 8],
            s[0][-5 - 8 : -5],
            s[3][18 : 18 + 8],
            s[3][-5 - 8 : -5],
        ]

    def _seqdiagseqstrings(self, tile):
        s = copy.copy(tile.strands)
        hp = s[3][0:13]
        s[3] = s[3][13:]
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16],
            s[1][16:24],
            s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            (s[2][:8] + "--" + s[2][8:16])[::-1],
            s[2][16:24][::-1],
            s[2][24:32],
            (s[2][32:40] + "--" + s[2][40:48]),
            hp[0:5] + "-" + hp[5:7] + "-" + hp[7:9],
            (hp[9:11] + "-" + hp[11:13] + "-")[::-1],
            (s[3][:5] + "--" + s[3][5:13])[::-1],
            (s[3][:-6:-1] + "--" + s[3][-6:-14:-1])[::-1],
        ]


class tile_daoe_3up_1h(tile_daoe_3up):
    _endtypes = ["hairpin", "DT", "TD", "TD"]


class tile_daoe_3up_3h(tile_daoe_3up):
    _endtypes = ["DT", "DT", "hairpin", "TD"]


class tile_daoe_3up_4h(tile_daoe_3up):
    _endtypes = ["DT", "DT", "TD", "hairpin"]


class tile_daoe_3up_2h(tile_daoe_3up):
    def _short_bound_full(self, tile):
        s = tile.strands
        return [s[0][5:-5], s[3][5:-18]]

    def _side_bound_regions(self, tile):
        s = tile.strands
        return [
            s[0][5 : 5 + 8],
            s[0][-5 - 8 : -5],
            s[3][5 : 5 + 8],
            s[3][-18 - 8 : -18],
        ]

    @property
    def rotations(self):
        return [
            (tile_daoe_5up_3h, (3, 2, 1, 0)),
            (tile_daoe_3up_1h, (1, 0, 3, 2)),
            (tile_daoe_5up_4h, (2, 3, 0, 1)),
        ]

    _endtypes = ["DT", "hairpin", "TD", "TD"]
    _endlocs = [(0, 21, None), (3, 21, None), (3, 0, 5), (0, 0, 5)]

    def _seqdiagseqstrings(self, tile):
        s = copy.copy(tile.strands)
        hp = s[3][26:]
        s[3] = s[3][:26]
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16],
            s[1][16:24],
            s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            (s[2][:8] + "--" + s[2][8:16])[::-1],
            s[2][16:24][::-1],
            s[2][24:32],
            (s[2][32:40] + "--" + s[2][40:48]),
            (s[3][:5] + "--" + s[3][5:13])[::-1],
            (s[3][:-6:-1] + "--" + s[3][-6:-14:-1])[::-1],
            "-" + hp[0:2] + "-" + hp[2:4],
            (hp[4:6] + "-" + hp[6:8] + "-" + hp[8:13])[::-1],
        ]


class tile_doublehoriz(object):
    _dirs = (0, 0, 1, 2, 2, 3)
    _abase = "tile_daoe_doublehoriz"
    _a_endlocs = ["northwest", "northeast", "east", "southeast", "southwest", "west"]
    singleends = ((0, None, 4, 5), (1, 2, 3, None))
    double = "h"


class tile_doublevert(object):
    _dirs = (0, 1, 1, 2, 3, 3)
    _abase = "tile_daoe_doublevert"
    _a_endlocs = ["north", "northeast", "southeast", "south", "southwest", "northwest"]
    singleends = ((0, 1, None, 5), (None, 2, 3, 4))
    double = "v"


class tile_daoe_doublehoriz(tile_doublehoriz, tile_daoe):
    ...


class tile_daoe_doublevert(tile_doublevert, tile_daoe):
    ...


class tile_daoe_doublehoriz_35up(tile_daoe_doublehoriz):
    _endtypes = ["DT", "TD", "TD", "DT", "TD", "TD"]
    _orient = ("5", "5")
    _endlocs = [
        (0, -5, None),
        (2, 0, 5),
        (5, 0, 5),
        (5, -5, None),
        (3, 0, 5),
        (0, 0, 5),
    ]

    @property
    def rotations(self):
        return [
            (tile_daoe_doublehoriz_35up, (3, 4, 5, 0, 1, 2)),
            (tile_daoe_doublevert_53up, (5, 4, 3, 2, 1, 0)),
            (tile_daoe_doublevert_53up, (2, 1, 0, 5, 4, 3)),
        ]

    def _short_bound_full(self, tile):
        s = tile.strands
        el = self._endlocs
        return [s[0][el[5][2] : el[0][1]], s[5][el[2][2] : el[3][1]]]

    def _side_bound_regions(self, tile):
        s = tile.strands
        el = self._endlocs

        def g(e):
            if e[2] is None or e[2] == len(s[e[0]]):
                return s[e[0]][-8 + e[1] : e[1]]
            elif e[1] == 0:
                return s[e[0]][e[2] : e[2] + 8]
            else:
                raise ValueError()

        return [g(e) for e in el]

    def _seqdiagseqstrings(self, tile):
        s = tile.strands
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16],
            s[1][16:24],
            s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            s[2][:5] + "--" + s[2][5:13],
            (
                s[2][13:21]
                + "--"
                + s[2][21:26]
                + "--"
                + s[2][26:34]
                + "--"
                + s[2][34:42]
            )[::-1],
            (s[2][42:50])[::-1],
            s[2][50:58],
            s[2][58:66] + "--" + s[2][66:74],
            (s[3][0:5] + "--" + s[3][5:13])[::-1],
            s[3][13:21] + "--" + s[3][21:26] + "--" + s[3][26:34] + "--" + s[3][34:42],
            s[3][42:50],
            s[3][50:58][::-1],
            (s[3][58:66] + "--" + s[3][66:74])[::-1],
            (s[4][0:8] + "--" + s[4][8:16])[::-1],
            (s[4][16:24])[::-1],
            s[4][24:32],
            s[4][32:40] + "--" + s[4][40:48],
            (s[5][0:5] + "--" + s[5][5:13])[::-1],
            s[5][13:21] + "--" + s[5][21:26],
        ]

    edotparen = "5.16(5.+8)16[16{8)+5.29(16]16}8(+5.29)16[16{8)+8(16]16}8(+5.16)5."


class tile_daoe_doublehoriz_53up(tile_daoe_doublehoriz):

    _endtypes = ["TD", "DT", "DT", "TD", "DT", "DT"]
    _orient = ("3", "3")
    _endlocs = [
        (0, 0, 5),
        (2, -5, None),
        (5, -5, None),
        (5, 0, 5),
        (3, -5, None),
        (0, -5, None),
    ]

    def _short_bound_full(self, tile):
        s = tile.strands
        el = self._endlocs
        return [s[0][el[5][2] : el[0][1]], s[5][el[2][2] : el[3][1]]]

    def _side_bound_regions(self, tile):
        s = tile.strands
        el = self._endlocs

        def g(e):
            if e[2] is None or e[2] == len(s[e[0]]):
                return s[e[0]][-8 + e[1] : e[1]]
            elif e[1] == 0:
                return s[e[0]][e[2] : e[2] + 8]
            else:
                raise ValueError()

        return [g(e) for e in el]

    def _seqdiagseqstrings(self, tile):
        s = tile.strands
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16],
            s[1][16:24],
            s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            s[2][:5] + "--" + s[2][5:13],
            (
                s[2][13:21]
                + "--"
                + s[2][21:26]
                + "--"
                + s[2][26:34]
                + "--"
                + s[2][34:42]
            )[::-1],
            (s[2][42:50])[::-1],
            s[2][50:58],
            s[2][58:66] + "--" + s[2][66:74],
            (s[3][0:5] + "--" + s[3][5:13])[::-1],
            s[3][13:21] + "--" + s[3][21:26] + "--" + s[3][26:34] + "--" + s[3][34:42],
            s[3][42:50],
            s[3][50:58][::-1],
            (s[3][58:66] + "--" + s[3][66:74])[::-1],
            (s[4][0:8] + "--" + s[4][8:16])[::-1],
            (s[4][16:24])[::-1],
            s[4][24:32],
            s[4][32:40] + "--" + s[4][40:48],
            (s[5][0:5] + "--" + s[5][5:13])[::-1],
            s[5][13:21] + "--" + s[5][21:26],
        ]


class tile_daoe_doublevert_35up(tile_daoe_doublevert):
    @property
    def rotations(self):
        return [
            (tile_daoe_doublevert_35up, (3, 4, 5, 0, 1, 2)),
            (tile_daoe_doublehoriz_53up, (5, 4, 3, 2, 1, 0)),
            (tile_daoe_doublehoriz_53up, (2, 1, 0, 5, 4, 3)),
        ]

    def _short_bound_full(self, tile):
        s = tile.strands
        el = self._endlocs
        return [s[0][el[5][2] : el[0][1]], s[5][el[2][2] : el[3][1]]]

    def _side_bound_regions(self, tile):
        s = tile.strands
        el = self._endlocs

        def g(e):
            if e[2] is None or e[2] == len(s[e[0]]):
                return s[e[0]][-8 + e[1] : e[1]]
            elif e[1] == 0:
                return s[e[0]][e[2] : e[2] + 8]
            else:
                raise ValueError()

        return [g(e) for e in el]

    _endtypes = ["DT", "DT", "TD", "DT", "DT", "TD"]
    _orient = ("3", "3")
    _endlocs = [
        (0, -5, None),
        (3, -5, None),
        (5, 0, 5),
        (5, -5, None),
        (2, -5, None),
        (0, 0, 5),
    ]


class tile_daoe_doublevert_53up(tile_daoe_doublevert):
    @property
    def rotations(self):  # FIXME: check this
        return [
            (tile_daoe_doublevert_53up, (3, 4, 5, 0, 1, 2)),
            (tile_daoe_doublehoriz_35up, (5, 4, 3, 2, 1, 0)),
            (tile_daoe_doublehoriz_35up, (2, 1, 0, 5, 4, 3)),
        ]

    def _short_bound_full(self, tile):
        raise NotImplementedError
        s = tile.strands
        el = self._endlocs
        return [s[0][el[5][2] : el[0][1]], s[5][el[2][2] : el[3][1]]]

    def _side_bound_regions(self, tile):
        raise NotImplementedError
        s = tile.strands
        el = self._endlocs

        def g(e):
            if e[2] is None or e[2] == len(s[e[0]]):
                return s[e[0]][-8 + e[1] : e[1]]
            elif e[1] == 0:
                return s[e[0]][e[2] : e[2] + 8]
            else:
                raise ValueError()

        return [g(e) for e in el]

    _endtypes = ["TD", "TD", "DT", "TD", "TD", "DT"]
    _orient = ("5", "5")
    _endlocs = [
        (0, 0, 5),
        (3, 0, 5),
        (5, -5, None),
        (5, 0, 5),
        (2, 0, 5),
        (0, -5, None),
    ]


class tile_daoe_doublehoriz_35up_1h2i(tile_daoe_doublehoriz_35up):
    edotparen = (
        "5.16(7(4.7)+8)16[16{8)+5(29(16]16}8(+5.29)16[16{8)5)+8(16]16}8(+5.16)5."
    )
    _endtypes = ["hairpin", "hairpin", "TD", "DT", "TD", "TD"]

    def _seqdiagseqstrings(self, tile):
        s = tile.strands
        return [
            s[0][:5] + "--" + s[0][5:13],
            (s[0][13:21] + "--" + s[0][21:26] + "-" + s[0][26:28] + "-" + s[0][28:30])[
                ::-1
            ],
            s[1][:8] + "--" + s[1][8:16],
            s[1][16:24],
            s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            s[2][:5] + "--" + s[2][5:13],
            (
                s[2][13:21]
                + "--"
                + s[2][21:26]
                + "--"
                + s[2][26:34]
                + "--"
                + s[2][34:42]
            )[::-1],
            (s[2][42:50])[::-1],
            s[2][50:58],
            s[2][58:66] + "--" + s[2][66:74],
            (s[3][0:5] + "--" + s[3][5:13])[::-1],
            s[3][13:21] + "--" + s[3][21:26] + "--" + s[3][26:34] + "--" + s[3][34:42],
            s[3][42:50],
            s[3][50:58][::-1],
            (s[3][58:66] + "--" + s[3][66:74] + "--" + s[3][74:79])[::-1],
            (s[4][0:8] + "--" + s[4][8:16])[::-1],
            (s[4][16:24])[::-1],
            s[4][24:32],
            s[4][32:40] + "--" + s[4][40:48],
            (s[5][0:5] + "--" + s[5][5:13])[::-1],
            s[5][13:21] + "--" + s[5][21:26],
            s[0][30:32] + "-" + s[0][32:34] + "-" + s[0][34:39],
        ]


class tile_daoe_doublehoriz_35up_4h5i(tile_daoe_doublehoriz_35up):
    edotparen = (
        "5.16(5.+8)16[16{8)+5.29(16]16}8(5(+5)29)16[16{8)+8(16]16}8(+5.16)7(4.7)"
    )
    _endtypes = ["DT", "TD", "TD", "hairpin", "blunt", "TD"]

    def _seqdiagseqstrings(self, tile):
        s = tile.strands
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16],
            s[1][16:24],
            s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            s[2][:5] + "--" + s[2][5:13],
            (
                s[2][13:21]
                + "--"
                + s[2][21:26]
                + "--"
                + s[2][26:34]
                + "--"
                + s[2][34:42]
            )[::-1],
            (s[2][42:50])[::-1],
            s[2][50:58],
            s[2][58:66] + "--" + s[2][66:74] + "--" + s[2][74:79],
            (s[3][0:5] + "--" + s[3][5:13])[::-1],
            s[3][13:21] + "--" + s[3][21:26] + "--" + s[3][26:34] + "--" + s[3][34:42],
            s[3][42:50],
            s[3][50:58][::-1],
            (s[3][58:66] + "--" + s[3][66:74])[::-1],
            (s[4][0:8] + "--" + s[4][8:16])[::-1],
            (s[4][16:24])[::-1],
            s[4][24:32],
            s[4][32:40] + "--" + s[4][40:48],
            (s[5][0:5] + "--" + s[5][5:13])[::-1],
            s[5][13:21] + "--" + s[5][21:26] + "-" + s[5][26:28] + "-" + s[5][28:30],
            (s[5][30:32] + "-" + s[5][32:34] + "-" + s[5][34:39])[::-1],
        ]


class tile_daoe_doublehoriz_35up_2h3h(tile_daoe_doublehoriz_35up):
    _endtypes = ["DT", "hairpin", "hairpin", "DT", "TD", "TD"]

    _endlocs = [
        (0, -5, None),
        (2, 0, 18),
        (5, 0, 18),
        (5, -5, None),
        (3, 0, 5),
        (0, 0, 5),
    ]

    edotparen = (
        "5.16(5.+8)16[16{8)+7(4.7)29(16]16}8(+5.29)16[16{8)+8(16]16}8(+7(4.23)5."
    )

    def _seqdiagseqstrings(self, tile):
        s = copy.copy(tile.strands)
        hp2 = s[2][0:13]
        s[2] = s[2][13:]
        hp3 = s[5][0:13]
        s[5] = s[5][13:]
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16],
            s[1][16:24],
            s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            (hp2[:5] + "-" + hp2[5:7] + "-" + hp2[7:9])[::-1],
            hp2[9:11] + "-" + hp2[11:13] + "-",
            s[2][:5] + "--" + s[2][5:13],
            (
                s[2][13:21]
                + "--"
                + s[2][21:26]
                + "--"
                + s[2][26:34]
                + "--"
                + s[2][34:42]
            )[::-1],
            (s[2][42:50])[::-1],
            s[2][50:58],
            s[2][58:66] + "--" + s[2][66:74],
            (s[3][0:5] + "--" + s[3][5:13])[::-1],
            s[3][13:21] + "--" + s[3][21:26] + "--" + s[3][26:34] + "--" + s[3][34:42],
            s[3][42:50],
            s[3][50:58][::-1],
            (s[3][58:66] + "--" + s[3][66:74])[::-1],
            (s[4][0:8] + "--" + s[4][8:16])[::-1],
            (s[4][16:24])[::-1],
            s[4][24:32],
            s[4][32:40] + "--" + s[4][40:48],
            hp3[:5] + "-" + hp3[5:7] + "-" + hp3[7:9],
            (hp3[9:11] + "-" + hp3[11:13] + "-")[::-1],
            (s[5][0:5] + "--" + s[5][5:13])[::-1],
            s[5][13:21] + "--" + s[5][21:26],
        ]


class tile_daoe_doublehoriz_53up_2h3h(tile_daoe_doublehoriz_53up):
    _endtypes = ["TD", "hairpin", "hairpin", "TD", "DT", "DT"]

    _endlocs = [
        (0, -5, None),
        (2, 0, 18),
        (5, 0, 18),
        (5, -5, None),
        (3, 0, 5),
        (0, 0, 5),
    ]

    edotparen = (
        "5.16(5.+8)16[16{8)+7(4.7)29(16]16}8(+5.29)16[16{8)+8(16]16}8(+7(4.23)5."
    )

    def _seqdiagseqstrings(self, tile):
        s = copy.copy(tile.strands)
        hp2 = s[2][0:13]
        s[2] = s[2][13:]
        hp3 = s[5][0:13]
        s[5] = s[5][13:]
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16],
            s[1][16:24],
            s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            (hp2[:5] + "-" + hp2[5:7] + "-" + hp2[7:9])[::-1],
            hp2[9:11] + "-" + hp2[11:13] + "-",
            s[2][:5] + "--" + s[2][5:13],
            (
                s[2][13:21]
                + "--"
                + s[2][21:26]
                + "--"
                + s[2][26:34]
                + "--"
                + s[2][34:42]
            )[::-1],
            (s[2][42:50])[::-1],
            s[2][50:58],
            s[2][58:66] + "--" + s[2][66:74],
            (s[3][0:5] + "--" + s[3][5:13])[::-1],
            s[3][13:21] + "--" + s[3][21:26] + "--" + s[3][26:34] + "--" + s[3][34:42],
            s[3][42:50],
            s[3][50:58][::-1],
            (s[3][58:66] + "--" + s[3][66:74])[::-1],
            (s[4][0:8] + "--" + s[4][8:16])[::-1],
            (s[4][16:24])[::-1],
            s[4][24:32],
            s[4][32:40] + "--" + s[4][40:48],
            hp3[:5] + "-" + hp3[5:7] + "-" + hp3[7:9],
            (hp3[9:11] + "-" + hp3[11:13] + "-")[::-1],
            (s[5][0:5] + "--" + s[5][5:13])[::-1],
            s[5][13:21] + "--" + s[5][21:26],
        ]


class tile_daoe_doublevert_35up_4h5h(tile_daoe_doublevert_35up):
    _endtypes = ["DT", "DT", "TD", "hairpin", "hairpin", "TD"]

    _endlocs = [
        (0, -5, None),
        (3, -5, None),
        (5, 0, 5),
        (5, -18, None),
        (2, -18, None),
        (0, 0, 5),
    ]

    edotparen = "5.16(5.+8)16[16{8)+8(16]16}8(5(16(7(4.7)+8)16[16{8)5)16)5.+8(16]16}8(+5.16)7(4.7)"

    def _seqdiagseqstrings(self, tile):
        s = copy.copy(tile.strands)
        return [
            s[0][:5] + "--" + s[0][5:13],
            s[0][:-6:-1] + "--" + s[0][-6:-14:-1],
            s[1][:8] + "--" + s[1][8:16],
            s[1][16:24],
            s[1][24:32][::-1],
            (s[1][32:40] + "--" + s[1][40:48])[::-1],
            (s[2][0:8] + "--" + s[2][8:16])[::-1],
            s[2][16:24][::-1],
            s[2][24:32],
            s[2][32:40] + "--" + s[2][40:48] + "--" + s[2][48:53] + "--" + s[2][53:61],
            (s[2][61:69] + "--" + s[2][69:74])[::-1],
            ("-" + s[2][74:76] + "-" + s[2][76:78])[::-1],
            s[2][78:80] + "-" + s[2][80:82] + "-" + s[2][82:87],
            s[3][0:8] + "--" + s[3][8:16],
            s[3][16:24],
            s[3][24:32][::-1],
            (
                s[3][32:40]
                + "--"
                + s[3][40:48]
                + "--"
                + s[3][48:53]
                + "--"
                + s[3][53:61]
            )[::-1],
            (s[3][61:69] + "--" + s[3][69:74]),
            (s[4][0:8] + "--" + s[4][8:16])[::-1],
            (s[4][16:24])[::-1],
            s[4][24:32],
            s[4][32:40] + "--" + s[4][40:48],
            (s[5][0:5] + "--" + s[5][5:13])[::-1],
            s[5][13:21] + "--" + s[5][21:26],
            "-" + s[5][26:28] + "-" + s[5][28:30],  # hp5
            (s[5][30:32] + "-" + s[5][32:34] + "-" + s[5][34:39])[::-1],  # hp5
        ]


class tile_daoe_doublevert_53up_4h5h(tile_daoe_doublevert_53up):
    _endtypes = ["TD", "TD", "DT", "hairpin", "hairpin", "DT"]


tilestructures = {
    "tile_daoe_5up": tile_daoe_5up,
    "tile_daoe_3up": tile_daoe_3up,
    "tile_daoe_5up_1h": tile_daoe_5up_1h,
    "tile_daoe_5up_2h": tile_daoe_5up_2h,
    "tile_daoe_5up_3h": tile_daoe_5up_3h,
    "tile_daoe_5up_3b": tile_daoe_5up_3b,
    "tile_daoe_5up_4h": tile_daoe_5up_4h,
    "tile_daoe_3up_1h": tile_daoe_3up_1h,
    "tile_daoe_3up_2h": tile_daoe_3up_2h,
    "tile_daoe_3up_3h": tile_daoe_3up_3h,
    "tile_daoe_3up_4h": tile_daoe_3up_4h,
    "tile_daoe_doublehoriz_53up": tile_daoe_doublehoriz_53up,
    "tile_daoe_doublehoriz_53up_2h3h": tile_daoe_doublehoriz_53up_2h3h,
    "tile_daoe_doublehoriz_35up": tile_daoe_doublehoriz_35up,
    "tile_daoe_doublehoriz_35up_2h3h": tile_daoe_doublehoriz_35up_2h3h,
    "tile_daoe_doublehoriz_35up_1h2i": tile_daoe_doublehoriz_35up_1h2i,
    "tile_daoe_doublehoriz_35up_4h5i": tile_daoe_doublehoriz_35up_4h5i,
    "tile_daoe_doublehoriz_35up_4h5b": tile_daoe_doublehoriz_35up_4h5i,
    "tile_daoe_doublevert_35up_4h5h": tile_daoe_doublevert_35up_4h5h,
    "tile_daoe_doublevert_53up_4h5h": tile_daoe_doublevert_53up_4h5h,
}


def getstructure(name, extra=None):
    if not name:
        return None
    if extra:
        name += "_" + extra
    return tilestructures[name]()


def compname(endname):
    if endname[-1] == "/":
        return endname[:-1]
    else:
        return endname + "/"


def order_pepper_strands(strandlist):
    # We're going to assume, for now, that they're already ordered by Pepper.
    # We'll see.
    return [strandseq for (strandname, strandseq) in strandlist]


def gettile(tset, tname):
    foundtiles = [x for x in tset["tiles"] if x["name"] == tname]
    assert len(foundtiles) == 1
    return foundtiles[0]
