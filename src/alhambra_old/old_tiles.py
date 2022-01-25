"""A module"""
import copy

from ruamel.yaml.comments import CommentedMap

from .ends import EndList
from .tilestructures import TileStructure, getstructure
from .util import NamedList


class Tile(dict):
    """A class representing one tile."""

    def __init__(self, val={}):
        dict.__init__(self, val)

        if "type" in self.keys():
            self["structure"] = self["type"]
            del self["type"]

        if "structure" in self.keys():
            self.structure = getstructure(
                self.get("structure", None), extra=self.get("extra", None)
            )
        if self.get("extra"):
            del self["extra"]

    def structure():
        doc = """The structure of the tile, set either as a string (to be looked up
with getstructure) or a TileStructure.  Get always returns the TileStructure."""

        def fget(self):
            return self._structure

        def fset(self, value):
            if value is None:
                del self["structure"]
            elif isinstance(value, TileStructure):
                self._structure = value
                self["structure"] = value.name
            else:
                self._structure = getstructure(value)
                self["structure"] = self._structure.name

        def fdel(self):
            self._structure = None
            del self["structure"]

        return locals()

    structure = property(**structure())

    def copy(self):
        return copy.deepcopy(self)

    @property
    def rotations(self):
        rl = self.structure.rotations
        tl = TileList()
        for ri, (structure, endorder) in enumerate(rl):
            t = copy.copy(self)
            if "input" in t.keys():
                del t["input"]
            t["ends"] = [t.ends[i] for i in endorder]
            t["rotation"] = ri
            t.structure = structure()
            t.structure._endtypes = [self.structure._endtypes[i] for i in endorder]
            tl.append(t)
        return tl

    @property
    def is_fake(self):
        """Is the tile fake?

        Returns
        -------
        bool
            True if the tile is fake, False otherwise.
        """
        return "fake" in self.keys()

    def named_rotations(self):
        rl = self.structure.rotations
        tl = TileList()
        for ri, (structure, endorder) in enumerate(rl):
            t = copy.copy(self)
            if "input" in t.keys():
                del t["input"]
            t["ends"] = [t.ends[i] for i in endorder]
            t["rotation"] = ri
            t["name"] += "_rot_{}".format(ri)
            t.structure = structure()
            t.structure._endtypes = [self.structure._endtypes[i] for i in endorder]
            tl.append(t)
        return tl

    def strands():
        doc = """Doc string"""

        def fget(self):
            return self.get("fullseqs", None)

        def fset(self, value):
            # Make sure there are the right number of ends.
            # FIXME: do a better job checking, raise a better error
            # if self.structure:
            #    assert len(self.structure._endlocs) == len(value)
            self["fullseqs"] = value

        def fdel(self):
            del self["fullseqs"]

        return locals()

    strands = property(**strands())

    def check_strands(self):
        self.structure.check_strands(self)

    def sequence_diagram(self):
        return self.structure.sequence_diagram(self)

    def abstract_diagram(self, tileset=None):
        return self.structure.abstract_diagram(self, tileset)

    def name():
        doc = """Doc string"""

        def fget(self):
            return self["name"]

        def fset(self, value):
            self["name"] = value

        return locals()

    name = property(**name())

    def ends():
        doc = """Doc string"""

        def fget(self):
            return self["ends"]

        def fset(self, value):
            # Make sure there are the right number of ends.
            # FIXME: do a better job checking, raise a better error
            if self.structure:
                assert len(self.structure._endlocs) == len(value)
            self["ends"] = value

        def fdel(self):
            del self["ends"]

        return locals()

    ends = property(**ends())

    @property
    def endlist(self):
        return self.structure.get_endlist(self)

    @property
    def orderableseqs(self):
        return self.structure.orderableseqs(self)

    def check_consistent(self):
        if self.structure:
            self.structure.check_consistent()
            if self.ends:
                assert self.structure.numends == len(self.ends)
            if self.strands:
                self.check_strands()

    def check_sequences(self):
        self.structure.check_strands(self.strands)

    def __deepcopy__(self, memo):
        c = self.__class__()
        memo[id(self)] = c
        for k, v in self.items():
            c[k] = copy.deepcopy(v)
        if self.structure:
            c._structure = copy.deepcopy(self._structure)
        return c


class TileList(NamedList[Tile]):
    """A list of `Tile` instances, taking into account tile names."""

    def __init__(self, val=[]):
        NamedList.__init__(self, val)
        for i, tile in enumerate(self):
            self[i] = Tile(tile)

    def check_consistent(self):
        NamedList.check_consistent(self)
        for tile in self:
            tile.check_consistent()

    def glues_from_tiles(self, fail_immediate: bool = True) -> EndList:
        """Extract sticky ends from the list of tiles.

        Parameters
        ----------

        fail_immediate : bool
            (default True) if True, immediately fail on a failure,
            with ValueError( tilename, exception ) from exception  If False, collect
            exceptions, then raise ValueError( "message", [(tilename, exception)] ,
            output ).

        """
        endlist = EndList()
        errors = []

        for tile in self:
            try:
                endlist.update(tile.endlist)
            except BaseException as e:
                if fail_immediate:
                    raise ValueError(tile.name) from e
                else:
                    errors.append((tile.name, e))

        if errors:
            raise ValueError("End list generation failed on:", errors, endlist)

        return endlist


from ruamel.yaml.representer import RoundTripRepresenter

RoundTripRepresenter.add_representer(TileList, RoundTripRepresenter.represent_list)
RoundTripRepresenter.add_representer(Tile, RoundTripRepresenter.represent_dict)
