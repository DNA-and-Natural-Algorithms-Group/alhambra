T = TypeVar("T")


class NamedList(CommentedSeq, Generic[T]):
    """A class for a list of dicts, where some dicts have a 'name' item
    which should be unique, but others might not.  Indexing works with
    either number or name.  Note that updating dicts may make the list
    inconsistent.

    """

    def __init__(self, x=[]):
        CommentedSeq.__init__(self, x)

    def __getitem__(self, i: Union[int, str]) -> T:
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

    def __setitem__(self, i: Union[str, int], v: T):
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

    def __delitem__(self, i: Union[str, int]):
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
