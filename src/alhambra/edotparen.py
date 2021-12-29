import collections
import re
import string

from . import seq as sq
from .seq import _VALID_NTS, revcomp

edp_closetoopen = {x: y for x, y in zip(string.ascii_lowercase, string.ascii_uppercase)}
edp_closetoopen.update({")": "(", "]": "[", "}": "{"})


def check_edotparen_consistency(expr):
    expr = expand_compact_edotparen(expr)
    expr = re.sub(r"\s+", "", expr)
    counts = collections.Counter()
    strand = 0
    strandloc = 0
    for s in expr:
        if s in edp_closetoopen.values():
            counts[s] += 1
        elif s in edp_closetoopen.keys():
            try:
                counts[edp_closetoopen[s]] -= 1
            except KeyError:
                raise ValueError("Opening not found", s, strand, strandloc)
        elif s == ".":
            pass
        elif s == "+":
            strand += 1
            strandloc = 0
            continue
        else:
            raise ValueError("Unknown char", s, strand, strandloc)
        strandloc += 1
    if max(counts.values()) > 0:
        raise ValueError(counts)


def check_edotparen_sequence(edotparen, sequence):
    expr = re.sub(r"\s+", "", expand_compact_edotparen(edotparen))
    seq = re.sub(r"\s+", "", sequence).lower()
    if len(expr) != len(seq):
        raise ValueError("Unequal lengths")
    stacks = {}
    strand = 0
    strandloc = 0
    for s, v in zip(expr, seq):
        if s in edp_closetoopen.values():
            if s not in stacks.keys():
                stacks[s] = []
            stacks[s].append(v)
        elif s in edp_closetoopen.keys():
            ss = edp_closetoopen[s]
            if ss not in stacks.keys():
                raise ValueError("Opening not found", s, strand, strandloc)
            vv = stacks[ss].pop()
            try:
                sq.merge(v, revcomp(vv))
            except sq.MergeConflictError:
                raise ValueError(
                    "{} != WC({}) at strand {} loc {} (both from 0)".format(
                        v, vv, strand, strandloc
                    ),
                    v,
                    vv,
                    strand,
                    strandloc,
                ) from None
        elif s == ".":
            assert v in _VALID_NTS
        elif s == "+":
            assert v == "+"
            strand += 1
            strandloc = 0
            continue
        else:
            raise ValueError("Unknown char", s, strand, strandloc)
        strandloc += 1
    if max(len(stack) for stack in stacks.values()) > 0:
        raise ValueError(stacks)


def expand_compact_edotparen(expr):
    return re.sub(
        r"(\d+)([\[\]\(\)\{\}A-Za-z\.])", lambda m: int(m.group(1)) * m.group(2), expr
    )


def prettify_edotparen(expr):
    # This is evil:
    return re.sub(
        r"(([\[\]\(\)\{\}A-Za-z\.])\2+)",
        lambda m: "{}{}".format(len(m.group(1)), m.group(2)),
        expr,
    )
