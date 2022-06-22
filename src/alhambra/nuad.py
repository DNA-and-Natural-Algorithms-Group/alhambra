from typing import (
    TYPE_CHECKING,
    Callable,
    Dict,
    Iterable,
    Literal,
    Mapping,
    Optional,
    Tuple,
)
import itertools
import logging

from alhambra.glues import GlueList, SSGlue
from alhambra.tiles import BaseSSTile
from alhambra.tiles import Tile, BaseSSTSingle

if TYPE_CHECKING:
    from alhambra.tilesets import TileSet
    import nuad.constraints as nc

log = logging.getLogger(__name__)


def load_nuad_design(
    tileset: "TileSet",
    design: "nc.Design",
    update_tiles: bool = True,
    inplace: bool = False,
) -> "TileSet":
    try:
        import nuad.constraints as nc
    except ImportError:
        raise ImportError("nuad must be installed for this function.")

    if not inplace:
        new_ts = tileset.copy()
    else:
        new_ts = tileset

    for domain in design.domains:
        new_ts.glues.add(
            SSGlue(domain.name, sequence=domain.sequence())
        )  # FIXME: determine domain type from pool

    if update_tiles:
        for t in new_ts.tiles:
            t.update_glues(new_ts.glues)

    for strand in design.strands:
        # FIXME: handle multi-strand tiles, etc
        # FIXME: handle non-update case
        tile: BaseSSTile = new_ts.tiles[strand.name]
        assert isinstance(tile, BaseSSTile)

        assert str(tile.sequence) == strand.sequence(
            delimiter="-"
        )  # FIXME: should not need str here.

    return new_ts


def tileset_to_nuad_design(
    tileset: "TileSet",
    groups: Literal["structure"]
    | Mapping[str, str]
    | Callable[[Tile], str] = "structure",
) -> "nc.Design":
    """
    From an Alhambra tileset, generate a Nuad Design with all of its domains.

    groups:
        Methods for setting the group of each strand.
    """
    try:
        import nuad.constraints as nc
    except ImportError:
        raise ImportError("nuad must be installed for this function.")

    # Determine group-setting function
    if groups == "structure":
        get_group = lambda tile: tile.structure
    elif isinstance(groups, Mapping):
        get_group = lambda tile: groups[tile]
    else:
        get_group = groups

    # Convert all tiles to strands.
    strands = []
    for tile in tileset.tiles:
        strand_group = get_group(tile)
        strands.append(
            nc.Strand(
                [domain.ident() for domain in tile.domains],
                name=tile.name,
                group=strand_group,
            )
        )

    # Nuad needs only non-complementary (ie, non-starred) domains.  Alhambra may contain
    # complementary domains that don't have a non-complementary counterpart.  So we'll need
    # to compile these.  FIXME: move out to other function?
    non_complementary_domains: GlueList[SSGlue] = GlueList()

    for domain in tileset.alldomains:
        # We only want domains that actually have sequences, and can be designed.
        # FIXME: should have more general sequence-containing class
        if not isinstance(domain, SSGlue):
            log.warning(f"Not adding domain {domain.name} to Nuad design: {domain}")
            continue

        if domain.is_complement:
            non_complementary_domains.add(domain.complement)
        else:
            non_complementary_domains.add(domain)

    pools: dict[tuple[str, int], nc.DomainPool] = {}

    des = nc.Design(strands=strands)

    ncdomains = []
    for ncdomain in des.domains:
        domain = non_complementary_domains[ncdomain.name]

        if domain.sequence.is_definite:
            ncdomain.set_fixed_sequence(domain.sequence.base_str)
        else:
            assert domain.sequence.is_null
            try:
                pool = pools[("SSGlue", domain.dna_length)]  # FIXME: determine type
            except KeyError:
                pool = nc.DomainPool(
                    "SSGlue", domain.dna_length
                )  # FIXME: determine type
                pools[("SSGlue", domain.dna_length)] = pool  # FIXME: determine type

            ncdomain.pool = pool
            ncdomains.append(ncdomain)

    # FIXME: do we need to do something with the pools?

    des.store_domain_pools()

    return des


def group_strand_pairs_by_groups_and_complementary_domains(
    design: "nc.Design", strands: "Optional[Iterable[nc.Strand]]" = None
) -> "Dict[Tuple[str, str, int], list[Tuple[nc.Strand, nc.Strand]]]":
    """
    Group pairs of strands by their groups (sorted) and number of complementary domains.
    """
    try:
        import nuad.constraints as nc
    except ImportError:
        raise ImportError("nuad must be installed for this function.")

    pairs: Dict[Tuple[str, str, int], list[Tuple[nc.Strand, nc.Strand]]] = {}

    if strands is None:
        strands = design.strands

    for strand1, strand2 in itertools.combinations_with_replacement(strands, 2):
        domains1_unstarred = strand1.unstarred_domains_set()
        domains2_unstarred = strand2.unstarred_domains_set()
        domains1_starred = strand1.starred_domains_set()
        domains2_starred = strand2.starred_domains_set()

        comp_domains = (domains1_unstarred & domains2_starred) | (
            domains2_unstarred & domains1_starred
        )
        comp_domain_names = [domain.name for domain in comp_domains]
        num_comp_domains = len(comp_domain_names)

        g1sorted, g2sorted = tuple(sorted([strand1.group, strand2.group]))

        key = (g1sorted, g2sorted, num_comp_domains)

        if key not in pairs:
            pairs[key] = []

        pairs[key].append((strand1, strand2))

    return pairs
