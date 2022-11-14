from __future__ import annotations
from typing import (
    TYPE_CHECKING,
    Callable,
    Dict,
    Iterable,
    Literal,
    Mapping,
    Optional,
    Tuple,
    cast,
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
        if strand.name in new_ts.tiles:
            tile: BaseSSTile = new_ts.tiles[strand.name]
            assert isinstance(tile, BaseSSTile)

            assert str(tile.sequence) == strand.sequence(
                delimiter="-"
            )  # FIXME: should not need str here.

    return new_ts


def update_nuad_design(
    tileset: "TileSet",
    old_design: "nc.Design",
    groups: Literal["structure"]
    | Mapping[str, str]
    | Callable[[Tile], str] = "structure",
) -> "nc.Design":
    """
    From an Alhambra tileset and an existing Nuad Design. IN-PLACE
    """
    try:
        import nuad.constraints as nc
    except ImportError:
        raise ImportError("nuad must be installed for this function.")

    # Determine group-setting function
    if groups == "structure":
        get_group = lambda tile: tile.structure
    elif isinstance(groups, Mapping):
        get_group = lambda tile: cast(Mapping, groups)[tile]
    else:
        get_group = groups

    def filter_strand(s: "nc.Strand"):
        if s.name not in tileset.tiles.asdict():
            log.info(f"Removing tile {s.name} from existing design.")
            return False
        return True

    old_design.compute_derived_fields()

    old_design.strands = [s for s in old_design.strands if filter_strand(s)]

    strands_by_name = {s.name: s for s in old_design.strands}

    pools = {(pool.name, pool.length): pool for pool in old_design.domain_pools()}
    new_ncdomains: dict[str, "nc.Domain"] = dict()

    for tile in tileset.tiles:
        if tile.ident() in strands_by_name:
            continue
        log.info(f"Adding new tile {tile.ident()} to existing design.")

        strand_group = get_group(tile)
        td: list[SSGlue] = tile.domains

        ncdomains = []
        for domain in td:
            try:
                ncdomains.append(old_design.domains_by_name[domain.basename()])
            except KeyError:
                try:
                    ncdomains.append(new_ncdomains[domain.basename()])
                except KeyError:
                    log.info(f"Adding domain {domain.basename()}.")
                    ncdomain = nc.Domain(domain.basename())

                    if domain.sequence.is_definite:
                        ncdomain.set_fixed_sequence(domain.sequence.base_str)
                    else:
                        assert domain.sequence.is_null
                        try:
                            pool = pools[
                                (f"SSGlue{domain.dna_length}", domain.dna_length)
                            ]  # FIXME: determine type
                        except KeyError:
                            pool = nc.DomainPool(
                                f"SSGlue{domain.dna_length}", domain.dna_length
                            )  # FIXME: determine type
                            pools[
                                (f"SSGlue{domain.dna_length}", domain.dna_length)
                            ] = pool  # FIXME: determine type
                    ncdomain.pool = pool
                    new_ncdomains[ncdomain.name] = ncdomain
                    ncdomains.append(ncdomain)

        old_design.strands.append(
            nc.Strand(
                domains=ncdomains,
                starred_domain_indices=[i for i, d in enumerate(td) if d.is_complement],
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

    old_design.compute_derived_fields()

    for ncdomain in old_design.domains:
        if not ncdomain.has_pool():
            log.warn(f"Domain {ncdomain.name} has no pool.")
        else:
            domain = non_complementary_domains[ncdomain.name]
            if ncdomain.pool.length != domain.dna_length:
                try:
                    ncdomain._pool = pools[
                        (f"SSGlue{domain.dna_length}", domain.dna_length)
                    ]
                except KeyError:
                    ncdomain._pool = nc.DomainPool(
                        f"SSGlue{domain.dna_length}", domain.dna_length
                    )  # FIXME: determine type
                    pools[
                        (f"SSGlue{domain.dna_length}", domain.dna_length)
                    ] = ncdomain._pool  # FIXME: determine type

    return old_design


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
        get_group = lambda tile: cast(Mapping, groups)[tile]
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
                pool = pools[
                    (f"SSGlue{domain.dna_length}", domain.dna_length)
                ]  # FIXME: determine type
            except KeyError:
                pool = nc.DomainPool(
                    "SSGlue", domain.dna_length
                )  # FIXME: determine type
                pools[
                    (f"SSGlue{domain.dna_length}", domain.dna_length)
                ] = pool  # FIXME: determine type

            ncdomain.pool = pool
            ncdomains.append(ncdomain)

    # FIXME: do we need to do something with the pools?

    des.store_domain_pools()

    return des


def group_strand_pairs_by_groups_and_complementary_domains(
    design: "nc.Design", strands: Optional["Optional[Iterable[nc.Strand]]"] = None
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
