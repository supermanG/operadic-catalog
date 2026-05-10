"""Catalog of operad generators per crystallographic point group.

For each of the 32 point groups, lists the generating operations (one
per irrep channel) and identifies relations: which generators are
obtained from others by composition or restriction.
"""

from __future__ import annotations

from dataclasses import dataclass

from operadic_catalog.operad import (
    PairingOperad,
    Operation,
    Color,
    is_subgroup,
    subgroups_of,
)
from operadic_catalog.sigrist_ueda import (
    POINT_GROUP_CHANNELS,
    ALL_POINT_GROUPS,
    CHANNEL_CYCLIC_ORDER,
)


@dataclass
class GeneratorRelation:
    """A relation between generators: target = restriction of source."""
    source_pg: str
    source_channel: str
    target_pg: str
    target_channel: str
    relation_type: str  # 'restriction', 'inherited', 'independent'


def catalog_generators(operad: PairingOperad | None = None) -> dict[str, list[Operation]]:
    """Return all generators organized by point group."""
    if operad is None:
        operad = PairingOperad()
    return {pg: operad.generators(pg) for pg in sorted(ALL_POINT_GROUPS)}


def independent_generators(operad: PairingOperad | None = None) -> dict[str, list[Operation]]:
    """Identify which generators per point group are independent
    (not obtainable by restriction from a supergroup).

    A generator (G, ch) is independent if there is no supergroup H > G
    such that (H, ch) restricts to (G, ch).
    """
    if operad is None:
        operad = PairingOperad()

    independent: dict[str, list[Operation]] = {}
    for pg in sorted(ALL_POINT_GROUPS):
        pg_gens = operad.generators(pg)
        ind = []
        for gen in pg_gens:
            is_inherited = False
            for super_pg in ALL_POINT_GROUPS:
                if super_pg == pg:
                    continue
                if not is_subgroup(pg, super_pg):
                    continue
                super_channels = POINT_GROUP_CHANNELS.get(super_pg, frozenset())
                if gen.channel in super_channels:
                    is_inherited = True
                    break
            if not is_inherited:
                ind.append(gen)
        independent[pg] = ind
    return independent


def generator_relations(operad: PairingOperad | None = None) -> list[GeneratorRelation]:
    """Enumerate all restriction relations between generators.

    For each pair (G, H) with H < G, and each channel ch present in
    both G and H, records that (H, ch) is the restriction of (G, ch).
    """
    if operad is None:
        operad = PairingOperad()

    relations: list[GeneratorRelation] = []
    for pg_g in sorted(ALL_POINT_GROUPS):
        channels_g = POINT_GROUP_CHANNELS.get(pg_g, frozenset())
        subs = subgroups_of(pg_g) - {pg_g}
        for pg_h in sorted(subs):
            channels_h = POINT_GROUP_CHANNELS.get(pg_h, frozenset())
            shared = channels_g & channels_h
            for ch in sorted(shared):
                relations.append(GeneratorRelation(
                    source_pg=pg_g,
                    source_channel=ch,
                    target_pg=pg_h,
                    target_channel=ch,
                    relation_type="restriction",
                ))
            lost = channels_g - channels_h - {"s"}
            for ch in sorted(lost):
                relations.append(GeneratorRelation(
                    source_pg=pg_g,
                    source_channel=ch,
                    target_pg=pg_h,
                    target_channel="(vanishes)",
                    relation_type="lost_under_restriction",
                ))
    return relations


def generator_summary() -> dict:
    """Summary of the generator catalog."""
    operad = PairingOperad()
    all_gens = catalog_generators(operad)
    indep = independent_generators(operad)
    rels = generator_relations(operad)

    total_gens = sum(len(v) for v in all_gens.values())
    total_indep = sum(len(v) for v in indep.values())
    restriction_rels = sum(1 for r in rels if r.relation_type == "restriction")
    lost_rels = sum(1 for r in rels if r.relation_type == "lost_under_restriction")

    return {
        "total_generators": total_gens,
        "total_independent": total_indep,
        "total_inherited": total_gens - total_indep,
        "restriction_relations": restriction_rels,
        "lost_under_restriction": lost_rels,
        "generators_per_pg": {pg: len(ops) for pg, ops in all_gens.items()},
        "independent_per_pg": {pg: len(ops) for pg, ops in indep.items()},
    }
