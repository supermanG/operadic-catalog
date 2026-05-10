"""Symmetry-reduction morphisms in the pairing operad.

For pairs (G, H) with H a subgroup of G, the inclusion H -> G induces
a restriction functor on operadic algebras. This module tabulates all
such morphisms among the 32 crystallographic point groups and builds
the morphism graph.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import networkx as nx

from operadic_catalog.operad import (
    PairingOperad,
    Operation,
    is_subgroup,
    subgroups_of,
    _MAXIMAL_SUBGROUPS,
)
from operadic_catalog.sigrist_ueda import (
    POINT_GROUP_CHANNELS,
    ALL_POINT_GROUPS,
)
from operadic_catalog.algebra import OperadAlgebra


@dataclass(frozen=True)
class RestrictionMorphism:
    """A morphism in the operad induced by subgroup inclusion H -> G."""
    source_pg: str
    target_pg: str
    preserved_channels: frozenset[str]
    lost_channels: frozenset[str]
    gained_channels: frozenset[str]

    @property
    def channel_retention_ratio(self) -> float:
        """Fraction of source channels preserved under restriction."""
        total = len(self.preserved_channels) + len(self.lost_channels)
        if total == 0:
            return 1.0
        return len(self.preserved_channels) / total

    def __repr__(self) -> str:
        return (
            f"Morph({self.source_pg} -> {self.target_pg}: "
            f"keep={set(self.preserved_channels)}, "
            f"lose={set(self.lost_channels)})"
        )


def compute_morphism(source_pg: str, target_pg: str) -> RestrictionMorphism | None:
    """Compute the restriction morphism from source to target.

    Returns None if target is not a subgroup of source.
    """
    if not is_subgroup(target_pg, source_pg):
        return None
    src_ch = POINT_GROUP_CHANNELS.get(source_pg, frozenset())
    tgt_ch = POINT_GROUP_CHANNELS.get(target_pg, frozenset())

    preserved = src_ch & tgt_ch
    lost = src_ch - tgt_ch
    gained = tgt_ch - src_ch

    return RestrictionMorphism(
        source_pg=source_pg,
        target_pg=target_pg,
        preserved_channels=frozenset(preserved),
        lost_channels=frozenset(lost),
        gained_channels=frozenset(gained),
    )


def all_morphisms(maximal_only: bool = False) -> list[RestrictionMorphism]:
    """Enumerate all restriction morphisms among the 32 point groups.

    If maximal_only=True, only returns morphisms between maximal
    subgroup pairs (the Hasse diagram edges).
    """
    morphisms: list[RestrictionMorphism] = []
    for pg_g in sorted(ALL_POINT_GROUPS):
        if maximal_only:
            targets = _MAXIMAL_SUBGROUPS.get(pg_g, set())
        else:
            targets = subgroups_of(pg_g) - {pg_g}
        for pg_h in sorted(targets):
            m = compute_morphism(pg_g, pg_h)
            if m is not None:
                morphisms.append(m)
    return morphisms


def build_morphism_graph(maximal_only: bool = True) -> nx.DiGraph:
    """Build the directed graph of restriction morphisms.

    Nodes: 32 crystallographic point groups.
    Edges: G -> H for each maximal subgroup relation (if maximal_only)
           or all subgroup relations.
    Edge attributes: preserved/lost/gained channels, retention ratio.
    """
    G = nx.DiGraph()

    for pg in sorted(ALL_POINT_GROUPS):
        channels = POINT_GROUP_CHANNELS.get(pg, frozenset())
        G.add_node(pg, n_channels=len(channels), channels=set(channels))

    for m in all_morphisms(maximal_only=maximal_only):
        G.add_edge(
            m.source_pg,
            m.target_pg,
            preserved=set(m.preserved_channels),
            lost=set(m.lost_channels),
            gained=set(m.gained_channels),
            retention=m.channel_retention_ratio,
        )

    return G


def restriction_path(source_pg: str, target_pg: str) -> list[str] | None:
    """Find the shortest restriction path from source to target.

    Returns the sequence of point groups along the path, or None
    if target is not a subgroup of source.
    """
    if not is_subgroup(target_pg, source_pg):
        return None
    G = build_morphism_graph(maximal_only=True)
    try:
        return nx.shortest_path(G, source_pg, target_pg)
    except nx.NetworkXNoPath:
        return None


def restrict_algebra(
    algebra: OperadAlgebra,
    target_pg: str,
    operad: PairingOperad | None = None,
) -> OperadAlgebra | None:
    """Restrict an O-algebra to a subgroup.

    The restricted algebra has the same cyclic weights but only the
    channels that survive restriction.
    """
    if not is_subgroup(target_pg, algebra.point_group):
        return None

    morphism = compute_morphism(algebra.point_group, target_pg)
    if morphism is None:
        return None

    new_channels = POINT_GROUP_CHANNELS.get(target_pg, frozenset())

    return OperadAlgebra(
        material_id=f"{algebra.material_id}|{target_pg}",
        formula=algebra.formula,
        point_group=target_pg,
        channels=new_channels,
        weights=dict(algebra.weights),
        metadata=dict(algebra.metadata),
    )


def morphism_summary() -> dict[str, Any]:
    """Summary statistics of the morphism structure."""
    maximal = all_morphisms(maximal_only=True)
    total = all_morphisms(maximal_only=False)
    G = build_morphism_graph(maximal_only=True)

    retention_values = [m.channel_retention_ratio for m in maximal]
    avg_retention = sum(retention_values) / len(retention_values) if retention_values else 0.0

    sources = set()
    sinks = set()
    for pg in ALL_POINT_GROUPS:
        if G.in_degree(pg) == 0:
            sources.add(pg)
        if G.out_degree(pg) == 0:
            sinks.add(pg)

    return {
        "n_maximal_morphisms": len(maximal),
        "n_total_morphisms": len(total),
        "avg_channel_retention": avg_retention,
        "sources": sorted(sources),
        "sinks": sorted(sinks),
        "longest_chain": nx.dag_longest_path_length(G) if nx.is_directed_acyclic_graph(G) else -1,
    }
