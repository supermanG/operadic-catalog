"""Operadic equivalence classes of materials.

Two materials are operadically equivalent if they define isomorphic
algebras over the pairing operad. In practice this means:
  1. Same point group.
  2. Same set of allowed channels (automatic from #1).
  3. Same cyclic weight profile up to a tolerance.

The "moduli quotient" groups the material pool into equivalence classes.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np

from operadic_catalog.algebra import OperadAlgebra, material_to_operad_algebra
from operadic_catalog.sigrist_ueda import CHANNEL_CYCLIC_ORDER


@dataclass
class EquivalenceClass:
    """An operadic equivalence class of materials."""
    class_id: int
    point_group: str
    channels: frozenset[str]
    representative: OperadAlgebra
    members: list[OperadAlgebra]
    weight_centroid: dict[tuple[int, int], float]

    @property
    def size(self) -> int:
        return len(self.members)

    def __repr__(self) -> str:
        return (
            f"Class({self.class_id}: {self.point_group}, "
            f"|channels|={len(self.channels)}, "
            f"n_members={self.size})"
        )


def _weight_vector(alg: OperadAlgebra, n_values: tuple[int, ...] = (3, 4, 5, 6)) -> np.ndarray:
    """Extract a flat weight vector for distance computation."""
    components = []
    for n in n_values:
        for k in range(n):
            components.append(alg.weights.get((n, k), 0.0))
    return np.array(components, dtype=float)


def _weight_distance(a: OperadAlgebra, b: OperadAlgebra) -> float:
    """L2 distance between cyclic weight vectors."""
    va = _weight_vector(a)
    vb = _weight_vector(b)
    return float(np.linalg.norm(va - vb))


def operad_equivalence_class(a: OperadAlgebra, b: OperadAlgebra,
                              weight_tol: float = 0.05) -> bool:
    """Test whether two O-algebras are operadically equivalent.

    Strict equivalence: same point group and channel set.
    Approximate equivalence: additionally, cyclic weight profiles
    are within tolerance.
    """
    if a.point_group != b.point_group:
        return False
    if a.channels != b.channels:
        return False
    if weight_tol <= 0:
        return a.operadic_fingerprint() == b.operadic_fingerprint()
    return _weight_distance(a, b) <= weight_tol


def enumerate_classes(
    algebras: Sequence[OperadAlgebra],
    weight_tol: float = 0.05,
    coarse: bool = False,
) -> list[EquivalenceClass]:
    """Partition a collection of O-algebras into equivalence classes.

    If coarse=True, uses only (point_group, channels) for grouping
    (ignoring weight data). This gives the "Sigrist-Ueda classes",
    the coarsest operadic partition.

    If coarse=False (default), further subdivides by cyclic weight
    profile using weight_tol.
    """
    if coarse:
        return _enumerate_coarse(algebras)
    return _enumerate_fine(algebras, weight_tol)


def _enumerate_coarse(algebras: Sequence[OperadAlgebra]) -> list[EquivalenceClass]:
    """Coarse partition: group by (point_group, channels) only."""
    buckets: dict[tuple[str, tuple[str, ...]], list[OperadAlgebra]] = defaultdict(list)
    for alg in algebras:
        key = (alg.point_group, alg.channel_signature)
        buckets[key].append(alg)

    classes = []
    for idx, ((pg, ch_sig), members) in enumerate(sorted(buckets.items())):
        centroid = _compute_centroid(members)
        classes.append(EquivalenceClass(
            class_id=idx,
            point_group=pg,
            channels=frozenset(ch_sig),
            representative=members[0],
            members=members,
            weight_centroid=centroid,
        ))
    return classes


def _enumerate_fine(
    algebras: Sequence[OperadAlgebra],
    weight_tol: float,
) -> list[EquivalenceClass]:
    """Fine partition: group by (point_group, channels), then subdivide
    by cyclic weight distance."""
    coarse = _enumerate_coarse(algebras)
    fine_classes: list[EquivalenceClass] = []
    class_id = 0

    for coarse_cls in coarse:
        if len(coarse_cls.members) <= 1 or weight_tol <= 0:
            coarse_cls.class_id = class_id
            fine_classes.append(coarse_cls)
            class_id += 1
            continue

        sub_clusters: list[list[OperadAlgebra]] = []
        for alg in coarse_cls.members:
            placed = False
            for cluster in sub_clusters:
                if _weight_distance(alg, cluster[0]) <= weight_tol:
                    cluster.append(alg)
                    placed = True
                    break
            if not placed:
                sub_clusters.append([alg])

        for cluster in sub_clusters:
            centroid = _compute_centroid(cluster)
            fine_classes.append(EquivalenceClass(
                class_id=class_id,
                point_group=coarse_cls.point_group,
                channels=coarse_cls.channels,
                representative=cluster[0],
                members=cluster,
                weight_centroid=centroid,
            ))
            class_id += 1

    return fine_classes


def _compute_centroid(
    members: list[OperadAlgebra],
    n_values: tuple[int, ...] = (3, 4, 5, 6),
) -> dict[tuple[int, int], float]:
    """Compute the centroid of cyclic weight vectors."""
    if not members:
        return {}
    vecs = np.array([_weight_vector(m, n_values) for m in members])
    mean_vec = vecs.mean(axis=0)
    centroid: dict[tuple[int, int], float] = {}
    idx = 0
    for n in n_values:
        for k in range(n):
            centroid[(n, k)] = float(mean_vec[idx])
            idx += 1
    return centroid


def class_summary(classes: list[EquivalenceClass]) -> dict[str, Any]:
    """Summary statistics for a set of equivalence classes."""
    sizes = [c.size for c in classes]
    pg_distribution: dict[str, int] = defaultdict(int)
    for c in classes:
        pg_distribution[c.point_group] += 1

    return {
        "n_classes": len(classes),
        "total_materials": sum(sizes),
        "largest_class": max(sizes) if sizes else 0,
        "smallest_class": min(sizes) if sizes else 0,
        "mean_class_size": float(np.mean(sizes)) if sizes else 0.0,
        "singleton_classes": sum(1 for s in sizes if s == 1),
        "pg_distribution": dict(pg_distribution),
    }
