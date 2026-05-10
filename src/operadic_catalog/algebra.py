"""Material -> O-algebra functor.

Each material M defines an algebra over the pairing operad O:
  A_M : O -> Vect_k
mapping each operation (G, chi_n) to the k-vector space of compatible
SC order parameters on M.

In practice, the algebra is specified by:
  1. The point group G of the material.
  2. The set of allowed pairing channels (from Sigrist-Ueda).
  3. The cyclic operad weights co_n{n}_w{k} from v13 (the "coordinates"
     of the algebra).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from operadic_catalog.operad import PairingOperad, Operation, Color
from operadic_catalog.sigrist_ueda import (
    POINT_GROUP_CHANNELS,
    CHANNEL_CYCLIC_ORDER,
    allowed_channels,
    spg_to_point_group,
)


@dataclass
class OperadAlgebra:
    """An algebra over the pairing operad, representing one material.

    Attributes:
        material_id: identifier (e.g. JARVIS jid)
        formula: chemical formula
        point_group: Schoenflies symbol
        channels: set of allowed pairing channels
        weights: dict mapping (cyclic_order, irrep_index) to weight
        metadata: additional material properties
    """
    material_id: str
    formula: str
    point_group: str
    channels: frozenset[str]
    weights: dict[tuple[int, int], float] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def channel_signature(self) -> tuple[str, ...]:
        """Sorted tuple of channels: the algebraic "type" of this material."""
        return tuple(sorted(self.channels))

    @property
    def cyclic_profile(self) -> dict[int, list[float]]:
        """Weights organized by cyclic order n -> [w0, w1, ..., w_{n-1}]."""
        profile: dict[int, list[float]] = {}
        for (n, k), w in sorted(self.weights.items()):
            if n not in profile:
                profile[n] = [0.0] * n
            if k < n:
                profile[n][k] = w
        return profile

    @property
    def dominant_channel(self) -> str | None:
        """The channel with the largest non-trivial weight."""
        best_ch = None
        best_w = 0.0
        for ch in self.channels:
            if ch == "s":
                continue
            n = CHANNEL_CYCLIC_ORDER.get(ch, 1)
            if n <= 1:
                continue
            for k in range(1, n):
                w = self.weights.get((n, k), 0.0)
                if w > best_w:
                    best_w = w
                    best_ch = ch
        return best_ch

    def evaluate(self, op: Operation) -> float:
        """Evaluate the algebra on an operation: return the weight.

        For a generating operation (G, ch), returns the sum of cyclic
        weights at the channel's cyclic order.
        """
        if op.point_group != self.point_group:
            return 0.0
        base_ch = op.channel.split("+")[0]
        if base_ch not in self.channels:
            return 0.0
        n = CHANNEL_CYCLIC_ORDER.get(base_ch, 1)
        if n <= 1:
            return 1.0
        return sum(
            self.weights.get((n, k), 0.0) for k in range(1, n)
        )

    def operadic_fingerprint(self) -> tuple:
        """A hashable fingerprint for operadic equivalence comparison.

        Two materials with the same fingerprint have isomorphic
        O-algebra structures.
        """
        return (
            self.point_group,
            self.channel_signature,
            tuple(sorted(
                ((n, k), round(w, 6))
                for (n, k), w in self.weights.items()
                if w > 1e-8
            )),
        )

    def __repr__(self) -> str:
        dom = self.dominant_channel or "s"
        return (
            f"Algebra({self.formula} [{self.point_group}], "
            f"channels={set(self.channels)}, dominant={dom})"
        )


def material_to_operad_algebra(row: dict[str, Any]) -> OperadAlgebra | None:
    """Convert a material data row to an O-algebra.

    Expected row keys:
      - 'jid' or 'material_id': identifier
      - 'formula': chemical formula
      - 'point_group' or 'spg': point group or spacegroup
      - 'co_n{n}_w{k}': cyclic operad weights (optional)
    """
    mat_id = row.get("jid") or row.get("material_id", "unknown")
    formula = row.get("formula", "?")

    pg = row.get("point_group")
    if pg is None:
        spg = row.get("spg")
        pg = spg_to_point_group(spg)
    if pg is None:
        return None

    channels = allowed_channels(pg)

    weights: dict[tuple[int, int], float] = {}
    for n in (3, 4, 5, 6):
        for k in range(n):
            key = f"co_n{n}_w{k}"
            val = row.get(key)
            if val is not None:
                try:
                    weights[(n, k)] = float(val)
                except (ValueError, TypeError):
                    pass

    metadata = {}
    for meta_key in ("Tc_central_K", "Tc_v11_K", "Tc_lit_K",
                     "z3", "anomaly_class", "h2_order",
                     "is_ductile", "is_magnetic", "has_strong_soc",
                     "synthesizable", "family"):
        if meta_key in row:
            metadata[meta_key] = row[meta_key]

    return OperadAlgebra(
        material_id=mat_id,
        formula=formula,
        point_group=pg,
        channels=channels,
        weights=weights,
        metadata=metadata,
    )


def free_algebra_generators(point_group: str) -> OperadAlgebra:
    """The free O-algebra on one generator for the given point group.

    This is the Sigrist-Ueda entry itself: all channels present with
    unit weight, no cyclic weight data. Represents the "universal"
    material of that symmetry type.
    """
    channels = allowed_channels(point_group)
    weights: dict[tuple[int, int], float] = {}
    for ch in channels:
        n = CHANNEL_CYCLIC_ORDER.get(ch, 1)
        if n > 1:
            for k in range(n):
                weights[(n, k)] = 1.0 / n
    return OperadAlgebra(
        material_id=f"free_{point_group}",
        formula=f"Free({point_group})",
        point_group=point_group,
        channels=channels,
        weights=weights,
    )
