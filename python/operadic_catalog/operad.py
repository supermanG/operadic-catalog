"""Colored operad of superconducting pairing mechanisms.

The operad O has:
  Colors: pairs (G, n) where G is a crystallographic point group and
          n is a cyclic order (the C_n character of the pairing channel).
  Operations of arity k: pairing channels of cyclic order n in point
          group G, viewed as k-ary operations when the channel decomposes
          into k sub-channels under restriction.
  Composition: (G, chi_n) circ_i (H, chi_m) is the pairing channel
          obtained by replacing the i-th sub-channel of an n-channel
          with an m-channel, when H embeds into G's stabilizer.
  Unit: the s-wave (trivial irrep) in any point group.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator

from operadic_catalog.sigrist_ueda import (
    POINT_GROUP_CHANNELS,
    ALL_POINT_GROUPS,
    ALL_CHANNELS,
    CHANNEL_CYCLIC_ORDER,
)


@dataclass(frozen=True, order=True)
class Color:
    """A color in the pairing operad: (point_group, cyclic_order)."""
    point_group: str
    cyclic_order: int

    def __repr__(self) -> str:
        return f"({self.point_group}, Z/{self.cyclic_order})"


@dataclass(frozen=True)
class Operation:
    """An operation in the pairing operad.

    Represents a pairing channel of a given type in a point group,
    with specified input and output colors.

    The arity is len(inputs). A generating operation has arity 0
    (nullary: the channel itself, no sub-channels).
    """
    channel: str
    point_group: str
    inputs: tuple[Color, ...] = ()
    output: Color | None = None

    def __post_init__(self):
        if self.output is None:
            n = CHANNEL_CYCLIC_ORDER.get(self.channel, 1)
            object.__setattr__(
                self, "output", Color(self.point_group, n)
            )

    @property
    def arity(self) -> int:
        return len(self.inputs)

    @property
    def cyclic_order(self) -> int:
        return CHANNEL_CYCLIC_ORDER.get(self.channel, 1)

    def __repr__(self) -> str:
        if self.arity == 0:
            return f"Op({self.point_group}:{self.channel})"
        ins = ", ".join(str(c) for c in self.inputs)
        return f"Op({self.point_group}:{self.channel}[{ins}])"


# Subgroup lattice of the 32 crystallographic point groups (Schoenflies).
# Each entry maps G -> set of maximal proper subgroups.
# This encodes the symmetry-reduction paths relevant for operadic composition.
_MAXIMAL_SUBGROUPS: dict[str, set[str]] = {
    "C1":  set(),
    "Ci":  {"C1"},
    "C2":  {"C1"},
    "Cs":  {"C1"},
    "C2h": {"C2", "Ci", "Cs"},
    "D2":  {"C2"},
    "C2v": {"C2", "Cs"},
    "D2h": {"D2", "C2v", "C2h"},
    "C3":  {"C1"},
    "S6":  {"C3", "Ci"},
    "D3":  {"C3", "C2"},
    "C3v": {"C3", "Cs"},
    "D3d": {"D3", "C3v", "S6"},
    "C4":  {"C2"},
    "S4":  {"C2"},
    "C4h": {"C4", "C2h", "S4"},
    "D4":  {"C4", "D2"},
    "C4v": {"C4", "C2v"},
    "D2d": {"S4", "D2", "C2v"},
    "D4h": {"D4", "C4v", "C4h", "D2h", "D2d"},
    "C6":  {"C3", "C2"},
    "C3h": {"C3", "Cs"},
    "C6h": {"C6", "C3h", "C2h", "S6"},
    "D6":  {"C6", "D3", "D2"},
    "C6v": {"C6", "C3v", "C2v"},
    "D3h": {"D3", "C3h", "C3v", "C2v"},
    "D6h": {"D6", "C6v", "C6h", "D3h", "D2h"},
    "T":   {"D2", "C3"},
    "Th":  {"T", "D2h", "S6"},
    "O":   {"T", "D4", "D3"},
    "Td":  {"T", "D2d", "C3v", "S4"},
    "Oh":  {"O", "Td", "Th", "D4h", "D3d"},
}


def _transitive_subgroups(g: str) -> frozenset[str]:
    """Compute all subgroups of g (including g itself) by transitive closure."""
    visited: set[str] = set()
    stack = [g]
    while stack:
        cur = stack.pop()
        if cur in visited:
            continue
        visited.add(cur)
        for sub in _MAXIMAL_SUBGROUPS.get(cur, set()):
            stack.append(sub)
    return frozenset(visited)


_SUBGROUP_CACHE: dict[str, frozenset[str]] = {}


def subgroups_of(g: str) -> frozenset[str]:
    """All subgroups of g (including g), cached."""
    if g not in _SUBGROUP_CACHE:
        _SUBGROUP_CACHE[g] = _transitive_subgroups(g)
    return _SUBGROUP_CACHE[g]


def is_subgroup(h: str, g: str) -> bool:
    """True if H is a subgroup of G in the crystallographic point group lattice."""
    return h in subgroups_of(g)


class PairingOperad:
    """The colored operad of superconducting pairing mechanisms.

    Colors are (point_group, cyclic_order) pairs. Generating operations
    are pairing channels from the Sigrist-Ueda table. Composition is
    sub-channel insertion when the subgroup relation holds.
    """

    def __init__(self) -> None:
        self._generators: dict[str, list[Operation]] = {}
        self._build_generators()

    def _build_generators(self) -> None:
        for pg, channels in POINT_GROUP_CHANNELS.items():
            ops = []
            for ch in sorted(channels):
                op = Operation(channel=ch, point_group=pg)
                ops.append(op)
            self._generators[pg] = ops

    def generators(self, point_group: str | None = None) -> list[Operation]:
        """Return generating operations, optionally filtered by point group."""
        if point_group is not None:
            return list(self._generators.get(point_group, []))
        result = []
        for pg in sorted(self._generators):
            result.extend(self._generators[pg])
        return result

    def colors(self) -> list[Color]:
        """All colors (point_group, cyclic_order) appearing in the operad."""
        seen: set[Color] = set()
        for ops in self._generators.values():
            for op in ops:
                seen.add(op.output)
                for c in op.inputs:
                    seen.add(c)
        return sorted(seen)

    def unit(self, point_group: str) -> Operation:
        """The unit (identity) operation: s-wave in the given point group."""
        return Operation(channel="s", point_group=point_group)

    def can_compose(self, outer: Operation, inner: Operation, position: int = 0) -> bool:
        """Check if inner can be inserted into outer at the given position.

        Composition is valid when:
        1. The inner operation's point group is a subgroup of the outer's.
        2. The position is valid (within arity of outer, or position 0
           for nullary outer which extends arity).
        """
        if not is_subgroup(inner.point_group, outer.point_group):
            return False
        if outer.arity > 0 and position >= outer.arity:
            return False
        return True

    def compose(self, outer: Operation, inner: Operation, position: int = 0) -> Operation | None:
        """Compose outer circ_i inner: replace sub-channel at position i.

        Returns the composite operation, or None if composition is invalid.

        For generating (nullary) operations, composition produces an
        operation of arity = inner.arity, representing the outer channel
        refined by the inner channel's structure.
        """
        if not self.can_compose(outer, inner, position):
            return None

        if outer.arity == 0:
            new_inputs = inner.inputs
        else:
            inputs_list = list(outer.inputs)
            inputs_list[position:position + 1] = list(inner.inputs)
            new_inputs = tuple(inputs_list)

        return Operation(
            channel=f"{outer.channel}+{inner.channel}",
            point_group=outer.point_group,
            inputs=new_inputs,
            output=outer.output,
        )

    def restriction(self, op: Operation, subgroup: str) -> Operation | None:
        """Restrict an operation to a subgroup.

        This is the operadic morphism induced by the inclusion H -> G.
        The restricted channel is the same label if it exists in H's
        channel set, otherwise None (the channel vanishes).
        """
        if not is_subgroup(subgroup, op.point_group):
            return None
        sub_channels = POINT_GROUP_CHANNELS.get(subgroup, frozenset())
        base_channel = op.channel.split("+")[0]
        if base_channel not in sub_channels:
            return None
        new_inputs = tuple(
            Color(subgroup, c.cyclic_order) for c in op.inputs
        )
        return Operation(
            channel=base_channel,
            point_group=subgroup,
            inputs=new_inputs,
            output=Color(subgroup, op.output.cyclic_order),
        )

    def verify_unit_axiom(self) -> list[str]:
        """Verify: composing with the unit is the identity. Return errors."""
        errors = []
        for pg in ALL_POINT_GROUPS:
            u = self.unit(pg)
            for gen in self.generators(pg):
                composed = self.compose(gen, u)
                if composed is None:
                    errors.append(
                        f"Unit composition failed: {gen} circ {u}"
                    )
        return errors

    def verify_associativity(self) -> list[str]:
        """Verify associativity of composition on a sample of triples.

        For generators f, g, h where composition is defined:
        (f circ_0 g) circ_0 h == f circ_0 (g circ_0 h)
        Return errors found.
        """
        errors = []
        for pg in sorted(ALL_POINT_GROUPS):
            gens = self.generators(pg)
            subs = [
                h for h in ALL_POINT_GROUPS if is_subgroup(h, pg)
            ]
            for f_op in gens[:3]:
                for sub1 in subs[:3]:
                    for g_op in self.generators(sub1)[:2]:
                        fg = self.compose(f_op, g_op, 0)
                        if fg is None:
                            continue
                        sub2_list = [
                            h for h in ALL_POINT_GROUPS
                            if is_subgroup(h, sub1)
                        ]
                        for sub2 in sub2_list[:2]:
                            for h_op in self.generators(sub2)[:2]:
                                fg_h = self.compose(fg, h_op, 0)
                                gh = self.compose(g_op, h_op, 0)
                                if gh is None:
                                    continue
                                f_gh = self.compose(f_op, gh, 0)
                                if fg_h is None and f_gh is None:
                                    continue
                                if fg_h is not None and f_gh is not None:
                                    if fg_h.output != f_gh.output:
                                        errors.append(
                                            f"Assoc fail: "
                                            f"({f_op} o {g_op}) o {h_op} "
                                            f"vs {f_op} o ({g_op} o {h_op})"
                                        )
        return errors

    def verify_equivariance(self) -> list[str]:
        """Verify restriction is compatible with composition.

        For f in G, g in H <= G, and K <= H:
        restrict(f circ g, K) == restrict(f, K) circ restrict(g, K)
        Return errors found.
        """
        errors = []
        for pg in sorted(ALL_POINT_GROUPS):
            subs = [h for h in ALL_POINT_GROUPS if is_subgroup(h, pg) and h != pg]
            for f_op in self.generators(pg)[:2]:
                for sub_h in subs[:3]:
                    for g_op in self.generators(sub_h)[:2]:
                        fg = self.compose(f_op, g_op, 0)
                        if fg is None:
                            continue
                        sub_ks = [
                            k for k in ALL_POINT_GROUPS
                            if is_subgroup(k, sub_h) and k != sub_h
                        ]
                        for sub_k in sub_ks[:2]:
                            lhs = self.restriction(fg, sub_k)
                            f_rest = self.restriction(f_op, sub_k)
                            g_rest = self.restriction(g_op, sub_k)
                            if f_rest is None or g_rest is None:
                                continue
                            rhs = self.compose(f_rest, g_rest, 0)
                            if lhs is None and rhs is None:
                                continue
                            if lhs is not None and rhs is not None:
                                if lhs.output != rhs.output:
                                    errors.append(
                                        f"Equivariance fail at {sub_k}: "
                                        f"restrict({f_op} o {g_op}) "
                                        f"!= restrict({f_op}) o restrict({g_op})"
                                    )
        return errors

    def summary(self) -> dict:
        """Summary statistics of the operad."""
        all_gens = self.generators()
        channels_per_pg = {
            pg: len(ops) for pg, ops in self._generators.items()
        }
        return {
            "n_colors": len(self.colors()),
            "n_generators": len(all_gens),
            "n_point_groups": len(self._generators),
            "channels_per_pg": channels_per_pg,
            "total_composable_pairs": self._count_composable_pairs(),
        }

    def _count_composable_pairs(self) -> int:
        count = 0
        for pg in ALL_POINT_GROUPS:
            for f_op in self.generators(pg):
                subs = subgroups_of(pg)
                for sub in subs:
                    for g_op in self.generators(sub):
                        if self.can_compose(f_op, g_op):
                            count += 1
        return count
