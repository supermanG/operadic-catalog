"""Tests for the pairing operad axioms: unit, associativity, equivariance."""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from operadic_catalog.operad import (
    PairingOperad,
    Operation,
    Color,
    is_subgroup,
    subgroups_of,
)
from operadic_catalog.sigrist_ueda import ALL_POINT_GROUPS, POINT_GROUP_CHANNELS


def test_all_32_point_groups_present():
    assert len(ALL_POINT_GROUPS) == 32


def test_every_pg_has_s_wave():
    for pg, channels in POINT_GROUP_CHANNELS.items():
        assert "s" in channels, f"{pg} missing s-wave"


def test_operad_construction():
    O = PairingOperad()
    assert O.summary()["n_point_groups"] == 32
    assert O.summary()["n_generators"] > 0


def test_unit_axiom():
    O = PairingOperad()
    errors = O.verify_unit_axiom()
    assert errors == [], f"Unit axiom violations: {errors}"


def test_associativity():
    O = PairingOperad()
    errors = O.verify_associativity()
    assert errors == [], f"Associativity violations: {errors}"


def test_equivariance():
    O = PairingOperad()
    errors = O.verify_equivariance()
    assert errors == [], f"Equivariance violations: {errors}"


def test_subgroup_reflexive():
    for pg in ALL_POINT_GROUPS:
        assert is_subgroup(pg, pg)


def test_subgroup_transitive():
    for g in ALL_POINT_GROUPS:
        for h in subgroups_of(g):
            for k in subgroups_of(h):
                assert is_subgroup(k, g), f"{k} <= {h} <= {g} but {k} not <= {g}"


def test_c1_is_universal_subgroup():
    for pg in ALL_POINT_GROUPS:
        assert is_subgroup("C1", pg), f"C1 should be subgroup of {pg}"


def test_oh_has_most_subgroups():
    sizes = {pg: len(subgroups_of(pg)) for pg in ALL_POINT_GROUPS}
    oh_size = sizes["Oh"]
    for pg, sz in sizes.items():
        assert sz <= oh_size, f"{pg} has {sz} subgroups, Oh has {oh_size}"


def test_composition_requires_subgroup():
    O = PairingOperad()
    d4h_d = Operation(channel="d", point_group="D4h")
    c3_z3 = Operation(channel="z3", point_group="C3")
    assert not O.can_compose(d4h_d, c3_z3), "C3 is not subgroup of D4h via D-type"


def test_composition_within_same_group():
    O = PairingOperad()
    oh_s = O.unit("Oh")
    oh_d = Operation(channel="d", point_group="Oh")
    assert O.can_compose(oh_d, oh_s)
    result = O.compose(oh_d, oh_s)
    assert result is not None
    assert result.output == oh_d.output


def test_restriction_preserves_channels():
    O = PairingOperad()
    d4h_d = Operation(channel="d", point_group="D4h")
    restricted = O.restriction(d4h_d, "D2h")
    assert restricted is not None
    assert restricted.channel == "d"
    assert restricted.point_group == "D2h"


def test_restriction_loses_channels():
    O = PairingOperad()
    d4h_z4 = Operation(channel="z4", point_group="D4h")
    restricted = O.restriction(d4h_z4, "D2h")
    assert restricted is None, "z4 should vanish when restricting D4h to D2h"


def test_colors():
    O = PairingOperad()
    colors = O.colors()
    assert len(colors) > 0
    for c in colors:
        assert c.point_group in ALL_POINT_GROUPS


if __name__ == "__main__":
    test_all_32_point_groups_present()
    test_every_pg_has_s_wave()
    test_operad_construction()
    test_unit_axiom()
    test_associativity()
    test_equivariance()
    test_subgroup_reflexive()
    test_subgroup_transitive()
    test_c1_is_universal_subgroup()
    test_oh_has_most_subgroups()
    test_composition_requires_subgroup()
    test_composition_within_same_group()
    test_restriction_preserves_channels()
    test_restriction_loses_channels()
    test_colors()
    print("All operad axiom tests passed.")
