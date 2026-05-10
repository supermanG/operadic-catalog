"""Tests for material -> O-algebra conversion."""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from operadic_catalog.algebra import (
    OperadAlgebra,
    material_to_operad_algebra,
    free_algebra_generators,
)
from operadic_catalog.operad import Operation
from operadic_catalog.sigrist_ueda import POINT_GROUP_CHANNELS


def _make_cuprate_row():
    """Synthetic cuprate-like material (D4h, dominant d_x2y2)."""
    return {
        "jid": "JVASP-TEST-001",
        "formula": "La2CuO4",
        "point_group": "D4h",
        "co_n3_w0": 0.4,
        "co_n3_w1": 0.3,
        "co_n3_w2": 0.3,
        "co_n4_w0": 0.2,
        "co_n4_w1": 0.1,
        "co_n4_w2": 0.6,
        "co_n4_w3": 0.1,
        "Tc_central_K": 38.0,
        "z3": 0.85,
        "anomaly_class": "Z2^3",
    }


def _make_heusler_row():
    """Synthetic Heusler-like material (Oh, dominant z3)."""
    return {
        "jid": "JVASP-TEST-002",
        "formula": "Co2MnSi",
        "point_group": "Oh",
        "co_n3_w0": 0.3,
        "co_n3_w1": 0.5,
        "co_n3_w2": 0.2,
        "Tc_central_K": 5.0,
    }


def test_material_to_algebra_cuprate():
    row = _make_cuprate_row()
    alg = material_to_operad_algebra(row)
    assert alg is not None
    assert alg.point_group == "D4h"
    assert alg.channels == POINT_GROUP_CHANNELS["D4h"]
    assert "d_x2y2" in alg.channels
    assert "z4" in alg.channels
    assert alg.weights[(4, 2)] == 0.6


def test_material_to_algebra_heusler():
    row = _make_heusler_row()
    alg = material_to_operad_algebra(row)
    assert alg is not None
    assert alg.point_group == "Oh"
    assert "f" in alg.channels
    assert "z3" in alg.channels


def test_material_to_algebra_from_spg():
    row = {"jid": "TEST-003", "formula": "NaCl", "spg": 225}
    alg = material_to_operad_algebra(row)
    assert alg is not None
    assert alg.point_group == "Oh"


def test_material_to_algebra_none_pg():
    row = {"jid": "TEST-004", "formula": "Unknown"}
    alg = material_to_operad_algebra(row)
    assert alg is None


def test_channel_signature():
    row = _make_cuprate_row()
    alg = material_to_operad_algebra(row)
    sig = alg.channel_signature
    assert sig == tuple(sorted(POINT_GROUP_CHANNELS["D4h"]))


def test_cyclic_profile():
    row = _make_cuprate_row()
    alg = material_to_operad_algebra(row)
    profile = alg.cyclic_profile
    assert 3 in profile
    assert 4 in profile
    assert len(profile[3]) == 3
    assert len(profile[4]) == 4


def test_dominant_channel_cuprate():
    row = _make_cuprate_row()
    alg = material_to_operad_algebra(row)
    dom = alg.dominant_channel
    assert dom is not None
    assert dom in alg.channels


def test_evaluate_on_allowed_channel():
    row = _make_cuprate_row()
    alg = material_to_operad_algebra(row)
    op = Operation(channel="z4", point_group="D4h")
    val = alg.evaluate(op)
    assert val > 0, f"z4 channel should have nonzero weight from co_n4 data, got {val}"


def test_evaluate_on_wrong_pg():
    row = _make_cuprate_row()
    alg = material_to_operad_algebra(row)
    op = Operation(channel="d", point_group="Oh")
    val = alg.evaluate(op)
    assert val == 0.0


def test_operadic_fingerprint_deterministic():
    row = _make_cuprate_row()
    a1 = material_to_operad_algebra(row)
    a2 = material_to_operad_algebra(row)
    assert a1.operadic_fingerprint() == a2.operadic_fingerprint()


def test_free_algebra():
    free_d4h = free_algebra_generators("D4h")
    assert free_d4h.point_group == "D4h"
    assert free_d4h.channels == POINT_GROUP_CHANNELS["D4h"]
    assert len(free_d4h.weights) > 0


def test_free_algebra_all_32():
    for pg in POINT_GROUP_CHANNELS:
        free = free_algebra_generators(pg)
        assert free.point_group == pg
        assert free.channels == POINT_GROUP_CHANNELS[pg]


if __name__ == "__main__":
    test_material_to_algebra_cuprate()
    test_material_to_algebra_heusler()
    test_material_to_algebra_from_spg()
    test_material_to_algebra_none_pg()
    test_channel_signature()
    test_cyclic_profile()
    test_dominant_channel_cuprate()
    test_evaluate_on_allowed_channel()
    test_evaluate_on_wrong_pg()
    test_operadic_fingerprint_deterministic()
    test_free_algebra()
    test_free_algebra_all_32()
    print("All algebra tests passed.")
