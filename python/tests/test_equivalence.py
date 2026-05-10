"""Tests for operadic equivalence classes."""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from operadic_catalog.algebra import (
    OperadAlgebra,
    material_to_operad_algebra,
)
from operadic_catalog.equivalence import (
    operad_equivalence_class,
    enumerate_classes,
    class_summary,
)
from operadic_catalog.morphisms import (
    compute_morphism,
    all_morphisms,
    build_morphism_graph,
    restriction_path,
    restrict_algebra,
    morphism_summary,
)
from operadic_catalog.sigrist_ueda import POINT_GROUP_CHANNELS


def _make_two_cuprates():
    """Two cuprate-like materials that should be equivalent."""
    row1 = {
        "jid": "CUP-001", "formula": "La2CuO4", "point_group": "D4h",
        "co_n4_w0": 0.2, "co_n4_w1": 0.1, "co_n4_w2": 0.6, "co_n4_w3": 0.1,
    }
    row2 = {
        "jid": "CUP-002", "formula": "YBa2Cu3O7", "point_group": "D4h",
        "co_n4_w0": 0.22, "co_n4_w1": 0.09, "co_n4_w2": 0.58, "co_n4_w3": 0.11,
    }
    return material_to_operad_algebra(row1), material_to_operad_algebra(row2)


def _make_cuprate_and_heusler():
    """A cuprate and a Heusler: should NOT be equivalent."""
    row1 = {
        "jid": "CUP-001", "formula": "La2CuO4", "point_group": "D4h",
        "co_n4_w2": 0.6,
    }
    row2 = {
        "jid": "HEU-001", "formula": "Co2MnSi", "point_group": "Oh",
        "co_n3_w1": 0.5,
    }
    return material_to_operad_algebra(row1), material_to_operad_algebra(row2)


def test_equivalent_cuprates():
    a, b = _make_two_cuprates()
    assert operad_equivalence_class(a, b, weight_tol=0.1)


def test_inequivalent_cuprate_heusler():
    a, b = _make_cuprate_and_heusler()
    assert not operad_equivalence_class(a, b)


def test_enumerate_coarse():
    rows = [
        {"jid": f"M{i}", "formula": f"M{i}", "point_group": "D4h"}
        for i in range(5)
    ] + [
        {"jid": f"N{i}", "formula": f"N{i}", "point_group": "Oh"}
        for i in range(3)
    ]
    algebras = [material_to_operad_algebra(r) for r in rows]
    algebras = [a for a in algebras if a is not None]
    classes = enumerate_classes(algebras, coarse=True)
    assert len(classes) == 2  # D4h and Oh


def test_enumerate_fine_splits():
    rows = [
        {"jid": "A1", "formula": "A1", "point_group": "D4h",
         "co_n4_w0": 0.9, "co_n4_w1": 0.0, "co_n4_w2": 0.1, "co_n4_w3": 0.0},
        {"jid": "A2", "formula": "A2", "point_group": "D4h",
         "co_n4_w0": 0.0, "co_n4_w1": 0.0, "co_n4_w2": 0.9, "co_n4_w3": 0.1},
    ]
    algebras = [material_to_operad_algebra(r) for r in rows]
    coarse = enumerate_classes(algebras, coarse=True)
    fine = enumerate_classes(algebras, weight_tol=0.05)
    assert len(coarse) == 1
    assert len(fine) == 2


def test_class_summary():
    rows = [
        {"jid": f"M{i}", "formula": f"M{i}", "point_group": "D4h"}
        for i in range(5)
    ]
    algebras = [material_to_operad_algebra(r) for r in rows]
    algebras = [a for a in algebras if a is not None]
    classes = enumerate_classes(algebras, coarse=True)
    summary = class_summary(classes)
    assert summary["n_classes"] == 1
    assert summary["total_materials"] == 5


def test_morphism_d4h_to_d2h():
    m = compute_morphism("D4h", "D2h")
    assert m is not None
    assert "d" in m.preserved_channels
    assert "z4" in m.lost_channels


def test_morphism_nonsubgroup():
    m = compute_morphism("D4h", "C3")
    assert m is None


def test_all_morphisms_nonempty():
    morphisms = all_morphisms(maximal_only=True)
    assert len(morphisms) > 0


def test_morphism_graph_is_dag():
    import networkx as nx
    G = build_morphism_graph(maximal_only=True)
    assert nx.is_directed_acyclic_graph(G)


def test_morphism_graph_has_32_nodes():
    G = build_morphism_graph()
    assert G.number_of_nodes() == 32


def test_restriction_path_d4h_to_c1():
    path = restriction_path("D4h", "C1")
    assert path is not None
    assert path[0] == "D4h"
    assert path[-1] == "C1"


def test_restrict_algebra():
    row = {
        "jid": "CUP-001", "formula": "La2CuO4", "point_group": "D4h",
        "co_n4_w2": 0.6,
    }
    alg = material_to_operad_algebra(row)
    restricted = restrict_algebra(alg, "D2h")
    assert restricted is not None
    assert restricted.point_group == "D2h"
    assert "z4" not in restricted.channels
    assert "d" in restricted.channels


def test_morphism_summary():
    summary = morphism_summary()
    assert summary["n_maximal_morphisms"] > 0
    assert summary["n_total_morphisms"] >= summary["n_maximal_morphisms"]
    assert "C1" in summary["sinks"]


if __name__ == "__main__":
    test_equivalent_cuprates()
    test_inequivalent_cuprate_heusler()
    test_enumerate_coarse()
    test_enumerate_fine_splits()
    test_class_summary()
    test_morphism_d4h_to_d2h()
    test_morphism_nonsubgroup()
    test_all_morphisms_nonempty()
    test_morphism_graph_is_dag()
    test_morphism_graph_has_32_nodes()
    test_restriction_path_d4h_to_c1()
    test_restrict_algebra()
    test_morphism_summary()
    print("All equivalence and morphism tests passed.")
