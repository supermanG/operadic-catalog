"""Tests for the generator catalog."""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from operadic_catalog.generators import (
    catalog_generators,
    independent_generators,
    generator_relations,
    generator_summary,
)
from operadic_catalog.operad import PairingOperad
from operadic_catalog.sigrist_ueda import ALL_POINT_GROUPS


def test_all_pgs_have_generators():
    gens = catalog_generators()
    for pg in ALL_POINT_GROUPS:
        assert pg in gens, f"Missing generators for {pg}"
        assert len(gens[pg]) > 0, f"Empty generators for {pg}"


def test_s_wave_is_always_generator():
    gens = catalog_generators()
    for pg, ops in gens.items():
        channels = {op.channel for op in ops}
        assert "s" in channels, f"{pg} missing s-wave generator"


def test_d4h_has_6_generators():
    gens = catalog_generators()
    assert len(gens["D4h"]) == 6  # s, p, d, d_x2y2, d_xy, z4


def test_c1_has_1_generator():
    gens = catalog_generators()
    assert len(gens["C1"]) == 1  # just s


def test_oh_has_5_generators():
    gens = catalog_generators()
    assert len(gens["Oh"]) == 5  # s, p, d, z3, f


def test_independent_subset_of_all():
    all_g = catalog_generators()
    indep = independent_generators()
    for pg in ALL_POINT_GROUPS:
        ind_channels = {op.channel for op in indep.get(pg, [])}
        all_channels = {op.channel for op in all_g.get(pg, [])}
        assert ind_channels <= all_channels


def test_relations_consistent():
    rels = generator_relations()
    assert len(rels) > 0
    for r in rels:
        assert r.source_pg in ALL_POINT_GROUPS
        assert r.target_pg in ALL_POINT_GROUPS
        assert r.relation_type in ("restriction", "lost_under_restriction")


def test_summary_counts():
    summary = generator_summary()
    assert summary["total_generators"] > 0
    assert summary["total_generators"] >= summary["total_independent"]
    assert summary["restriction_relations"] > 0


if __name__ == "__main__":
    test_all_pgs_have_generators()
    test_s_wave_is_always_generator()
    test_d4h_has_6_generators()
    test_c1_has_1_generator()
    test_oh_has_5_generators()
    test_independent_subset_of_all()
    test_relations_consistent()
    test_summary_counts()
    print("All generator tests passed.")
