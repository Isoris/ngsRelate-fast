"""
Tests for adaptive.scheduler.confidence — Phase 4.

Exhaustive coverage of the 7×7 prior×per-chrom transition table — 49
cases — per IMPLEMENTATION_PLAN.md §2 Phase 4.
"""
from __future__ import annotations
import pytest

from adaptive.scheduler.edge_class import EdgeClass
from adaptive.scheduler.confidence import (
    EscalationDecision,
    INTERESTING_TRANSITIONS,
    REASON_AT_BOUNDARY,
    REASON_INTERESTING_DISAGREEMENT,
    REASON_NONE,
    chrom_class_from_res_row,
    distance_to_nearest_boundary,
    is_at_boundary,
    is_interesting_transition,
    should_escalate,
)
from adaptive.scheduler.config import (
    KING_THETA_DUPLICATE,
    KING_THETA_FIRST_DEGREE,
    KING_THETA_SECOND_DEG,
    KING_THETA_THIRD_DEG,
    SchedulerConfig,
    EscalationConfig,
)


# ---- distance_to_nearest_boundary ----------------------------------------

@pytest.mark.parametrize("theta, expected_dist_below", [
    (KING_THETA_DUPLICATE,       1e-12),
    (KING_THETA_FIRST_DEGREE,    1e-12),
    (KING_THETA_SECOND_DEG,      1e-12),
    (KING_THETA_THIRD_DEG,       1e-12),
])
def test_distance_zero_at_threshold(theta, expected_dist_below):
    assert distance_to_nearest_boundary(theta) <= expected_dist_below


def test_distance_in_middle_of_band():
    # midway between 0.0442 and 0.0884 → distance to nearest = 0.0221
    theta = (KING_THETA_THIRD_DEG + KING_THETA_SECOND_DEG) / 2
    d = distance_to_nearest_boundary(theta)
    assert abs(d - 0.0221) < 1e-6


def test_is_at_boundary():
    assert is_at_boundary(KING_THETA_THIRD_DEG, 0.001) is True
    assert is_at_boundary(KING_THETA_THIRD_DEG + 0.0005, 0.001) is True
    assert is_at_boundary(0.20, 0.001) is False


# ---- is_interesting_transition: exhaustive 7×7 table ---------------------

ALL_CLASSES = list(EdgeClass)


def test_interesting_transitions_table_is_49_pairs():
    pairs = [(p, c) for p in ALL_CLASSES for c in ALL_CLASSES]
    assert len(pairs) == 49


@pytest.mark.parametrize("prior,chrom", [
    (p, c) for p in ALL_CLASSES for c in ALL_CLASSES
])
def test_interesting_transition_each_cell(prior, chrom):
    """Validate every cell against the documented rule:
       - AMBIGUOUS_FIRST_DEGREE prior: interesting iff chrom != AMBIGUOUS_FIRST_DEGREE
       - All others: interesting iff (prior, chrom) in INTERESTING_TRANSITIONS
    """
    result = is_interesting_transition(prior, chrom)
    if prior == EdgeClass.AMBIGUOUS_FIRST_DEGREE:
        expected = (chrom != EdgeClass.AMBIGUOUS_FIRST_DEGREE)
    else:
        expected = (prior, chrom) in INTERESTING_TRANSITIONS
    assert result is expected


def test_interesting_transitions_includes_po_to_unrelated():
    # The PO-on-chrom-X violation case — Tier-1 Gate 6 depends on this
    # being treated as interesting (IMPLEMENTATION_PLAN.md §2 Phase 8).
    assert is_interesting_transition(
        EdgeClass.PARENT_OFFSPRING, EdgeClass.UNRELATED) is True


# ---- chrom_class_from_res_row --------------------------------------------

def test_chrom_class_from_res_row_basic():
    row = {"theta": "0.40", "IBS0": "0.000"}
    assert chrom_class_from_res_row(row) is EdgeClass.DUPLICATE_OR_CLONE


def test_chrom_class_from_res_row_po_vs_fs():
    assert chrom_class_from_res_row({"theta": "0.25", "IBS0": "0.001"}) \
        is EdgeClass.PARENT_OFFSPRING
    assert chrom_class_from_res_row({"theta": "0.25", "IBS0": "0.020"}) \
        is EdgeClass.FULL_SIBLING


def test_chrom_class_from_res_row_missing_ibs0_first_degree():
    # No IBS0 in row → first-degree falls back to ambiguous
    assert chrom_class_from_res_row({"theta": "0.25"}) \
        is EdgeClass.AMBIGUOUS_FIRST_DEGREE


# ---- should_escalate -----------------------------------------------------

def _conf(boundary_epsilon: float = 0.005) -> SchedulerConfig:
    return SchedulerConfig(
        escalation=EscalationConfig(boundary_epsilon_theta=boundary_epsilon))


def test_should_escalate_interesting_disagreement_wins():
    # Prior PO, per-chrom unrelated → interesting AND near a boundary at 0.0442
    row = {"theta": "0.045", "IBS0": "0.20"}
    decision = should_escalate(row, EdgeClass.PARENT_OFFSPRING, config=_conf(0.005))
    assert decision.should_escalate is True
    assert decision.reason == REASON_INTERESTING_DISAGREEMENT


def test_should_escalate_boundary_only():
    # Prior unrelated, per-chrom 3rd → not interesting (UNRELATED→THIRD_DEG not in set)
    # but theta exactly at 3rd-deg boundary → boundary
    row = {"theta": f"{KING_THETA_THIRD_DEG:.6f}", "IBS0": "0.20"}
    decision = should_escalate(row, EdgeClass.UNRELATED, config=_conf(0.001))
    assert decision.should_escalate is True
    assert decision.reason == REASON_AT_BOUNDARY


def test_should_escalate_none():
    # Prior unrelated, per-chrom unrelated, theta in mid-band → no escalation
    row = {"theta": "0.001", "IBS0": "0.30"}
    decision = should_escalate(row, EdgeClass.UNRELATED, config=_conf(0.005))
    assert decision.should_escalate is False
    assert decision.reason == REASON_NONE


def test_should_escalate_ambiguous_prior_any_concrete_chrom():
    row = {"theta": "0.25", "IBS0": "0.001"}  # → PO
    decision = should_escalate(row, EdgeClass.AMBIGUOUS_FIRST_DEGREE,
                               config=_conf(0.001))
    assert decision.should_escalate is True
    assert decision.reason == REASON_INTERESTING_DISAGREEMENT


def test_escalation_decision_dataclass():
    d = EscalationDecision(True, REASON_AT_BOUNDARY, 0.001, EdgeClass.UNRELATED)
    assert d.as_tuple() == (True, REASON_AT_BOUNDARY)
