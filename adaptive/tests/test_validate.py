"""
Tests for adaptive.scheduler.validate — Phase 8.

The validator's own correctness gets tested by feeding it hand-crafted
baseline + adaptive .res files where we know the verdict in advance.
"""
from __future__ import annotations
import json
import pytest
from pathlib import Path

from adaptive.scheduler.edge_class import EdgeClass
from adaptive.scheduler.validate import (
    GateVerdict,
    _load_pairs,
    gate1_edge_class_precision,
    gate4_stage2_semantic_compatibility,
    gate5_reproducibility,
    gate6_interesting_disagreement_preservation,
    validate_one_chromosome,
)


HEADER = "a\tb\tida\tidb\tnSites\ttheta\tIBS0\tKING"


def _write_res(path, rows):
    with open(path, "w") as fh:
        fh.write(HEADER + "\n")
        for i, (ida, idb, theta, ibs0) in enumerate(rows):
            fh.write("\t".join([
                str(i), str(i + 1), ida, idb, "5000",
                f"{theta:.6f}", f"{ibs0:.6f}", f"{theta * 1.9:.6f}"
            ]) + "\n")
    return path


# ---- _load_pairs ----------------------------------------------------------

def test_load_pairs_canonicalizes(tmp_path):
    p = _write_res(tmp_path / "x.res", [
        ("Z9", "A1", 0.01, 0.20),
        ("B2", "C3", 0.40, 0.001),
    ])
    pairs = _load_pairs(p)
    assert ("A1", "Z9") in pairs
    assert ("B2", "C3") in pairs


# ---- Gate 1 ---------------------------------------------------------------

def test_gate1_passes_when_all_agree(tmp_path):
    rows = [("S1", "S2", 0.40, 0.000),
            ("S1", "S3", 0.01, 0.200)]
    baseline = _load_pairs(_write_res(tmp_path / "base.res", rows))
    adaptive = _load_pairs(_write_res(tmp_path / "adap.res", rows))
    v = gate1_edge_class_precision(baseline, adaptive)
    assert v.passed
    assert v.detail["fraction_agree"] == 1.0


def test_gate1_fails_below_min(tmp_path):
    base_rows = [("S1", "S2", 0.40, 0.000),  # duplicate
                 ("S1", "S3", 0.40, 0.000)]
    ad_rows   = [("S1", "S2", 0.01, 0.200),  # unrelated — disagrees
                 ("S1", "S3", 0.01, 0.200)]  # unrelated — disagrees
    baseline = _load_pairs(_write_res(tmp_path / "b.res", base_rows))
    adaptive = _load_pairs(_write_res(tmp_path / "a.res", ad_rows))
    v = gate1_edge_class_precision(baseline, adaptive)
    assert v.passed is False
    assert v.detail["fraction_agree"] == 0.0


# ---- Gate 4 (tier 1) ------------------------------------------------------

def test_gate4_pass_on_perfect_agreement(tmp_path):
    rows = [("S1", "S2", 0.25, 0.001),   # PO
            ("S1", "S3", 0.25, 0.020),   # FS
            ("S1", "S4", 0.01, 0.200)]
    baseline = _load_pairs(_write_res(tmp_path / "b.res", rows))
    adaptive = _load_pairs(_write_res(tmp_path / "a.res", rows))
    v = gate4_stage2_semantic_compatibility(baseline, adaptive)
    assert v.passed is True
    assert v.detail["po_set_equality"] is True
    assert v.detail["fs_jaccard"] == 1.0


def test_gate4_fails_if_po_lost(tmp_path):
    base = [("S1", "S2", 0.25, 0.001),   # PO
            ("S3", "S4", 0.25, 0.001)]   # PO
    ad   = [("S1", "S2", 0.25, 0.001),   # PO
            ("S3", "S4", 0.10, 0.05)]    # 2nd — PO LOST
    baseline = _load_pairs(_write_res(tmp_path / "b.res", base))
    adaptive = _load_pairs(_write_res(tmp_path / "a.res", ad))
    v = gate4_stage2_semantic_compatibility(baseline, adaptive)
    assert v.passed is False
    assert ["S3", "S4"] in v.detail["po_missing_in_adaptive"]


def test_gate4_fails_if_theta_drift_exceeds_tol(tmp_path):
    base = [("S1", "S2", 0.10, 0.05)]    # 2nd
    ad   = [("S1", "S2", 0.20, 0.05)]    # 2nd, but Δθ=0.10 > 0.05
    baseline = _load_pairs(_write_res(tmp_path / "b.res", base))
    adaptive = _load_pairs(_write_res(tmp_path / "a.res", ad))
    v = gate4_stage2_semantic_compatibility(baseline, adaptive)
    assert v.passed is False
    assert v.detail["delta_theta_violations_count"] == 1


# ---- Gate 5 ---------------------------------------------------------------

def test_gate5_bit_identical(tmp_path):
    p1 = tmp_path / "x.res"; p1.write_text("HEADER\nROW\n")
    p2 = tmp_path / "y.res"; p2.write_text("HEADER\nROW\n")
    v = gate5_reproducibility(p1, p2)
    assert v.passed


def test_gate5_drift_detected(tmp_path):
    p1 = tmp_path / "x.res"; p1.write_text("HEADER\nROW_A\n")
    p2 = tmp_path / "y.res"; p2.write_text("HEADER\nROW_B\n")
    v = gate5_reproducibility(p1, p2)
    assert v.passed is False


def test_gate5_skipped_without_pair():
    v = gate5_reproducibility(None, None)
    assert v.passed is True
    assert v.detail.get("skipped") is True


# ---- Gate 6 (tier 1) ------------------------------------------------------

def test_gate6_zero_loss_of_po_to_unrelated(tmp_path):
    # Prior PO for (S1, S2). Baseline per-chrom: unrelated (the PO-violation).
    base = [("S1", "S2", 0.01, 0.20)]
    # Adaptive: also unrelated → interesting transition preserved.
    ad   = [("S1", "S2", 0.01, 0.20)]
    baseline = _load_pairs(_write_res(tmp_path / "b.res", base))
    adaptive = _load_pairs(_write_res(tmp_path / "a.res", ad))
    prior_map = {("S1", "S2"): EdgeClass.PARENT_OFFSPRING}
    v = gate6_interesting_disagreement_preservation(baseline, adaptive, prior_map)
    assert v.passed is True
    assert v.detail["interesting_recall"] == 1.0


def test_gate6_fails_if_po_violation_lost(tmp_path):
    # Baseline shows PO→unrelated (interesting). Adaptive papers over it.
    base = [("S1", "S2", 0.01, 0.20)]   # per-chrom: unrelated
    ad   = [("S1", "S2", 0.25, 0.001)]  # per-chrom: PO (PO-violation MASKED)
    baseline = _load_pairs(_write_res(tmp_path / "b.res", base))
    adaptive = _load_pairs(_write_res(tmp_path / "a.res", ad))
    prior_map = {("S1", "S2"): EdgeClass.PARENT_OFFSPRING}
    v = gate6_interesting_disagreement_preservation(baseline, adaptive, prior_map)
    assert v.passed is False
    assert ["S1", "S2"] in v.detail["po_to_unrelated_lost"]


# ---- Top-level integration ------------------------------------------------

def test_validate_one_chromosome_end_to_end(tmp_path):
    rows = [("S1", "S2", 0.40, 0.000),
            ("S1", "S3", 0.01, 0.200)]
    base = _write_res(tmp_path / "b.res", rows)
    adap = _write_res(tmp_path / "a.res", rows)
    report = validate_one_chromosome(
        baseline_res_path=base, adaptive_res_path=adap)
    assert report["overall_pass"] is True
    assert report["tier_1_pass"] is True
    assert report["tier_2_pass"] is True


def test_validate_one_chromosome_emits_json_serialisable(tmp_path):
    rows = [("S1", "S2", 0.25, 0.001)]
    base = _write_res(tmp_path / "b.res", rows)
    adap = _write_res(tmp_path / "a.res", rows)
    report = validate_one_chromosome(
        baseline_res_path=base, adaptive_res_path=adap)
    json.dumps(report)  # must not raise


def test_validate_cli_writes_report(tmp_path):
    rows = [("S1", "S2", 0.40, 0.000),
            ("S1", "S3", 0.01, 0.200)]
    base = _write_res(tmp_path / "b.res", rows)
    adap = _write_res(tmp_path / "a.res", rows)
    out  = tmp_path / "report.json"
    from adaptive.scheduler.validate import main as validate_main
    rc = validate_main([
        "--baseline-res", str(base),
        "--adaptive-res", str(adap),
        "--report-json", str(out),
    ])
    assert rc == 0
    assert out.exists()
    doc = json.loads(out.read_text())
    assert doc["overall_pass"] is True
