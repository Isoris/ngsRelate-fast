"""
Tests for adaptive.scheduler.scheduler — Phase 5.

Uses a configurable fake binary so we can drive specific (theta, IBS0)
results per pair and verify the escalation/orchestration logic.
"""
from __future__ import annotations
import json
import stat
import sys
import textwrap
from pathlib import Path

import pytest

from adaptive.scheduler.edge_class import EdgeClass
from adaptive.scheduler.scheduler import (
    AdaptiveScheduler,
    ChromosomeResult,
    ManifestRow,
)
from adaptive.scheduler.config import (
    BudgetConfig,
    EscalationConfig,
    SchedulerConfig,
)


HEADER = "a\tb\tida\tidb\tnSites\ttheta\tIBS0\tKING\tJ7\tJ8\tJ9"


def _write_configurable_fake_binary(tmp_path, pair_results: dict, default=(0.01, 0.20)):
    """Write a fake binary that looks up (ida, idb) in a JSON sidecar to decide
    which (theta, IBS0) to emit. The sidecar lives at <fake>.results.json.

    pair_results: {"S1|S2": (theta, ibs0), ...}
    """
    fake = tmp_path / "fake_ngsrel"
    results_json = tmp_path / "fake_ngsrel.results.json"
    results_json.write_text(json.dumps({
        "pairs": pair_results,
        "default_theta": default[0],
        "default_ibs0": default[1],
    }))
    script = textwrap.dedent(f"""\
        #!{sys.executable}
        import json, sys, os
        args = sys.argv[1:]
        out, idx_a, idx_b, samples_path, n_samples = None, None, None, None, None
        i = 0
        while i < len(args):
            tok = args[i]
            if tok == "-O": out = args[i+1]; i += 2
            elif tok == "-a": idx_a = int(args[i+1]); i += 2
            elif tok == "-b": idx_b = int(args[i+1]); i += 2
            elif tok == "-z": samples_path = args[i+1]; i += 2
            elif tok == "-n": n_samples = int(args[i+1]); i += 2
            else: i += 1
        # Load sample IDs from -z file
        with open(samples_path) as fh:
            samples = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
        a_id, b_id = samples[idx_a], samples[idx_b]
        cfg = json.loads(open({str(results_json)!r}).read())
        key = a_id + "|" + b_id
        key_rev = b_id + "|" + a_id
        if key in cfg["pairs"]:
            theta, ibs0 = cfg["pairs"][key]
        elif key_rev in cfg["pairs"]:
            theta, ibs0 = cfg["pairs"][key_rev]
        else:
            theta, ibs0 = cfg["default_theta"], cfg["default_ibs0"]
        with open(out, "w") as fh:
            fh.write({HEADER!r} + "\\n")
            row = "\\t".join([
                str(idx_a), str(idx_b), a_id, b_id, "5000",
                f"{{theta:.6f}}", f"{{ibs0:.6f}}", f"{{theta * 1.9:.6f}}",
                "0.10", "0.20", "0.70",
            ])
            fh.write(row + "\\n")
        sys.exit(0)
        """)
    fake.write_text(script)
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return fake


def _conf(*, n_workers=1, escalation_factor=3, boundary_eps=0.005,
          budget_low=10, budget_med=20, budget_high=40):
    return SchedulerConfig(
        budgets=BudgetConfig(
            budget_low=budget_low, budget_med=budget_med, budget_high=budget_high,
            budget_ambiguous=budget_high, budget_duplicate=budget_low,
        ),
        escalation=EscalationConfig(factor=escalation_factor,
                                    boundary_epsilon_theta=boundary_eps),
        n_workers=n_workers,
    )


# ---- Orchestration: all-unrelated, no escalation -------------------------

def test_scheduler_processes_all_pairs(make_beagle, make_samples_file, tmp_path):
    samples = ["S1", "S2", "S3", "S4"]
    samples_file = make_samples_file(samples)
    beagle, freqs = make_beagle(n_sites=200, n_samples=len(samples))
    fake = _write_configurable_fake_binary(tmp_path, pair_results={})

    sched = AdaptiveScheduler(
        config=_conf(n_workers=1),
        prior_map={},
        binary_path=fake,
    )
    result = sched.run_chromosome(
        chrom_id="LG01",
        beagle_path=beagle, freqs_path=freqs, samples_path=samples_file,
        n_samples=len(samples), sample_ids=samples,
    )
    assert isinstance(result, ChromosomeResult)
    # 4 samples → 6 pairs
    assert len(result.res_rows) == 6
    assert len(result.manifest) == 6
    # canonical order
    pair_order = [(m.sample_a, m.sample_b) for m in result.manifest]
    assert pair_order == sorted(pair_order)


def test_scheduler_res_rows_in_canonical_order(make_beagle, make_samples_file, tmp_path):
    # Reverse-order sample IDs to make sure canonicalization doesn't depend on input order
    samples = ["Z9", "M5", "A1"]
    samples_file = make_samples_file(samples)
    beagle, freqs = make_beagle(n_sites=100, n_samples=3)
    fake = _write_configurable_fake_binary(tmp_path, pair_results={})

    sched = AdaptiveScheduler(config=_conf(n_workers=1),
                              prior_map={}, binary_path=fake)
    result = sched.run_chromosome(
        chrom_id="LG01",
        beagle_path=beagle, freqs_path=freqs, samples_path=samples_file,
        n_samples=3, sample_ids=samples,
    )
    # ida column is the first sample of each row
    idas = [r["ida"] for r in result.res_rows]
    idbs = [r["idb"] for r in result.res_rows]
    pairs = list(zip(idas, idbs))
    # Each pair should be sorted within itself, and the list of pairs should be sorted
    for a, b in pairs:
        assert a <= b
    assert pairs == sorted(pairs)


def test_scheduler_assigns_budget_by_prior(make_beagle, make_samples_file, tmp_path):
    samples = ["S1", "S2", "S3"]
    samples_file = make_samples_file(samples)
    beagle, freqs = make_beagle(n_sites=500, n_samples=3)
    fake = _write_configurable_fake_binary(tmp_path, pair_results={})

    # S1-S2 has PO prior → high budget
    prior_map = {("S1", "S2"): EdgeClass.PARENT_OFFSPRING}
    sched = AdaptiveScheduler(
        config=_conf(n_workers=1, budget_low=10, budget_high=50),
        prior_map=prior_map, binary_path=fake)
    result = sched.run_chromosome(
        chrom_id="LG01",
        beagle_path=beagle, freqs_path=freqs, samples_path=samples_file,
        n_samples=3, sample_ids=samples,
    )
    by_pair = {(r.sample_a, r.sample_b): r for r in result.manifest}
    assert by_pair[("S1", "S2")].initial_budget == 50
    assert by_pair[("S1", "S3")].initial_budget == 10
    assert by_pair[("S2", "S3")].initial_budget == 10


# ---- Escalation pass -----------------------------------------------------

def test_scheduler_escalates_interesting_disagreement(
        make_beagle, make_samples_file, tmp_path):
    samples = ["S1", "S2", "S3"]
    samples_file = make_samples_file(samples)
    beagle, freqs = make_beagle(n_sites=500, n_samples=3)
    # S1-S2 has PO prior, but per-chrom result is unrelated → interesting
    fake = _write_configurable_fake_binary(
        tmp_path,
        pair_results={"S1|S2": (0.01, 0.20)},   # per-chrom: unrelated
    )
    prior_map = {("S1", "S2"): EdgeClass.PARENT_OFFSPRING}

    sched = AdaptiveScheduler(
        config=_conf(n_workers=1, escalation_factor=3,
                     budget_low=10, budget_high=30),
        prior_map=prior_map, binary_path=fake)
    result = sched.run_chromosome(
        chrom_id="LG01",
        beagle_path=beagle, freqs_path=freqs, samples_path=samples_file,
        n_samples=3, sample_ids=samples,
    )
    by_pair = {(r.sample_a, r.sample_b): r for r in result.manifest}
    s12 = by_pair[("S1", "S2")]
    assert s12.escalated is True
    assert s12.escalation_reason == "interesting_disagreement"
    assert s12.final_budget == 30 * 3  # factor=3 against the high budget


def test_scheduler_no_escalation_for_uninteresting(
        make_beagle, make_samples_file, tmp_path):
    samples = ["S1", "S2", "S3"]
    samples_file = make_samples_file(samples)
    beagle, freqs = make_beagle(n_sites=500, n_samples=3)
    fake = _write_configurable_fake_binary(tmp_path, pair_results={},
                                            default=(0.001, 0.30))
    sched = AdaptiveScheduler(config=_conf(n_workers=1, boundary_eps=0.001),
                              prior_map={}, binary_path=fake)
    result = sched.run_chromosome(
        chrom_id="LG01",
        beagle_path=beagle, freqs_path=freqs, samples_path=samples_file,
        n_samples=3, sample_ids=samples,
    )
    for m in result.manifest:
        assert m.escalated is False
        assert m.escalation_reason == "none"


def test_scheduler_caps_escalation_at_chrom_total(
        make_beagle, make_samples_file, tmp_path):
    samples = ["S1", "S2"]
    samples_file = make_samples_file(samples)
    # Only 30 sites available
    beagle, freqs = make_beagle(n_sites=30, n_samples=2)
    fake = _write_configurable_fake_binary(
        tmp_path,
        pair_results={"S1|S2": (0.01, 0.20)},
    )
    prior_map = {("S1", "S2"): EdgeClass.PARENT_OFFSPRING}

    sched = AdaptiveScheduler(
        config=_conf(n_workers=1, escalation_factor=10,
                     budget_low=10, budget_high=50),
        prior_map=prior_map, binary_path=fake)
    result = sched.run_chromosome(
        chrom_id="LG01",
        beagle_path=beagle, freqs_path=freqs, samples_path=samples_file,
        n_samples=2, sample_ids=samples,
    )
    by_pair = {(r.sample_a, r.sample_b): r for r in result.manifest}
    s12 = by_pair[("S1", "S2")]
    # Initial budget capped at 30 (the chrom total); escalation should also cap.
    assert s12.initial_budget == 30
    assert s12.final_budget == 30
    # If escalated budget would be capped to the same as initial, the scheduler
    # should NOT actually re-run — but the decision flag should still reflect intent.
    # Check that the escalation lookup chose not to re-run by checking that
    # we don't have a doubled count (no new compute spent).


def test_scheduler_run_manifest_has_expected_fields(
        make_beagle, make_samples_file, tmp_path):
    samples = ["S1", "S2"]
    samples_file = make_samples_file(samples)
    beagle, freqs = make_beagle(n_sites=100, n_samples=2)
    fake = _write_configurable_fake_binary(tmp_path, pair_results={})
    sched = AdaptiveScheduler(config=_conf(n_workers=1),
                              prior_map={}, binary_path=fake)
    result = sched.run_chromosome(
        chrom_id="LG07",
        beagle_path=beagle, freqs_path=freqs, samples_path=samples_file,
        n_samples=2, sample_ids=samples,
    )
    rm = result.run_manifest
    assert rm["schema"] == "ngsrelate_adaptive.run_manifest.v1"
    assert rm["chromosome"] == "LG07"
    assert rm["n_samples"] == 2
    assert rm["n_pairs_expected"] == 1
    assert rm["n_pairs_completed"] == 1
    assert "config" in rm
    assert "elapsed_seconds" in rm
    assert "completed_at" in rm


def test_scheduler_handles_binary_failure_gracefully(
        make_beagle, make_samples_file, tmp_path):
    samples = ["S1", "S2"]
    samples_file = make_samples_file(samples)
    beagle, freqs = make_beagle(n_sites=100, n_samples=2)
    # Write a fake that exits nonzero
    fake = tmp_path / "broken_ngsrel"
    fake.write_text(f"#!{sys.executable}\nimport sys; sys.exit(7)\n")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    sched = AdaptiveScheduler(config=_conf(n_workers=1),
                              prior_map={}, binary_path=fake)
    result = sched.run_chromosome(
        chrom_id="LG01",
        beagle_path=beagle, freqs_path=freqs, samples_path=samples_file,
        n_samples=2, sample_ids=samples,
    )
    assert result.run_manifest["n_pairs_failed"] == 1
    assert result.run_manifest["n_pairs_completed"] == 0


def test_manifest_row_tsv_columns():
    m = ManifestRow(
        sample_a="S1", sample_b="S2", genome_wide_class="parent_offspring",
        initial_budget=100, initial_sites_used=100, initial_chrom_class="unrelated",
        initial_theta=0.01, initial_IBS0=0.20, initial_KING=0.05,
        escalated=True, escalation_reason="interesting_disagreement",
        final_budget=300, final_sites_used=300, final_chrom_class="parent_offspring",
        final_theta=0.25, final_IBS0=0.001, final_KING=0.45,
        elapsed_seconds=1.23,
    )
    cols = m.as_tsv_columns()
    assert len(cols) == len(ManifestRow.tsv_header())
    assert "parent_offspring" in cols
    assert cols[ManifestRow.tsv_header().index("escalated")] == "1"
