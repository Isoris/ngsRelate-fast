"""
Tests for adaptive.scheduler.output — Phase 6.
"""
from __future__ import annotations
import json
import pytest

from adaptive.scheduler.scheduler import (
    ChromosomeResult,
    ManifestRow,
    _config_snapshot,
)
from adaptive.scheduler.config import SchedulerConfig
from adaptive.scheduler.output import (
    write_res,
    write_manifest,
    write_run_manifest,
    write_all,
)


def _mk_result(*, chrom_id="LG01", with_rows=True):
    header = ["a", "b", "ida", "idb", "nSites", "theta", "IBS0", "KING"]
    rows = []
    manifest = []
    if with_rows:
        rows = [
            {"a": "0", "b": "1", "ida": "S1", "idb": "S2", "nSites": "5000",
             "theta": "0.250000", "IBS0": "0.001000", "KING": "0.475000"},
            {"a": "0", "b": "2", "ida": "S1", "idb": "S3", "nSites": "3000",
             "theta": "0.010000", "IBS0": "0.200000", "KING": "0.019000"},
        ]
        manifest = [
            ManifestRow(
                sample_a="S1", sample_b="S2", genome_wide_class="parent_offspring",
                initial_budget=30, initial_sites_used=30, initial_chrom_class="parent_offspring",
                initial_theta=0.25, initial_IBS0=0.001, initial_KING=0.475,
                escalated=False, escalation_reason="none",
                final_budget=30, final_sites_used=30, final_chrom_class="parent_offspring",
                final_theta=0.25, final_IBS0=0.001, final_KING=0.475,
                elapsed_seconds=0.05,
            ),
            ManifestRow(
                sample_a="S1", sample_b="S3", genome_wide_class="unrelated",
                initial_budget=10, initial_sites_used=10, initial_chrom_class="unrelated",
                initial_theta=0.01, initial_IBS0=0.20, initial_KING=0.019,
                escalated=False, escalation_reason="none",
                final_budget=10, final_sites_used=10, final_chrom_class="unrelated",
                final_theta=0.01, final_IBS0=0.20, final_KING=0.019,
                elapsed_seconds=0.03,
            ),
        ]
    return ChromosomeResult(
        chrom_id=chrom_id,
        res_header=header if with_rows else [],
        res_rows=rows,
        manifest=manifest,
        run_manifest={
            "schema": "ngsrelate_adaptive.run_manifest.v1",
            "chromosome": chrom_id,
            "n_samples": 3,
            "n_pairs_expected": 3,
            "n_pairs_completed": len(rows),
            "n_pairs_escalated": 0,
            "n_pairs_failed": 3 - len(rows),
            "chrom_total_sites": 1000,
            "binary_path": "/path/to/ngsRelate-fast",
            "config": _config_snapshot(SchedulerConfig()),
            "budgets_built": [10, 30],
            "elapsed_seconds": 1.23,
            "completed_at": "2026-05-21T00:00:00+00:00",
            "failed_pairs": [],
        },
    )


def test_write_res_round_trips(tmp_path):
    result = _mk_result()
    out = write_res(tmp_path / "relatedness.res", result)
    text = out.read_text()
    lines = text.splitlines()
    assert lines[0] == "a\tb\tida\tidb\tnSites\ttheta\tIBS0\tKING"
    assert len(lines) == 1 + len(result.res_rows)
    # Verbatim string values preserved
    assert "0.250000" in lines[1]


def test_write_res_canonical_order(tmp_path):
    result = _mk_result()
    out = write_res(tmp_path / "relatedness.res", result)
    lines = out.read_text().splitlines()
    # ida, idb columns from the rows: (S1, S2), (S1, S3) — already canonical
    for line in lines[1:]:
        parts = line.split("\t")
        assert parts[2] <= parts[3]


def test_write_res_empty_raises(tmp_path):
    result = _mk_result(with_rows=False)
    with pytest.raises(ValueError, match="empty"):
        write_res(tmp_path / "x.res", result)


def test_write_manifest(tmp_path):
    result = _mk_result()
    out = write_manifest(tmp_path / "x.adaptive_manifest.tsv", result)
    text = out.read_text()
    lines = text.splitlines()
    header = lines[0].split("\t")
    assert "genome_wide_class" in header
    assert "initial_budget" in header
    assert "final_chrom_class" in header
    assert "escalation_reason" in header
    assert len(lines) == 1 + len(result.manifest)


def test_write_run_manifest(tmp_path):
    result = _mk_result()
    anchor = tmp_path / "genomewide.res"
    anchor.write_text("a\tb\tida\tidb\ttheta\n0\t1\tS1\tS2\t0.25\n")
    res_path = tmp_path / "relatedness.res"
    write_res(res_path, result)
    out = write_run_manifest(
        tmp_path / "x.adaptive_run_manifest.json",
        result,
        anchor_path=anchor,
        sample_set_match=True,
        res_path=res_path,
    )
    doc = json.loads(out.read_text())
    assert doc["schema"] == "ngsrelate_adaptive.run_manifest.v1"
    assert doc["anchor"]["genome_wide_res_path"].endswith("genomewide.res")
    assert doc["anchor"]["genome_wide_res_sha256"] is not None
    assert doc["anchor"]["sample_set_match"] is True
    assert "outputs" in doc and "res" in doc["outputs"]
    assert "stats" in doc and "mean_sites_per_pair" in doc["stats"]


def test_run_manifest_validates_against_schema(tmp_path):
    """If jsonschema is available, the run manifest must validate."""
    jsonschema = pytest.importorskip("jsonschema")
    from pathlib import Path
    schema_path = Path(__file__).resolve().parents[2] / "contracts" / \
                  "ngsrelate_adaptive.run_manifest.v1.schema.json"
    schema = json.loads(schema_path.read_text())

    result = _mk_result()
    anchor = tmp_path / "genomewide.res"
    anchor.write_text("a\tb\tida\tidb\ttheta\n0\t1\tS1\tS2\t0.25\n")
    out = write_run_manifest(
        tmp_path / "x.adaptive_run_manifest.json",
        result, anchor_path=anchor, sample_set_match=True,
    )
    doc = json.loads(out.read_text())
    jsonschema.validate(doc, schema)


def test_write_all_creates_three_files(tmp_path):
    result = _mk_result()
    anchor = tmp_path / "genomewide.res"
    anchor.write_text("a\tb\ttheta\n0\t1\t0.25\n")
    paths = write_all(
        output_dir=tmp_path / "out",
        result=result,
        anchor_path=anchor,
    )
    assert paths["res"].exists()
    assert paths["manifest"].exists()
    assert paths["run_manifest"].exists()
    assert paths["run_manifest"].name.endswith(".adaptive_run_manifest.json")
