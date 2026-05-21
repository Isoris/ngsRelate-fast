"""
output.py — Phase 6: three writers for one chromosome's adaptive run.

1. `write_res(path, result)` — the per-chrom .res file. Same 23-column
   schema as the binary's output. ngsPedigree Stage 2 consumes this
   byte-for-byte. Rows are written in canonical (alphabetical) pair
   order; column values are written back verbatim from the binary's
   per-pair output strings so per-column number formatting is preserved.

2. `write_manifest(path, result)` — `.adaptive_manifest.tsv` audit
   sidecar. NOT consumed by ngsPedigree. Per-pair: prior class, budget,
   sites used, escalation reason, before/after theta etc.

3. `write_run_manifest(path, result, anchor_info, ...)` — JSON file
   matching ngsrelate_adaptive.run_manifest.v1.schema.json. Adds the
   prior-source anchor info, output file refs, and stats.
"""

from __future__ import annotations
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Any, Dict, Optional

from .scheduler import ChromosomeResult, ManifestRow


def _sha256_of(path: Optional[Path]) -> Optional[str]:
    if path is None or not Path(path).exists():
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _file_ref(path: Optional[Path], *, with_hash: bool = True) -> Optional[Dict[str, Any]]:
    if path is None:
        return None
    p = Path(path)
    if not p.exists():
        return {"path": str(p), "missing": True}
    out: Dict[str, Any] = {"path": str(p.resolve()),
                           "size_bytes": p.stat().st_size}
    if with_hash:
        out["sha256"] = _sha256_of(p)
    return out


# ---------------------------------------------------------------------------
# 1. .res writer — sacred schema, byte-for-byte ngsPedigree-compatible.
# ---------------------------------------------------------------------------

def write_res(path, result: ChromosomeResult) -> Path:
    """Write the per-chrom .res file in canonical row order.

    The header and per-column string values come straight from the binary's
    output — we do not reformat any number, so column format matches the
    binary's printf format exactly. This is what makes Gate 5 (per-laptop
    vs per-LANTA reproducibility) tractable.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not result.res_header:
        raise ValueError(
            f"refusing to write empty .res for {result.chrom_id}: "
            f"no completed pairs ({result.run_manifest.get('n_pairs_failed', 0)} failed)")
    with open(path, "w") as fh:
        fh.write("\t".join(result.res_header) + "\n")
        for row in result.res_rows:
            fh.write("\t".join(row.get(col, "") for col in result.res_header) + "\n")
    return path


# ---------------------------------------------------------------------------
# 2. .adaptive_manifest.tsv audit sidecar.
# ---------------------------------------------------------------------------

def write_manifest(path, result: ChromosomeResult) -> Path:
    """Write the per-pair audit TSV. NOT a Stage-2 input — informational."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(ManifestRow.tsv_header()) + "\n")
        for row in result.manifest:
            fh.write("\t".join(row.as_tsv_columns()) + "\n")
    return path


# ---------------------------------------------------------------------------
# 3. .adaptive_run_manifest.json — extends the contract pattern.
# ---------------------------------------------------------------------------

def write_run_manifest(
    path,
    result: ChromosomeResult,
    *,
    anchor_path,
    anchor_input_contract_path: Optional[Path] = None,
    sample_set_match: Optional[bool] = None,
    allow_mismatched_anchor: bool = False,
    res_path: Optional[Path] = None,
    adaptive_manifest_path: Optional[Path] = None,
) -> Path:
    """Write the per-chrom run manifest JSON.

    `anchor_path` is required — the adaptive scheduler refuses to run without
    a genome-wide anchor (CLARIFICATION_NOTE §2), so the manifest always
    records it.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    sites_used = [r.final_sites_used or r.initial_sites_used for r in result.manifest if r.final_sites_used or r.initial_sites_used]
    stats = {
        "elapsed_seconds": result.run_manifest.get("elapsed_seconds"),
        "mean_sites_per_pair":   (sum(sites_used) / len(sites_used)) if sites_used else None,
        "median_sites_per_pair": median(sites_used) if sites_used else None,
    }

    doc = dict(result.run_manifest)
    doc.update({
        "schema": "ngsrelate_adaptive.run_manifest.v1",
        "anchor": {
            "genome_wide_res_path": str(Path(anchor_path).resolve()),
            "genome_wide_res_sha256": _sha256_of(Path(anchor_path)),
            "genome_wide_input_contract_sha256": (
                _sha256_of(Path(anchor_input_contract_path))
                if anchor_input_contract_path else None),
            "sample_set_match": sample_set_match,
            "allow_mismatched_anchor": allow_mismatched_anchor,
        },
        "outputs": {k: v for k, v in {
            "res":               _file_ref(res_path),
            "adaptive_manifest": _file_ref(adaptive_manifest_path),
        }.items() if v is not None},
        "stats": stats,
        "completed_at": doc.get("completed_at", datetime.now(tz=timezone.utc).isoformat()),
    })

    with open(path, "w") as fh:
        json.dump(doc, fh, indent=2, sort_keys=False)
    return path


# ---------------------------------------------------------------------------
# Convenience: write all three at once.
# ---------------------------------------------------------------------------

def write_all(
    *,
    output_dir,
    result: ChromosomeResult,
    res_basename: str = "relatedness.res",
    anchor_path,
    anchor_input_contract_path: Optional[Path] = None,
    sample_set_match: Optional[bool] = None,
    allow_mismatched_anchor: bool = False,
) -> Dict[str, Path]:
    """Write .res + .adaptive_manifest.tsv + .adaptive_run_manifest.json."""
    output_dir = Path(output_dir)
    res_path      = output_dir / res_basename
    manifest_path = output_dir / f"{res_basename}.adaptive_manifest.tsv"
    runman_path   = output_dir / f"{res_basename}.adaptive_run_manifest.json"

    write_res(res_path, result)
    write_manifest(manifest_path, result)
    write_run_manifest(
        runman_path, result,
        anchor_path=anchor_path,
        anchor_input_contract_path=anchor_input_contract_path,
        sample_set_match=sample_set_match,
        allow_mismatched_anchor=allow_mismatched_anchor,
        res_path=res_path,
        adaptive_manifest_path=manifest_path,
    )
    return {"res": res_path, "manifest": manifest_path, "run_manifest": runman_path}
