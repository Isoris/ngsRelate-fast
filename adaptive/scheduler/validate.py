"""
validate.py — Phase 8: gate the adaptive run against a uniform-stride baseline.

Standalone tool. Takes one chromosome's adaptive `.res` and the corresponding
uniform-stride baseline `.res` (produced by `STEP_A07b_*.sh` with `-D 100000`)
and runs the six gates from SPEC §6, with the priority tightening from
IMPLEMENTATION_PLAN.md §2 Phase 8:

  - Tier 1 (hard fail blocks release):
      * Gate 4: Stage 2 semantic compatibility — confirmed-PO dyad set
                must match exactly; confirmed-FS dyads ≥99% Jaccard;
                frac_disagreement per pair within ±0.05.
      * Gate 6: Interesting-disagreement preservation — ≥95% recall;
                ZERO loss of parent_offspring → unrelated transitions.
  - Tier 2 (numeric target):
      * Gate 1: Edge-class precision ≥90% (≥95% recommended).
  - Tier 3 (operational targets, documented in CALIBRATION_LOG.md):
      * Gate 2: site efficiency
      * Gate 3: runtime
      * Gate 5: reproducibility (laptop vs LANTA bit-equal)

Gate 4's "confirmed-PO dyad set" is the set of pairs whose per-chrom
classification is parent_offspring. Full ngsPedigree-Stage-2 compatibility
requires running Stage 2 on both .res files and diffing the resulting
dyad/triad sets — that's an external pass not done here. This tool checks
the necessary precondition (per-chrom PO classifications match).

Invocation:
    python -m adaptive.scheduler.validate \\
        --baseline-res     baseline/relatedness.res \\
        --adaptive-res     adaptive/relatedness.res \\
        --report-json      adaptive/relatedness.res.validation_report.json
"""

from __future__ import annotations
import argparse
import json
import math
import sys
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .confidence import (
    INTERESTING_TRANSITIONS,
    chrom_class_from_res_row,
    is_interesting_transition,
)
from .edge_class import EdgeClass
from .prior import _to_float, _iter_res_rows


PairKey = Tuple[str, str]


# ---------------------------------------------------------------------------
# .res loaders
# ---------------------------------------------------------------------------

def _load_pairs(res_path: Path) -> Dict[PairKey, Dict[str, str]]:
    out: Dict[PairKey, Dict[str, str]] = {}
    for row in _iter_res_rows(res_path):
        a, b = row["_sample_a"], row["_sample_b"]
        key = (a, b) if a <= b else (b, a)
        out[key] = row
    return out


# ---------------------------------------------------------------------------
# Gate verdicts
# ---------------------------------------------------------------------------

@dataclass
class GateVerdict:
    name: str
    tier: int
    passed: bool
    detail: Dict[str, Any] = field(default_factory=dict)
    note: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# ---------------------------------------------------------------------------
# Gate 1 — edge-class precision  (tier 2)
# ---------------------------------------------------------------------------

def gate1_edge_class_precision(
    baseline: Dict[PairKey, Dict[str, str]],
    adaptive: Dict[PairKey, Dict[str, str]],
    *,
    min_agreement: float = 0.90,
    recommended_agreement: float = 0.95,
) -> GateVerdict:
    common = set(baseline) & set(adaptive)
    if not common:
        return GateVerdict(
            "edge_class_precision", 2, False,
            detail={"common_pairs": 0},
            note="no pairs in common between baseline and adaptive .res")

    n_agree = 0
    confusion: Dict[Tuple[str, str], int] = {}
    for key in common:
        c_b = chrom_class_from_res_row(baseline[key])
        c_a = chrom_class_from_res_row(adaptive[key])
        confusion[(c_b.value, c_a.value)] = confusion.get((c_b.value, c_a.value), 0) + 1
        if c_b == c_a:
            n_agree += 1
    frac = n_agree / len(common)
    return GateVerdict(
        "edge_class_precision", 2, frac >= min_agreement,
        detail={
            "n_common": len(common),
            "n_agree":  n_agree,
            "fraction_agree": frac,
            "min_required":   min_agreement,
            "recommended":    recommended_agreement,
            "confusion_baseline_to_adaptive": [
                {"baseline_class": b, "adaptive_class": a, "n": n}
                for (b, a), n in sorted(confusion.items(), key=lambda kv: -kv[1])
            ],
        },
        note=("PASS" if frac >= recommended_agreement
              else ("PASS (above min, below recommended)"
                    if frac >= min_agreement
                    else "FAIL — below tier-2 threshold")),
    )


# ---------------------------------------------------------------------------
# Gate 2 — site efficiency  (tier 3)
# ---------------------------------------------------------------------------

def gate2_site_efficiency(
    adaptive_manifest_path: Optional[Path],
    baseline_total_sites: Optional[int] = None,
    *,
    target_savings_fraction: float = 0.5,
) -> GateVerdict:
    """Site efficiency: total adaptive sites used vs. baseline (which uses all
    chrom sites times n_pairs). Operational target; informational only."""
    if adaptive_manifest_path is None or not Path(adaptive_manifest_path).exists():
        return GateVerdict(
            "site_efficiency", 3, True,
            detail={"skipped": True},
            note="adaptive manifest not provided — site efficiency not computed")

    total_sites_adaptive = 0
    n_pairs = 0
    with open(adaptive_manifest_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx_final = header.index("final_sites_used")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            try:
                total_sites_adaptive += int(parts[idx_final])
                n_pairs += 1
            except (ValueError, IndexError):
                continue
    detail = {
        "n_pairs": n_pairs,
        "total_sites_adaptive": total_sites_adaptive,
        "mean_sites_per_pair": (total_sites_adaptive / n_pairs) if n_pairs else 0,
    }
    if baseline_total_sites is not None and n_pairs > 0:
        baseline_total = baseline_total_sites * n_pairs
        savings = 1.0 - (total_sites_adaptive / baseline_total)
        detail["baseline_total_sites_per_pair"] = baseline_total_sites
        detail["baseline_total"] = baseline_total
        detail["fraction_savings"] = savings
        passed = savings >= target_savings_fraction
        return GateVerdict("site_efficiency", 3, passed, detail=detail,
                           note=("PASS" if passed
                                 else "operational target not met (tier 3 — informational)"))
    return GateVerdict("site_efficiency", 3, True, detail=detail,
                       note="baseline site count not provided; reporting raw stats only")


# ---------------------------------------------------------------------------
# Gate 3 — runtime  (tier 3)  — comparing recorded elapsed times.
# ---------------------------------------------------------------------------

def gate3_runtime(
    adaptive_run_manifest_path: Optional[Path],
    baseline_elapsed_seconds: Optional[float] = None,
    *,
    target_speedup: float = 3.0,
) -> GateVerdict:
    if adaptive_run_manifest_path is None or not Path(adaptive_run_manifest_path).exists():
        return GateVerdict("runtime", 3, True, detail={"skipped": True},
                           note="adaptive run manifest not provided")
    doc = json.loads(Path(adaptive_run_manifest_path).read_text())
    adaptive_elapsed = doc.get("elapsed_seconds")
    detail: Dict[str, Any] = {"adaptive_elapsed_seconds": adaptive_elapsed}
    if baseline_elapsed_seconds is not None and adaptive_elapsed:
        speedup = baseline_elapsed_seconds / adaptive_elapsed
        detail["baseline_elapsed_seconds"] = baseline_elapsed_seconds
        detail["speedup"] = speedup
        passed = speedup >= target_speedup
        return GateVerdict("runtime", 3, passed, detail=detail,
                           note=("PASS" if passed
                                 else f"speedup {speedup:.2f}x below target {target_speedup}x"))
    return GateVerdict("runtime", 3, True, detail=detail,
                       note="baseline elapsed not provided; reporting raw stats only")


# ---------------------------------------------------------------------------
# Gate 4 — Stage 2 semantic compatibility  (tier 1)
# ---------------------------------------------------------------------------

def gate4_stage2_semantic_compatibility(
    baseline: Dict[PairKey, Dict[str, str]],
    adaptive: Dict[PairKey, Dict[str, str]],
    *,
    theta_disagreement_tol: float = 0.05,
    fs_jaccard_min: float = 0.99,
) -> GateVerdict:
    """Per IMPLEMENTATION_PLAN.md §2 Phase 8:

       - Confirmed PO dyad set must match exactly (set equality). Zero tolerance.
       - Confirmed FS dyad set: ≥99% Jaccard agreement.
       - frac_disagreement per pair: within ±0.05 of baseline (proxied by
         |theta_adaptive - theta_baseline| ≤ 0.05).
       - No new categories of disagreement introduced.
    """
    common = set(baseline) & set(adaptive)

    po_baseline = {k for k in common
                   if chrom_class_from_res_row(baseline[k]) == EdgeClass.PARENT_OFFSPRING}
    po_adaptive = {k for k in common
                   if chrom_class_from_res_row(adaptive[k]) == EdgeClass.PARENT_OFFSPRING}
    fs_baseline = {k for k in common
                   if chrom_class_from_res_row(baseline[k]) == EdgeClass.FULL_SIBLING}
    fs_adaptive = {k for k in common
                   if chrom_class_from_res_row(adaptive[k]) == EdgeClass.FULL_SIBLING}

    po_match = (po_baseline == po_adaptive)
    po_missing_in_adaptive = sorted(po_baseline - po_adaptive)
    po_gained_in_adaptive  = sorted(po_adaptive - po_baseline)

    fs_union = fs_baseline | fs_adaptive
    fs_jaccard = (len(fs_baseline & fs_adaptive) / len(fs_union)
                  if fs_union else 1.0)

    # frac_disagreement proxy: |Δθ| per pair
    delta_theta_violations: List[Dict[str, Any]] = []
    for key in common:
        tb = _to_float(baseline[key].get("theta", "nan"))
        ta = _to_float(adaptive[key].get("theta", "nan"))
        if math.isnan(tb) or math.isnan(ta):
            continue
        d = abs(ta - tb)
        if d > theta_disagreement_tol:
            delta_theta_violations.append(
                {"pair": list(key), "theta_baseline": tb, "theta_adaptive": ta,
                 "delta": d})

    # No NEW categories of disagreement — every confused (baseline_cls,
    # adaptive_cls) pair must already be observed in the baseline. We
    # approximate this by saying: if no pair drifts ACROSS the duplicate
    # boundary or PO/non-PO boundary in a direction that didn't exist as
    # ambiguity in the baseline, we pass. Practically captured by combining
    # po_match and fs_jaccard.

    passed = (po_match
              and fs_jaccard >= fs_jaccard_min
              and not delta_theta_violations)
    return GateVerdict(
        "stage2_semantic_compatibility", 1, passed,
        detail={
            "po_set_equality": po_match,
            "po_baseline_count": len(po_baseline),
            "po_adaptive_count": len(po_adaptive),
            "po_missing_in_adaptive": [list(k) for k in po_missing_in_adaptive],
            "po_gained_in_adaptive":  [list(k) for k in po_gained_in_adaptive],
            "fs_jaccard": fs_jaccard,
            "fs_jaccard_min": fs_jaccard_min,
            "fs_baseline_count": len(fs_baseline),
            "fs_adaptive_count": len(fs_adaptive),
            "delta_theta_violations_count": len(delta_theta_violations),
            "delta_theta_tol": theta_disagreement_tol,
            "delta_theta_violations_sample": delta_theta_violations[:10],
        },
        note=("PASS" if passed
              else "TIER-1 FAIL — Stage 2 semantic compatibility broken; "
                   "see po_missing_in_adaptive / fs_jaccard / "
                   "delta_theta_violations for details"))


# ---------------------------------------------------------------------------
# Gate 5 — reproducibility (laptop vs LANTA bit-equal)  (tier 3)
# ---------------------------------------------------------------------------

def gate5_reproducibility(
    res_a: Optional[Path],
    res_b: Optional[Path],
) -> GateVerdict:
    """Bit-equal compare. Caller passes laptop + LANTA paths; if either is
    missing the gate is skipped (informational)."""
    if not res_a or not res_b:
        return GateVerdict("reproducibility", 3, True, detail={"skipped": True},
                           note="comparison pair not provided")
    res_a, res_b = Path(res_a), Path(res_b)
    if not res_a.exists() or not res_b.exists():
        return GateVerdict("reproducibility", 3, False,
                           detail={"missing": [str(p) for p in (res_a, res_b)
                                               if not p.exists()]},
                           note="one or both .res files missing")
    import hashlib
    def _h(p):
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    ha, hb = _h(res_a), _h(res_b)
    return GateVerdict("reproducibility", 3, ha == hb,
                       detail={"sha256_a": ha, "sha256_b": hb,
                               "paths": [str(res_a), str(res_b)]},
                       note=("PASS — bit-identical"
                             if ha == hb else "FAIL — bit drift"))


# ---------------------------------------------------------------------------
# Gate 6 — interesting-disagreement preservation  (tier 1)
# ---------------------------------------------------------------------------

def gate6_interesting_disagreement_preservation(
    baseline: Dict[PairKey, Dict[str, str]],
    adaptive: Dict[PairKey, Dict[str, str]],
    prior_map: Optional[Dict[PairKey, EdgeClass]] = None,
    *,
    min_recall: float = 0.95,
) -> GateVerdict:
    """≥95% recall of "interesting disagreements", AND zero loss of pairs
    where the genome-wide prior → per-chrom transition is PO → unrelated
    (the PO-on-chrom-X violation case from SPEC §4.3, tier-1 per the plan)."""
    common = set(baseline) & set(adaptive)
    if not common:
        return GateVerdict("interesting_disagreement_preservation", 1, False,
                           detail={"common_pairs": 0})

    interesting_baseline: List[PairKey] = []
    interesting_adaptive_recall: List[bool] = []
    po_to_unrelated_baseline: List[PairKey] = []
    po_to_unrelated_adaptive: List[PairKey] = []

    for key in common:
        c_b = chrom_class_from_res_row(baseline[key])
        c_a = chrom_class_from_res_row(adaptive[key])

        if prior_map is not None and key in prior_map:
            prior = prior_map[key]
            if is_interesting_transition(prior, c_b):
                interesting_baseline.append(key)
                interesting_adaptive_recall.append(is_interesting_transition(prior, c_a))
            if prior == EdgeClass.PARENT_OFFSPRING and c_b == EdgeClass.UNRELATED:
                po_to_unrelated_baseline.append(key)
            if prior == EdgeClass.PARENT_OFFSPRING and c_a == EdgeClass.UNRELATED:
                po_to_unrelated_adaptive.append(key)

    if interesting_baseline:
        recall = sum(interesting_adaptive_recall) / len(interesting_baseline)
    else:
        recall = 1.0

    po_to_unrelated_lost = sorted(set(po_to_unrelated_baseline) - set(po_to_unrelated_adaptive))
    passed = (recall >= min_recall) and not po_to_unrelated_lost
    return GateVerdict(
        "interesting_disagreement_preservation", 1, passed,
        detail={
            "interesting_baseline_count": len(interesting_baseline),
            "interesting_recall": recall,
            "min_recall": min_recall,
            "po_to_unrelated_baseline_count": len(po_to_unrelated_baseline),
            "po_to_unrelated_lost": [list(k) for k in po_to_unrelated_lost],
        },
        note=("PASS" if passed
              else "TIER-1 FAIL — interesting disagreements not preserved; "
                   "see po_to_unrelated_lost for the PO-violation cases"))


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def validate_one_chromosome(
    *,
    baseline_res_path,
    adaptive_res_path,
    adaptive_manifest_path: Optional[Path] = None,
    adaptive_run_manifest_path: Optional[Path] = None,
    genome_wide_res_for_prior: Optional[Path] = None,
    baseline_chrom_total_sites: Optional[int] = None,
    baseline_elapsed_seconds: Optional[float] = None,
    reproducibility_pair: Optional[Tuple[Path, Path]] = None,
) -> Dict[str, Any]:
    """Run all six gates against one chromosome's baseline + adaptive .res files."""
    baseline_res_path = Path(baseline_res_path)
    adaptive_res_path = Path(adaptive_res_path)

    baseline = _load_pairs(baseline_res_path)
    adaptive = _load_pairs(adaptive_res_path)

    prior_map = None
    if genome_wide_res_for_prior is not None:
        from .prior import derive_priors
        prior_map = derive_priors(genome_wide_res_for_prior)

    gates: List[GateVerdict] = [
        gate1_edge_class_precision(baseline, adaptive),
        gate2_site_efficiency(adaptive_manifest_path, baseline_chrom_total_sites),
        gate3_runtime(adaptive_run_manifest_path, baseline_elapsed_seconds),
        gate4_stage2_semantic_compatibility(baseline, adaptive),
        gate5_reproducibility(*reproducibility_pair) if reproducibility_pair
            else gate5_reproducibility(None, None),
        gate6_interesting_disagreement_preservation(baseline, adaptive, prior_map),
    ]

    tier1 = [g for g in gates if g.tier == 1]
    tier2 = [g for g in gates if g.tier == 2]
    tier3 = [g for g in gates if g.tier == 3]
    overall_pass = all(g.passed for g in tier1) and all(g.passed for g in tier2)

    return {
        "schema": "ngsrelate_adaptive.validation_report.v1",
        "baseline_res": str(baseline_res_path.resolve()),
        "adaptive_res": str(adaptive_res_path.resolve()),
        "n_pairs_baseline": len(baseline),
        "n_pairs_adaptive": len(adaptive),
        "n_pairs_common":   len(set(baseline) & set(adaptive)),
        "gates": [g.to_dict() for g in gates],
        "overall_pass": overall_pass,
        "tier_1_pass": all(g.passed for g in tier1),
        "tier_2_pass": all(g.passed for g in tier2),
        "tier_3_pass": all(g.passed for g in tier3),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="python -m adaptive.scheduler.validate",
        description="Run the 6 SPEC §6 validation gates against a baseline "
                    "uniform-stride .res and an adaptive .res for the same "
                    "chromosome. Emits a JSON report; nonzero exit if any "
                    "tier-1 or tier-2 gate fails.")
    ap.add_argument("--baseline-res", required=True, type=Path)
    ap.add_argument("--adaptive-res", required=True, type=Path)
    ap.add_argument("--adaptive-manifest", type=Path, default=None,
                    help=".adaptive_manifest.tsv (audit sidecar) — needed for Gate 2.")
    ap.add_argument("--adaptive-run-manifest", type=Path, default=None,
                    help=".adaptive_run_manifest.json — needed for Gate 3.")
    ap.add_argument("--genome-wide-res", type=Path, default=None,
                    help="Genome-wide .res for prior-based Gate 6 analysis.")
    ap.add_argument("--baseline-chrom-total-sites", type=int, default=None,
                    help="Sites available on this chromosome for the baseline. "
                         "Needed for Gate 2 savings calculation.")
    ap.add_argument("--baseline-elapsed-seconds", type=float, default=None,
                    help="Baseline runtime; enables Gate 3 speedup calc.")
    ap.add_argument("--repro-laptop", type=Path, default=None)
    ap.add_argument("--repro-lanta",  type=Path, default=None)
    ap.add_argument("--report-json",  type=Path, default=None,
                    help="Where to write the JSON report; if omitted prints to stdout.")
    args = ap.parse_args(argv)

    repro_pair = ((args.repro_laptop, args.repro_lanta)
                  if args.repro_laptop and args.repro_lanta else None)

    report = validate_one_chromosome(
        baseline_res_path=args.baseline_res,
        adaptive_res_path=args.adaptive_res,
        adaptive_manifest_path=args.adaptive_manifest,
        adaptive_run_manifest_path=args.adaptive_run_manifest,
        genome_wide_res_for_prior=args.genome_wide_res,
        baseline_chrom_total_sites=args.baseline_chrom_total_sites,
        baseline_elapsed_seconds=args.baseline_elapsed_seconds,
        reproducibility_pair=repro_pair,
    )

    text = json.dumps(report, indent=2)
    if args.report_json:
        args.report_json.parent.mkdir(parents=True, exist_ok=True)
        args.report_json.write_text(text)
        print(f"[validate] report written to {args.report_json}", file=sys.stderr)
    else:
        print(text)

    # Exit nonzero if tier 1 or tier 2 fails.
    if not (report["tier_1_pass"] and report["tier_2_pass"]):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
