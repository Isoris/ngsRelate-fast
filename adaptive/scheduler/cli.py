"""
cli.py — Phase 7: single-chromosome entry point for the adaptive scheduler.

Invocation:

    python -m adaptive.scheduler \\
        --binary           bin/ngsRelate-fast \\
        --beagle           per_chrom.beagle.gz \\
        --freqs            per_chrom.freqs \\
        --samples          samples.txt \\
        --genome-wide-res  genomewide.res \\
        --chrom            C_gar_LG12 \\
        --out-dir          out/LG12/ \\
        --n-samples        226 \\
        --workers          8

Writes three files into `--out-dir`:
    relatedness.res                            (Stage-2-compatible)
    relatedness.res.adaptive_manifest.tsv      (audit sidecar)
    relatedness.res.adaptive_run_manifest.json (run metadata)

The genome-wide .res is the **prior source** per CLARIFICATION_NOTE §1.
The scheduler refuses to start if it's missing (CLARIFICATION_NOTE §2),
and warns/refuses on a stale anchor (CLARIFICATION_NOTE §3) unless
`--allow-mismatched-anchor` is passed.
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path
from typing import List, Optional

from .config import (
    BudgetConfig,
    EscalationConfig,
    InversionBoostConfig,
    SchedulerConfig,
)
from .output import write_all
from .prior import derive_priors
from .runner import load_sample_index
from .scheduler import AdaptiveScheduler


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m adaptive.scheduler",
        description="Adaptive per-pair site-budgeting for one chromosome of "
                    "ngsRelate-fast. Reads a genome-wide .res for the prior, "
                    "runs each pair at its budgeted site count, escalates "
                    "borderline/interesting pairs, writes a Stage-2-compatible "
                    ".res plus audit + run manifests. See "
                    "adaptive/docs/IMPLEMENTATION_PLAN.md.")
    # Required.
    p.add_argument("--binary",          required=True, type=Path)
    p.add_argument("--beagle",          required=True, type=Path)
    p.add_argument("--freqs",           required=True, type=Path)
    p.add_argument("--samples",         required=True, type=Path)
    p.add_argument("--genome-wide-res", required=True, type=Path,
                   help="Genome-wide .res from a prior ngsRelate-fast run. "
                        "Required prior source (CLARIFICATION_NOTE §2).")
    p.add_argument("--chrom",           required=True,
                   help="Chromosome label (used in manifests and tempdir names).")
    p.add_argument("--out-dir",         required=True, type=Path)
    p.add_argument("--n-samples",       required=True, type=int)
    # Concurrency.
    p.add_argument("--workers", type=int, default=None,
                   help="Worker pool size. Default: os.cpu_count().")
    p.add_argument("--threads-per-pair", type=int, default=1,
                   help="Threads passed to ngsRelate-fast per pair (-p). "
                        "Default 1 — parallelism is across pairs, not within.")
    p.add_argument("--per-pair-timeout", type=float, default=600.0,
                   help="Seconds before killing a single pair's binary call.")

    # Budget knobs (overrides for calibration).
    p.add_argument("--budget-low",       type=int, default=None)
    p.add_argument("--budget-med",       type=int, default=None)
    p.add_argument("--budget-high",      type=int, default=None)
    p.add_argument("--budget-ambiguous", type=int, default=None)
    p.add_argument("--budget-duplicate", type=int, default=None)

    # Escalation knobs.
    p.add_argument("--escalation-factor",  type=int, default=None,
                   help="SPEC §4.3 literal: new_budget = factor * initial. Default 3.")
    p.add_argument("--softer-escalation", action="store_true",
                   help="Use SPEC §7 OQ2 softer rule: max(B*2, B+5000).")
    p.add_argument("--boundary-epsilon", type=float, default=None,
                   help="θ distance from any KING boundary to trigger "
                        "at_boundary escalation. Default 0.01.")

    # Inversion boost.
    p.add_argument("--boost-inversion-chroms", default=None,
                   help="Comma-separated list of chromosome IDs whose "
                        "BUDGET_HIGH should be doubled (SPEC §7 OQ8).")

    # Anchor staleness.
    p.add_argument("--allow-mismatched-anchor", action="store_true",
                   help="Escape hatch (CLARIFICATION_NOTE §3). Default strict.")
    return p


def _build_config(args) -> SchedulerConfig:
    base = SchedulerConfig()
    budgets = BudgetConfig(
        budget_low       = args.budget_low       if args.budget_low       is not None else base.budgets.budget_low,
        budget_med       = args.budget_med       if args.budget_med       is not None else base.budgets.budget_med,
        budget_high      = args.budget_high      if args.budget_high      is not None else base.budgets.budget_high,
        budget_ambiguous = args.budget_ambiguous if args.budget_ambiguous is not None else base.budgets.budget_ambiguous,
        budget_duplicate = args.budget_duplicate if args.budget_duplicate is not None else base.budgets.budget_duplicate,
        budget_floor_fraction_of_available = base.budgets.budget_floor_fraction_of_available,
    )
    escalation = EscalationConfig(
        factor = args.escalation_factor if args.escalation_factor is not None else base.escalation.factor,
        softer_rule = args.softer_escalation,
        boundary_epsilon_theta = (args.boundary_epsilon
                                  if args.boundary_epsilon is not None
                                  else base.escalation.boundary_epsilon_theta),
    )
    boosted = tuple(c.strip() for c in args.boost_inversion_chroms.split(",")
                    if c.strip()) if args.boost_inversion_chroms else ()
    inv = InversionBoostConfig(enabled=bool(boosted), boosted_chroms=boosted)
    return SchedulerConfig(
        budgets=budgets,
        escalation=escalation,
        inversion_boost=inv,
        n_workers=args.workers,
        threads_per_pair=args.threads_per_pair,
        allow_mismatched_anchor=args.allow_mismatched_anchor,
    )


def _check_anchor(args, anchor_sample_set: set) -> bool:
    """Compare the per-chrom samples file to the genome-wide anchor's sample
    set. Returns True iff they match (set equality).

    CLARIFICATION_NOTE §3: refuses to run on mismatch unless --allow-mismatched-anchor.
    """
    per_chrom_samples = set(load_sample_index(args.samples).keys())
    match = (per_chrom_samples == anchor_sample_set)
    if not match and not args.allow_mismatched_anchor:
        only_chrom  = sorted(per_chrom_samples - anchor_sample_set)
        only_anchor = sorted(anchor_sample_set - per_chrom_samples)
        msg = (f"\nAnchor sample set ({len(anchor_sample_set)} IDs) does not match "
               f"per-chrom samples file ({len(per_chrom_samples)} IDs).\n")
        if only_chrom:
            msg += f"  Only in per-chrom: {only_chrom[:5]}{' ...' if len(only_chrom) > 5 else ''}\n"
        if only_anchor:
            msg += f"  Only in anchor:    {only_anchor[:5]}{' ...' if len(only_anchor) > 5 else ''}\n"
        msg += ("Pass --allow-mismatched-anchor to override (CLARIFICATION_NOTE §3). "
                "The default is strict because a stale anchor silently corrupts "
                "the per-pair prior — see DOWNSTREAM_CONSUMERS.md.")
        raise SystemExit(msg)
    return match


def _samples_in_genome_wide_res(path: Path) -> set:
    """Read sample IDs from a genome-wide .res (ida/idb columns)."""
    out = set()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "ida" not in header or "idb" not in header:
            # No ida/idb means no sample IDs to anchor against — return empty
            # and the caller will treat sample_set_match as None (unknown).
            return set()
        ai, bi = header.index("ida"), header.index("idb")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > max(ai, bi):
                out.add(parts[ai])
                out.add(parts[bi])
    return out


def main(argv: Optional[List[str]] = None) -> int:
    args = _build_argparser().parse_args(argv)

    # ---- Preflight: anchor exists, has rows, sample-set check ---------------
    if not args.genome_wide_res.exists():
        print(f"ERROR: --genome-wide-res not found: {args.genome_wide_res}",
              file=sys.stderr)
        print("The adaptive scheduler refuses to run without a genome-wide anchor "
              "(CLARIFICATION_NOTE §2). Run ngsRelate-fast genome-wide first, "
              "then point --genome-wide-res at its .res output.", file=sys.stderr)
        return 2

    print(f"[adaptive] Deriving priors from {args.genome_wide_res}",
          file=sys.stderr)
    prior_map = derive_priors(args.genome_wide_res)
    if not prior_map:
        print(f"ERROR: anchor .res yielded zero pairs: {args.genome_wide_res}",
              file=sys.stderr)
        return 2

    anchor_samples = _samples_in_genome_wide_res(args.genome_wide_res)
    sample_set_match = (None if not anchor_samples
                        else _check_anchor(args, anchor_samples))
    if sample_set_match is False and args.allow_mismatched_anchor:
        print("[adaptive] WARNING: anchor sample set mismatch, "
              "--allow-mismatched-anchor was passed — proceeding.",
              file=sys.stderr)

    # ---- Build scheduler + run ---------------------------------------------
    config = _build_config(args)
    sched = AdaptiveScheduler(
        config=config,
        prior_map=prior_map,
        binary_path=args.binary,
        per_pair_timeout=args.per_pair_timeout,
    )

    print(f"[adaptive] Running {args.chrom} "
          f"(n_samples={args.n_samples}, workers={config.n_workers or 'cpu_count()'})",
          file=sys.stderr)
    result = sched.run_chromosome(
        chrom_id=args.chrom,
        beagle_path=args.beagle,
        freqs_path=args.freqs,
        samples_path=args.samples,
        n_samples=args.n_samples,
    )

    # ---- Write outputs ------------------------------------------------------
    paths = write_all(
        output_dir=args.out_dir,
        result=result,
        anchor_path=args.genome_wide_res,
        sample_set_match=sample_set_match,
        allow_mismatched_anchor=args.allow_mismatched_anchor,
    )
    rm = result.run_manifest
    print(f"[adaptive] {args.chrom}: {rm['n_pairs_completed']}/"
          f"{rm['n_pairs_expected']} pairs, "
          f"{rm['n_pairs_escalated']} escalated, "
          f"{rm['n_pairs_failed']} failed, "
          f"{rm['elapsed_seconds']:.1f}s elapsed", file=sys.stderr)
    print(f"[adaptive] .res:           {paths['res']}", file=sys.stderr)
    print(f"[adaptive] audit manifest: {paths['manifest']}", file=sys.stderr)
    print(f"[adaptive] run manifest:   {paths['run_manifest']}", file=sys.stderr)

    return 0 if rm["n_pairs_failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
