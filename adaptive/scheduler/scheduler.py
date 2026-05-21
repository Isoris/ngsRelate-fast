"""
scheduler.py — Phase 5: the per-chromosome adaptive scheduler.

The scheduler:

1. Builds a BeagleSubsetCache for the chromosome.
2. For every pair, looks up the prior class and the corresponding budget.
3. Submits all pairs to a multiprocessing pool.
4. Collects results, scores each for escalation.
5. Re-runs escalated pairs at a larger budget.
6. Assembles the final per-chrom .res rows in canonical order (sorted by
   (sample_a, sample_b) with both alphabetized) and returns them along
   with the audit manifest and the run manifest.

The scheduler does NOT write files — that's Phase 6 (output.py). Keeping
write logic out of the scheduler makes it easier to test against in-memory
expectations and to drive from both the CLI and from external callers.

Concurrency model: a single multiprocessing.Pool, one task per pair.
Workers are short-lived subprocesses that each fork an ngsRelate-fast
binary invocation with `-a` and `-b` set. The BEAGLE subset cache is
built BEFORE the pool launches so all workers see consistent, already-
materialised inputs. No shared writable state.
"""

from __future__ import annotations
import multiprocessing as mp
import os
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .edge_class import EdgeClass
from .runner import (
    PairResult,
    RunnerError,
    load_sample_index,
    run_pair_on_chrom,
)
from .subset import BeagleSubsetCache
from .confidence import (
    EscalationDecision,
    REASON_NONE,
    should_escalate,
)
from . import config as cfg


PairKey = Tuple[str, str]


# ---------------------------------------------------------------------------
# Output containers (consumed by output.py)
# ---------------------------------------------------------------------------

@dataclass
class ManifestRow:
    """One row of the audit sidecar (.adaptive_manifest.tsv)."""
    sample_a: str
    sample_b: str
    genome_wide_class: str
    initial_budget: int
    initial_sites_used: int
    initial_chrom_class: str
    initial_theta: float
    initial_IBS0: float
    initial_KING: float
    escalated: bool
    escalation_reason: str
    final_budget: int
    final_sites_used: int
    final_chrom_class: str
    final_theta: float
    final_IBS0: float
    final_KING: float
    elapsed_seconds: float

    def as_tsv_columns(self) -> List[str]:
        return [
            self.sample_a, self.sample_b, self.genome_wide_class,
            str(self.initial_budget), str(self.initial_sites_used),
            self.initial_chrom_class,
            _fmt(self.initial_theta), _fmt(self.initial_IBS0), _fmt(self.initial_KING),
            "1" if self.escalated else "0",
            self.escalation_reason,
            str(self.final_budget), str(self.final_sites_used),
            self.final_chrom_class,
            _fmt(self.final_theta), _fmt(self.final_IBS0), _fmt(self.final_KING),
            f"{self.elapsed_seconds:.3f}",
        ]

    @staticmethod
    def tsv_header() -> List[str]:
        return [
            "sample_a", "sample_b", "genome_wide_class",
            "initial_budget", "initial_sites_used", "initial_chrom_class",
            "initial_theta", "initial_IBS0", "initial_KING",
            "escalated", "escalation_reason",
            "final_budget", "final_sites_used", "final_chrom_class",
            "final_theta", "final_IBS0", "final_KING",
            "elapsed_seconds",
        ]


@dataclass
class ChromosomeResult:
    chrom_id: str
    res_header: List[str]
    res_rows: List[Dict[str, str]]    # canonical-ordered per-chrom .res rows
    manifest: List[ManifestRow]       # 1 row per pair, canonical order
    run_manifest: Dict[str, Any]


def _fmt(x: float) -> str:
    if x != x:           # NaN
        return "nan"
    return f"{x:.6f}"


def _safe_float(s: Optional[str]) -> float:
    if s is None:
        return float("nan")
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


# ---------------------------------------------------------------------------
# Worker function — must be at module top level for multiprocessing pickling.
# ---------------------------------------------------------------------------

def _pool_worker(task: Tuple) -> Tuple[str, PairKey, Any]:
    """Run a single pair. Returns ('ok', pair, PairResult) or ('err', pair, str).

    Task tuple: (binary_path, beagle_path, freqs_path, samples_path,
                 sample_a, sample_b, n_samples, threads, timeout)
    """
    (binary_path, beagle_path, freqs_path, samples_path,
     sample_a, sample_b, n_samples, threads, timeout) = task
    try:
        result = run_pair_on_chrom(
            binary_path=binary_path,
            beagle_path=beagle_path,
            freqs_path=freqs_path,
            samples_path=samples_path,
            sample_a=sample_a,
            sample_b=sample_b,
            n_samples=n_samples,
            threads=threads,
            timeout=timeout,
        )
        return ("ok", (sample_a, sample_b), result)
    except Exception as e:   # capture and report; do not crash the pool
        return ("err", (sample_a, sample_b),
                f"{type(e).__name__}: {e}")


# ---------------------------------------------------------------------------
# Scheduler
# ---------------------------------------------------------------------------

class AdaptiveScheduler:
    """Per-chromosome adaptive scheduler.

    Construct once with a config + prior map + binary path; call
    `run_chromosome(...)` for each chromosome.
    """

    def __init__(
        self,
        *,
        config: cfg.SchedulerConfig,
        prior_map: Dict[PairKey, EdgeClass],
        binary_path,
        per_pair_timeout: Optional[float] = None,
    ):
        self.config = config
        self.prior_map = prior_map
        self.binary_path = Path(binary_path).resolve()
        if not self.binary_path.exists():
            raise FileNotFoundError(f"binary not found: {self.binary_path}")
        self.per_pair_timeout = per_pair_timeout

    # ---- Budget computation ------------------------------------------------

    def _budget_for(self, prior_class: EdgeClass, chrom_id: str,
                    n_total: int) -> int:
        b = self.config.budgets.budget_for(prior_class)

        # SPEC §7 OQ8 — boost on inversion chromosomes.
        if (self.config.inversion_boost.enabled
                and chrom_id in self.config.inversion_boost.boosted_chroms
                and prior_class in (EdgeClass.PARENT_OFFSPRING,
                                    EdgeClass.FULL_SIBLING,
                                    EdgeClass.AMBIGUOUS_FIRST_DEGREE)):
            b *= 2

        # Cap at available sites — there's no way to ask for more than exist.
        return max(1, min(b, n_total))

    # ---- Pair listing ------------------------------------------------------

    def _enumerate_pairs(self, sample_ids: List[str]) -> List[Tuple[str, str, EdgeClass]]:
        """All (a, b) pairs with a < b (alphabetical), plus prior class.

        Iterates over a sorted copy of sample_ids so each emitted (a, b)
        is internally canonical AND the overall list is in canonical
        order. This is the same order the .res file will be written in
        (Phase 6) and the prior_map is keyed in.

        Unknown pairs default to UNRELATED (low budget) — the safe default.
        """
        sorted_ids = sorted(sample_ids)
        out: List[Tuple[str, str, EdgeClass]] = []
        for i in range(len(sorted_ids)):
            for j in range(i + 1, len(sorted_ids)):
                a, b = sorted_ids[i], sorted_ids[j]   # a <= b guaranteed
                prior = self.prior_map.get((a, b), EdgeClass.UNRELATED)
                out.append((a, b, prior))
        return out

    # ---- Pool dispatch -----------------------------------------------------

    def _dispatch(
        self,
        tasks: List[Tuple],
    ) -> Dict[PairKey, Any]:
        """Run tasks through a pool (or serially if n_workers<=1).
        Returns dict pair_key -> PairResult (or RuntimeError on failure)."""
        results: Dict[PairKey, Any] = {}
        n_workers = self.config.n_workers or os.cpu_count() or 1

        if n_workers <= 1 or len(tasks) == 1:
            for t in tasks:
                status, key, payload = _pool_worker(t)
                results[key] = payload if status == "ok" else RuntimeError(payload)
            return results

        with mp.Pool(processes=n_workers) as pool:
            for status, key, payload in pool.imap_unordered(_pool_worker, tasks,
                                                            chunksize=64):
                results[key] = payload if status == "ok" else RuntimeError(payload)
        return results

    # ---- Per-chromosome entry point ---------------------------------------

    def run_chromosome(
        self,
        *,
        chrom_id: str,
        beagle_path,
        freqs_path,
        samples_path,
        n_samples: int,
        sample_ids: Optional[List[str]] = None,
        cache_tmpdir: Optional[Path] = None,
    ) -> ChromosomeResult:
        """Run the adaptive scheduler against one chromosome.

        Returns a ChromosomeResult ready to be written by output.py.
        Does not touch the filesystem outside of cache_tmpdir.
        """
        beagle_path = Path(beagle_path)
        freqs_path  = Path(freqs_path)
        samples_path = Path(samples_path)

        if sample_ids is None:
            idx_map = load_sample_index(samples_path)
            sample_ids = sorted(idx_map, key=idx_map.get)

        wall_t0 = time.time()

        with BeagleSubsetCache(beagle_path, freqs_path,
                               tmpdir=cache_tmpdir,
                               chrom_label=chrom_id) as cache:
            n_total = cache.n_total
            pairs = self._enumerate_pairs(sample_ids)

            # ---- INITIAL PASS -------------------------------------------------
            initial_budget_for: Dict[PairKey, int] = {}
            initial_tasks: List[Tuple] = []
            for a, b, prior in pairs:
                budget = self._budget_for(prior, chrom_id, n_total)
                initial_budget_for[(a, b)] = budget
                sub_b, sub_f = cache.get(budget)
                initial_tasks.append((
                    self.binary_path, sub_b, sub_f, samples_path,
                    a, b, n_samples, self.config.threads_per_pair,
                    self.per_pair_timeout,
                ))

            initial_results = self._dispatch(initial_tasks)

            # ---- ESCALATION DECISION -----------------------------------------
            escalation_decisions: Dict[PairKey, EscalationDecision] = {}
            escalated_keys: List[PairKey] = []
            for a, b, prior in pairs:
                key = (a, b)
                result = initial_results.get(key)
                if isinstance(result, Exception) or result is None:
                    escalation_decisions[key] = EscalationDecision(
                        False, "initial_run_failed", float("nan"),
                        EdgeClass.UNRELATED)
                    continue
                decision = should_escalate(result.row, prior, config=self.config)
                escalation_decisions[key] = decision
                if decision.should_escalate:
                    escalated_keys.append(key)

            # ---- ESCALATION PASS ----------------------------------------------
            final_results: Dict[PairKey, PairResult] = {}
            final_budget_for: Dict[PairKey, int] = dict(initial_budget_for)
            for a, b, prior in pairs:
                key = (a, b)
                result = initial_results.get(key)
                if not isinstance(result, Exception) and result is not None:
                    final_results[key] = result

            escalation_tasks: List[Tuple] = []
            escalation_lookup: Dict[PairKey, int] = {}   # pair → escalated budget
            for key in escalated_keys:
                a, b = key
                init_b = initial_budget_for[key]
                new_b = min(self.config.escalation.escalated_budget(init_b), n_total)
                if new_b <= init_b:
                    # Already at the chromosome's ceiling — no point re-running.
                    continue
                sub_b, sub_f = cache.get(new_b)
                escalation_lookup[key] = new_b
                final_budget_for[key] = new_b
                escalation_tasks.append((
                    self.binary_path, sub_b, sub_f, samples_path,
                    a, b, n_samples, self.config.threads_per_pair,
                    self.per_pair_timeout,
                ))

            if escalation_tasks:
                escalation_results = self._dispatch(escalation_tasks)
                for key, payload in escalation_results.items():
                    if not isinstance(payload, Exception):
                        final_results[key] = payload
                    else:
                        # Keep the initial result; flag in the manifest below.
                        pass

            # ---- ASSEMBLE OUTPUTS ---------------------------------------------
            # Canonical order: sort pairs lexicographically.
            sorted_pair_keys = sorted(final_results.keys())

            res_header: List[str] = []
            res_rows: List[Dict[str, str]] = []
            for key in sorted_pair_keys:
                r = final_results[key]
                if not res_header:
                    res_header = list(r.header)
                res_rows.append(r.row)

            manifest_rows: List[ManifestRow] = []
            for a, b, prior in pairs:
                key = (a, b)
                init = initial_results.get(key)
                final = final_results.get(key)
                decision = escalation_decisions.get(
                    key,
                    EscalationDecision(False, REASON_NONE, float("nan"),
                                       EdgeClass.UNRELATED))
                init_row = init.row if isinstance(init, PairResult) else {}
                final_row = final.row if final is not None else {}
                manifest_rows.append(ManifestRow(
                    sample_a=a, sample_b=b,
                    genome_wide_class=prior.value,
                    initial_budget=initial_budget_for[key],
                    initial_sites_used=getattr(init, "sites_used", 0) if isinstance(init, PairResult) else 0,
                    initial_chrom_class=(
                        decision.chrom_class.value if isinstance(init, PairResult)
                        else "initial_run_failed"),
                    initial_theta=_safe_float(init_row.get("theta")),
                    initial_IBS0=_safe_float(init_row.get("IBS0")),
                    initial_KING=_safe_float(init_row.get("KING")),
                    escalated=key in escalation_lookup,
                    escalation_reason=decision.reason,
                    final_budget=final_budget_for[key],
                    final_sites_used=getattr(final, "sites_used", 0) if final else 0,
                    final_chrom_class=(
                        # Recompute on final row in case the escalated run
                        # changed the per-chrom class.
                        _recompute_chrom_class(final.row).value
                        if final else "missing"),
                    final_theta=_safe_float(final_row.get("theta")),
                    final_IBS0=_safe_float(final_row.get("IBS0")),
                    final_KING=_safe_float(final_row.get("KING")),
                    elapsed_seconds=(
                        (init.elapsed_seconds if isinstance(init, PairResult) else 0.0)
                        + (final.elapsed_seconds if (final and final is not init) else 0.0)),
                ))
            manifest_rows.sort(key=lambda r: (r.sample_a, r.sample_b))

            n_pairs_expected = len(sample_ids) * (len(sample_ids) - 1) // 2
            n_pairs_actual = len(res_rows)
            failed_pairs = [
                {"pair": list(k), "error": str(v)}
                for k, v in initial_results.items()
                if isinstance(v, Exception)
            ]

            run_manifest = {
                "schema": "ngsrelate_adaptive.run_manifest.v1",
                "chromosome": chrom_id,
                "n_samples": len(sample_ids),
                "n_pairs_expected": n_pairs_expected,
                "n_pairs_completed": n_pairs_actual,
                "n_pairs_escalated": len(escalation_lookup),
                "n_pairs_failed": len(failed_pairs),
                "chrom_total_sites": n_total,
                "binary_path": str(self.binary_path),
                "config": _config_snapshot(self.config),
                "budgets_built": cache.budgets_built(),
                "elapsed_seconds": round(time.time() - wall_t0, 3),
                "completed_at": datetime.now(tz=timezone.utc).isoformat(),
                "failed_pairs": failed_pairs,
            }

        return ChromosomeResult(
            chrom_id=chrom_id,
            res_header=res_header,
            res_rows=res_rows,
            manifest=manifest_rows,
            run_manifest=run_manifest,
        )


def _recompute_chrom_class(row: Dict[str, str]) -> EdgeClass:
    # Lazy import to avoid circular dependency at module load time.
    from .confidence import chrom_class_from_res_row
    return chrom_class_from_res_row(row)


def _config_snapshot(c: cfg.SchedulerConfig) -> Dict[str, Any]:
    """Serialize the config to a dict for the run manifest."""
    return {
        "budgets": {
            "low":       c.budgets.budget_low,
            "med":       c.budgets.budget_med,
            "high":      c.budgets.budget_high,
            "ambiguous": c.budgets.budget_ambiguous,
            "duplicate": c.budgets.budget_duplicate,
            "floor_fraction_of_available": c.budgets.budget_floor_fraction_of_available,
        },
        "escalation": {
            "factor": c.escalation.factor,
            "softer_rule": c.escalation.softer_rule,
            "softer_multiplier": c.escalation.softer_multiplier,
            "softer_floor_addition": c.escalation.softer_floor_addition,
            "boundary_epsilon_theta": c.escalation.boundary_epsilon_theta,
        },
        "inversion_boost": {
            "enabled": c.inversion_boost.enabled,
            "boosted_chroms": list(c.inversion_boost.boosted_chroms),
        },
        "ambiguous_first_degree_ibs0_band": list(c.ambiguous_first_degree_ibs0_band),
        "n_workers": c.n_workers,
        "threads_per_pair": c.threads_per_pair,
        "allow_mismatched_anchor": c.allow_mismatched_anchor,
    }
