"""
config.py — defaults for the adaptive scheduler.

**Calibration status:** every value marked `# TODO calibrate (SPEC §7 OQX)`
is a placeholder. Defaults are LOCKED only after Phase 9 (calibration)
runs the validation tool on LG12, all six gates pass, and the locked
values are recorded in adaptive/docs/CALIBRATION_LOG.md.

Never quietly change one of these defaults during implementation. Surface
the question — that is the point of the open-question handling in
IMPLEMENTATION_PLAN.md §3.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional

from .edge_class import EdgeClass


# ---------------------------------------------------------------------------
# KING-robust thresholds on the genome-wide kinship coefficient θ
# (Manichaikul et al. 2010). Per CLARIFICATION_NOTE §1 — these are the
# exact class boundaries used to derive the per-pair prior class from the
# genome-wide .res. Do NOT relax these without a methodological discussion.
# ---------------------------------------------------------------------------

KING_THETA_DUPLICATE    = 0.354
KING_THETA_FIRST_DEGREE = 0.177
KING_THETA_SECOND_DEG   = 0.0884
KING_THETA_THIRD_DEG    = 0.0442

# IBS0 split inside the first-degree band (CLARIFICATION_NOTE §1.1).
IBS0_PO_MAX = 0.008

# Ambiguous first-degree band: IBS0 within ±20% of the PO/FS threshold.
AMBIGUOUS_FIRST_DEGREE_IBS0_BAND = (0.0064, 0.0096)


# ---------------------------------------------------------------------------
# Budget defaults per prior class.
# Initial values chosen for the first calibration pass; revisit on LG12.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BudgetConfig:
    # Per-pair site budget per prior class.
    budget_low:       int = 3_000    # TODO calibrate (SPEC §7 OQ1) — unrelated / third_degree
    budget_med:       int = 10_000   # TODO calibrate (SPEC §7 OQ1) — second_degree
    budget_high:      int = 30_000   # TODO calibrate (SPEC §7 OQ1) — first_degree (PO / FS)
    budget_ambiguous: int = 30_000   # TODO calibrate (SPEC §7 OQ1) — ambiguous_first_degree
    budget_duplicate: int = 3_000    # TODO calibrate (SPEC §7 OQ1) — duplicate_or_clone (already confirmed)

    # Floor: never use fewer than this fraction of available sites on chrom
    # (SPEC §7 OQ3). Default 0.5 means "if a budget tier requests <50% of
    # the chromosome's available sites, that's fine; if it requests more,
    # we cap at the available sites."  This knob exists mainly so small
    # chromosomes don't get accidentally over-budgeted.
    budget_floor_fraction_of_available: float = 0.5  # TODO calibrate (SPEC §7 OQ3)

    def budget_for(self, prior_class: EdgeClass) -> int:
        if prior_class == EdgeClass.DUPLICATE_OR_CLONE:
            return self.budget_duplicate
        if prior_class == EdgeClass.PARENT_OFFSPRING:
            return self.budget_high
        if prior_class == EdgeClass.FULL_SIBLING:
            return self.budget_high
        if prior_class == EdgeClass.AMBIGUOUS_FIRST_DEGREE:
            return self.budget_ambiguous
        if prior_class == EdgeClass.SECOND_DEGREE:
            return self.budget_med
        if prior_class == EdgeClass.THIRD_DEGREE:
            return self.budget_low
        if prior_class == EdgeClass.UNRELATED:
            return self.budget_low
        raise ValueError(f"unknown prior class: {prior_class!r}")


# ---------------------------------------------------------------------------
# Escalation rules (SPEC §4.3).
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EscalationConfig:
    # SPEC §4.3 literal: escalation budget = factor * initial_budget.
    factor: int = 3   # TODO calibrate (SPEC §7 OQ2)

    # SPEC §7 OQ2 softer alternative: max(B*2, B+5000).
    # Selected by setting softer_rule=True.
    softer_rule: bool = False
    softer_multiplier: int = 2
    softer_floor_addition: int = 5_000

    # Per-class boundary epsilon: a pair is "at the boundary" if its
    # per-chrom θ lands within ±epsilon of the nearest KING threshold.
    # Initial values per spec, marked TODO — refine from histograms in Phase 9.
    boundary_epsilon_theta: float = 0.01  # TODO calibrate (SPEC §7 OQ4)

    def escalated_budget(self, initial_budget: int) -> int:
        if self.softer_rule:
            return max(initial_budget * self.softer_multiplier,
                       initial_budget + self.softer_floor_addition)
        return initial_budget * self.factor


# ---------------------------------------------------------------------------
# Inversion-aware boost (SPEC §7 OQ8).
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class InversionBoostConfig:
    enabled: bool = False
    # List of chromosome IDs (matching BEAGLE site-id chromosomes) whose
    # BUDGET_HIGH gets doubled. Populated on the CLI side from
    # --boost-inversion-chroms or kept empty.
    boosted_chroms: tuple = ()


# ---------------------------------------------------------------------------
# Top-level scheduler config.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SchedulerConfig:
    budgets: BudgetConfig = field(default_factory=BudgetConfig)
    escalation: EscalationConfig = field(default_factory=EscalationConfig)
    inversion_boost: InversionBoostConfig = field(default_factory=InversionBoostConfig)
    ambiguous_first_degree_ibs0_band: tuple = AMBIGUOUS_FIRST_DEGREE_IBS0_BAND

    # Worker pool size for per-pair runs. None → use os.cpu_count().
    n_workers: Optional[int] = None

    # Per-worker thread count for the ngsRelate-fast binary itself.
    # Default 1 — the parallelism is across pairs, not within one pair.
    threads_per_pair: int = 1

    # Strict anchor staleness check (CLARIFICATION_NOTE §3).
    allow_mismatched_anchor: bool = False
