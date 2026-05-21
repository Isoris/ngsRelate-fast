"""
confidence.py — Phase 4: decide whether a per-chrom result needs escalation.

Two responsibilities:

1. `chrom_class_from_res_row(row)` — apply the same KING thresholds as the
   genome-wide prior derivation, but to a per-chrom result row. This is
   the per-chrom classification used for "interesting disagreement"
   comparisons against the prior.

2. `should_escalate(chrom_class, theta, prior_class, config)` — implement
   SPEC §4.3: a pair is escalated to a larger budget if it is at a class
   boundary OR if the prior-vs-per-chrom transition is in
   INTERESTING_TRANSITIONS.

INTERESTING_TRANSITIONS encodes the SPEC §4.3 rule as data. Changing this
set is a methodological choice (it shifts the bar for re-running pairs at
a larger budget) — it is NOT a bug fix. Discuss before touching.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

from .edge_class import EdgeClass
from .prior import _classify_one, _to_float
from . import config as cfg


# ---- Interesting transitions table (SPEC §4.3) ----------------------------

INTERESTING_TRANSITIONS = frozenset([
    (EdgeClass.UNRELATED,        EdgeClass.SECOND_DEGREE),
    (EdgeClass.UNRELATED,        EdgeClass.FULL_SIBLING),
    (EdgeClass.UNRELATED,        EdgeClass.PARENT_OFFSPRING),
    (EdgeClass.FULL_SIBLING,     EdgeClass.UNRELATED),
    (EdgeClass.FULL_SIBLING,     EdgeClass.THIRD_DEGREE),
    (EdgeClass.PARENT_OFFSPRING, EdgeClass.UNRELATED),
    (EdgeClass.PARENT_OFFSPRING, EdgeClass.THIRD_DEGREE),
    (EdgeClass.PARENT_OFFSPRING, EdgeClass.SECOND_DEGREE),
    (EdgeClass.PARENT_OFFSPRING, EdgeClass.FULL_SIBLING),
])
# AMBIGUOUS_FIRST_DEGREE prior is handled specially: any per-chrom class
# except itself counts as interesting (see is_interesting_transition).


# ---- Escalation reasons ---------------------------------------------------

REASON_NONE                    = "none"
REASON_AT_BOUNDARY             = "at_boundary"
REASON_INTERESTING_DISAGREEMENT = "interesting_disagreement"


# ---- KING thresholds used for boundary distance ---------------------------

_KING_THRESHOLDS = (
    cfg.KING_THETA_DUPLICATE,
    cfg.KING_THETA_FIRST_DEGREE,
    cfg.KING_THETA_SECOND_DEG,
    cfg.KING_THETA_THIRD_DEG,
)


def distance_to_nearest_boundary(theta: float) -> float:
    """Min absolute distance from `theta` to any KING-class boundary."""
    return min(abs(theta - t) for t in _KING_THRESHOLDS)


def is_at_boundary(theta: float, epsilon: float) -> bool:
    return distance_to_nearest_boundary(theta) <= epsilon


def is_interesting_transition(
    prior_class: EdgeClass,
    chrom_class: EdgeClass,
) -> bool:
    """True if (prior → per-chrom) is an interesting disagreement per SPEC §4.3."""
    if prior_class == EdgeClass.AMBIGUOUS_FIRST_DEGREE:
        return chrom_class != EdgeClass.AMBIGUOUS_FIRST_DEGREE
    return (prior_class, chrom_class) in INTERESTING_TRANSITIONS


# ---- Public API ----------------------------------------------------------

def chrom_class_from_res_row(
    row: Dict[str, str],
    *,
    ibs0_ambiguous_band: tuple = cfg.AMBIGUOUS_FIRST_DEGREE_IBS0_BAND,
) -> EdgeClass:
    """Apply KING thresholds to a per-chrom .res row, returning EdgeClass.

    Mirrors prior._classify_one but consumes the dict-shaped row returned
    by the binary. Missing IBS0 in the first-degree band falls back to
    AMBIGUOUS_FIRST_DEGREE (same rule as the prior path).
    """
    theta = _to_float(row.get("theta", "nan"))
    ibs0  = _to_float(row.get("IBS0",  "nan")) if "IBS0" in row else None
    return _classify_one(theta, ibs0, ibs0_ambiguous_band)


@dataclass(frozen=True)
class EscalationDecision:
    should_escalate: bool
    reason: str            # one of REASON_*
    boundary_distance: float
    chrom_class: EdgeClass

    def as_tuple(self) -> Tuple[bool, str]:
        return self.should_escalate, self.reason


def should_escalate(
    chrom_row: Dict[str, str],
    prior_class: EdgeClass,
    *,
    config: Optional[cfg.SchedulerConfig] = None,
) -> EscalationDecision:
    """Decide whether a per-chrom result needs a larger-budget re-run.

    Order of precedence: interesting_disagreement > at_boundary > none.
    A pair flagged as interesting_disagreement is also reported even if it
    happens to land near a boundary — the disagreement is the louder signal.
    """
    conf = config or cfg.SchedulerConfig()
    epsilon = conf.escalation.boundary_epsilon_theta

    chrom_class = chrom_class_from_res_row(chrom_row)
    theta = _to_float(chrom_row.get("theta", "nan"))
    dist  = distance_to_nearest_boundary(theta) if theta == theta else float("inf")

    interesting = is_interesting_transition(prior_class, chrom_class)
    at_boundary = dist <= epsilon

    if interesting:
        return EscalationDecision(True, REASON_INTERESTING_DISAGREEMENT,
                                  dist, chrom_class)
    if at_boundary:
        return EscalationDecision(True, REASON_AT_BOUNDARY, dist, chrom_class)
    return EscalationDecision(False, REASON_NONE, dist, chrom_class)
