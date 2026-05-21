"""
prior.py — Phase 1: derive per-pair prior class from the genome-wide .res.

Per CLARIFICATION_NOTE §1, the per-pair prior used by the adaptive
scheduler is computed directly from the genome-wide ngsRelate-fast .res
using KING-robust thresholds. It is NOT downstream of ngsPedigree Stage 1.

Public API:

    derive_priors(res_path, *, config=None) -> dict[(sample_a, sample_b), EdgeClass]

Pair keys are tuples of sample IDs sorted lexicographically — this is the
canonical ordering used throughout the package.
"""

from __future__ import annotations
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

from .edge_class import EdgeClass
from . import config as cfg


PairKey = Tuple[str, str]


class PriorDerivationError(Exception):
    """Raised when the genome-wide .res cannot be parsed into priors."""


def _canonical_pair(a: str, b: str) -> PairKey:
    return (a, b) if a <= b else (b, a)


def _classify_one(theta: float, ibs0: Optional[float],
                  ibs0_ambiguous_band: tuple) -> EdgeClass:
    """Apply the KING thresholds + IBS0 PO/FS split (clarification note §1)."""
    if theta is None or math.isnan(theta):
        # Missing θ → treat as unrelated; the scheduler will run with
        # BUDGET_LOW and the per-chrom result will speak for itself.
        return EdgeClass.UNRELATED

    if theta >= cfg.KING_THETA_DUPLICATE:
        return EdgeClass.DUPLICATE_OR_CLONE

    if theta >= cfg.KING_THETA_FIRST_DEGREE:
        # First-degree band: split PO vs FS using IBS0 (clarification §1.1).
        if ibs0 is None or math.isnan(ibs0):
            return EdgeClass.AMBIGUOUS_FIRST_DEGREE
        lo, hi = ibs0_ambiguous_band
        if lo <= ibs0 <= hi:
            return EdgeClass.AMBIGUOUS_FIRST_DEGREE
        return EdgeClass.PARENT_OFFSPRING if ibs0 <= cfg.IBS0_PO_MAX else EdgeClass.FULL_SIBLING

    if theta >= cfg.KING_THETA_SECOND_DEG:
        return EdgeClass.SECOND_DEGREE

    if theta >= cfg.KING_THETA_THIRD_DEG:
        return EdgeClass.THIRD_DEGREE

    return EdgeClass.UNRELATED


def _to_float(s: str) -> float:
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


def _iter_res_rows(res_path: Path) -> Iterable[dict]:
    """Yield each row as a dict keyed by header column name.

    The .res may have either (a, b, ida, idb, ...) when -z was used, or
    just (a, b, ...) otherwise. Sample IDs come from ida/idb; if those
    are absent we fall back to the integer indices in a/b, which is
    enough to derive priors but means downstream pair-key matching
    against per-chrom runs must use the same indexing.
    """
    with open(res_path) as fh:
        header_line = fh.readline()
        if not header_line:
            raise PriorDerivationError(
                f"genome-wide .res is empty: {res_path}")
        header = header_line.rstrip("\n").split("\t")
        has_ids = "ida" in header and "idb" in header
        if "theta" not in header:
            raise PriorDerivationError(
                f"genome-wide .res missing 'theta' column: {res_path} "
                f"(header={header})")
        for lineno, line in enumerate(fh, start=2):
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                raise PriorDerivationError(
                    f"malformed row at line {lineno} of {res_path}: "
                    f"got {len(parts)} columns, expected {len(header)}")
            row = dict(zip(header, parts))
            if has_ids:
                row["_sample_a"] = row["ida"]
                row["_sample_b"] = row["idb"]
            else:
                # Fall back to integer indices as string sample IDs.
                row["_sample_a"] = row["a"]
                row["_sample_b"] = row["b"]
            yield row


def derive_priors(
    res_path,
    *,
    config: Optional[cfg.SchedulerConfig] = None,
) -> Dict[PairKey, EdgeClass]:
    """Read the genome-wide .res and return a per-pair prior class map.

    Args:
        res_path: path to the genome-wide ngsRelate-fast .res file.
        config:   optional SchedulerConfig (for the ambiguous-first-degree
                  IBS0 band). Uses defaults if not supplied.

    Returns:
        dict keyed by sorted-tuple (sample_a, sample_b), value is EdgeClass.

    Raises:
        PriorDerivationError: malformed input.
        FileNotFoundError:    res_path does not exist.
    """
    res_path = Path(res_path)
    if not res_path.exists():
        raise FileNotFoundError(f"genome-wide .res not found: {res_path}")

    conf = config or cfg.SchedulerConfig()
    band = conf.ambiguous_first_degree_ibs0_band

    out: Dict[PairKey, EdgeClass] = {}
    for row in _iter_res_rows(res_path):
        a, b = row["_sample_a"], row["_sample_b"]
        theta = _to_float(row.get("theta", "nan"))
        ibs0  = _to_float(row.get("IBS0", "nan")) if "IBS0" in row else None
        klass = _classify_one(theta, ibs0, band)
        out[_canonical_pair(a, b)] = klass

    return out


def _main(argv=None):
    """CLI: print class counts for a genome-wide .res. Independently useful."""
    import argparse
    from collections import Counter

    ap = argparse.ArgumentParser(
        prog="python -m adaptive.scheduler.prior",
        description="Derive KING-class priors from a genome-wide ngsRelate-fast .res "
                    "and print per-class counts.")
    ap.add_argument("res", help="path to genome-wide .res file")
    args = ap.parse_args(argv)

    try:
        priors = derive_priors(args.res)
    except (PriorDerivationError, FileNotFoundError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    counts = Counter(priors.values())
    print(f"# {len(priors)} pairs from {args.res}")
    for klass in EdgeClass:
        n = counts.get(klass, 0)
        print(f"{klass.value:<25s}\t{n}")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
