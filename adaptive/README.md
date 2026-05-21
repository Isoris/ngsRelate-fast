# adaptive/ — Per-Pair Adaptive Site-Budgeting for Per-Chromosome ngsRelate

**Status:** documentation-only. Implementation has not begun.

This subdirectory will eventually hold a Python package that wraps
`ngsRelate-fast` to do **per-pair adaptive site-budgeting** for the
per-chromosome ngsRelate runs that feed ngsPedigree Stage 2. It is an
optimization on top of the existing per-chrom workflow (`scripts/STEP_A07b_*`),
not a replacement for the genome-wide ngsRelate run.

## Current contents

- `docs/SPEC_v0.1_CLARIFICATION_NOTE.md` — audit-equivalent note clarifying
  that the per-pair prior is derived from the genome-wide `.res` directly
  (via KING thresholds), not from ngsPedigree Stage 1 output. Resolves the
  ordering ambiguity in the v0.1 spec.
- `docs/IMPLEMENTATION_PLAN.md` — the build plan: 9 phases, what to test,
  what NOT to silently decide, calibration workflow, definition of done.
- `docs/DOWNSTREAM_CONSUMERS.md` — the chain of work that depends on the
  per-chromosome `.res` (Stage 2 → Stage 3 → HPP → KBC → manuscript). Read
  this before relaxing any format-preservation constraints in the plan.

## Originating spec

`SPEC_NGSRELATE_ADAPTIVE_v0_1.md` (stored separately, not in this repo —
it's the working spec marked SPEC ONLY). Read the spec, then the
clarification note, then the implementation plan, in that order.

## What this is NOT

- Not a replacement for `ngsRelate-fast` (the binary).
- Not a replacement for the genome-wide ngsRelate run.
- Not a change to ngsPedigree (which consumes per-chrom `.res` unchanged).
- Not a change to the `.res` file format (sacred — 23 columns, identical
  to upstream).

## When implementation starts

Follow `docs/IMPLEMENTATION_PLAN.md` phase by phase. Do not skip Phase 8
(validation) or Phase 9 (calibration). Defaults in any future
`scheduler/config.py` must be marked `# TODO calibrate` until the
calibration log says otherwise.
