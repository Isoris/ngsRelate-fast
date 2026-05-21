# CALIBRATION_LOG.md — Phase 9 record + locked defaults

**Status:** awaiting first real-binary calibration run on LG12.

This log records the calibration procedure (per IMPLEMENTATION_PLAN.md
§2 Phase 9), the values tested at each iteration, and — once Tier-1
gates 4 and 6 and Tier-2 gate 1 pass — the locked defaults that
`scheduler/config.py` will be updated to.

Until the calibration log declares values "LOCKED", every default in
`config.py` is marked `# TODO calibrate (SPEC §7 OQX)` and is a
placeholder, not a recommendation.

---

## Procedure (canonical, do not deviate)

1. Pick LG12 (moderate size; carries first-degree pairs per SPEC §7 OQ1).
2. Run uniform-stride baseline:
   ```bash
   bash scripts/STEP_A07b_fast_relatedness_per_chrom.sh
   # one task on LG12 with -D 100000
   ```
3. Run adaptive with current `config.py` defaults:
   ```bash
   python -m adaptive.scheduler --chrom C_gar_LG12 ...
   ```
4. Diff:
   ```bash
   python -m adaptive.scheduler.validate \
       --baseline-res          .../baseline/relatedness.res \
       --adaptive-res          .../adaptive/relatedness.res \
       --adaptive-manifest     .../adaptive/relatedness.res.adaptive_manifest.tsv \
       --adaptive-run-manifest .../adaptive/relatedness.res.adaptive_run_manifest.json \
       --genome-wide-res       .../genomewide.res \
       --report-json           .../validation_report.json
   ```
5. Inspect `validation_report.json`. If Gate 4 or Gate 6 (tier-1) fail,
   first investigate **which pairs** drifted and **what direction** — a
   PO-set difference of even one dyad must be explained, not knob-tuned
   away. See DOWNSTREAM_CONSUMERS.md for why.
6. If Gate 1 (tier-2) is below 90%, adjust `BUDGET_LOW` upward first
   (most pairs are unrelated and use this); only adjust `BUDGET_HIGH`
   if the failing pairs are concentrated in the first-degree band.
7. If escalation rate exceeds ~10% of pairs, tighten
   `boundary_epsilon_theta` (smaller window). If escalation is missing
   pairs that drifted across boundaries, widen it.
8. Iterate until Tier-1 gates 4 and 6 and Tier-2 gate 1 pass on LG12.
9. Apply the locked config to all 28 chromosomes. Re-run validation per
   chromosome.
10. Record the final values in the "Locked defaults" section below and
    update `scheduler/config.py` to match (replacing the
    `# TODO calibrate` comments with `# LOCKED <date>`).

---

## Iteration log

### Iteration 0 — Initial placeholder values (no run yet)

| Knob                     | Value     | Source                                 |
|--------------------------|-----------|----------------------------------------|
| `budget_low`             | 3,000     | TODO calibrate (SPEC §7 OQ1)           |
| `budget_med`             | 10,000    | TODO calibrate (SPEC §7 OQ1)           |
| `budget_high`            | 30,000    | TODO calibrate (SPEC §7 OQ1)           |
| `budget_ambiguous`       | 30,000    | TODO calibrate (SPEC §7 OQ1)           |
| `budget_duplicate`       | 3,000     | TODO calibrate (SPEC §7 OQ1)           |
| `escalation.factor`      | 3         | TODO calibrate (SPEC §7 OQ2; literal)  |
| `escalation.softer_rule` | false     | Default; flag exists per SPEC §7 OQ2   |
| `boundary_epsilon_theta` | 0.01      | TODO calibrate (SPEC §7 OQ4)           |
| `floor_fraction_of_available` | 0.5  | TODO calibrate (SPEC §7 OQ3)           |

These are placeholders for the first real-binary run. **Do not cite
these values as the calibrated defaults.**

### Iteration N (template; copy when running)

- **Date:**
- **LG12 baseline elapsed:**
- **LG12 adaptive elapsed:**
- **Gate 1 fraction agree:**
- **Gate 4 PO set equality:**
- **Gate 4 FS Jaccard:**
- **Gate 6 interesting recall:**
- **Gate 6 PO→unrelated lost:**
- **Knobs changed since prior iteration:**
- **Decision:**

---

## Locked defaults

**Status: NOT LOCKED.** Awaiting Phase 9 completion.

Once locked, this section will list each knob, its locked value, the
LG12 iteration that produced it, and the validation report path that
attests to it. After locking, `scheduler/config.py` will be updated to
match, the `# TODO calibrate` comments will be replaced with
`# LOCKED <date>`, and `CHANGELOG.md` will get a `v1.0.0-adaptive` entry.

---

## Per-chromosome sweep (post-lock)

A table populated after the locked config is applied to all 28
chromosomes. Each row records the per-chrom validation report path and
the verdicts for the six gates. Any tier-1 gate failure on any
chromosome blocks release and triggers re-calibration.

| Chromosome | Elapsed (s) | Gate 1 | Gate 4 | Gate 6 | Report |
|------------|-------------|--------|--------|--------|--------|
| LG01       |             |        |        |        |        |
| LG02       |             |        |        |        |        |
| ...        |             |        |        |        |        |
| LG28       |             |        |        |        |        |
