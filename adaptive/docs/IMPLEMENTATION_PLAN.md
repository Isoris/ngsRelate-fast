# NGSRELATE_ADAPTIVE — Implementation Plan for Claude Code

**Plan version:** 1.0
**Spec basis:** `SPEC_NGSRELATE_ADAPTIVE_v0_1.md` + `SPEC_v0.1_CLARIFICATION_NOTE.md`
**Target repo:** `ngsRelate-fast/`, new subdirectory `adaptive/`
**Target consumer:** ngsPedigree Stage 2 (no changes to it; this is purely an upstream optimization)
**Status when started:** implementation begins after this plan is read and approved

---

## 0. Read these first (in this order)

1. `CLAUDE.md` (repo root) — orientation, what is load-bearing
2. `README.md` — what `ngsRelate-fast` is and isn't
3. `SPEC_NGSRELATE_ADAPTIVE_v0_1.md` (uploaded separately) — the actual spec
4. `adaptive/docs/SPEC_v0.1_CLARIFICATION_NOTE.md` — fixes the prior-source ambiguity
5. `contracts/ngsrelate_fast.input.v1.schema.json` and `…output.v1.schema.json` — the JSON contract pattern this work extends
6. **This file.**

Then start.

---

## 1. What you are building

A Python package `adaptive/` that wraps the `ngsRelate-fast` binary and implements per-pair adaptive site-budgeting for **per-chromosome** ngsRelate runs.

**Inputs to the package:**
- Per-chromosome BEAGLE files (one per chromosome, already exists from upstream pipeline)
- The **genome-wide `.res`** produced by `ngsRelate-fast` (this is the prior source — see clarification note §1)
- Sample list (one ID per line)

**Outputs from the package:**
- Per-chromosome `.res` files in the **exact format ngsPedigree Stage 2 expects** (23 columns, identical to non-adaptive `STEP_A07b` output)
- Per-chromosome `.adaptive_manifest.tsv` (audit sidecar, NOT consumed by ngsPedigree)
- Per-chromosome `.adaptive_run_manifest.json` (run metadata, extends the existing `output.v1` contract pattern)

**Critical constraint:** the `.res` schema is sacred. ngsPedigree Stage 2 consumes it byte-for-byte. If you find yourself wanting to add columns, **stop** — they go in the sidecar manifest instead.

**Why this constraint is non-negotiable — downstream blast radius:**

The per-chromosome `.res` does not just feed ngsPedigree Stage 2. Stage 2 outputs (confirmed PO dyads/triads) feed Stage 3 (inheritance maps: Gold/Silver/Bronze segments, `parental_hap_inherited`). Stage 3 outputs feed HPP, which projects MODULE_CONSERVATION-scored variants (SnpEff, SIFT4G, VESM, splice) onto inheritance segments. HPP outputs feed KBC cross-checks against arrangement assignments and the manuscript figures.

A semantic drift in Stage 2 caused by `.res` format change or borderline-pair reclassification does not just "make ngsPedigree complain." It silently propagates six steps downstream and the failure surfaces as wrong damaging-variant projections in HPP outputs that nobody notices until the manuscript figures look off. See `adaptive/docs/DOWNSTREAM_CONSUMERS.md` for the full chain.

This is why Gate 4 (Stage 2 semantic compatibility, SPEC §6) is the most consequential of the validation gates — see §2 Phase 8 below for the tightened acceptance criteria.

---

## 2. Build order (do not reorder)

This order is chosen so each phase is independently testable and produces something useful even if the next phase is delayed.

### Phase 1 — Genome-wide prior derivation (~½ day)

**Module:** `adaptive/scheduler/prior.py`

A pure function `derive_priors(genome_wide_res_path) -> dict[(sample_a, sample_b), EdgeClass]`.

- Read the genome-wide `.res` (use existing `scripts/contract_io.py` if `.input.json`/`.output.json` sidecars exist; fall back to plain TSV read if they don't)
- For each pair, apply the KING thresholds from clarification note §1
- Apply the IBS0 split for the first-degree band per clarification note §1.1
- Return a dict keyed by sorted pair tuple

**Tests** (in `adaptive/tests/test_prior.py`):
- Hand-crafted `.res` rows hitting each class boundary
- Hand-crafted ambiguous PO/FS pair (IBS0 near 0.008)
- Empty input → empty output, no crash
- Malformed `.res` → clear error message, not silent miscount

**Deliverable:** a CLI `python -m adaptive.scheduler.prior <res>` that prints class counts. This is independently useful as a sanity check on any genome-wide `.res`.

### Phase 2 — Per-pair BEAGLE subsetting + cache (~½ day)

**Module:** `adaptive/scheduler/subset.py`

A class `BeagleSubsetCache(chrom_beagle_path, tmpdir)`:
- On `.get(budget)` returns a path to a temporary BEAGLE containing exactly `budget` sites sampled by deterministic stride from `chrom_beagle_path`
- Caches by `budget` — second call with same budget returns cached path
- Cleanup method removes all temp files

**Important details:**
- Deterministic stride matches `ngsRelate-fast` `-D` convention (see SPEC §4.2): `stride = floor(N_total / budget)`, kept indices `[0, stride, 2*stride, ...]`
- If `budget >= N_total`, return the input path itself (no subsetting needed)
- Memory-map the source BEAGLE; do not load the whole file into Python memory
- Preserve BEAGLE gzip compression in the output

**Tests** (`adaptive/tests/test_subset.py`):
- Tiny BEAGLE (~100 sites), budget=10, verify stride correctness
- Budget > N_total → returns original path
- Same budget called twice → returns same cached path, no second write
- Cleanup removes all temp files

**Deliverable:** standalone module, importable, with unit tests passing.

### Phase 3 — Pair runner (single pair, no parallelism) (~1 day)

**Module:** `adaptive/scheduler/runner.py`

Function `run_pair_on_chrom(binary_path, beagle_path, freqs_path, sample_a, sample_b, n_samples, threads=1) -> ngsRelateResRow`:

- Invokes `ngsRelate-fast` with `-a <idx_a> -b <idx_b>` to run just that one pair
- Passes `-D 0` to disable the binary's built-in downsampling (the scheduler is doing the subsetting itself — see SPEC §7 OQ7)
- Captures the `.res` row, returns it as a dict matching the 23-column schema
- Captures stderr for diagnostics

**Open question to verify in this phase:** does our fork preserve the upstream `-a`/`-b` flags? It should — we didn't touch that path — but the very first thing this module does on the first real run is verify the flags work. Add a startup sanity check.

**Tests** (`adaptive/tests/test_runner.py`):
- Mock the binary (use a fake script that emits a canned `.res` row) and verify parsing
- Real-binary test marked `@pytest.mark.requires_binary`, skipped in CI without binary

**Deliverable:** working single-pair runner against a real binary (manual test on LANTA or laptop with a tiny BEAGLE).

### Phase 4 — Confidence scoring (~½ day)

**Module:** `adaptive/scheduler/confidence.py`

Two functions:

- `chrom_class_from_res_row(row) -> EdgeClass` — same KING thresholds as Phase 1, applied to per-chrom result
- `should_escalate(chrom_row, prior_class, config) -> tuple[bool, str]` — implements SPEC §4.3 escalation triggers; returns `(should_escalate, reason_string)` where `reason ∈ {at_boundary, interesting_disagreement, none}`

**Implementation of "interesting disagreement"** — make this a precise predicate, not prose:

```python
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
# AMBIGUOUS_FIRST_DEGREE prior → any per-chrom class except itself is interesting
```

This is the SPEC §4.3 "interesting disagreement" rule converted to data. Document in the module docstring that any change here is a methodological choice, not a bug fix.

**Tests:** boundary-distance correctness, transition lookup correctness, exhaustive coverage of the 7×7 transition table (49 cases).

### Phase 5 — Scheduler (the main loop) (~1 day)

**Module:** `adaptive/scheduler/scheduler.py`

Class `AdaptiveScheduler(config, prior_map, binary_path, n_threads)`:

- Method `run_chromosome(chrom_id, beagle_path, freqs_path, samples) -> tuple[ResFile, ManifestRows, RunManifest]`
- Internally: builds a `BeagleSubsetCache`, computes the budget for every pair, runs them through a multiprocessing pool, collects results, applies escalation, runs escalations through the pool again, assembles outputs

**Concurrency model:** `multiprocessing.Pool` over pairs. The cache is built once before the pool starts, then the pool workers each open the cached subset BEAGLEs read-only. No shared writable state.

**Order of operations per chromosome:**

1. Build `BeagleSubsetCache(beagle_path)`
2. For each pair: look up prior class → look up budget → get cached subset path
3. Submit all `(pair, subset_path)` tuples to the pool
4. Collect results into a dict keyed by pair
5. For each result: compute confidence, decide escalation
6. For escalated pairs: get larger-budget subset path, re-run through pool
7. Final pass: assemble `.res` in canonical pair order (lexicographic sample IDs)
8. Write `.res`, `.adaptive_manifest.tsv`, `.adaptive_run_manifest.json`
9. Cleanup cache

### Phase 6 — Output writers + manifest (~½ day)

**Module:** `adaptive/scheduler/output.py`

Three writer functions, one per output file. The `.res` writer must produce a file **bit-identical** to what `STEP_A07b_relatedness_per_chrom.sh` produces given the same site set. The `.adaptive_run_manifest.json` reuses the schema pattern from `contracts/`; add a new schema file `contracts/ngsrelate_adaptive.run_manifest.v1.schema.json`.

**Critical:** the `.res` row order must be canonical (sorted by `(sample_a, sample_b)` with both alphabetized) regardless of the order the pool returned results in. Stage 5 of SPEC §6 (Reproducibility gate) depends on this.

### Phase 7 — CLI + driver scripts (~½ day)

**Module:** `adaptive/scheduler/cli.py`

`python -m adaptive.scheduler <args>` — single-chromosome entry point. Flags match SPEC §5.1.

**Scripts:**
- `adaptive/scripts/RUN_NGSRELATE_ADAPTIVE_LOCAL.sh` — laptop driver, 28-chrom loop
- `adaptive/scripts/RUN_NGSRELATE_ADAPTIVE_LANTA.sh` — same but with LANTA paths and 32 threads

Both invoke the same Python CLI. Differ only in path constants and thread count.

### Phase 8 — Validation (~1 day) — DO NOT SKIP

**Module:** `adaptive/scheduler/validate.py`

A standalone tool that takes one chromosome's adaptive `.res` and the corresponding uniform-stride baseline `.res` and runs all six gates from SPEC §6. Outputs a JSON report.

This is the **calibration tool**. The open questions in SPEC §7 (budget defaults, escalation factor, boundary epsilon) cannot be answered without this tool. Phase 8 happens before any default in `config.py` is "locked."

**Gate priority (tightened):** Gates 4 and 6 are tier-1 (hard fails block release); Gate 1 is tier-2 (numeric target); Gates 2, 3, 5 are tier-3 (operational targets).

- **Tier 1 — Gate 4 (Stage 2 semantic compatibility).** A Stage 2 drift cascades through Stage 3 → HPP → KBC and ends up wrong in manuscript figures. Acceptance: run actual ngsPedigree Stage 2 on both adaptive and baseline `.res` and diff:
  - **Confirmed PO dyad set:** must match exactly (set equality). One missing PO dyad means a missing Stage 3 inheritance map, which is a missing HPP projection. Zero tolerance.
  - **Confirmed FS dyad set:** ≥99% Jaccard agreement. <1% drift acceptable only if all drifted pairs are flagged for manual review in both runs.
  - **`frac_disagreement` per pair:** within ±0.05 of baseline (already in SPEC §6).
  - **No new categories of disagreement introduced** (already in SPEC §6).
- **Tier 1 — Gate 6 (interesting-disagreement preservation).** These are the inversion-detection signal. Acceptance: ≥95% recall vs uniform-stride baseline as the spec says, BUT add: zero loss of pairs where the genome-wide → per-chrom transition is `parent_offspring → unrelated` (the PO-on-chrom-X violation case, SPEC §4.3).
- **Tier 2 — Gate 1 (edge-class precision).** ≥90% per spec. Drift below 90% on any chromosome is a hard fail; drift between 90–95% on >5 of 28 chromosomes triggers re-calibration before release.
- **Tier 3 — Gates 2, 3, 5.** Operational targets. Misses are documented in `CALIBRATION_LOG.md` but do not block release if Tiers 1 and 2 pass.

**Calibration workflow** (run in implementation Phase 9, below):
1. Pick LG12 (moderate size, has first-degree pairs per spec §7 OQ1)
2. Run uniform-stride baseline (`-D 100000`) to get reference `.res`
3. Run adaptive with current defaults
4. Run validation tool → see which gates pass, which fail
5. Adjust defaults in `config.py`, re-run
6. Repeat until **Tier-1 Gates 4 and 6 pass and Tier-2 Gate 1 passes on LG12**
7. Apply locked defaults to all 28 chromosomes, re-run validation per chrom
8. Document the locked values + calibration log in `adaptive/docs/CALIBRATION_LOG.md`

### Phase 9 — Calibration (variable, ~1 day if defaults are close, ~3 days if not)

See above. This is where SPEC §7 open questions OQ1, OQ2, OQ4 get answered empirically.

---

## 3. What NOT to silently decide

The spec lists 9 open questions in §7. **None of them should be quietly resolved by picking a value during implementation.** Each gets explicit handling:

| Open Q | Spec text | Implementation handling |
|---|---|---|
| OQ1 | budget defaults | Initial values in `config.py`, marked `# TODO calibrate (SPEC §7 OQ1)`. Locked only after Phase 9. |
| OQ2 | escalation factor too aggressive? | Implement BOTH the literal `factor=3` and the softer `max(B*2, B+5000)` rule. Default to literal. `EscalationConfig.softer_rule` flag selects. |
| OQ3 | per-chromosome budget floor | Implement `BUDGET_FLOOR_FRACTION_OF_AVAILABLE = 0.5` as a config knob. Default 0.5 per spec. Document chrom site-pool size in run manifest. |
| OQ4 | boundary epsilon calibration | Initial values from spec, marked TODO. Phase 9 produces histograms of pair-distance-to-boundary; epsilon set to capture top ~10%. |
| OQ5 | anchor staleness | Implemented per clarification note §3 — sample-set match check + content hashes. Default strict; `--allow-mismatched-anchor` escape hatch. |
| OQ6 | no genome-wide anchor | Resolved by clarification note §2 — genome-wide is a hard prerequisite. Scheduler refuses to run without it. |
| OQ7 | `-D` flag interaction | Scheduler always passes `-D 0` to the binary. Document in help text. |
| OQ8 | inversion-zone awareness | Add `--boost-inversion-chroms` flag that doubles `BUDGET_HIGH` for chromosomes in a config list. Default off. |
| OQ9 | biomod packaging | Defer per spec. |

If during implementation you find yourself wanting to "just pick a reasonable default" for one of these, **stop and surface it.** The point of the open-question handling is that defaults are deliberate, not accidental.

---

## 4. Naming compliance (hard rules)

Per SPEC §9 and per project-wide naming conventions:

- Edge classes: `parent_offspring`, `full_sibling`, `second_degree`, `third_degree`, `unrelated`, `duplicate_or_clone`, `ambiguous_first_degree` — **never** bare `PO`, `FS`, `2nd`, etc.
- "POD-compatible," never "POD found"
- "POD_candidate_variant," never "POD variant"
- `HWE_FIS`, never bare `FIS`
- `arrangement_FST_like`, never bare `FST`
- "adaptive scheduler" or "NGSRELATE_ADAPTIVE" — the module name; **do not** call it "fast ngsRelate v2" or anything that suggests it replaces the fork's binary
- Scheduler output is "per-chromosome adaptive `.res`" — clarify "adaptive" when it might be confused with the uniform-stride per-chrom `.res`

---

## 5. Tests are not optional

Coverage requirements:

- Every module in `adaptive/scheduler/` has a corresponding test file in `adaptive/tests/`
- Phase 1, 2, 4 modules: ≥90% line coverage (pure functions, easy to test)
- Phase 3, 5, 6 modules: integration tests with mocked binary + tiny synthetic BEAGLE
- Phase 8 validation tool itself has tests (a test that runs the validator on hand-crafted adaptive + baseline `.res` files and checks the gate verdicts)

CI in `.github/workflows/test.yml` runs `pytest adaptive/tests/` on every push. Tests requiring the real binary are marked `@pytest.mark.requires_binary` and skipped in CI.

---

## 6. Repo integration

The adaptive code lives in `adaptive/` under the existing `ngsRelate-fast` repo. It does **not** modify:

- The patch in `patches/`
- The binary itself
- The existing `scripts/STEP_A07_*` SLURM scripts
- The contracts in `contracts/` (the adaptive run manifest is a new schema, not a modification of the existing two)
- `contract_io.py` (the adaptive scheduler may import it, but does not change it)

The adaptive code adds:

```
adaptive/
├── README.md                       what this subdirectory is
├── docs/
│   ├── SPEC_v0.1_CLARIFICATION_NOTE.md   (already exists)
│   ├── IMPLEMENTATION_PLAN.md             (this file)
│   └── CALIBRATION_LOG.md                 (populated during Phase 9)
├── scheduler/
│   ├── __init__.py
│   ├── config.py                   defaults, marked TODO where calibration pending
│   ├── prior.py                    Phase 1
│   ├── subset.py                   Phase 2
│   ├── runner.py                   Phase 3
│   ├── confidence.py               Phase 4
│   ├── scheduler.py                Phase 5
│   ├── output.py                   Phase 6
│   ├── validate.py                 Phase 8
│   └── cli.py                      Phase 7
├── scripts/
│   ├── RUN_NGSRELATE_ADAPTIVE_LOCAL.sh
│   └── RUN_NGSRELATE_ADAPTIVE_LANTA.sh
└── tests/
    ├── test_prior.py
    ├── test_subset.py
    ├── test_runner.py
    ├── test_confidence.py
    ├── test_scheduler.py
    ├── test_output.py
    └── test_validate.py
```

Update repo-root files:

- `CLAUDE.md`: add a "Adaptive scheduling" section pointing at `adaptive/docs/IMPLEMENTATION_PLAN.md`
- `CHANGELOG.md`: add `[Unreleased]` entry "feat(adaptive): per-pair adaptive site budgeting for per-chrom runs"
- `README.md`: add a one-paragraph "Adaptive per-pair scheduling" section near the bottom with a link to `adaptive/README.md`
- `Makefile`: add `make adaptive-test` target → `pytest adaptive/tests/`

---

## 7. Definition of done

Phases 1–8 complete and:

- All tests passing in CI
- One real chromosome (LG12) processed end-to-end on a laptop and on LANTA, both producing identical `.res` output (Phase 8 Gate 5 reproducibility)
- Validation tool reports all 6 gates passing on LG12
- Defaults in `config.py` locked, calibration log written
- `adaptive/README.md` written for end users (not a re-spec, just how to run it)
- Methods paragraph drafted (SPEC §9) and stored in `adaptive/docs/METHODS_DRAFT.md`

When all of the above is true, mark the work `v1.0.0-adaptive` in `CHANGELOG.md` and tag the git commit.

---

## 8. What to ask before starting

Before you write any code, surface any of these if they don't have an obvious answer:

1. Does our fork's `ngsRelate-fast` binary actually accept `-a` and `-b` flags? (We didn't touch that code path; the answer should be yes, but the very first thing Phase 3 does is verify.)
2. What is the exact 23-column `.res` schema upstream produces? Do we have a header line cached anywhere in the repo to compare against?
3. Is there a representative tiny BEAGLE in the repo for testing, or does Phase 2/3 need to construct synthetic data?
4. Where will the genome-wide `.res` live in practice? Path will be configured in the SLURM scripts; what's the convention on LANTA?

These are not blockers — they're things to confirm with the user before sinking effort into a wrong assumption.

---

**End of implementation plan. Read in conjunction with SPEC v0.1 + clarification note. Begin at Phase 1.**
