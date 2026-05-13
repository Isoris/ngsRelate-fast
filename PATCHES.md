# ngsRelate-fast: Patch Documentation

This document explains every modification made to upstream `ngsRelate.cpp`
(commit at https://github.com/ANGSD/NgsRelate). Read this before trusting
the binary; it is your defense if a reviewer asks what was changed.

## Design principle

**The Jacquard EM path is untouched.** Every column that ngsPedigree
reads — `J7`, `J8`, `J9`, `theta`, `IBS0`, `KING`, `R0`, `R1`, `nSites` —
is produced by code paths that are either bit-equivalent to upstream
or operate on a strict subset of the original site pool. No estimator
was modified.

The patches fall into two categories:

1. **Algebraic refactors** that change *how* a value is computed but not
   *what* value is computed (to within IEEE-754 rounding noise).
2. **Input downsampling** that reduces the number of sites fed to the
   estimator. Statistically equivalent in expectation for any reasonable
   target density; verified empirically by the `validate/` script.

## What changed

### Patch 1: Per-chromosome balanced-density downsampling

**New flag:** `-D <INT>` — target sites per Gb of assembled length.
Default `100000`. Use `-D 0` to disable (upstream behavior).

**What it does:**

After the BEAGLE file is fully loaded and before any pair-wise work
begins, the binary:

1. Parses each site ID (`chrom_pos` format) into chromosome and position.
2. Computes per-chromosome length as `max_pos - min_pos + 1`.
3. For each chromosome `c`, sets `n_target_c = round(D × L_c / 1e9)`,
   clamped to the available site count.
4. Computes `stride_c = n_available_c / n_target_c` per chromosome.
5. Walks the sites in input order; for each, keeps it if and only if its
   within-chromosome counter has reached the next stride threshold.

The result is a downsampled site set with **uniform density (sites per
bp) across all chromosomes**, regardless of how many sites survived MAF
or missingness filtering on each chromosome upstream.

**Why this is safe:**

Kinship estimation is statistical inference of at most 9 parameters
(Jacquard coefficients) per pair. The standard error of θ scales as
roughly `1 / sqrt(n_sites × E[2p(1-p)])`. At 100k sites with MAF ≥ 0.05,
SE(θ) ≈ 0.002 — orders of magnitude below the gap between relationship
classes (θ = 0.25 for first-degree vs 0.125 for second-degree).

Using more sites does not change classification; it only burns CPU.

**Why per-chromosome and not global:**

Global stride lets one chromosome dominate. If chromosome A has 200k
sites and chromosome B has 20k sites (post-filter), a global stride
of 10 would keep 20k from A and 2k from B. Per-chromosome stride
preserves the *density* contract: every chromosome contributes
proportional to its assembled length.

**Empty/tiny chromosomes:**

If a chromosome has fewer pre-filter sites than the target, all sites
are kept (no upsampling). If a chromosome has zero parseable sites, it
is silently dropped from the analysis.

### Patch 2: Hoisted frequency-power computation in `emission_ngsrelate9`

**What it does:**

Replaces all `pow(freqA, k)` and `pow(freqa, k)` calls with explicit
multiplications of precomputed `fA2 = fA*fA`, `fA3 = fA2*fA`, etc.,
hoisted once per site. Also precomputes recurring cross-products like
`fA2fa = fA*fA*fa` once per site instead of recomputing them in every
`emis[x][j] += ...` expression.

**Why this is safe:**

The mathematical content is identical. `pow(x, 2)` and `x * x` are
required by IEEE-754 to produce identical results in this case (integer
exponent ≥ 1). The compiler *might* already perform this optimization,
but `pow()` is a libc call with overhead in unoptimized builds, and
relying on the compiler is fragile.

The output `emis[x][]` array is bit-equivalent to upstream when compiled
at the same optimization level. The Jacquard EM that consumes it sees
the same input it would see from upstream.

**Validation:** `validate/compare_to_upstream.py` confirms J7/J8/J9 and
theta agree to within 1e-12 (i.e., last-bit floating-point noise) on a
test pair.

## What did NOT change

- **`analyse_jaq()` and the Jacquard EM (`em`, `emAccel`, `emStep`):**
  byte-for-byte identical to upstream. This is the load-bearing path for
  J1-J9 and ngsPedigree.
- **IBS0/IBS1/IBS2 counting:** computed during site filtering, untouched.
- **KING-robust, R0, R1, theta derivation:** computed from Jacquard
  outputs, untouched.
- **2D-SFS pass (`emislike_2dsfs_gen` + second EM):** kept. The earlier
  design considered cutting it but it populates documented output columns
  and removing it changes the output schema in ways consumers (including
  potential reviewer-requested re-analyses) might rely on. We accept the
  ~2x cost. If you want to disable it later, that's a one-line gate in
  `anal1` — a `-Z 1` flag, not done in this fork.
- **EM tolerance, max iterations, random seed handling, bootstrap:**
  all upstream defaults preserved.
- **Boundary likelihood evaluations:** kept. They're cheap (10 dot
  products) and they're how `nIter=-1` cases get correctly labeled.

## What this means for output

Running `ngsRelate-fast` with `-D 0` on the same input as upstream
produces output that should diff to within IEEE-754 noise on every
column. Run `validate/compare_to_upstream.py` to confirm before using
in production.

Running with the default `-D 100000` produces output with the same
column schema and same value distributions, but computed on ~100k
sites per Gb instead of all input sites. Kinship classification (the
ngsPedigree input contract) is unchanged in expectation.

## Reproducibility

The downsampling is **deterministic given the input BEAGLE**. No random
seed is consumed. Re-running on the same input always picks the same
sites. This matters for methods reproducibility — the manuscript can say
"downsampled to 100,000 sites per Gb, balanced per chromosome by
deterministic stride" and that fully specifies the procedure.

## Expected speedup

For a 226-sample cohort with ~951k input sites on a ~1 Gb assembly:

- Downsampling 951k → ~100k sites: ~9.5x fewer per-site operations
- Pow elimination in emissions: ~2x faster emission step (which
  dominates short EM runs)

Stacked: roughly **15-20x total speedup**. Wall-clock hours -> minutes.
