# Changelog

All notable changes to ngsRelate-fast will be documented in this file.
Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- **feat(adaptive): per-pair adaptive site budgeting for per-chromosome
  runs.** New `adaptive/` Python package that wraps `ngsRelate-fast` to
  budget sites per pair based on the pair's KING-class prior (derived
  from the genome-wide `.res`). Most pairs (unrelated, third-degree)
  use a small budget and run fast; first-degree and ambiguous pairs use
  a larger budget; borderline / interesting-disagreement pairs are
  re-run at an escalated budget. Output is a Stage-2-compatible per-
  chromosome `.res` plus an audit `.adaptive_manifest.tsv` and a
  `.adaptive_run_manifest.json` (new schema:
  `ngsrelate_adaptive.run_manifest.v1`).
- `contracts/ngsrelate_adaptive.run_manifest.v1.schema.json` — schema
  for the new run manifest; sits alongside the existing input/output
  contracts.
- `adaptive/docs/SPEC_v0.1_CLARIFICATION_NOTE.md`,
  `adaptive/docs/IMPLEMENTATION_PLAN.md`,
  `adaptive/docs/DOWNSTREAM_CONSUMERS.md`,
  `adaptive/docs/CALIBRATION_LOG.md` (stub),
  `adaptive/docs/METHODS_DRAFT.md` (draft).
- `make adaptive-test` Makefile target.

### Not changed

- `patches/01_ngsRelate_fast.patch` — untouched.
- `contracts/ngsrelate_fast.input.v1.schema.json`,
  `contracts/ngsrelate_fast.output.v1.schema.json` — untouched.
- `scripts/contract_io.py`, `scripts/ngsrelate_fast_run.py`,
  `scripts/STEP_A07*.sh` — untouched. The adaptive scheduler is a
  new, additive entry point.
- The `.res` 23-column schema — preserved byte-for-byte (the adaptive
  writer copies binary output verbatim).



## [1.0.0] - 2026-05-13

Initial release. Pre-validation; do not use in production until
`validate/run_validation.sh` Stage 1 passes on at least one chromosome.

### Added

- Fork of upstream ANGSD/NgsRelate with two changes:
  - **Per-chromosome balanced-density downsampling** via new `-D` flag
    (default 100,000 sites per Gb). Density target applied independently
    to each chromosome by deterministic stride.
  - **Algebraic refactor in `emission_ngsrelate9`** — replaces `pow()`
    calls with precomputed power products. Mathematically identical to
    upstream.
- JSON contract system (modeled on ANGSD `.arg` files):
  - `ngsrelate_fast.input.v1` schema — written before run
  - `ngsrelate_fast.output.v1` schema — written after successful run
  - `scripts/ngsrelate_fast_run.py` — wrapper that generates contracts
  - `scripts/contract_io.py` — consumer library for downstream tools
- Two-stage validation pipeline:
  - Stage 1: identity mode (`-D 0`), expects 1e-6 agreement with upstream
  - Stage 2: downsampling mode, expects 0 first-degree class disagreements
- SLURM submission scripts for whole-genome and per-chromosome runs
- `build.sh` that clones upstream, applies patch, compiles both binaries

### Not changed

- Jacquard EM (`analyse_jaq`, `em`, `emAccel`, `emStep`) — byte-for-byte
  identical to upstream
- 2D-SFS pass and all output columns it populates
- IBS counting, KING-robust, R0, R1 derivation
- Output schema: 23 columns, identical column names and order

[Unreleased]: https://github.com/USER/ngsRelate-fast/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/USER/ngsRelate-fast/releases/tag/v1.0.0
