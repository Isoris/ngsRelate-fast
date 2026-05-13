# CLAUDE.md — orientation for Claude Code sessions

This is `ngsRelate-fast`, a surgical fork of [ANGSD/NgsRelate](https://github.com/ANGSD/NgsRelate).
Read this before making changes; it explains what is load-bearing and
what is safe to touch.

## What this repo is

A patched fork of ngsRelate v2 that runs 15-20× faster on whole-genome
cohorts. Two changes from upstream:

1. **Per-chromosome balanced-density downsampling** (`-D` flag, default
   100k sites/Gb). Implementation in `patches/01_ngsRelate_fast.patch`.
2. **Algebraic refactor in `emission_ngsrelate9`** — replaces `pow()`
   calls with explicit multiplications. Mathematically identical.

The Jacquard EM is **untouched** — that's the load-bearing path that
ngsPedigree and downstream tools depend on.

## What is load-bearing (do NOT modify)

- **`patches/01_ngsRelate_fast.patch`** — the actual scientific change.
  Touch this only with explicit user discussion; any edit invalidates
  validation runs and the manuscript methods text.
- **`contracts/*.schema.json`** — JSON Schema definitions. Adding fields
  requires a version bump (`v1` → `v2`) and an entry in `CHANGELOG.md`.
  Never silently change schema semantics.
- **`scripts/contract_io.py`** — consumed by downstream pipelines
  (ngsPedigree). Changing the `Run` dataclass shape is a breaking change.
- **The .res output schema** — 23 columns, identical to upstream. The
  patch must never change column count/order/names.

## What is safe to touch

- `README.md`, `PATCHES.md`, `CHANGELOG.md` — documentation
- `scripts/STEP_A07*.sh` — SLURM wrappers, project-specific paths
- `validate/*` — testing/validation infrastructure
- `build.sh` — build pipeline; may need fixups as upstream evolves
- `scripts/ngsrelate_fast_run.py` — wrapper; can add features as long as
  it still writes valid `.input.json` and `.output.json` per the schemas

## Repo layout

```
.
├── README.md             # GitHub landing page
├── PATCHES.md            # per-patch audit trail
├── CLAUDE.md             # THIS FILE — Claude Code orientation
├── CHANGELOG.md          # version history
├── CONTRIBUTING.md       # contribution rules
├── LICENSE               # GPL-3.0
├── Makefile              # shortcuts: make build / make validate / make clean
├── build.sh              # clone upstream + apply patch + compile
├── patches/              # the actual algorithmic change
├── contracts/            # JSON schemas + filled examples
├── scripts/              # SLURM submission + contract wrapper + consumer lib
├── validate/             # bit-equivalence + classification-agreement tests
├── build/                # [gitignored] upstream + patched source clones
└── bin/                  # [gitignored] compiled binaries
```

## Working environment

- Target HPC: **LANTA** (SLURM, account `lt200308`)
- Conda env: `assembly` (has g++, htslib, python3)
- Base path on LANTA: `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/`
- Suggested location for repo: `${BASE}/tools/ngsRelate-fast/`

## Common tasks

### Build
```bash
make build      # = bash build.sh
```

### Validate against canonical upstream
```bash
make validate BEAGLE=... FREQS=... N=226 SAMPLES=...
```

### Run on cohort (whole genome)
```bash
sbatch scripts/STEP_A07_fast_relatedness.sh
```

### Inspect a finished run's contracts
```bash
python3 scripts/contract_io.py /path/to/relatedness.res
```

## Status (current)

- Patch written, not yet compiled on LANTA
- Validation pipeline written, not yet executed
- Contracts schema at v1
- **DO NOT use in production until validation Stage 1 (identity mode)
  passes on at least one chromosome.**

## Working principles for this repo

1. **Don't touch the patch without explicit discussion.** The whole point
   of the validation pipeline is to verify a specific patch.
2. **Schema changes are breaking changes.** Bump versions, update
   `KNOWN_*_SCHEMAS` in `contract_io.py`, document in `CHANGELOG.md`.
3. **Defaults exist for a reason.** The default `-D 100000` was chosen
   from sampling-variance arguments. Don't change defaults silently;
   discuss first.
4. **Failure modes leave evidence.** The `.input.json` is written before
   the run for exactly this reason — a failed job still leaves a trail
   of what was attempted.
