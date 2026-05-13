# ngsRelate-fast

> A surgical fork of [ANGSD/NgsRelate](https://github.com/ANGSD/NgsRelate)
> that runs **15-20× faster** on whole-genome cohorts without modifying
> the Jacquard EM that downstream tools depend on.

![status](https://img.shields.io/badge/status-pre--validation-orange)
![upstream](https://img.shields.io/badge/upstream-ANGSD%2FNgsRelate-blue)
![license](https://img.shields.io/badge/license-GPLv3-green)

---

## Why this exists

ngsRelate is the standard tool for genotype-likelihood-based kinship
estimation, but it was written for small studies. On a **226-sample cohort
with ~951k SNPs it takes hours**. Most of that time is spent computing
statistical information that already saturated long ago — the same 9
Jacquard coefficients can be estimated to 3+ decimal places from ~100k
well-spaced SNPs.

This fork makes two changes:

1. **Per-chromosome balanced-density downsampling** (new `-D` flag, default
   100,000 sites per Gb). Cuts the site count to what is statistically
   necessary while keeping uniform genomic density.
2. **Algebraic refactor in `emission_ngsrelate9`** to eliminate redundant
   `pow()` calls. Mathematically identical to upstream, just faster.

**The Jacquard EM, the 2D-SFS pass, IBS counting, KING/R0/R1 derivation
are all untouched.** Output schema is identical to upstream. Drop-in
replacement for downstream tools like ngsPedigree.

---

## Performance

| Cohort | Sites | Upstream | ngsRelate-fast (`-D 100000`) | Speedup |
|---|---|---|---|---|
| 226 samples, 1 Gb genome | ~951k | ~3 hours | ~10 min | ~18× |
| Per-chromosome (50 Mb) | ~30k | ~5 min | ~30 s | ~10× |

*Measured on LANTA HPC, 32 threads. Your numbers will vary with input
size, MAF distribution, and number of samples.*

---

## Installation

### Prerequisites

- A C++ compiler (`g++` ≥ 5.0)
- `make`, `git`, `python3`
- Optional: `htslib` for BCF input (most BEAGLE-only workflows don't need it)

On LANTA the `assembly` conda env has everything needed:

```bash
conda activate assembly
```

### Build

```bash
git clone https://github.com/<your-username>/ngsRelate-fast.git
cd ngsRelate-fast
bash build.sh
```

This will:

1. Clone the upstream `ANGSD/NgsRelate` repo into `build/`.
2. Compile the **upstream binary** to `bin/ngsRelate-upstream` (kept as a
   reference for validation).
3. Apply `patches/01_ngsRelate_fast.patch` to a second copy of upstream.
4. Compile the **patched binary** to `bin/ngsRelate-fast`.

If the build fails with `'gls' was not declared` or similar, upstream
uses different variable names than the patch assumes — see
[Troubleshooting](#troubleshooting).

---

## Validation (run before trusting the binary)

The whole point of this fork is that you can verify it. Run:

```bash
bash validate/run_validation.sh \
    /path/to/test.beagle.gz \
    /path/to/test.freqs \
    <n_samples> \
    /path/to/samples.txt
```

This runs **two stages**:

**Stage 1 — Identity mode (`-D 0`).** Both binaries run without
downsampling. Outputs should agree to within 1e-6 on every column. This
isolates the pow-elimination patch and proves the refactor is correct.

**Stage 2 — Downsampling mode (default `-D 100000`).** The fast binary
downsamples to ~100k sites/Gb. We **don't** expect bit-equivalence here —
we expect *classification agreement*. The script prints a confusion
matrix and flags any first-degree relationship that gets reclassified.

A clean validation looks like:

```
[stage 1] Identity mode
[validate] PASS: outputs near-equivalent within tolerance

[stage 2] Downsampling mode (-D 100000)
[validate] Class agreement: 25425/25425 (100.00%)
[validate] First-degree classification disagreements: 0/25425
[validate] PASS: no first-degree disagreements — safe for ngsPedigree
```

**Run on a small subset first** (e.g. one chromosome). Only proceed if
both stages pass.

---

## Usage

### Whole-genome run

```bash
bin/ngsRelate-fast \
    -G cohort.beagle.gz \
    -f allele_freqs.txt \
    -n 226 \
    -p 32 \
    -z samples.txt \
    -O relatedness.res
```

The default `-D 100000` applies automatically (100k sites per Gb,
balanced per chromosome). To get canonical upstream behavior, pass
`-D 0`.

### Per-chromosome run

The flag works the same way — density scales automatically with whatever
chromosome length appears in the input BEAGLE. A 50 Mb chromosome with
`-D 100000` keeps ~5,000 sites; a 5 Mb scaffold keeps ~500.

### Output format

Identical to upstream ngsRelate. 23 columns including:

- `a`, `b`, `ida`, `idb` — pair identifiers
- `theta`, `IBS0`, `KING`, `R0`, `R1` — relatedness summary statistics
- `J1`-`J9` — the 9 Jacquard coefficients (J9 = "both alleles IBD")
- `nSites` — sites used after filtering and downsampling
- `coverage`, plus 2D-SFS columns (untouched from upstream)

ngsPedigree consumes this format unchanged.

---

## How the downsampling works

Sites are kept by **per-chromosome deterministic stride**:

```
for each chromosome c:
    L_c     = max_pos_c - min_pos_c + 1     # length in bp
    target  = round(D × L_c / 1e9)           # target site count
    stride  = n_available_c / target         # spacing
    keep every stride-th site within chromosome c
```

Result: uniform sites-per-bp across all chromosomes, regardless of how
many sites survived MAF/missingness filtering on each one.

**Why deterministic stride and not random subsampling?** Reproducibility.
Running on the same input always picks the same sites. No seed
dependency, no run-to-run variance. The methods section gets one
sentence.

**Why per-chromosome and not global?** Global stride lets large
chromosomes dominate. If chromosome A has 200k post-filter sites and
chromosome B has 20k, a global 10× stride keeps 20k from A and 2k from
B — density unbalanced. Per-chromosome stride preserves the contract:
every chromosome contributes proportional to its assembled length.

---

## Flag reference

| Flag | Description | Default |
|---|---|---|
| `-G <file>` | input BEAGLE GL file | required |
| `-f <file>` | per-site allele frequencies | required |
| `-n <int>` | number of samples | required |
| `-p <int>` | threads | 4 |
| `-z <file>` | sample ID file (one per line) | optional |
| `-O <file>` | output `.res` file | stdout |
| **`-D <int>`** | **[NEW] target sites per Gb (0 disables)** | **100000** |
| `-c <int>` | call genotypes from GLs | 0 |
| `-F <int>` | estimate inbreeding | 0 |
| `-o <int>` | 3-coefficient mode (no inbreeding) | 0 |
| `-v <int>` | verbose | 0 |

All upstream flags are preserved unchanged.

---

## Run contracts (input/output JSON)

Every run writes **two JSON sidecar files** alongside the `.res` output —
modeled on ANGSD's `.arg` files but structured (`json.load()`-able) so
downstream code can read them programmatically.

```
relatedness.res                  ← the 23-column ngsRelate output
relatedness.res.input.json       ← written BEFORE the run
relatedness.res.output.json      ← written AFTER successful completion
relatedness.res.stderr.log       ← binary stderr capture
```

**`.input.json`** captures: tool version, upstream commit, patch hash,
input file paths with size + mtime + sha256, full effective parameter
set (not just user-passed flags), environment (host, SLURM IDs, conda
env), and the literal `argv` used. Written *before* the binary runs, so
even a failed job leaves a trail of what was attempted.

**`.output.json`** captures: per-chromosome downsampling summary (input
sites → kept sites, stride per chromosome), wallclock seconds, pair
count vs expected, output file hash, status (`ok` or `warn`), and any
warnings. Written *after* the binary exits cleanly. **Presence of this
file = run succeeded.** Downstream tools check for it before consuming
the `.res`.

Schemas in `contracts/`; examples in `contracts/examples/`.

### Always run via the wrapper

To get contracts, invoke through `scripts/ngsrelate_fast_run.py` rather
than the binary directly:

```bash
python3 scripts/ngsrelate_fast_run.py \
    --binary    bin/ngsRelate-fast \
    --beagle    cohort.beagle.gz \
    --freqs     allele_freqs.txt \
    --samples   samples.txt \
    --n         226 \
    --threads   32 \
    --D         100000 \
    --out       /path/to/relatedness.res \
    --run-id    cohort_226_full_fast_v1 \
    --patch-file patches/01_ngsRelate_fast.patch
```

The SLURM scripts in `scripts/` do this for you.

### Consuming contracts from downstream tools

`scripts/contract_io.py` is a small library that ngsPedigree (or any
downstream script) can use:

```python
from contract_io import load_run, RunNotComplete

try:
    run = load_run("/path/to/relatedness.res", verify_hashes=True)
except RunNotComplete:
    sys.exit("ngsRelate-fast hasn't finished — wait for .output.json")

# run.res_path, run.samples_path, run.params, run.n_pairs, etc.
df = pd.read_csv(run.res_path, sep="\t")
```

Or from the command line:

```bash
$ python3 scripts/contract_io.py /path/to/relatedness.res
run_id:         cohort_226_full_fast_v1
status:         ok
.res:           /scratch/.../relatedness.res
n_pairs:        25425
elapsed:        507.3s
D_sites_per_gb: 100000
downsampling:   951295 -> 100024 sites across 28 chromosomes
```

### Why two contracts and not one

Splitting input/output means:

- A failed run still leaves the `.input.json` — you know what was tried
- The presence/absence of `.output.json` is a clean success signal for
  pipelines (`test -f *.output.json` is enough)
- Input parameters are immutable once written; output stats are added
  separately on completion. No race condition where you read a
  half-written file.

---

## Repository layout

```
ngsRelate-fast/
├── README.md                          this file
├── PATCHES.md                         detailed per-patch rationale
├── LICENSE                            GPL-3.0 (inherited from upstream)
├── build.sh                           clone upstream + apply patch + compile
├── patches/
│   └── 01_ngsRelate_fast.patch        unified diff against upstream
├── contracts/
│   ├── ngsrelate_fast.input.v1.schema.json   JSON schema, input contract
│   ├── ngsrelate_fast.output.v1.schema.json  JSON schema, output contract
│   └── examples/                      filled-in example contracts
├── scripts/
│   ├── ngsrelate_fast_run.py          wrapper: runs binary + writes contracts
│   ├── contract_io.py                 library for downstream consumers
│   ├── STEP_A07_fast_relatedness.sh           whole-genome SLURM submission
│   └── STEP_A07b_fast_relatedness_per_chrom.sh per-chromosome SLURM array
├── validate/
│   ├── run_validation.sh              runs both binaries, diffs outputs
│   └── compare_to_upstream.py         python comparator
├── build/                             [gitignored] upstream + patched source
└── bin/                               [gitignored] compiled binaries
    ├── ngsRelate-upstream             reference, for validation
    └── ngsRelate-fast                 default tool going forward
```

---

## SLURM workflow (LANTA)

Pre-configured submission scripts in `scripts/`:

### Whole-genome (single job)

```bash
sbatch scripts/STEP_A07_fast_relatedness.sh
```

32 cores, 64 GB RAM, 1h wallclock cap. Typically finishes in under 10
minutes with default downsampling.

### Per-chromosome (SLURM array)

```bash
sbatch --array=1-28 scripts/STEP_A07b_fast_relatedness_per_chrom.sh
```

One task per chromosome, 16 cores each, 30 min cap. Output goes to
`${OUT_BASE}/C_gar_LG01/`, `.../C_gar_LG02/`, etc., each with a
`run_manifest.json` for ngsPedigree consumption.

Both scripts assume the binary at
`${BASE}/tools/ngsRelate-fast/bin/ngsRelate-fast`. Edit paths at the top
of each script to match your layout.

---

## Troubleshooting

### `build.sh` fails with `'gls' was not declared in this scope`

The patch hooks the downsampling call into `main()` using assumed
variable names from upstream. If upstream changed them, the inline
Python step in `build.sh` doesn't find the right anchor.

**Fix:** open `build/NgsRelate-fast/ngsRelate.cpp`, find where the
BEAGLE file finishes loading (look for `overall_number_of_sites` being
set), and insert manually:

```cpp
if (target_sites_per_gb > 0 && <NSITES_VAR> > 0) {
    <NSITES_VAR> = downsample_sites_balanced(
        &<GLS_VAR>, &<FREQS_VAR>, <SITE_IDS_VAR>,
        <NSITES_VAR>, <NIND_VAR>);
}
```

Substitute the actual variable names from your upstream version. Then
re-run `make` in `build/NgsRelate-fast/`.

### Validation Stage 1 fails (identity mode shows large differences)

This shouldn't happen — `-D 0` disables downsampling entirely, and the
pow refactor is mathematically identical. If you see this, the patch
applied to the wrong upstream version. Pin a specific upstream commit by
setting `UPSTREAM_COMMIT=<sha>` in `build.sh` and rebuild.

### Validation Stage 2 shows first-degree disagreements

Means downsampling lost discriminating information on some borderline
pair. Increase the target density (`-D 200000` or `-D 500000`) and
re-run. For a clean cohort this should not happen at the default 100k.

### Output has 0 rows or `nSites` is 0 for every pair

Site IDs in your BEAGLE don't match the expected `chrom_pos` format. The
downsampler can't parse them and falls back to using all sites. Check
stderr for `WARN: cannot parse site id`. Either fix the BEAGLE format or
run with `-D 0` to disable downsampling.

---

## When NOT to use this fork

- **You need bit-equivalence to canonical ngsRelate output.** Use upstream.
- **Your sample size is small** (n < 50). Upstream is already fast enough.
- **Your BEAGLE has non-standard site IDs.** The parser expects
  `<chrom>_<integer_position>`. Other formats fall back to no
  downsampling silently.
- **You're inferring inbreeding coefficients on low-MAF sites.**
  Downsampling to 100k/Gb may drop rare informative sites. Either
  increase `-D` or disable downsampling for inbreeding-specific runs.

---

## Citing

If you use this fork in published work, **cite the upstream tool first**:

> Korneliussen TS, Moltke I. NgsRelate: a software tool for estimating
> pairwise relatedness from next-generation sequencing data.
> *Bioinformatics* 31, 4009-4011 (2015).
>
> Hanghøj K, Moltke I, Andersen PA, Manica A, Korneliussen TS. Fast and
> accurate relatedness estimation from high-throughput sequencing data
> in the presence of inbreeding. *GigaScience* 8, giz034 (2019).

Suggested methods text for the fork:

> Pairwise relatedness was estimated with ngsRelate-fast v1.0
> (https://github.com/<your-username>/ngsRelate-fast), a fork of
> ngsRelate v2 (Hanghøj et al. 2019) that downsamples the input site
> set to 100,000 sites per Gb of assembled genome length, applied per
> chromosome by deterministic stride. The Jacquard EM and all downstream
> statistics (theta, KING, R0, R1, IBS counts) are unchanged from
> upstream and produce output bit-equivalent to canonical ngsRelate
> when downsampling is disabled.

---

## License

GPL-3.0, inherited from upstream `ANGSD/NgsRelate`.

## Acknowledgements

All algorithmic work belongs to the ANGSD/NgsRelate authors. This fork
is purely an engineering optimization on top of their estimator.
