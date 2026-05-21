# METHODS_DRAFT.md — manuscript-ready methods paragraph

**Status:** DRAFT. Final wording locked when calibration completes
(IMPLEMENTATION_PLAN.md §7 Definition of done).

---

## Per-chromosome adaptive site-budgeting for relatedness estimation

To reduce per-chromosome ngsRelate runtime without compromising pairwise
classification used downstream by ngsPedigree Stage 2, we developed an
adaptive scheduler that runs each sample pair at a site budget tailored
to the pair's *prior* relationship class. The prior class for pair
*(i, j)* is derived directly from the genome-wide ngsRelate-fast `.res`
using KING-robust thresholds on the kinship coefficient *θ*
(Manichaikul *et al.* 2010): θ ≥ 0.354 → *duplicate*; 0.177 ≤ θ < 0.354
→ *first degree* (split into *parent-offspring* if IBS0 ≤ 0.008,
otherwise *full sibling*; ambiguous when IBS0 ∈ [0.0064, 0.0096]);
0.0884 ≤ θ < 0.177 → *second degree*; 0.0442 ≤ θ < 0.0884 → *third
degree*; otherwise *unrelated*.

Per-class site budgets — `BUDGET_LOW`, `BUDGET_MED`, `BUDGET_HIGH` —
were calibrated on linkage group LG12 against a uniform-stride
baseline (`-D 100000`) until the per-chrom edge-class agreement against
that baseline reached ≥95%, the set of pairs classified as
parent-offspring matched exactly, and the recall of "interesting
disagreements" (per-pair prior-vs-per-chrom transitions implying a
pedigree or inversion-driven discrepancy; see SPEC §4.3) reached
≥95% with zero loss of *parent-offspring → unrelated* transitions.

The scheduler escalates a pair to a larger site budget when its
per-chromosome θ lands within a calibrated window of a KING boundary or
when the prior-vs-per-chrom class transition is in a pre-defined set of
nine biologically interesting disagreements (e.g., *parent-offspring →
unrelated* on a single chromosome, which signals either a haplotype-
inheritance violation in an inversion region or a sample swap). For
each chromosome, the scheduler pre-thins BEAGLE input files to a small
number of distinct budget tiers (cached on disk), then dispatches one
ngsRelate-fast subprocess per pair with the `-a`/`-b` flags. Pairs are
written to the output `.res` in canonical (alphabetical sample-ID)
order; the file is byte-compatible with the per-chromosome output
consumed by ngsPedigree Stage 2.

Validation against the uniform-stride baseline on linkage group LG12
confirmed all six gates of the validation protocol (SPEC §6); per-pair
budget allocations, escalation events, and the prior source `.res`
hash are recorded in a per-chromosome run manifest
(`*.adaptive_run_manifest.json`) for reproducibility.

---

## Citations to insert

- Manichaikul, A. *et al.* (2010). Robust relationship inference in
  genome-wide association studies. *Bioinformatics*, 26(22), 2867–2873.
- (your existing ngsRelate / ngsRelate-fast citations)
- (ngsPedigree citation when published)
