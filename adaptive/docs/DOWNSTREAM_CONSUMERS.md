# Downstream Consumers of the Per-Chromosome `.res`

**Purpose:** capture the full chain of work that depends, directly or transitively, on the per-chromosome ngsRelate `.res` files. This is the blast radius for any format change, semantic drift, or borderline-pair reclassification introduced by the adaptive scheduler.

**Audience:** anyone tempted to relax the "`.res` schema is sacred" constraint in `IMPLEMENTATION_PLAN.md` §1. Read this first.

---

## The chain

```
Per-chromosome .res files  (the adaptive scheduler's output)
      ↓
ngsPedigree Stage 2  (per-chrom QC; produces confirmed PO/FS dyads + triads)
      ↓
ngsPedigree Stage 3  (inheritance maps: Gold/Silver/Bronze segments,
                      parental_hap_inherited)
      ↓
HPP  (projects MODULE_CONSERVATION variant scores onto inheritance segments
      — this is "haplotype-projected pathogenicity")
      ↓
KBC cross-checks  (projection consistency vs arrangement assignments inside
                   known inversion intervals)
      ↓
Manuscript figures + tables (MS_Inversions_North_african_catfish v20)
```

Six hops from the per-chrom `.res` to manuscript output. Each hop is silent — a drift introduced at the top does not announce itself; it surfaces as wrong-looking figures at the bottom, attributed to whatever the most recent change was, not to the per-chrom `.res` six steps back.

---

## What each consumer needs from the `.res`

### ngsPedigree Stage 2 — direct consumer

Reads: full 23-column `.res` per chromosome.
Critical columns: `theta`, `IBS0`, `KING`, `J7`, `J8`, `J9`.
Failure mode if these drift: PO/FS classification differs from genome-wide call, pair gets flagged for review when it shouldn't or escapes review when it should.

### ngsPedigree Stage 3 — indirect (via Stage 2)

Reads: confirmed PO dyad and triad lists from Stage 2 + BEAGLE GLs of those pairs.
Failure mode if Stage 2 drifts: missing dyad → missing inheritance map. The map is unrecoverable downstream; HPP just sees a gap.

### HPP — indirect (via Stage 3)

Reads: Stage 3 inheritance maps (Gold/Silver/Bronze segments, `parental_hap_inherited`) + `variant_master_scored.tsv` (from MODULE_CONSERVATION) + joint multisample VCF.
Failure mode: projects damaging variants onto the wrong haplotype segments. The output is a per-sample damaging-burden estimate that is silently miscounted.

### KBC cross-checks — indirect (via HPP)

Reads: HPP outputs + `kbc_variant_arrangement_assignments.tsv` (table B) + PCAngsd K=3 karyotype calls + per-sample ROH BEDs.
Failure mode: cross-check passes spuriously, or fails spuriously, because the underlying HPP projection was already wrong.

### Manuscript figures — indirect (via KBC)

Failure mode: figures look fine but numbers are off. Reviewer flags inconsistency. Worst case: figure is regenerated multiple times trying to fix the surface issue without finding the upstream cause.

---

## Inputs HPP needs that do NOT come from `.res`/BEAGLE

For completeness — these are HPP's other dependencies, not affected by the adaptive scheduler but listed so the full input contract is clear:

| Input | Source module | Why HPP needs it |
|---|---|---|
| `variant_master_scored.tsv` | MODULE_CONSERVATION STEP 16 (catfish-variant-analysis) | the projection targets — SnpEff, SIFT4G, VESM_650M, splice scores |
| Joint multisample VCF | MODULE_CONSERVATION STEP 03 | per-parent genotypes at every variant; BEAGLE GLs at thin-500 are too sparse |
| Reference FASTA `fClaHyb_Gar_LG.fa` | existing | variant normalisation + sanity checks |
| Sample metadata | sample sheet | sample IDs, batch, cohort flag |
| Parent-het phase | open (HANDOFF #1) | which parent-het site phases to hap-1 vs hap-2 |

Optional but recommended downstream:

| Input | What it adds |
|---|---|
| KBC table B (`kbc_variant_arrangement_assignments.tsv`) | cross-check for variants inside inversion intervals (table E) |
| ROH BEDs per sample | confirms unambiguous projection in ROH-resident segments |
| PCAngsd K=3 karyotype calls | translates `hap_copy` → arrangement for KBC cross-check |

**None of these are inputs to or outputs of the adaptive scheduler.** They appear here so anyone reading this document has the complete picture of what HPP consumes, and understands that a per-chrom `.res` drift is the *only* failure mode the adaptive scheduler can introduce into this chain.

---

## Implication for the adaptive scheduler's design

The "adaptive scheduling is safe because the `.res` schema is preserved" claim only holds if Gate 4 (Stage 2 semantic compatibility) actually holds — not just in the abstract but for the **specific pairs that feed Stage 3 → HPP**.

Concretely:

- A confirmed-PO pair lost in Stage 2 = missing Stage 3 inheritance map = missing HPP projection = manuscript figure with one less sample.
- A confirmed-PO pair *gained* in Stage 2 (false positive caused by per-chrom drift) = inheritance map computed on a non-PO pair = HPP projection that's biologically meaningless.

Both are bad. Both are silent. Neither is caught by Gate 1 (edge-class precision) or Gate 2 (site efficiency) — both of those operate at the per-pair level, not at the Stage 2 set-membership level.

**This is why `IMPLEMENTATION_PLAN.md` §2 Phase 8 elevates Gate 4 to tier-1** with the specific requirement: *confirmed PO dyad set must match exactly (set equality)*. Anything less than set equality on PO dyads compromises the entire downstream chain documented above.

---

## What this document does NOT cover

- Specs for HPP, MODULE_CONSERVATION, KBC, PCAngsd, or parent-het phase — those live in their own specs and repos.
- The genome-wide `.res` chain — see the clarification note in `SPEC_v0.1_CLARIFICATION_NOTE.md` for how the genome-wide `.res` feeds the adaptive scheduler as prior.
- Implementation details of the adaptive scheduler itself — see `IMPLEMENTATION_PLAN.md`.

---

**End of downstream consumers map.**
