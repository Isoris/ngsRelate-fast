# SPEC_NGSRELATE_ADAPTIVE v0.1 — Clarification Note

**Note date:** 2026-05-13
**Spec clarified:** `SPEC_NGSRELATE_ADAPTIVE_v0_1.md` (Status: SPEC ONLY — awaiting audit)
**Status:** AUDIT NOTE — to be folded into v0.2 of the spec, or kept as a permanent companion document
**Outcome:** Spec is implementable as-is. This note resolves one ambiguity about the prior source and one open question (§7 OQ6). All other open questions stand as calibration tasks for the implementation phase.

---

## 1. Headline clarification — prior source is the genome-wide `.res`, not ngsPedigree Stage 1

The v0.1 spec refers to the per-pair prior as the "genome-wide Stage 1 classification" (§3 Inputs table; §4.1 Step 1). This phrasing implies that **ngsPedigree must have already run** before the adaptive per-chromosome scheduler can run, which would make the adaptive scheduler downstream of ngsPedigree rather than a peer of it.

This is not the intended architecture. The intended architecture is:

```
   ngsRelate-fast genome-wide  ──→  genome-wide .res  ──┬──→  ngsPedigree Stage 1
                                                        │
                                                        └──→  KING-class derivation
                                                                   │
                                                                   ↓  (prior for adaptive)
   ngsRelate-fast per-chrom (adaptive)  ──→  per-chrom .res
                                                                   ↓
                                                        ngsPedigree Stage 2
```

Both ngsRelate runs (genome-wide and per-chromosome) are **siblings**, both upstream of ngsPedigree. ngsPedigree Stage 1 consumes the genome-wide `.res`; Stage 2 consumes the per-chromosome `.res` files. The adaptive scheduler operates between the two ngsRelate runs, **not** after ngsPedigree.

**Clarification:** The per-pair prior class used by the adaptive scheduler (§4.1) is derived **directly from the genome-wide `.res`** using KING-robust thresholds (Manichaikul et al. 2010), not from ngsPedigree output. Concretely, for each pair `(i, j)` with genome-wide kinship coefficient `θ`:

| Threshold on θ | Class assignment |
|---|---|
| `θ ≥ 0.354` | `duplicate_or_clone` |
| `0.177 ≤ θ < 0.354` | `first_degree` (PO or FS; see §1.1 below for splitting) |
| `0.0884 ≤ θ < 0.177` | `second_degree` |
| `0.0442 ≤ θ < 0.0884` | `third_degree` |
| `θ < 0.0442` | `unrelated` |

These are the exact KING-robust class boundaries; they require no separate inference step and are reproducible from the genome-wide `.res` alone.

### 1.1 PO vs FS split within `first_degree`

The KING θ band `[0.177, 0.354)` lumps PO and FS together. For the adaptive scheduler's budget allocation (§4.1 table), splitting them matters because:

- The spec assigns both `full_sibling` and `parent_offspring` to `BUDGET_HIGH` — so for budget purposes the split doesn't change behavior.
- But the "interesting disagreement" rule (§4.3) **does** distinguish them: a genome-wide PO → per-chrom non-PO is "interesting" in a way that genome-wide FS → per-chrom non-FS is not necessarily.

**Recommendation:** Use `IBS0` to split first-degree pairs at the standard threshold `IBS0 ≤ 0.008` → `parent_offspring`; otherwise → `full_sibling`. This is the same threshold ngsPedigree Stage 1 uses (per spec §2 "Stage 2's `ibs0_po_max = 0.008`"), so the prior is consistent with what Stage 1 would produce on the same data.

`ambiguous_first_degree` (the spec's seventh class) is reserved for pairs where the IBS0 split itself is uncertain — operationally, IBS0 within ±20% of the 0.008 threshold (i.e., IBS0 ∈ [0.0064, 0.0096]). This is the only class that requires a tunable parameter; document it as `ambiguous_first_degree_ibs0_band` in the run manifest.

## 2. This resolves §7 OQ6 ("no genome-wide anchor")

The original open question reads:

> First-run case where genome-wide ngsRelateFast hasn't been done yet. In that case, fall back to BUDGET_MED for all pairs (uniform-ish), produce per-chrom .res, then a downstream pass can recompute genome-wide from per-chrom output. Out of scope for v0.1 but worth noting.

With the clarification above, OQ6 dissolves: **the genome-wide ngsRelate-fast run is a hard prerequisite, not optional**. The adaptive scheduler refuses to run without a genome-wide `.res`. The "first-run case" is simply "run ngsRelate-fast genome-wide first" — which is already the standard first step of the pipeline and takes ~10 minutes with the `-D 100000` default.

**Operationally:** the adaptive scheduler's preflight check reads the genome-wide `.res` and refuses to start if it's missing, empty, or fewer rows than expected. No fallback to `BUDGET_MED`-for-all; that would be a silent failure mode.

## 3. Staleness check (§7 OQ5) — concrete proposal

The spec's OQ5 asks how to detect a stale genome-wide anchor. Concrete proposal:

The adaptive scheduler's run manifest records:

- `anchor_genome_wide_res_path` — path to the genome-wide `.res` used as prior
- `anchor_genome_wide_res_sha256` — content hash of that file
- `anchor_genome_wide_input_contract_sha256` — if the genome-wide run was produced by `ngsrelate_fast_run.py` (which it should be), copy the `.input.json` hash too

**Preflight rule:** the per-chrom BEAGLE's sample list must match (set equality) the genome-wide `.res`'s sample list. If not, refuse to run unless `--allow-mismatched-anchor` is passed. The flag exists only as an escape hatch; the default is strict.

This is stronger than v0.1's "BEAGLE input panel match" check because it uses the contract system from `ngsRelate-fast` directly — no separate panel-version tracking needed.

## 4. What does NOT change

All other v0.1 spec content stands:

- §1 Motivation — unchanged
- §2 Scope — unchanged
- §3 Inputs — replace "ngsPedigree Stage 1 classification" with "genome-wide `.res` via KING-class derivation"
- §4.1 Step 1 — replace "Load the genome-wide Stage 1 classification" with "Compute per-pair KING-class from genome-wide `.res` per §1 of this note"
- §4.2 – §4.6 — unchanged
- §5 Local execution path — unchanged
- §6 Validation gates — unchanged, but add: **Gate 0: Genome-wide anchor exists and passes preflight checks (§3 of this note).** No adaptive run starts without it.
- §7 Open questions — OQ5 (concrete proposal in §3 above); OQ6 (resolved, see §2 above); OQ1–4, OQ7–9 remain open and are calibration tasks for the implementation phase
- §8 Deliverables checklist — unchanged
- §9 Manuscript integration — unchanged, but the methods text "based on the pair's genome-wide classification" should be "based on the pair's KING-robust class derived from the genome-wide ngsRelate-fast output" to be precise

## 5. Implications for naming

The spec uses "genome-wide Stage 1 classification" and "genome-wide class" interchangeably. After this clarification, prefer:

- **"genome-wide KING-class"** or **"genome-wide prior class"** — derived directly from the genome-wide `.res`
- NOT "Stage 1 classification" — that's ngsPedigree output, which is downstream

The audit sidecar column `genome_wide_class` (§4.5 manifest) stays the same name but its definition is "genome-wide KING-class per the rules in this clarification note."

---

**End of clarification note. This document folds into SPEC v0.2 when v0.2 is written; until then it is read in conjunction with v0.1.**
