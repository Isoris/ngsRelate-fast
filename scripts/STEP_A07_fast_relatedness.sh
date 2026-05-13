#!/usr/bin/env bash
###############################################################################
# STEP_A07_fast_relatedness.sh
#
# Run ngsRelate-fast on the full 226-sample C. gariepinus hatchery cohort.
# Uses ngsrelate_fast_run.py wrapper, which generates input/output JSON
# contracts alongside the .res output. Downstream tools (ngsPedigree) check
# the .output.json before consuming the .res.
#
# Outputs (in OUT_DIR):
#   - relatedness.res                 23-column ngsRelate output
#   - relatedness.res.input.json      input contract (written BEFORE run)
#   - relatedness.res.output.json     output contract (written AFTER run)
#   - relatedness.res.stderr.log      binary stderr capture
###############################################################################
#SBATCH --job-name=ngsRelate_fast_226
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --output=logs/ngsRelate_fast_226.%j.out
#SBATCH --error=logs/ngsRelate_fast_226.%j.err

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================
BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
TOOL_DIR="${BASE}/tools/ngsRelate-fast"

NGSRELATE_FAST="${TOOL_DIR}/bin/ngsRelate-fast"
WRAPPER="${TOOL_DIR}/scripts/ngsrelate_fast_run.py"
PATCH_FILE="${TOOL_DIR}/patches/01_ngsRelate_fast.patch"

BEAGLE="${BASE}/popstruct_thin/04_beagle_byRF_majmin/catfish.wholegenome.byRF.thin_500.beagle.gz"
FREQS="${BASE}/popstruct_thin/05_ngsrelate/freqs.aligned_to_beagle.txt"
SAMPLES="${BASE}/popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv"
N_SAMPLES=226
DSAMPLE=100000

RUN_ID="cohort_226_full_fast_v1"
OUT_DIR="${BASE}/popstruct_thin/05_ngsrelate/${RUN_ID}"
mkdir -p "${OUT_DIR}" logs

# ============================================================================
# Pre-flight checks
# ============================================================================
for f in "${NGSRELATE_FAST}" "${WRAPPER}" "${BEAGLE}" "${FREQS}" "${SAMPLES}"; do
    if [[ ! -e "${f}" ]]; then
        echo "[A07-fast] ERROR: missing: ${f}" >&2
        exit 1
    fi
done

cp "${SAMPLES}" "${OUT_DIR}/samples.txt"

# ============================================================================
# Run via wrapper (writes contracts, captures stderr, validates output)
# ============================================================================
echo "[A07-fast] Run ID: ${RUN_ID}"
echo "[A07-fast] Output: ${OUT_DIR}/relatedness.res"

python3 "${WRAPPER}" \
    --binary     "${NGSRELATE_FAST}" \
    --beagle     "${BEAGLE}" \
    --freqs      "${FREQS}" \
    --samples    "${OUT_DIR}/samples.txt" \
    --n          "${N_SAMPLES}" \
    --threads    "${SLURM_CPUS_PER_TASK}" \
    --D          "${DSAMPLE}" \
    --out        "${OUT_DIR}/relatedness.res" \
    --run-id     "${RUN_ID}" \
    --patch-file "${PATCH_FILE}"

# ============================================================================
# Confirm the contracts wrote OK
# ============================================================================
if [[ ! -f "${OUT_DIR}/relatedness.res.output.json" ]]; then
    echo "[A07-fast] ERROR: output contract missing — run failed" >&2
    exit 2
fi

echo "[A07-fast] DONE"
echo "[A07-fast] Inspect contracts:"
echo "  python3 ${TOOL_DIR}/scripts/contract_io.py ${OUT_DIR}/relatedness.res"
echo "[A07-fast] Next: ngsPedigree Stage 1"
