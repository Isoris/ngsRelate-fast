#!/usr/bin/env bash
###############################################################################
# STEP_A07b_fast_relatedness_per_chrom.sh
#
# Per-chromosome ngsRelate-fast for ngsPedigree Stage 2. SLURM array job:
# one task per chromosome. Uses ngsrelate_fast_run.py wrapper for contracts.
#
# Submit:
#   sbatch --array=1-28 STEP_A07b_fast_relatedness_per_chrom.sh
#
# Each task produces relatedness.res + .input.json + .output.json + .stderr.log
###############################################################################
#SBATCH --job-name=ngsRelate_fast_perchrom
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --output=logs/ngsRelate_fast_perchrom_%A_%a.out
#SBATCH --error=logs/ngsRelate_fast_perchrom_%A_%a.err

set -euo pipefail

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
TOOL_DIR="${BASE}/tools/ngsRelate-fast"

NGSRELATE_FAST="${TOOL_DIR}/bin/ngsRelate-fast"
WRAPPER="${TOOL_DIR}/scripts/ngsrelate_fast_run.py"
PATCH_FILE="${TOOL_DIR}/patches/01_ngsRelate_fast.patch"

BEAGLE_DIR="${BASE}/popstruct_thin/04_beagle_byRF_majmin"
FREQS_DIR="${BASE}/popstruct_thin/05_ngsrelate/per_chrom_freqs"
SAMPLES="${BASE}/popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv"
N_SAMPLES=226
DSAMPLE=100000

RUN_ID="cohort_226_perchrom_fast_v1"
OUT_BASE="${BASE}/popstruct_thin/05_ngsrelate/${RUN_ID}"
mkdir -p "${OUT_BASE}" logs

CHR_ID=$(printf "C_gar_LG%02d" "${SLURM_ARRAY_TASK_ID}")
BEAGLE="${BEAGLE_DIR}/catfish.${CHR_ID}.byRF.thin_500.beagle.gz"
FREQS="${FREQS_DIR}/freqs.${CHR_ID}.txt"
OUT_DIR="${OUT_BASE}/${CHR_ID}"
mkdir -p "${OUT_DIR}"

echo "[A07b-fast] Task ${SLURM_ARRAY_TASK_ID} -> ${CHR_ID}"

for f in "${NGSRELATE_FAST}" "${WRAPPER}" "${BEAGLE}" "${FREQS}" "${SAMPLES}"; do
    if [[ ! -e "${f}" ]]; then
        echo "[A07b-fast] ERROR: missing: ${f}" >&2
        exit 1
    fi
done

cp "${SAMPLES}" "${OUT_DIR}/samples.txt"

python3 "${WRAPPER}" \
    --binary     "${NGSRELATE_FAST}" \
    --beagle     "${BEAGLE}" \
    --freqs      "${FREQS}" \
    --samples    "${OUT_DIR}/samples.txt" \
    --n          "${N_SAMPLES}" \
    --threads    "${SLURM_CPUS_PER_TASK}" \
    --D          "${DSAMPLE}" \
    --out        "${OUT_DIR}/relatedness.res" \
    --run-id     "${RUN_ID}.${CHR_ID}" \
    --patch-file "${PATCH_FILE}"

if [[ ! -f "${OUT_DIR}/relatedness.res.output.json" ]]; then
    echo "[A07b-fast] ERROR: output contract missing for ${CHR_ID}" >&2
    exit 2
fi

echo "[A07b-fast] ${CHR_ID}: DONE"
