#!/usr/bin/env bash
###############################################################################
# RUN_NGSRELATE_ADAPTIVE_LANTA.sh
#
# LANTA SLURM driver for the adaptive per-pair scheduler. One task per
# chromosome via SBATCH array.
#
# Submit:
#   sbatch --array=1-28 adaptive/scripts/RUN_NGSRELATE_ADAPTIVE_LANTA.sh
#
# Differs from the LOCAL script only in: SLURM headers + LANTA paths + worker
# count. The Python entrypoint is the same module.
###############################################################################
#SBATCH --job-name=ngsRelate_adaptive_perchrom
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:30:00
#SBATCH --mem=48G
#SBATCH --output=logs/ngsRelate_adaptive_perchrom_%A_%a.out
#SBATCH --error=logs/ngsRelate_adaptive_perchrom_%A_%a.err

set -euo pipefail

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
TOOL_DIR="${BASE}/tools/ngsRelate-fast"
NGSRELATE_FAST="${TOOL_DIR}/bin/ngsRelate-fast"

BEAGLE_DIR="${BASE}/popstruct_thin/04_beagle_byRF_majmin"
FREQS_DIR="${BASE}/popstruct_thin/05_ngsrelate/per_chrom_freqs"
SAMPLES="${BASE}/popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv"
GENOME_WIDE_RES="${BASE}/popstruct_thin/05_ngsrelate/cohort_226_full_fast_v1/relatedness.res"
N_SAMPLES=226

RUN_ID="cohort_226_perchrom_adaptive_v1"
OUT_BASE="${BASE}/popstruct_thin/05_ngsrelate/${RUN_ID}"
mkdir -p "${OUT_BASE}" logs

CHR_ID=$(printf "C_gar_LG%02d" "${SLURM_ARRAY_TASK_ID}")
BEAGLE="${BEAGLE_DIR}/catfish.${CHR_ID}.byRF.thin_500.beagle.gz"
FREQS="${FREQS_DIR}/freqs.${CHR_ID}.txt"
OUT_DIR="${OUT_BASE}/${CHR_ID}"
mkdir -p "${OUT_DIR}"

echo "[adaptive-LANTA] Task ${SLURM_ARRAY_TASK_ID} -> ${CHR_ID}"
for f in "${NGSRELATE_FAST}" "${BEAGLE}" "${FREQS}" "${SAMPLES}" "${GENOME_WIDE_RES}"; do
    if [[ ! -e "${f}" ]]; then
        echo "[adaptive-LANTA] ERROR: missing: ${f}" >&2
        exit 1
    fi
done

cd "${TOOL_DIR}"
python3 -m adaptive.scheduler \
    --binary           "${NGSRELATE_FAST}" \
    --beagle           "${BEAGLE}" \
    --freqs            "${FREQS}" \
    --samples          "${SAMPLES}" \
    --genome-wide-res  "${GENOME_WIDE_RES}" \
    --chrom            "${CHR_ID}" \
    --out-dir          "${OUT_DIR}" \
    --n-samples        "${N_SAMPLES}" \
    --workers          "${SLURM_CPUS_PER_TASK}" \
    --threads-per-pair 1

echo "[adaptive-LANTA] ${CHR_ID}: done"
