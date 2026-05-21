#!/usr/bin/env bash
###############################################################################
# RUN_NGSRELATE_ADAPTIVE_LOCAL.sh
#
# Laptop driver: loop over 28 chromosomes, run the adaptive scheduler on each.
# Same Python entrypoint as the LANTA script; just paths + worker count differ.
#
# Usage:
#   bash adaptive/scripts/RUN_NGSRELATE_ADAPTIVE_LOCAL.sh
#
# Requires that you've already produced a genome-wide .res (the prior source).
# If you haven't:
#   1. Build the binary:    make build
#   2. Run genome-wide:     python scripts/ngsrelate_fast_run.py --D 100000 ...
#   3. Then run this script.
#
# Tune WORKERS to your laptop's CPU count (use os.cpu_count() if unsure).
###############################################################################
set -euo pipefail

# ---- Edit these paths for your local setup ---------------------------------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
NGSRELATE_FAST="${REPO_ROOT}/bin/ngsRelate-fast"

# Per-chrom inputs (28 BEAGLEs, 28 freqs).
BEAGLE_DIR="${BEAGLE_DIR:-${HOME}/data/popstruct_thin/04_beagle_byRF_majmin}"
FREQS_DIR="${FREQS_DIR:-${HOME}/data/popstruct_thin/05_ngsrelate/per_chrom_freqs}"
SAMPLES="${SAMPLES:-${HOME}/data/popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv}"
GENOME_WIDE_RES="${GENOME_WIDE_RES:-${HOME}/data/popstruct_thin/05_ngsrelate/cohort_226_full_fast_v1/relatedness.res}"
N_SAMPLES="${N_SAMPLES:-226}"
WORKERS="${WORKERS:-8}"

OUT_BASE="${OUT_BASE:-${HOME}/data/popstruct_thin/05_ngsrelate/cohort_226_perchrom_adaptive_v1}"
mkdir -p "${OUT_BASE}"

# ---- Sanity check ---------------------------------------------------------
for f in "${NGSRELATE_FAST}" "${SAMPLES}" "${GENOME_WIDE_RES}"; do
    if [[ ! -e "${f}" ]]; then
        echo "[adaptive-local] ERROR: missing: ${f}" >&2
        exit 1
    fi
done

# ---- Loop over chromosomes ------------------------------------------------
for i in $(seq 1 28); do
    CHR_ID=$(printf "C_gar_LG%02d" "${i}")
    BEAGLE="${BEAGLE_DIR}/catfish.${CHR_ID}.byRF.thin_500.beagle.gz"
    FREQS="${FREQS_DIR}/freqs.${CHR_ID}.txt"
    OUT_DIR="${OUT_BASE}/${CHR_ID}"
    mkdir -p "${OUT_DIR}"

    if [[ ! -e "${BEAGLE}" || ! -e "${FREQS}" ]]; then
        echo "[adaptive-local] WARN: skipping ${CHR_ID}, missing input" >&2
        continue
    fi

    echo "[adaptive-local] ${CHR_ID} starting (workers=${WORKERS})"
    python3 -m adaptive.scheduler \
        --binary           "${NGSRELATE_FAST}" \
        --beagle           "${BEAGLE}" \
        --freqs            "${FREQS}" \
        --samples          "${SAMPLES}" \
        --genome-wide-res  "${GENOME_WIDE_RES}" \
        --chrom            "${CHR_ID}" \
        --out-dir          "${OUT_DIR}" \
        --n-samples        "${N_SAMPLES}" \
        --workers          "${WORKERS}" \
        --threads-per-pair 1
    echo "[adaptive-local] ${CHR_ID} done"
done

echo "[adaptive-local] All 28 chromosomes complete."
