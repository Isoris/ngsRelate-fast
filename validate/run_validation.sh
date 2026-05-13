#!/usr/bin/env bash
###############################################################################
# run_validation.sh — Run both binaries on the same input and diff.
#
# Usage:
#   bash run_validation.sh <beagle.gz> <freqs> <n_samples> [samples_file]
#
# Produces two .res files in validate/output/ and runs the Python comparator
# in both identity mode (-D 0) and downsampling mode (-D 100000).
#
# RECOMMENDED FIRST RUN: use a SMALL subset (e.g. one chromosome) to confirm
# identity mode passes. Only then trust the fast binary for full-cohort use.
###############################################################################
set -euo pipefail

BEAGLE="${1:?usage: $0 <beagle.gz> <freqs> <n_samples> [samples_file]}"
FREQS="${2:?usage}"
N="${3:?usage}"
SAMPLES="${4:-}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${REPO_ROOT}/validate/output"
mkdir -p "${OUT}"

UPSTREAM="${REPO_ROOT}/bin/ngsRelate-upstream"
FAST="${REPO_ROOT}/bin/ngsRelate-fast"

if [[ ! -x "${UPSTREAM}" || ! -x "${FAST}" ]]; then
    echo "ERROR: binaries missing. Run build.sh first." >&2
    exit 1
fi

Z_FLAG=""
if [[ -n "${SAMPLES}" ]]; then
    Z_FLAG="-z ${SAMPLES}"
fi

# ----------------------------------------------------------------------------
# Stage 1: IDENTITY MODE — both binaries with downsampling disabled
# This isolates the pow-elimination patch
# ----------------------------------------------------------------------------
echo "[stage 1] Identity mode (-D 0 on fast, upstream as-is)"
"${UPSTREAM}" -G "${BEAGLE}" -f "${FREQS}" -n "${N}" -p 8 ${Z_FLAG} \
    -O "${OUT}/upstream.res"
"${FAST}"     -G "${BEAGLE}" -f "${FREQS}" -n "${N}" -p 8 ${Z_FLAG} -D 0 \
    -O "${OUT}/fast_D0.res"

python3 "${REPO_ROOT}/validate/compare_to_upstream.py" \
    --upstream-res "${OUT}/upstream.res" \
    --fast-res     "${OUT}/fast_D0.res" \
    --mode identity

# ----------------------------------------------------------------------------
# Stage 2: DOWNSAMPLING MODE — fast with default -D 100000
# This tests that downsampling preserves classification
# ----------------------------------------------------------------------------
echo ""
echo "[stage 2] Downsampling mode (fast with -D 100000)"
"${FAST}" -G "${BEAGLE}" -f "${FREQS}" -n "${N}" -p 8 ${Z_FLAG} -D 100000 \
    -O "${OUT}/fast_D100k.res"

python3 "${REPO_ROOT}/validate/compare_to_upstream.py" \
    --upstream-res "${OUT}/upstream.res" \
    --fast-res     "${OUT}/fast_D100k.res" \
    --mode downsampling

echo ""
echo "[done] Validation outputs in: ${OUT}"
