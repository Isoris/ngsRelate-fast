#!/usr/bin/env bash
###############################################################################
# build.sh — Clone upstream ngsRelate, apply patch, compile.
#
# Produces:
#   ./bin/ngsRelate-upstream    (canonical reference, for validation diffs)
#   ./bin/ngsRelate-fast        (patched binary, default tool going forward)
#
# Run from the repo root:
#   bash build.sh
#
# On LANTA: activate the "assembly" conda env first for htslib if you need
# BCF support; for BEAGLE-only workflows the system toolchain suffices.
###############################################################################
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${REPO_ROOT}/build"
BIN_DIR="${REPO_ROOT}/bin"
UPSTREAM_URL="https://github.com/ANGSD/NgsRelate.git"
UPSTREAM_COMMIT=""   # leave empty for HEAD; pin once you've validated once

mkdir -p "${BUILD_DIR}" "${BIN_DIR}"

# ----------------------------------------------------------------------------
# 1. Get a pristine upstream clone
# ----------------------------------------------------------------------------
if [[ ! -d "${BUILD_DIR}/NgsRelate-upstream" ]]; then
    echo "[build] Cloning upstream..."
    git clone "${UPSTREAM_URL}" "${BUILD_DIR}/NgsRelate-upstream"
fi

cd "${BUILD_DIR}/NgsRelate-upstream"
git fetch
git reset --hard
git clean -fdx
if [[ -n "${UPSTREAM_COMMIT}" ]]; then
    git checkout "${UPSTREAM_COMMIT}"
fi
echo "[build] Upstream at $(git rev-parse --short HEAD)"
cd "${REPO_ROOT}"

# ----------------------------------------------------------------------------
# 2. Build the unpatched reference binary
# ----------------------------------------------------------------------------
echo "[build] Building upstream reference binary..."
cd "${BUILD_DIR}/NgsRelate-upstream"
make clean >/dev/null 2>&1 || true
make
cp ngsRelate "${BIN_DIR}/ngsRelate-upstream"
cd "${REPO_ROOT}"

# ----------------------------------------------------------------------------
# 3. Make a patched copy
# ----------------------------------------------------------------------------
PATCHED_DIR="${BUILD_DIR}/NgsRelate-fast"
rm -rf "${PATCHED_DIR}"
cp -r "${BUILD_DIR}/NgsRelate-upstream" "${PATCHED_DIR}"

cd "${PATCHED_DIR}"
cp ngsRelate.cpp ngsRelate.cpp.orig
patch -p0 < "${REPO_ROOT}/patches/01_ngsRelate_fast.patch"

# ----------------------------------------------------------------------------
# 4. Add the -D flag parser to the getopt loop in main()
# ----------------------------------------------------------------------------
# The patch adds the -D help text and the downsampling function, but the
# getopt() call string in main() needs "D:" appended and a case branch.
# Done here with sed to keep the patch readable.
python3 - <<'PYEOF'
import re, sys, pathlib
p = pathlib.Path("ngsRelate.cpp")
src = p.read_text()

# Find the getopts string in main() and add D:
m = re.search(r'getopt\(argc,\s*argv,\s*"([^"]+)"', src)
if not m:
    sys.exit("ERROR: could not locate getopt() in main()")
opts = m.group(1)
if "D:" not in opts:
    new_opts = opts + "D:"
    src = src.replace(m.group(0), m.group(0).replace(opts, new_opts))

# Add the case branch. We insert it next to an existing case as anchor.
anchor = "case 'p': num_threads = atoi(optarg); break;"
case_D = "case 'D': target_sites_per_gb = atoi(optarg); break;"
if case_D not in src:
    if anchor not in src:
        sys.exit("ERROR: could not find anchor case 'p' in main()")
    src = src.replace(anchor, anchor + "\n      " + case_D)

# Hook the downsampler call right before the pair work begins. The original
# main has a line like "fprintf(stderr,\"\\t-> nsites:\" ..." after data load;
# we add the call before it. The exact location depends on upstream version;
# this regex finds the first "FINISHED=0" or pair dispatch line.
hook_anchor = "FINISHED=0;"
hook_call = """
  // ===== ngsRelate-fast: downsample BEFORE pair work =====
  if (target_sites_per_gb > 0 && overall_number_of_sites > 0) {
      overall_number_of_sites = downsample_sites_balanced(
          &gls, &freqs, site_ids, overall_number_of_sites, nind);
  }
  // =======================================================
"""
# NOTE: gls / freqs / site_ids / nind names are conventional in upstream;
# verify against the actual main() variable names after first build attempt.
# If your upstream version uses different names, edit this block.
if "downsample_sites_balanced" not in src.replace(case_D, ""):
    if hook_anchor in src:
        src = src.replace(hook_anchor, hook_call + "\n  " + hook_anchor, 1)
    else:
        sys.stderr.write("WARN: could not find FINISHED=0 anchor; you'll need to\n")
        sys.stderr.write("      manually insert the downsample_sites_balanced() call\n")
        sys.stderr.write("      in main() after data load.\n")

p.write_text(src)
print("[build] Patched main() (getopt + case + hook)")
PYEOF

# ----------------------------------------------------------------------------
# 5. Build the patched binary
# ----------------------------------------------------------------------------
echo "[build] Building patched binary..."
make clean >/dev/null 2>&1 || true
make
cp ngsRelate "${BIN_DIR}/ngsRelate-fast"

cd "${REPO_ROOT}"
echo ""
echo "[build] DONE"
echo "  Upstream:  ${BIN_DIR}/ngsRelate-upstream"
echo "  Patched:   ${BIN_DIR}/ngsRelate-fast"
echo ""
echo "Next: run validation"
echo "  bash validate/run_validation.sh <your.beagle.gz> <your.freqs> <n_samples>"
