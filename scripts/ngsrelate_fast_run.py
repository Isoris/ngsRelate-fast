#!/usr/bin/env python3
"""
ngsrelate_fast_run.py — Run ngsRelate-fast with full input/output JSON contracts.

Writes <out>.input.json BEFORE invoking the binary (records inputs, params,
environment, exact argv). Runs the binary. Parses stderr for downsampling
stats. Writes <out>.output.json AFTER successful completion (records outputs,
stats, warnings).

Presence of <out>.output.json signals to downstream tools (ngsPedigree, etc.)
that the run completed successfully. Absence = either still running or failed.

Usage:
    python ngsrelate_fast_run.py \\
        --binary    /path/to/bin/ngsRelate-fast \\
        --beagle    cohort.beagle.gz \\
        --freqs     allele_freqs.txt \\
        --samples   samples.txt \\
        --n         226 \\
        --threads   32 \\
        --D         100000 \\
        --out       /path/to/relatedness.res \\
        --run-id    cohort_226_full_fast_v1

Produces alongside the .res:
    relatedness.res.input.json
    relatedness.res.output.json
    relatedness.res.stderr.log
"""

import argparse
import hashlib
import json
import os
import re
import socket
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

def sha256_of(path, chunk=1 << 20):
    """Stream-hash a file. Returns None if path missing."""
    if not path or not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            buf = fh.read(chunk)
            if not buf:
                break
            h.update(buf)
    return h.hexdigest()


def file_ref(path, with_hash=True):
    """Build a {path, size_bytes, mtime, sha256} block."""
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        return {"path": str(p), "missing": True}
    st = p.stat()
    ref = {
        "path": str(p.resolve()),
        "size_bytes": st.st_size,
        "mtime": datetime.fromtimestamp(st.st_mtime, tz=timezone.utc).isoformat(),
    }
    if with_hash:
        ref["sha256"] = sha256_of(str(p))
    return ref


def iso_now():
    return datetime.now(tz=timezone.utc).isoformat()


# ----------------------------------------------------------------------------
# Stderr parsing: extracts downsampling summary printed by the binary.
# Expected lines (from src/ngsRelate.cpp downsample_sites_balanced):
#   [ngsRelate-fast] Downsampling: target density 100000 sites/Gb (1.000e-04 sites/bp)
#   [ngsRelate-fast] Detected N chromosome(s) in input
#   [ngsRelate-fast]   C_gar_LG01: 35421 sites in 52.34 Mb -> keep 5234 (stride 6.77)
#   [ngsRelate-fast] Total: 951295 -> 100024 sites (genome length 1.003 Gb)
# ----------------------------------------------------------------------------

RE_TARGET   = re.compile(r"target density (\d+) sites/Gb")
RE_PERCHROM = re.compile(
    r"^\[ngsRelate-fast\]\s+(\S+):\s+(\d+) sites in ([\d.]+) Mb"
    r"\s+->\s+keep\s+(\d+)\s+\(stride\s+([\d.]+)\)"
)
RE_TOTAL    = re.compile(
    r"Total:\s+(\d+)\s+->\s+(\d+) sites\s+\(genome length\s+([\d.]+)\s+Gb\)"
)
RE_DISABLED = re.compile(r"-D 0:\s+downsampling disabled")


def parse_downsampling(stderr_text, target_D):
    """Parse the downsampling summary block from binary stderr."""
    out = {
        "enabled": target_D > 0,
        "target_sites_per_gb": target_D,
        "total_input_sites": None,
        "total_kept_sites": None,
        "total_genome_bp": None,
        "per_chromosome": [],
    }
    if RE_DISABLED.search(stderr_text):
        out["enabled"] = False
        return out

    for line in stderr_text.splitlines():
        m = RE_PERCHROM.search(line)
        if m:
            out["per_chromosome"].append({
                "chromosome": m.group(1),
                "length_bp":  int(float(m.group(3)) * 1e6),
                "input_sites": int(m.group(2)),
                "kept_sites":  int(m.group(4)),
                "stride":      float(m.group(5)),
            })
            continue
        m = RE_TOTAL.search(line)
        if m:
            out["total_input_sites"] = int(m.group(1))
            out["total_kept_sites"]  = int(m.group(2))
            out["total_genome_bp"]   = int(float(m.group(3)) * 1e9)
    return out


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary",  required=True, help="Path to ngsRelate-fast binary")
    ap.add_argument("--beagle",  required=True)
    ap.add_argument("--freqs",   required=True)
    ap.add_argument("--samples", default=None, help="Optional sample ID file (-z)")
    ap.add_argument("--n",       type=int, required=True, help="Number of samples")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--D",       type=int, default=100000, help="Sites per Gb; 0 disables")
    ap.add_argument("--out",     required=True, help="Output .res path")
    ap.add_argument("--run-id",  required=True)
    ap.add_argument("--tool-version",     default="v1.0")
    ap.add_argument("--upstream-commit",  default=None)
    ap.add_argument("--patch-file",       default=None,
                    help="Path to the .patch file (for sha256 in contract)")
    ap.add_argument("--extra", nargs="*", default=[],
                    help="Extra flags to pass through, e.g. --extra -F 1")
    args = ap.parse_args()

    out_res     = Path(args.out).resolve()
    out_dir     = out_res.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    input_json_path  = out_dir / f"{out_res.name}.input.json"
    output_json_path = out_dir / f"{out_res.name}.output.json"
    stderr_path      = out_dir / f"{out_res.name}.stderr.log"

    # ------- Build argv -------
    argv = [
        args.binary,
        "-G", args.beagle,
        "-f", args.freqs,
        "-n", str(args.n),
        "-p", str(args.threads),
        "-D", str(args.D),
        "-O", str(out_res),
    ]
    if args.samples:
        argv.extend(["-z", args.samples])
    argv.extend(args.extra)

    # ========================================================================
    # STAGE 1 — Write input contract BEFORE running.
    # If the binary dies, this file still tells you what was attempted.
    # ========================================================================
    input_contract = {
        "schema": "ngsrelate_fast.input.v1",
        "run_id": args.run_id,
        "tool": {
            "name": "ngsRelate-fast",
            "version": args.tool_version,
            "upstream_commit": args.upstream_commit,
            "patch_sha256": sha256_of(args.patch_file) if args.patch_file else None,
            "binary_path": str(Path(args.binary).resolve()),
            "binary_sha256": sha256_of(args.binary),
        },
        "inputs": {
            "beagle":  file_ref(args.beagle),
            "freqs":   file_ref(args.freqs),
            "samples": file_ref(args.samples) if args.samples else None,
        },
        "params": {
            "D_sites_per_gb": args.D,
            "n_samples": args.n,
            "threads": args.threads,
            "call_genotypes": 0,
            "estimate_inbreeding": 0,
            "three_coef_mode": 0,
            "verbose": 0,
        },
        "environment": {
            "hostname": socket.gethostname(),
            "submitted_at": iso_now(),
            "user": os.environ.get("USER"),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
            "conda_env": os.environ.get("CONDA_DEFAULT_ENV"),
            "working_dir": str(Path.cwd()),
        },
        "invocation": {
            "argv": argv,
            "stdout_path": None,
            "stderr_path": str(stderr_path),
        },
    }
    with open(input_json_path, "w") as fh:
        json.dump(input_contract, fh, indent=2)
    print(f"[ngsrelate_fast_run] Wrote input contract: {input_json_path}",
          file=sys.stderr)

    # ========================================================================
    # STAGE 2 — Run the binary.
    # ========================================================================
    print(f"[ngsrelate_fast_run] Invoking: {' '.join(argv)}", file=sys.stderr)
    start = time.time()
    with open(stderr_path, "w") as ferr:
        rc = subprocess.call(argv, stderr=ferr)
    elapsed = time.time() - start

    if rc != 0:
        print(f"[ngsrelate_fast_run] Binary exited with code {rc} after "
              f"{elapsed:.1f}s. NOT writing output contract.", file=sys.stderr)
        print(f"[ngsrelate_fast_run] See {stderr_path} for details.",
              file=sys.stderr)
        sys.exit(rc)

    # ========================================================================
    # STAGE 3 — Parse stderr, validate output, write output contract.
    # ========================================================================
    stderr_text = stderr_path.read_text()
    downsampling = parse_downsampling(stderr_text, args.D)

    warnings = []
    if not out_res.exists() or out_res.stat().st_size == 0:
        print(f"[ngsrelate_fast_run] ERROR: .res missing or empty: {out_res}",
              file=sys.stderr)
        sys.exit(2)

    with open(out_res) as fh:
        next(fh, None)  # skip header
        n_pairs = sum(1 for _ in fh)
    n_expected = args.n * (args.n - 1) // 2
    if n_pairs != n_expected:
        warnings.append(f"row count {n_pairs} != expected {n_expected}")

    output_contract = {
        "schema": "ngsrelate_fast.output.v1",
        "run_id": args.run_id,
        "input_contract": str(input_json_path),
        "status": "warn" if warnings else "ok",
        "warnings": warnings,
        "outputs": {
            "res":              file_ref(str(out_res)),
            "stderr_log":       file_ref(str(stderr_path), with_hash=False),
            "samples_sidecar":  file_ref(args.samples) if args.samples else None,
        },
        "downsampling": downsampling,
        "stats": {
            "elapsed_seconds": round(elapsed, 2),
            "n_pairs": n_pairs,
            "n_pairs_expected": n_expected,
            "peak_memory_mb": None,  # populated by SLURM script via /usr/bin/time
            "completed_at": iso_now(),
        },
    }
    with open(output_json_path, "w") as fh:
        json.dump(output_contract, fh, indent=2)
    print(f"[ngsrelate_fast_run] Wrote output contract: {output_json_path}",
          file=sys.stderr)
    print(f"[ngsrelate_fast_run] DONE in {elapsed:.1f}s, {n_pairs} pairs, "
          f"status={output_contract['status']}", file=sys.stderr)


if __name__ == "__main__":
    main()
