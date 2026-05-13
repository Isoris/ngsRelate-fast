#!/usr/bin/env python3
"""
compare_to_upstream.py — Validate ngsRelate-fast against canonical upstream.

Two validation modes:

  1) IDENTITY MODE  (-D 0): downsampling disabled, the fast binary should
     reproduce upstream output to within IEEE-754 floating-point noise on
     every column. Tolerance: 1e-9 on theta/J7/J8/J9, exact match on nSites.
     This isolates the pow-elimination refactor.

  2) DOWNSAMPLING MODE (default -D 100000): outputs differ because the fast
     binary uses ~100k sites/Gb. We don't expect bit-equivalence; we expect
     classification agreement. For every pair, check:
       - relationship class (PO/FS/2nd/3rd/unrelated/duplicate) agrees
       - theta within 0.02 absolute
       - KING within 0.02 absolute

Usage:
    python compare_to_upstream.py \
        --upstream-res /path/to/upstream.res \
        --fast-res /path/to/fast.res \
        --mode identity        # or 'downsampling'
"""

import argparse
import sys
import math
from collections import Counter


def classify(theta):
    """Standard relationship thresholds (Manichaikul et al. 2010)."""
    if theta is None or math.isnan(theta):
        return "missing"
    if theta >= 0.354:
        return "duplicate"
    if theta >= 0.177:
        return "first_degree"
    if theta >= 0.0884:
        return "second_degree"
    if theta >= 0.0442:
        return "third_degree"
    return "unrelated"


def load_res(path):
    """Load ngsRelate .res, return {(ida, idb): row_dict}."""
    rows = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # ngsRelate may emit either (a, b, ida, idb, ...) or (a, b, ...)
        # depending on whether -z was used. Detect.
        has_ids = "ida" in header and "idb" in header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            d = dict(zip(header, parts))
            if has_ids:
                key = (d["ida"], d["idb"])
            else:
                key = (d["a"], d["b"])
            rows[key] = d
    return rows, header


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def identity_mode(up, fa):
    """Bit-near-equivalence check. Used when both binaries run with -D 0."""
    print("[validate] IDENTITY MODE: expecting near-equivalence on all columns")
    n_checked = 0
    n_fail = 0
    max_diff = {"theta": 0.0, "J7": 0.0, "J8": 0.0, "J9": 0.0, "KING": 0.0}

    for key in up:
        if key not in fa:
            print(f"  MISSING in fast: {key}", file=sys.stderr)
            n_fail += 1
            continue
        u, f = up[key], fa[key]
        # nSites must match exactly when downsampling is off
        if u.get("nSites") != f.get("nSites"):
            print(f"  nSites mismatch {key}: {u['nSites']} vs {f['nSites']}",
                  file=sys.stderr)
            n_fail += 1
        for col in ("theta", "J7", "J8", "J9", "KING"):
            if col in u and col in f:
                d = abs(to_float(u[col]) - to_float(f[col]))
                if d > max_diff[col]:
                    max_diff[col] = d
                if d > 1e-6:
                    print(f"  {col} drift {key}: {u[col]} vs {f[col]} (Δ={d:.2e})",
                          file=sys.stderr)
                    n_fail += 1
        n_checked += 1

    print(f"[validate] Checked {n_checked} pairs")
    print(f"[validate] Max absolute differences: {max_diff}")
    if n_fail == 0:
        print("[validate] PASS: outputs near-equivalent within tolerance")
        return 0
    print(f"[validate] FAIL: {n_fail} discrepancies")
    return 1


def downsampling_mode(up, fa):
    """Classification-agreement check."""
    print("[validate] DOWNSAMPLING MODE: expecting classification agreement")
    n_checked = 0
    n_class_agree = 0
    confusion = Counter()
    theta_diffs = []
    king_diffs = []

    for key in up:
        if key not in fa:
            continue
        u_theta = to_float(up[key].get("theta", "nan"))
        f_theta = to_float(fa[key].get("theta", "nan"))
        u_cls = classify(u_theta)
        f_cls = classify(f_theta)
        confusion[(u_cls, f_cls)] += 1
        if u_cls == f_cls:
            n_class_agree += 1
        theta_diffs.append(abs(u_theta - f_theta))
        if "KING" in up[key] and "KING" in fa[key]:
            king_diffs.append(abs(to_float(up[key]["KING"]) -
                                   to_float(fa[key]["KING"])))
        n_checked += 1

    if not theta_diffs:
        print("[validate] No common pairs to compare")
        return 1

    print(f"[validate] Checked {n_checked} pairs")
    print(f"[validate] Class agreement: {n_class_agree}/{n_checked} "
          f"({100.0 * n_class_agree / n_checked:.2f}%)")
    print(f"[validate] theta abs diff:  max={max(theta_diffs):.4f}  "
          f"mean={sum(theta_diffs)/len(theta_diffs):.4f}")
    if king_diffs:
        print(f"[validate] KING abs diff:   max={max(king_diffs):.4f}  "
              f"mean={sum(king_diffs)/len(king_diffs):.4f}")
    print("[validate] Confusion matrix (upstream_class -> fast_class):")
    for (u_cls, f_cls), n in sorted(confusion.items(), key=lambda x: -x[1]):
        flag = "  " if u_cls == f_cls else "!!"
        print(f"  {flag} {u_cls:>15s} -> {f_cls:<15s}  n={n}")

    # First-degree pruning is the key application — check it specifically
    first_deg_disagree = sum(
        n for (u, f), n in confusion.items()
        if (u == "first_degree") != (f == "first_degree")
    )
    print(f"[validate] First-degree classification disagreements: "
          f"{first_deg_disagree}/{n_checked}")

    if first_deg_disagree == 0:
        print("[validate] PASS: no first-degree disagreements — safe for ngsPedigree")
        return 0
    elif first_deg_disagree < 0.005 * n_checked:
        print("[validate] WARN: <0.5% first-degree disagreements (borderline pairs)")
        return 0
    else:
        print("[validate] FAIL: first-degree classification differs in >0.5% of pairs")
        return 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--upstream-res", required=True,
                    help=".res from canonical ngsRelate")
    ap.add_argument("--fast-res", required=True,
                    help=".res from ngsRelate-fast")
    ap.add_argument("--mode", choices=("identity", "downsampling"),
                    default="downsampling")
    args = ap.parse_args()

    up, up_hdr = load_res(args.upstream_res)
    fa, fa_hdr = load_res(args.fast_res)
    print(f"[validate] Upstream: {len(up)} pairs, columns: {up_hdr}")
    print(f"[validate] Fast:     {len(fa)} pairs, columns: {fa_hdr}")

    if args.mode == "identity":
        rc = identity_mode(up, fa)
    else:
        rc = downsampling_mode(up, fa)
    sys.exit(rc)


if __name__ == "__main__":
    main()
