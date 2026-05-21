"""
runner.py — Phase 3: invoke ngsRelate-fast for a single pair.

The adaptive scheduler runs the binary one pair at a time so each pair
can use its own budgeted BEAGLE subset. Upstream ngsRelate v2 supports
this via the `-a` and `-b` flags (0-indexed individuals into the BEAGLE);
when both are set the binary computes only that pair.

The first time the scheduler calls into this module on a real run it
should invoke `probe_binary_single_pair_support()` against a tiny BEAGLE
to confirm the fork still honors the `-a`/`-b` flags. The fork doesn't
touch that path (we only added `-D` and the algebraic refactor), but
"should" is not the same as "does" — see IMPLEMENTATION_PLAN.md §2 Phase 3.

The runner always passes `-D 0` to the binary because the scheduler is
doing the subsetting itself; letting the binary downsample again would
double-thin the input (SPEC §7 OQ7).
"""

from __future__ import annotations
import os
import shlex
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional


class RunnerError(Exception):
    """ngsRelate-fast binary failed for a single pair invocation."""


class BinaryDoesNotSupportSinglePair(Exception):
    """Probe found the binary ignoring -a/-b. Fork compatibility broken."""


@dataclass
class PairResult:
    """One pair's per-chrom output. `row` is the full .res row as a dict
    keyed by column name; `header` is the column order from the binary."""
    sample_a: str
    sample_b: str
    header: List[str]
    row: Dict[str, str]
    stderr: str
    elapsed_seconds: float
    sites_used: int   # what we asked the binary to use (= len(BEAGLE body))


def load_sample_index(samples_path) -> Dict[str, int]:
    """Read a one-ID-per-line samples file into {sample_id: 0-based index}."""
    out: Dict[str, int] = {}
    with open(samples_path) as fh:
        for i, line in enumerate(fh):
            sid = line.strip()
            if not sid or sid.startswith("#"):
                continue
            if sid in out:
                raise ValueError(
                    f"duplicate sample ID '{sid}' at line {i + 1} of {samples_path}")
            out[sid] = i
    if not out:
        raise ValueError(f"samples file is empty: {samples_path}")
    return out


def _count_beagle_body_rows(beagle_path: Path) -> int:
    """Count BEAGLE body rows; used to record sites_used in the result."""
    import gzip
    opener = gzip.open if str(beagle_path).endswith(".gz") else open
    with opener(beagle_path, "rt") as fh:
        next(fh, None)
        return sum(1 for _ in fh)


def _build_argv(
    binary_path: Path,
    beagle_path: Path,
    freqs_path: Path,
    samples_path: Path,
    out_res: Path,
    n_samples: int,
    idx_a: int,
    idx_b: int,
    threads: int,
) -> List[str]:
    return [
        str(binary_path),
        "-G", str(beagle_path),
        "-f", str(freqs_path),
        "-z", str(samples_path),
        "-n", str(n_samples),
        "-p", str(threads),
        "-D", "0",
        "-a", str(idx_a),
        "-b", str(idx_b),
        "-O", str(out_res),
    ]


def _parse_single_pair_res(res_path: Path) -> tuple[List[str], Dict[str, str]]:
    with open(res_path) as fh:
        header_line = fh.readline()
        body_line   = fh.readline()
    if not header_line:
        raise RunnerError(f".res from binary is empty: {res_path}")
    if not body_line:
        raise RunnerError(
            f".res from binary has no data row (binary may have run all pairs "
            f"and emitted nothing for the requested single pair): {res_path}")
    header = header_line.rstrip("\n").split("\t")
    body   = body_line.rstrip("\n").split("\t")
    if len(header) != len(body):
        raise RunnerError(
            f"binary .res row has {len(body)} columns but header has {len(header)} "
            f"at {res_path}")
    return header, dict(zip(header, body))


def run_pair_on_chrom(
    *,
    binary_path,
    beagle_path,
    freqs_path,
    samples_path,
    sample_a: str,
    sample_b: str,
    n_samples: int,
    threads: int = 1,
    tmpdir: Optional[Path] = None,
    timeout: Optional[float] = None,
    sample_index: Optional[Dict[str, int]] = None,
) -> PairResult:
    """Run ngsRelate-fast for a single (sample_a, sample_b) pair.

    Args:
        binary_path:   path to ngsRelate-fast.
        beagle_path:   per-pair BEAGLE (usually from BeagleSubsetCache).
        freqs_path:    matching freqs file.
        samples_path:  one-ID-per-line file matching the BEAGLE columns.
        sample_a:      first sample ID (must appear in samples_path).
        sample_b:      second sample ID (must appear in samples_path).
        n_samples:     total number of samples in the BEAGLE (the `-n` flag).
        threads:       threads to give the binary itself (`-p`).
        tmpdir:        optional directory for the temporary .res output.
        timeout:       seconds; passed through to subprocess.run.
        sample_index:  precomputed {id: idx} to skip re-reading samples_path
                       for every pair in a tight loop.

    Returns:
        PairResult with the parsed .res row + diagnostics.

    Raises:
        RunnerError:                  binary returned non-zero, or produced
                                      malformed .res output.
        FileNotFoundError:            binary/beagle/freqs/samples missing.
        ValueError:                   sample IDs not in samples_path.
    """
    binary_path  = Path(binary_path)
    beagle_path  = Path(beagle_path)
    freqs_path   = Path(freqs_path)
    samples_path = Path(samples_path)

    for p, label in [(binary_path, "binary"),
                     (beagle_path, "BEAGLE"),
                     (freqs_path,  "freqs"),
                     (samples_path, "samples")]:
        if not p.exists():
            raise FileNotFoundError(f"{label} not found: {p}")

    idx_map = sample_index or load_sample_index(samples_path)
    if sample_a not in idx_map:
        raise ValueError(f"sample_a '{sample_a}' not in {samples_path}")
    if sample_b not in idx_map:
        raise ValueError(f"sample_b '{sample_b}' not in {samples_path}")
    idx_a, idx_b = idx_map[sample_a], idx_map[sample_b]
    if idx_a == idx_b:
        raise ValueError(f"sample_a and sample_b are the same: {sample_a!r}")

    # Temp .res output. Use a real tempfile so concurrent workers don't collide.
    tmp_kw = {"prefix": f"ngsrel_pair_{idx_a}_{idx_b}_", "suffix": ".res"}
    if tmpdir is not None:
        Path(tmpdir).mkdir(parents=True, exist_ok=True)
        tmp_kw["dir"] = str(tmpdir)
    fd, out_res_str = tempfile.mkstemp(**tmp_kw)
    os.close(fd)
    out_res = Path(out_res_str)

    argv = _build_argv(binary_path, beagle_path, freqs_path, samples_path,
                       out_res, n_samples, idx_a, idx_b, threads)

    import time
    t0 = time.time()
    try:
        proc = subprocess.run(
            argv,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as e:
        out_res.unlink(missing_ok=True)
        raise RunnerError(
            f"binary timed out after {timeout}s on pair ({sample_a}, {sample_b}): "
            f"argv={shlex.join(argv)}") from e
    elapsed = time.time() - t0

    stderr_text = proc.stderr.decode("utf-8", errors="replace")
    if proc.returncode != 0:
        out_res.unlink(missing_ok=True)
        raise RunnerError(
            f"binary exited with code {proc.returncode} for pair "
            f"({sample_a}, {sample_b})\nargv: {shlex.join(argv)}\nstderr:\n{stderr_text}")

    try:
        header, row = _parse_single_pair_res(out_res)
    finally:
        out_res.unlink(missing_ok=True)

    sites_used = _count_beagle_body_rows(beagle_path)

    return PairResult(
        sample_a=sample_a,
        sample_b=sample_b,
        header=header,
        row=row,
        stderr=stderr_text,
        elapsed_seconds=elapsed,
        sites_used=sites_used,
    )


def probe_binary_single_pair_support(
    *,
    binary_path,
    beagle_path,
    freqs_path,
    samples_path,
    n_samples: int,
    timeout: float = 30.0,
) -> None:
    """Sanity check: confirm the binary honors `-a`/`-b` for single-pair runs.

    Picks (0, 1) and runs the binary. If the output .res contains more than
    one row, the binary is not honoring the flags and the scheduler must
    refuse to run. This is meant to be called once before launching the
    worker pool.

    Raises:
        BinaryDoesNotSupportSinglePair: -a/-b flags ignored.
        RunnerError:                    other binary failure.
    """
    # Read first two sample IDs from samples_path to call into run_pair_on_chrom.
    idx_map = load_sample_index(samples_path)
    if len(idx_map) < 2:
        raise ValueError("probe requires at least 2 samples in samples_path")
    by_idx = sorted(idx_map.items(), key=lambda kv: kv[1])
    sample_a, sample_b = by_idx[0][0], by_idx[1][0]

    res = run_pair_on_chrom(
        binary_path=binary_path,
        beagle_path=beagle_path,
        freqs_path=freqs_path,
        samples_path=samples_path,
        sample_a=sample_a,
        sample_b=sample_b,
        n_samples=n_samples,
        threads=1,
        timeout=timeout,
        sample_index=idx_map,
    )

    # _parse_single_pair_res reads only the first body row. Re-read to count
    # all rows: if the binary ignored -a/-b, the .res would have n*(n-1)/2 rows
    # — but we deleted the file. So we use a different probe path: invoke the
    # binary again with -O pointing to a path we keep, and count rows ourselves.
    import tempfile, os, shlex, subprocess
    fd, probe_out_str = tempfile.mkstemp(prefix="ngsrel_probe_", suffix=".res")
    os.close(fd)
    probe_out = Path(probe_out_str)
    argv = _build_argv(
        Path(binary_path), Path(beagle_path), Path(freqs_path),
        Path(samples_path), probe_out, n_samples,
        idx_map[sample_a], idx_map[sample_b], 1,
    )
    proc = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          timeout=timeout, check=False)
    if proc.returncode != 0:
        probe_out.unlink(missing_ok=True)
        raise RunnerError(
            f"probe invocation failed: {shlex.join(argv)}\n"
            f"stderr:\n{proc.stderr.decode('utf-8', errors='replace')}")
    try:
        with open(probe_out) as fh:
            n_rows = sum(1 for _ in fh) - 1   # excluding header
    finally:
        probe_out.unlink(missing_ok=True)

    if n_rows > 1:
        raise BinaryDoesNotSupportSinglePair(
            f"binary emitted {n_rows} pair rows for a single -a/-b request; "
            f"the fork no longer honors single-pair invocation. The adaptive "
            f"scheduler cannot run safely against this binary.")
