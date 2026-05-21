"""
Tests for adaptive.scheduler.runner — Phase 3.

Real-binary tests are marked @pytest.mark.requires_binary and skipped
when NGSRELATE_FAST_BIN is not set. The default tests use a fake binary
(a Python script) that emits canned .res output.
"""
from __future__ import annotations
import os
import stat
import sys
import textwrap
import pytest
from pathlib import Path

from adaptive.scheduler.runner import (
    PairResult,
    RunnerError,
    load_sample_index,
    run_pair_on_chrom,
    _build_argv,
    _parse_single_pair_res,
)


HEADER_LINE = "a\tb\tida\tidb\tnSites\ttheta\tIBS0\tKING\tJ7\tJ8\tJ9"
GOOD_ROW    = "0\t1\tS1\tS2\t5000\t0.25\t0.001\t0.45\t0.10\t0.20\t0.70"


def _write_fake_binary(tmp_path, *, n_body_rows: int = 1,
                       header: str = HEADER_LINE,
                       row: str = GOOD_ROW,
                       exit_code: int = 0,
                       extra_body: str = "",
                       stderr_text: str = "[fake] ok"):
    """Write a fake binary that reads -O <out> from argv, writes a canned .res,
    prints stderr, and exits."""
    fake = tmp_path / "fake_ngsrelate"
    body_lines = "\n".join([row] * n_body_rows)
    if extra_body:
        body_lines = body_lines + "\n" + extra_body
    script = textwrap.dedent(f"""\
        #!{sys.executable}
        import sys
        args = sys.argv[1:]
        out = None
        i = 0
        while i < len(args):
            if args[i] == "-O":
                out = args[i + 1]
                i += 2
            else:
                i += 1
        sys.stderr.write({stderr_text!r} + "\\n")
        if out:
            with open(out, "w") as fh:
                fh.write({header!r} + "\\n")
                if {n_body_rows!r}:
                    fh.write({body_lines!r} + "\\n")
        sys.exit({exit_code!r})
        """)
    fake.write_text(script)
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return fake


# ---- load_sample_index ----------------------------------------------------

def test_load_sample_index(tmp_path):
    p = tmp_path / "samples.txt"
    p.write_text("S1\nS2\nS3\n")
    assert load_sample_index(p) == {"S1": 0, "S2": 1, "S3": 2}


def test_load_sample_index_skips_blank_and_comments(tmp_path):
    p = tmp_path / "samples.txt"
    p.write_text("S1\n\n# comment\nS2\nS3\n")
    idx = load_sample_index(p)
    assert idx == {"S1": 0, "S2": 3, "S3": 4}


def test_load_sample_index_duplicate_raises(tmp_path):
    p = tmp_path / "samples.txt"
    p.write_text("S1\nS2\nS1\n")
    with pytest.raises(ValueError, match="duplicate"):
        load_sample_index(p)


def test_load_sample_index_empty_raises(tmp_path):
    p = tmp_path / "samples.txt"
    p.write_text("")
    with pytest.raises(ValueError, match="empty"):
        load_sample_index(p)


# ---- _build_argv passes -D 0 ----------------------------------------------

def test_build_argv_passes_D_zero(tmp_path):
    argv = _build_argv(
        binary_path=tmp_path / "bin",
        beagle_path=tmp_path / "b.gz",
        freqs_path=tmp_path / "f",
        samples_path=tmp_path / "s",
        out_res=tmp_path / "o",
        n_samples=10,
        idx_a=0, idx_b=1, threads=1,
    )
    # Must contain "-D" followed by "0".
    i = argv.index("-D")
    assert argv[i + 1] == "0"
    assert "-a" in argv and argv[argv.index("-a") + 1] == "0"
    assert "-b" in argv and argv[argv.index("-b") + 1] == "1"


# ---- _parse_single_pair_res ----------------------------------------------

def test_parse_single_pair_res(tmp_path):
    p = tmp_path / "x.res"
    p.write_text(HEADER_LINE + "\n" + GOOD_ROW + "\n")
    header, row = _parse_single_pair_res(p)
    assert "theta" in row
    assert row["theta"] == "0.25"


def test_parse_single_pair_res_empty_raises(tmp_path):
    p = tmp_path / "empty.res"
    p.write_text("")
    with pytest.raises(RunnerError, match="empty"):
        _parse_single_pair_res(p)


def test_parse_single_pair_res_no_body_raises(tmp_path):
    p = tmp_path / "hdr.res"
    p.write_text(HEADER_LINE + "\n")
    with pytest.raises(RunnerError, match="no data row"):
        _parse_single_pair_res(p)


# ---- run_pair_on_chrom against fake binary --------------------------------

def test_run_pair_on_chrom_fake_binary(tmp_path, make_beagle, make_samples_file):
    beagle, freqs = make_beagle(n_sites=30, n_samples=3)
    samples = make_samples_file(["S1", "S2", "S3"])
    fake = _write_fake_binary(tmp_path)
    result = run_pair_on_chrom(
        binary_path=fake,
        beagle_path=beagle,
        freqs_path=freqs,
        samples_path=samples,
        sample_a="S1", sample_b="S2",
        n_samples=3, threads=1,
    )
    assert isinstance(result, PairResult)
    assert result.sample_a == "S1"
    assert result.sample_b == "S2"
    assert result.row["theta"] == "0.25"
    assert result.sites_used == 30


def test_run_pair_on_chrom_binary_nonzero(tmp_path, make_beagle, make_samples_file):
    beagle, freqs = make_beagle(n_sites=10, n_samples=3)
    samples = make_samples_file(["S1", "S2", "S3"])
    fake = _write_fake_binary(tmp_path, exit_code=2,
                              stderr_text="[fake] some error")
    with pytest.raises(RunnerError, match="exited with code 2"):
        run_pair_on_chrom(
            binary_path=fake,
            beagle_path=beagle,
            freqs_path=freqs,
            samples_path=samples,
            sample_a="S1", sample_b="S2",
            n_samples=3,
        )


def test_run_pair_on_chrom_rejects_unknown_sample(tmp_path, make_beagle, make_samples_file):
    beagle, freqs = make_beagle(n_sites=10, n_samples=3)
    samples = make_samples_file(["S1", "S2", "S3"])
    fake = _write_fake_binary(tmp_path)
    with pytest.raises(ValueError, match="sample_a"):
        run_pair_on_chrom(
            binary_path=fake,
            beagle_path=beagle, freqs_path=freqs, samples_path=samples,
            sample_a="ZZZ", sample_b="S2", n_samples=3,
        )


def test_run_pair_on_chrom_rejects_same_sample(tmp_path, make_beagle, make_samples_file):
    beagle, freqs = make_beagle(n_sites=10, n_samples=3)
    samples = make_samples_file(["S1", "S2", "S3"])
    fake = _write_fake_binary(tmp_path)
    with pytest.raises(ValueError, match="same"):
        run_pair_on_chrom(
            binary_path=fake,
            beagle_path=beagle, freqs_path=freqs, samples_path=samples,
            sample_a="S1", sample_b="S1", n_samples=3,
        )


def test_run_pair_on_chrom_missing_inputs_raise(tmp_path, make_samples_file):
    samples = make_samples_file(["S1", "S2"])
    fake = _write_fake_binary(tmp_path)
    with pytest.raises(FileNotFoundError, match="BEAGLE"):
        run_pair_on_chrom(
            binary_path=fake,
            beagle_path=tmp_path / "nope.gz",
            freqs_path=tmp_path / "nope.freqs",
            samples_path=samples,
            sample_a="S1", sample_b="S2", n_samples=2,
        )
