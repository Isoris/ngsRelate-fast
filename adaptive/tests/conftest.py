"""
conftest.py — shared pytest fixtures and markers for adaptive/.

Tests requiring the real ngsRelate-fast binary are marked
`@pytest.mark.requires_binary`. They are skipped in CI when the binary
isn't built. Set NGSRELATE_FAST_BIN to a binary path to enable them
locally.
"""

from __future__ import annotations
import gzip
import os
from pathlib import Path

import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "requires_binary: skip unless NGSRELATE_FAST_BIN is set to a real binary path")


def pytest_collection_modifyitems(config, items):
    bin_path = os.environ.get("NGSRELATE_FAST_BIN")
    skip = pytest.mark.skip(reason="NGSRELATE_FAST_BIN not set")
    has_binary = bin_path and Path(bin_path).is_file() and os.access(bin_path, os.X_OK)
    for item in items:
        if "requires_binary" in item.keywords and not has_binary:
            item.add_marker(skip)


# ----------------------------------------------------------------------------
# Synthetic BEAGLE + freqs builder — used everywhere we need a tiny input.
# ----------------------------------------------------------------------------

def _build_beagle_rows(chrom: str, n_sites: int, n_samples: int, seed: int = 0):
    """Yield BEAGLE TSV body rows (one per site) with deterministic GLs."""
    import random
    rng = random.Random(seed)
    for i in range(n_sites):
        pos = (i + 1) * 1000
        marker = f"{chrom}_{pos}"
        a1, a2 = "A", "C"
        gls = []
        for _ in range(n_samples):
            x = rng.random()
            y = rng.random() * (1 - x)
            z = 1 - x - y
            gls.extend([f"{x:.6f}", f"{y:.6f}", f"{z:.6f}"])
        yield "\t".join([marker, a1, a2, *gls])


def _build_beagle_header(n_samples: int) -> str:
    cols = ["marker", "allele1", "allele2"]
    for s in range(n_samples):
        sample = f"Ind{s}"
        cols.extend([sample, sample, sample])
    return "\t".join(cols)


@pytest.fixture
def make_beagle(tmp_path):
    """Factory: write a gzipped BEAGLE + matching freqs.

    Usage:
        beagle_path, freqs_path = make_beagle(chrom="LG01", n_sites=100, n_samples=4)
    """
    def _make(*, chrom="C_gar_LG01", n_sites=100, n_samples=4, seed=0,
              suffix=""):
        beagle = tmp_path / f"test{suffix}.{chrom}.beagle.gz"
        freqs  = tmp_path / f"test{suffix}.{chrom}.freqs"
        with gzip.open(beagle, "wt") as fh:
            fh.write(_build_beagle_header(n_samples) + "\n")
            for row in _build_beagle_rows(chrom, n_sites, n_samples, seed=seed):
                fh.write(row + "\n")
        # freqs: one allele frequency per site, in [0.05, 0.95]
        import random
        rng = random.Random(seed + 1)
        with open(freqs, "w") as fh:
            for _ in range(n_sites):
                fh.write(f"{0.05 + 0.9 * rng.random():.6f}\n")
        return beagle, freqs
    return _make


@pytest.fixture
def make_samples_file(tmp_path):
    """Factory: write a samples file (one ID per line)."""
    def _make(sample_ids, name="samples.txt"):
        p = tmp_path / name
        with open(p, "w") as fh:
            for s in sample_ids:
                fh.write(s + "\n")
        return p
    return _make


@pytest.fixture
def make_genome_wide_res(tmp_path):
    """Factory: write a synthetic genome-wide .res with chosen (sample_a, sample_b, theta, IBS0) rows.

    rows = [(ida, idb, theta, ibs0), ...]
    """
    def _make(rows, *, name="genomewide.res"):
        p = tmp_path / name
        header = ["a", "b", "ida", "idb", "nSites", "theta", "IBS0", "KING"]
        with open(p, "w") as fh:
            fh.write("\t".join(header) + "\n")
            for i, (ida, idb, theta, ibs0) in enumerate(rows):
                fh.write("\t".join([
                    str(i), str(i + 1), ida, idb, "10000",
                    f"{theta:.6f}", f"{ibs0:.6f}", f"{theta * 1.9:.6f}",
                ]) + "\n")
        return p
    return _make
