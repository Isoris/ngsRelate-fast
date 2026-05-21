"""
subset.py — Phase 2: per-pair BEAGLE + freqs subsetting cache.

The adaptive scheduler hands the binary a thinned BEAGLE per pair so each
pair runs against exactly its budgeted number of sites. To avoid writing
the same subset N times, we cache by budget per chromosome.

BEAGLE and freqs MUST be subsetted with the same site indices, otherwise
the binary's per-site freq lookup desynchronises and silently produces
garbage. The cache returns a (beagle_path, freqs_path) tuple to enforce
this — never call `.get(budget)` and pass the original freqs file
alongside a subsetted BEAGLE.

The stride logic matches the convention used by the binary's `-D` flag
(see patches/01_ngsRelate_fast.patch, downsample_sites_balanced): for a
target of B sites out of N total, `stride = N/B` and we keep site i when
`i >= next_pick`, then advance next_pick by stride. This means the
scheduler's subset for budget B is bit-identical to what the binary would
produce with `-D` set to deliver B sites on the same input — so the
scheduler can safely pass `-D 0` to the binary (no double-downsampling).
"""

from __future__ import annotations
import gzip
import os
import shutil
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def select_indices(n_total: int, budget: int) -> List[int]:
    """Indices to keep when subsetting `n_total` sites down to `budget`.

    Matches the binary's downsample_sites_balanced stride convention.
    If `budget >= n_total`, returns `list(range(n_total))`. If `budget <= 0`,
    raises ValueError — a zero budget has no defensible interpretation here.
    """
    if budget <= 0:
        raise ValueError(f"budget must be > 0, got {budget}")
    if budget >= n_total:
        return list(range(n_total))

    stride = n_total / budget
    next_pick = 0.0
    out: List[int] = []
    for i in range(n_total):
        if float(i) >= next_pick:
            out.append(i)
            next_pick += stride
            if len(out) == budget:
                break
    return out


def _open_text(path):
    """Open .gz or plain text transparently for reading."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _count_beagle_sites(beagle_path: Path) -> int:
    """Count BEAGLE body rows (header excluded)."""
    n = 0
    with _open_text(beagle_path) as fh:
        next(fh, None)   # header
        for _ in fh:
            n += 1
    return n


def _count_freqs_rows(freqs_path: Path) -> int:
    n = 0
    with open(freqs_path) as fh:
        for _ in fh:
            n += 1
    return n


class BeagleSubsetCache:
    """Cache thinned BEAGLE + freqs copies, one pair per budget tier.

    Construct once per chromosome before launching the worker pool. The
    workers each open returned paths read-only; no shared writable state.

    Call `.cleanup()` (or use as a context manager) to remove temp files.
    """

    def __init__(
        self,
        beagle_path,
        freqs_path,
        tmpdir: Optional[Path] = None,
        *,
        chrom_label: str = "chrom",
    ):
        self.beagle_path = Path(beagle_path).resolve()
        self.freqs_path  = Path(freqs_path).resolve()
        if not self.beagle_path.exists():
            raise FileNotFoundError(f"BEAGLE not found: {self.beagle_path}")
        if not self.freqs_path.exists():
            raise FileNotFoundError(f"freqs not found: {self.freqs_path}")
        self.chrom_label = chrom_label

        self._owns_tmpdir = tmpdir is None
        self.tmpdir = Path(tmpdir) if tmpdir else Path(
            tempfile.mkdtemp(prefix=f"ngsrelate_adaptive_{chrom_label}_"))
        self.tmpdir.mkdir(parents=True, exist_ok=True)

        # Site counts. We require BEAGLE and freqs to align.
        self._n_total: Optional[int] = None
        self._cache: Dict[int, Tuple[Path, Path]] = {}

    # -- public API -----------------------------------------------------------

    @property
    def n_total(self) -> int:
        if self._n_total is None:
            n_b = _count_beagle_sites(self.beagle_path)
            n_f = _count_freqs_rows(self.freqs_path)
            if n_b != n_f:
                raise ValueError(
                    f"BEAGLE/freqs site count mismatch for {self.chrom_label}: "
                    f"BEAGLE has {n_b} sites, freqs has {n_f} rows. These must align.")
            self._n_total = n_b
        return self._n_total

    def get(self, budget: int) -> Tuple[Path, Path]:
        """Return (beagle_subset, freqs_subset) for the requested budget."""
        if budget <= 0:
            raise ValueError(f"budget must be > 0, got {budget}")

        n_total = self.n_total
        # No subsetting needed: hand back the originals.
        if budget >= n_total:
            return self.beagle_path, self.freqs_path

        if budget in self._cache:
            return self._cache[budget]

        kept = set(select_indices(n_total, budget))
        out_beagle = self.tmpdir / f"subset_b{budget}.beagle.gz"
        out_freqs  = self.tmpdir / f"subset_b{budget}.freqs"

        # Stream BEAGLE: copy header, copy rows whose 0-based index is in `kept`.
        with _open_text(self.beagle_path) as fin, \
             gzip.open(out_beagle, "wt") as fout:
            header = next(fin, None)
            if header is None:
                raise ValueError(f"BEAGLE has no header: {self.beagle_path}")
            fout.write(header)
            for i, line in enumerate(fin):
                if i in kept:
                    fout.write(line)

        # Stream freqs.
        with open(self.freqs_path) as fin, open(out_freqs, "w") as fout:
            for i, line in enumerate(fin):
                if i in kept:
                    fout.write(line)

        self._cache[budget] = (out_beagle, out_freqs)
        return out_beagle, out_freqs

    def budgets_built(self) -> List[int]:
        return sorted(self._cache.keys())

    def cleanup(self):
        if self._owns_tmpdir and self.tmpdir.exists():
            shutil.rmtree(self.tmpdir, ignore_errors=True)
        self._cache.clear()

    # -- context manager ------------------------------------------------------

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.cleanup()
        return False
