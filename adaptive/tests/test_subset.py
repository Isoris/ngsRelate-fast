"""
Tests for adaptive.scheduler.subset — Phase 2.
"""
from __future__ import annotations
import gzip
import pytest

from adaptive.scheduler.subset import (
    BeagleSubsetCache,
    select_indices,
    _count_beagle_sites,
)


# ---- select_indices stride correctness ------------------------------------

def test_select_indices_evenly_divisible():
    # N=100, B=10, stride=10 → indices 0,10,20,...,90
    assert select_indices(100, 10) == [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]


def test_select_indices_budget_equals_total():
    assert select_indices(50, 50) == list(range(50))


def test_select_indices_budget_greater_than_total():
    assert select_indices(40, 100) == list(range(40))


def test_select_indices_uneven():
    # N=99, B=10, stride=9.9 — keeps exactly 10 indices, evenly spread
    out = select_indices(99, 10)
    assert len(out) == 10
    assert out[0] == 0
    # First skip should be at i=1..9, next pick at i=10 (since 10 >= 9.9)
    assert out[1] == 10
    # Last index should be ≤ 99 - 1
    assert out[-1] < 99


def test_select_indices_budget_one():
    assert select_indices(1000, 1) == [0]


def test_select_indices_budget_zero_raises():
    with pytest.raises(ValueError):
        select_indices(100, 0)


def test_select_indices_negative_budget_raises():
    with pytest.raises(ValueError):
        select_indices(100, -5)


# ---- BeagleSubsetCache end-to-end ----------------------------------------

def _read_markers(path):
    """Read markers (column 0) from a gzipped BEAGLE body."""
    with gzip.open(path, "rt") as fh:
        next(fh)  # header
        return [line.split("\t", 1)[0] for line in fh]


def test_cache_returns_original_when_budget_exceeds_total(make_beagle):
    beagle, freqs = make_beagle(n_sites=100, n_samples=3)
    with BeagleSubsetCache(beagle, freqs) as cache:
        b, f = cache.get(1000)
        # Returned paths should be the originals — no subsetting.
        assert b == beagle.resolve()
        assert f == freqs.resolve()


def test_cache_subsets_to_exact_budget(make_beagle):
    beagle, freqs = make_beagle(n_sites=100, n_samples=3)
    with BeagleSubsetCache(beagle, freqs) as cache:
        sub_b, sub_f = cache.get(10)
        markers = _read_markers(sub_b)
        assert len(markers) == 10
        with open(sub_f) as fh:
            freq_lines = fh.readlines()
        assert len(freq_lines) == 10


def test_cache_is_deterministic(make_beagle):
    beagle, freqs = make_beagle(n_sites=200, n_samples=3)
    with BeagleSubsetCache(beagle, freqs) as cache:
        sub_b, _ = cache.get(20)
        markers = _read_markers(sub_b)
    with BeagleSubsetCache(beagle, freqs) as cache2:
        sub_b2, _ = cache2.get(20)
        markers2 = _read_markers(sub_b2)
    assert markers == markers2


def test_cache_second_call_returns_cached_paths(make_beagle):
    beagle, freqs = make_beagle(n_sites=100, n_samples=3)
    with BeagleSubsetCache(beagle, freqs) as cache:
        first_b, first_f = cache.get(10)
        first_b_mtime = first_b.stat().st_mtime_ns
        # second call must return the SAME path and not rewrite
        second_b, second_f = cache.get(10)
        assert second_b == first_b
        assert second_f == first_f
        assert second_b.stat().st_mtime_ns == first_b_mtime
        assert cache.budgets_built() == [10]


def test_cache_subsets_beagle_and_freqs_with_identical_indices(make_beagle):
    """The whole point of the cache: BEAGLE and freqs are aligned by site."""
    beagle, freqs = make_beagle(n_sites=100, n_samples=2)
    # Read original freqs once.
    with open(freqs) as fh:
        original_freqs = [line.strip() for line in fh]
    with BeagleSubsetCache(beagle, freqs) as cache:
        sub_b, sub_f = cache.get(10)
        markers = _read_markers(sub_b)
        with open(sub_f) as fh:
            sub_freqs = [line.strip() for line in fh]
    # Indices kept should be [0, 10, 20, ..., 90].
    assert markers[0].endswith("_1000")
    assert markers[1].endswith("_11000")
    assert sub_freqs[0] == original_freqs[0]
    assert sub_freqs[1] == original_freqs[10]


def test_cache_cleanup_removes_files(make_beagle, tmp_path):
    beagle, freqs = make_beagle(n_sites=100, n_samples=3)
    cache = BeagleSubsetCache(beagle, freqs)
    cache.get(10)
    cache.get(20)
    tmpdir = cache.tmpdir
    assert tmpdir.exists()
    files = list(tmpdir.iterdir())
    assert len(files) >= 4   # 2 budgets × (beagle + freqs)
    cache.cleanup()
    assert not tmpdir.exists()


def test_cache_rejects_misaligned_inputs(tmp_path, make_beagle):
    beagle, freqs = make_beagle(n_sites=100, n_samples=2)
    # write a freqs file with the wrong row count
    bad_freqs = tmp_path / "bad.freqs"
    bad_freqs.write_text("0.5\n" * 99)
    cache = BeagleSubsetCache(beagle, bad_freqs)
    with pytest.raises(ValueError, match="mismatch"):
        cache.n_total
    cache.cleanup()


def test_count_beagle_sites(make_beagle):
    beagle, _ = make_beagle(n_sites=42, n_samples=3)
    assert _count_beagle_sites(beagle) == 42
