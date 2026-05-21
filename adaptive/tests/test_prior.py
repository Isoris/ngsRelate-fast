"""
Tests for adaptive.scheduler.prior — Phase 1.
"""
from __future__ import annotations
import pytest

from adaptive.scheduler.prior import (
    derive_priors,
    _classify_one,
    PriorDerivationError,
)
from adaptive.scheduler.edge_class import EdgeClass
from adaptive.scheduler.config import (
    AMBIGUOUS_FIRST_DEGREE_IBS0_BAND,
    KING_THETA_DUPLICATE,
    KING_THETA_FIRST_DEGREE,
    KING_THETA_SECOND_DEG,
    KING_THETA_THIRD_DEG,
)


BAND = AMBIGUOUS_FIRST_DEGREE_IBS0_BAND


# ---- _classify_one boundary table -----------------------------------------

@pytest.mark.parametrize("theta, ibs0, expected", [
    # Duplicate boundary
    (0.50,                       0.000, EdgeClass.DUPLICATE_OR_CLONE),
    (KING_THETA_DUPLICATE,       0.000, EdgeClass.DUPLICATE_OR_CLONE),
    (KING_THETA_DUPLICATE - 1e-6, 0.001, EdgeClass.PARENT_OFFSPRING),

    # First-degree band: PO vs FS via IBS0
    (0.25,                       0.001, EdgeClass.PARENT_OFFSPRING),
    (0.25,                       0.020, EdgeClass.FULL_SIBLING),
    (KING_THETA_FIRST_DEGREE,    0.001, EdgeClass.PARENT_OFFSPRING),

    # Ambiguous first-degree band
    (0.25,                       0.008, EdgeClass.AMBIGUOUS_FIRST_DEGREE),
    (0.25,                       0.0070, EdgeClass.AMBIGUOUS_FIRST_DEGREE),
    (0.25,                       0.0090, EdgeClass.AMBIGUOUS_FIRST_DEGREE),

    # Second / third degree bands
    (KING_THETA_SECOND_DEG,      0.05,  EdgeClass.SECOND_DEGREE),
    (KING_THETA_FIRST_DEGREE - 1e-6, 0.05, EdgeClass.SECOND_DEGREE),
    (KING_THETA_THIRD_DEG,       0.10,  EdgeClass.THIRD_DEGREE),
    (KING_THETA_SECOND_DEG - 1e-6, 0.10, EdgeClass.THIRD_DEGREE),

    # Unrelated
    (0.0,                        0.20,  EdgeClass.UNRELATED),
    (KING_THETA_THIRD_DEG - 1e-6, 0.20, EdgeClass.UNRELATED),
    (-0.05,                      0.30,  EdgeClass.UNRELATED),
])
def test_classify_boundaries(theta, ibs0, expected):
    assert _classify_one(theta, ibs0, BAND) is expected


def test_classify_first_degree_no_ibs0_is_ambiguous():
    assert _classify_one(0.25, None, BAND) is EdgeClass.AMBIGUOUS_FIRST_DEGREE


def test_classify_nan_theta_is_unrelated():
    assert _classify_one(float("nan"), 0.10, BAND) is EdgeClass.UNRELATED


# ---- derive_priors end-to-end ---------------------------------------------

def test_derive_priors_full_table(make_genome_wide_res):
    res = make_genome_wide_res([
        ("S1", "S2", 0.40,  0.000),  # duplicate
        ("S1", "S3", 0.25,  0.001),  # PO
        ("S1", "S4", 0.25,  0.020),  # FS
        ("S2", "S3", 0.25,  0.008),  # ambiguous
        ("S3", "S4", 0.10,  0.05),   # 2nd
        ("S4", "S5", 0.06,  0.10),   # 3rd
        ("S5", "S6", 0.01,  0.20),   # unrelated
    ])
    priors = derive_priors(res)
    assert priors[("S1", "S2")] is EdgeClass.DUPLICATE_OR_CLONE
    assert priors[("S1", "S3")] is EdgeClass.PARENT_OFFSPRING
    assert priors[("S1", "S4")] is EdgeClass.FULL_SIBLING
    assert priors[("S2", "S3")] is EdgeClass.AMBIGUOUS_FIRST_DEGREE
    assert priors[("S3", "S4")] is EdgeClass.SECOND_DEGREE
    assert priors[("S4", "S5")] is EdgeClass.THIRD_DEGREE
    assert priors[("S5", "S6")] is EdgeClass.UNRELATED


def test_derive_priors_canonicalizes_pair_keys(make_genome_wide_res):
    # write with (idb < ida); the prior map must still key on sorted tuple.
    res = make_genome_wide_res([
        ("Z9", "A1", 0.01, 0.2),
    ])
    priors = derive_priors(res)
    assert ("A1", "Z9") in priors
    assert ("Z9", "A1") not in priors


def test_derive_priors_ambiguous_pair_at_threshold(make_genome_wide_res):
    res = make_genome_wide_res([
        ("S1", "S2", 0.20, 0.0080),
        ("S3", "S4", 0.20, 0.0064),
        ("S5", "S6", 0.20, 0.0096),
    ])
    priors = derive_priors(res)
    assert priors[("S1", "S2")] is EdgeClass.AMBIGUOUS_FIRST_DEGREE
    assert priors[("S3", "S4")] is EdgeClass.AMBIGUOUS_FIRST_DEGREE
    assert priors[("S5", "S6")] is EdgeClass.AMBIGUOUS_FIRST_DEGREE


def test_derive_priors_empty_res_raises(tmp_path):
    p = tmp_path / "empty.res"
    p.write_text("")
    with pytest.raises(PriorDerivationError, match="empty"):
        derive_priors(p)


def test_derive_priors_missing_theta_column_raises(tmp_path):
    p = tmp_path / "no_theta.res"
    p.write_text("a\tb\tKING\n0\t1\t0.4\n")
    with pytest.raises(PriorDerivationError, match="theta"):
        derive_priors(p)


def test_derive_priors_malformed_row_raises(tmp_path):
    p = tmp_path / "bad.res"
    p.write_text("a\tb\tida\tidb\tnSites\ttheta\tIBS0\tKING\n"
                 "0\t1\tS1\tS2\t10000\t0.25\n")  # too few columns
    with pytest.raises(PriorDerivationError, match="malformed"):
        derive_priors(p)


def test_derive_priors_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        derive_priors(tmp_path / "no_such_file.res")


def test_derive_priors_falls_back_to_integer_indices(tmp_path):
    p = tmp_path / "no_ids.res"
    p.write_text("a\tb\ttheta\tIBS0\n"
                 "0\t1\t0.40\t0.000\n"
                 "0\t2\t0.01\t0.200\n")
    priors = derive_priors(p)
    assert priors[("0", "1")] is EdgeClass.DUPLICATE_OR_CLONE
    assert priors[("0", "2")] is EdgeClass.UNRELATED


def test_cli_prints_class_counts(make_genome_wide_res, capsys):
    res = make_genome_wide_res([
        ("S1", "S2", 0.40, 0.000),
        ("S1", "S3", 0.01, 0.200),
        ("S2", "S3", 0.01, 0.200),
    ])
    from adaptive.scheduler.prior import _main
    rc = _main([str(res)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "duplicate_or_clone" in out
    assert "unrelated" in out
