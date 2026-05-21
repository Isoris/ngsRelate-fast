"""adaptive — per-pair adaptive site-budgeting for per-chromosome ngsRelate-fast runs.

Marker module that promotes `adaptive/` from a PEP 420 namespace package
to a regular package. This is what lets pytest's rootpath walk up from
`adaptive/tests/test_X.py` all the way to the repo root, putting the
repo root on sys.path so absolute imports like
`from adaptive.scheduler.prior import derive_priors` work whether
pytest is invoked as `pytest` (cwd NOT on sys.path) or as
`python -m pytest` (cwd on sys.path).

See adaptive/README.md for usage; this file is intentionally empty
otherwise.
"""
