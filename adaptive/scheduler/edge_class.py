"""
edge_class.py — the seven relationship classes used by the adaptive scheduler.

These names are the project-wide canonical naming per SPEC §9 and
IMPLEMENTATION_PLAN.md §4. **Never use bare abbreviations** (PO, FS, 2nd,
etc.) anywhere in this package — the string values here are the only form
that should appear in code, logs, manifests, or schemas.
"""

from __future__ import annotations
from enum import Enum


class EdgeClass(str, Enum):
    DUPLICATE_OR_CLONE      = "duplicate_or_clone"
    PARENT_OFFSPRING        = "parent_offspring"
    FULL_SIBLING            = "full_sibling"
    AMBIGUOUS_FIRST_DEGREE  = "ambiguous_first_degree"
    SECOND_DEGREE           = "second_degree"
    THIRD_DEGREE            = "third_degree"
    UNRELATED               = "unrelated"

    def __str__(self) -> str:
        return self.value
