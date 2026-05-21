"""
adaptive.scheduler — per-pair adaptive site-budgeting for per-chromosome
ngsRelate-fast runs.

Public API:

    from adaptive.scheduler.prior      import derive_priors
    from adaptive.scheduler.subset     import BeagleSubsetCache
    from adaptive.scheduler.runner     import run_pair_on_chrom
    from adaptive.scheduler.confidence import chrom_class_from_res_row, should_escalate
    from adaptive.scheduler.scheduler  import AdaptiveScheduler
    from adaptive.scheduler.output     import write_res, write_manifest, write_run_manifest

See adaptive/docs/IMPLEMENTATION_PLAN.md for the build order and what each
module is responsible for. See adaptive/docs/SPEC_v0.1_CLARIFICATION_NOTE.md
for the prior-source clarification (genome-wide .res, not ngsPedigree Stage 1).
"""

from .edge_class import EdgeClass

__all__ = ["EdgeClass"]
