"""
contract_io.py — Load and verify ngsRelate-fast input/output contracts.

For downstream tools (ngsPedigree, manuscript figure scripts, etc.) to
consume the .res output safely. Standard pattern:

    from contract_io import load_run, RunNotComplete

    try:
        run = load_run("/path/to/relatedness.res")
    except RunNotComplete:
        sys.exit("Upstream ngsRelate-fast run not finished or failed")

    # run.res_path, run.samples_path, run.params, etc.
    df = pd.read_csv(run.res_path, sep="\\t")

Schema versions are checked. Unknown/future schemas raise SchemaMismatch.
"""

from __future__ import annotations
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


KNOWN_INPUT_SCHEMAS  = {"ngsrelate_fast.input.v1"}
KNOWN_OUTPUT_SCHEMAS = {"ngsrelate_fast.output.v1"}


class RunNotComplete(Exception):
    """Output contract is missing — run hasn't finished or failed."""


class SchemaMismatch(Exception):
    """JSON schema version not recognized."""


class ContractCorrupt(Exception):
    """Contract file exists but can't be parsed or fails self-checks."""


@dataclass
class Run:
    """Resolved, verified view of one ngsRelate-fast run."""
    run_id: str
    res_path: Path
    samples_path: Optional[Path]
    params: dict
    n_pairs: int
    elapsed_seconds: float
    downsampling: dict
    status: str
    warnings: list
    input_contract: dict
    output_contract: dict


def _sidecar_paths(res_path):
    p = Path(res_path).resolve()
    return (
        p.parent / f"{p.name}.input.json",
        p.parent / f"{p.name}.output.json",
    )


def load_run(res_path, verify_hashes: bool = False) -> Run:
    """
    Load and verify the contracts beside a .res file.

    Args:
        res_path: path to the .res file
        verify_hashes: if True, recompute sha256 of the .res and check it
                       matches the output contract (catches post-run edits)

    Raises:
        RunNotComplete: output contract missing
        ContractCorrupt: files exist but inconsistent or malformed
        SchemaMismatch: unknown schema version
    """
    res_path = Path(res_path).resolve()
    in_json, out_json = _sidecar_paths(res_path)

    if not in_json.exists():
        raise ContractCorrupt(f"input contract missing: {in_json}")
    if not out_json.exists():
        raise RunNotComplete(f"output contract missing: {out_json}")

    try:
        inp  = json.loads(in_json.read_text())
        outp = json.loads(out_json.read_text())
    except json.JSONDecodeError as e:
        raise ContractCorrupt(f"malformed JSON: {e}")

    if inp.get("schema") not in KNOWN_INPUT_SCHEMAS:
        raise SchemaMismatch(f"unknown input schema: {inp.get('schema')}")
    if outp.get("schema") not in KNOWN_OUTPUT_SCHEMAS:
        raise SchemaMismatch(f"unknown output schema: {outp.get('schema')}")

    if inp.get("run_id") != outp.get("run_id"):
        raise ContractCorrupt(
            f"run_id mismatch: input={inp.get('run_id')} output={outp.get('run_id')}")

    if not res_path.exists():
        raise ContractCorrupt(f".res file missing: {res_path}")

    if verify_hashes:
        expected = outp.get("outputs", {}).get("res", {}).get("sha256")
        if expected:
            h = hashlib.sha256()
            with open(res_path, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
            actual = h.hexdigest()
            if actual != expected:
                raise ContractCorrupt(
                    f".res hash mismatch: expected {expected[:12]}..., got {actual[:12]}...")

    samples_ref = inp.get("inputs", {}).get("samples")
    samples_path = Path(samples_ref["path"]) if samples_ref else None

    return Run(
        run_id=inp["run_id"],
        res_path=res_path,
        samples_path=samples_path,
        params=inp["params"],
        n_pairs=outp["stats"]["n_pairs"],
        elapsed_seconds=outp["stats"]["elapsed_seconds"],
        downsampling=outp.get("downsampling", {}),
        status=outp["status"],
        warnings=outp.get("warnings", []),
        input_contract=inp,
        output_contract=outp,
    )


if __name__ == "__main__":
    # CLI: quick contract dump for debugging
    import sys
    if len(sys.argv) != 2:
        print("usage: contract_io.py <path/to/relatedness.res>", file=sys.stderr)
        sys.exit(1)
    try:
        run = load_run(sys.argv[1], verify_hashes=True)
    except (RunNotComplete, ContractCorrupt, SchemaMismatch) as e:
        print(f"FAIL: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(2)
    print(f"run_id:         {run.run_id}")
    print(f"status:         {run.status}")
    print(f".res:           {run.res_path}")
    print(f"n_pairs:        {run.n_pairs}")
    print(f"elapsed:        {run.elapsed_seconds:.1f}s")
    print(f"D_sites_per_gb: {run.params['D_sites_per_gb']}")
    if run.downsampling.get("enabled"):
        print(f"downsampling:   {run.downsampling['total_input_sites']} -> "
              f"{run.downsampling['total_kept_sites']} sites "
              f"across {len(run.downsampling['per_chromosome'])} chromosomes")
    if run.warnings:
        print(f"warnings:       {run.warnings}")
