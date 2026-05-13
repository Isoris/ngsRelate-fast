# Contributing to ngsRelate-fast

This fork has a narrow remit: make ngsRelate faster on whole-genome
cohorts **without changing the Jacquard EM**. Contributions outside
that scope should go upstream to [ANGSD/NgsRelate](https://github.com/ANGSD/NgsRelate).

## Rules

### Hard constraints

1. **The Jacquard EM is sacred.** Do not modify `analyse_jaq`, `em`,
   `emAccel`, `emStep`, or any function consumed by them. These produce
   the J1-J9 columns that downstream tools (ngsPedigree) depend on.
2. **The `.res` output schema is sacred.** 23 columns, identical names
   and order to upstream. No additions, no removals, no reorderings.
3. **Validation must pass before merge.** Both stages of
   `validate/run_validation.sh` on at least one test cohort.
4. **Schema changes are major version bumps.** Adding a field to
   `*.input.v1.schema.json` requires a new file `*.input.v2.schema.json`
   and a corresponding update to `KNOWN_INPUT_SCHEMAS` in `contract_io.py`.
   Document in `CHANGELOG.md`.

### Soft constraints

- Keep the patch small and focused. If `patches/01_ngsRelate_fast.patch`
  grows beyond ~300 lines, split it into multiple numbered patches.
- Match upstream code style (K&R-ish C++, 2-space indent).
- New flags need: help text in `print_info`, a `case` in `getopt`, an
  entry in the README flag table, and an entry in the input contract
  `params` block.

## Workflow

```bash
git checkout -b feature/short-description
# make changes
make build
make validate BEAGLE=... FREQS=... N=... SAMPLES=...
git commit -m "feat: <what changed> + <why>"
git push -u origin feature/short-description
# open PR
```

## Commit message convention

Prefix with one of: `feat:`, `fix:`, `docs:`, `perf:`, `refactor:`,
`test:`, `build:`, `chore:`. Reference issues by `#NN`. Keep the
subject line under 72 chars.

Examples:
- `feat: add --min-sites-per-chrom floor for tiny scaffolds`
- `fix: handle empty BEAGLE without segfault`
- `docs: clarify -D 0 semantics in README`
- `perf: cache freq powers outside per-pair loop`
