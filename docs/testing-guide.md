# drVizer Testing Guide

## Test command

Use the DRS conda environment for repository validation:

```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest -q
```

Run from the repository root:

```text
/datf/hanxi/software/drVizer/repo
```

## Test layout

Tests live in `tests/`.

Important coverage areas:

| File | Coverage area |
| --- | --- |
| `tests/test_bed_parser_acceleration.py` | BED parser behavior, Python/Cython parity, region filtering, BED3 support. |
| `tests/test_bam_parallel.py` | BAM coverage aggregation, worker caps, int64 accumulators, fallback behavior. |
| `tests/test_build_parallel_preparation.py` | Deferred track preparation and ordering. |
| `tests/test_gtf_parser_acceleration.py` | GTF parser fast path behavior and fallbacks. |
| `tests/test_projection_acceleration.py` | Transcript-to-genomic projection behavior. |
| `tests/test_split_by_transcript_api.py` | Public API validation and split transcript data preparation. |
| `tests/test_split_by_transcript_visualizer.py` | Split transcript rendering behavior. |

## Test data style

Most tests construct small temporary GTF, BED, and BAM-like fixtures at runtime. Prefer minimal fixtures that show one behavior rather than copying large real data into the repository.

Good patterns:

- Use `tmp_path` for temporary files.
- Use monkeypatching for external behavior such as multiprocessing pools or pysam calls.
- Assert public API behavior when changing `DrViz` workflow behavior.
- Assert lower-level parser behavior when changing parser semantics.
- Keep Python and Cython path expectations aligned when both paths exist.

## TDD expectations

For bug fixes and behavior changes:

1. Add or update a test that fails for the current behavior.
2. Run the targeted test and confirm the failure is the expected one.
3. Make the smallest code change that passes the test.
4. Run the targeted test again.
5. Run the full pytest suite.

## Targeted commands

Examples:

```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest -q tests/test_bed_parser_acceleration.py
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest -q tests/test_bam_parallel.py
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest -q tests/test_split_by_transcript_api.py
```

## Cython extension validation

The repository has Cython sources for acceleration:

- `src/drvizer/_cython_projection.pyx`
- `src/drvizer/_cython_gtf.pyx`
- `src/drvizer/_cython_bed.pyx`

If a `.pyx` file changes, rebuild extensions from the repository root before relying on Cython-path tests:

```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python setup.py build_ext --inplace
```

Tests should either exercise the available compiled fast path or skip fast-path checks when the extension is unavailable.

## BAM testing notes

BAM tests should avoid depending on large external files. For parser-level behavior, use monkeypatched `pysam.AlignmentFile`, fake reads, and fake pool objects where possible.

When real BAM data is needed for local validation, keep it outside the formal repository and document the path in workflow material rather than adding it to `tests/`.

## Regression checklist

Before reporting code changes complete:

- Targeted tests for changed behavior pass.
- Full pytest suite passes.
- Cython source and compiled behavior are in sync when `.pyx` files changed.
- New tests cover both public API behavior and low-level parser semantics when both are affected.
- Workflow-only validation files remain outside the source repository.
