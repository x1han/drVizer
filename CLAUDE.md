# CLAUDE.md

This file provides repository guidance to Claude Code when working in `/datf/hanxi/software/drVizer/source`.

## Repository boundary

This directory is the **formal GitHub-linked software repository** for drVizer.

- Push, pull, branch, tag, release, and worktree actions belong here.
- This is the correct layer for `source/` worktrees and subagent-driven software development.
- Local-only plans, specs, and migration notes belong in `/datf/hanxi/software/drVizer/workflow`, not in this repository.

## Dual-git operating model

The surrounding workspace uses two git layers:

1. **Outer workspace git** at `/datf/hanxi/software/drVizer/.git`
   - Local-only history for workspace structure and `workflow/`
   - Never the default target for software pushes

2. **Inner repository git** at `/datf/hanxi/software/drVizer/source/.git`
   - Authoritative history for drVizer software
   - Owns the GitHub remote and all software worktrees

## Dual-CLAUDE operating model

This inner `CLAUDE.md` governs the **software repository layer**.

- It owns software-repository authority, implementation guidance, and repo-layer worktree rules.
- The outer `/datf/hanxi/software/drVizer/CLAUDE.md` governs the **workspace layer**.
- If the layered operating model changes, update both files together so their boundaries stay aligned.

## Surgical operating rules

Use these rules before merge, push, reset, branch deletion, or worktree removal:

1. Verify the current layer with:
   - `pwd`
   - `git rev-parse --show-toplevel`
   - `git rev-parse --git-common-dir`
   - `git branch --show-current`
2. Do not act until those values match the intended layer.
3. Do not use error-swallowing command chains to decide critical git state.
4. Push, pull, branch, tag, release, and software worktree actions belong to this repository layer.
5. If both git layers change, commit this repository first, then the outer workspace layer.

### Commit rules

- If you changed code, tracked tests, packaging, or public docs in this repository, commit here.
- If you only changed local plans/specs/workspace structure, commit the outer workspace git instead.
- If both layers changed, commit this repository first, then commit the outer workspace.

## Project overview

`drVizer` is a Python library and CLI for visualizing gene transcript structures from GTF annotations, optionally overlaying BED annotations and BAM-derived coverage tracks.

The repository exposes both:
- the installed CLI entry point in `src/drvizer/cli.py`
- a package API under `src/drvizer/`, including the chainable notebook-friendly API in `src/drvizer/api.py`
- the repository-local wrapper `drvizer_cli.py`

The codebase is small and focused: most important behavior lives in a few modules under `src/drvizer/`.

## Common commands

### Environment setup
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pip install -r requirements.txt
```

If you need BAM support, `src/drvizer/bam_parser.py` imports `pysam`, but `requirements.txt` currently only lists pandas/matplotlib/numpy. BAM functionality may require installing `pysam` separately.

### Run the CLI
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python drvizer_cli.py --gtf genes.gtf --gene ENSG00000136997 --output gene_plot.png
/datf/hanxi/software/miniconda3/envs/DRS/bin/python drvizer_cli.py --gtf genes.gtf --bed repeats.bed --gene TP53 --output merged_plot.png
```

### Tests
Use the DRS environment for all test and build validation:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest
```

Legacy local scripts may exist outside this repository root, but formal repository validation should move toward tracked pytest coverage.

## Architecture

### Execution paths
- `drvizer_cli.py` is the repository-local CLI wrapper.
- `src/drvizer/cli.py` is the installed console entry point.
- `src/drvizer/__init__.py` re-exports the main public API.
- `src/drvizer/api.py` provides the higher-level chainable API (`DrViz`, `PreparedDataSource`, `ReusableParser`) for notebook and interactive usage.

### Core data flow
1. Parse annotation data:
   - `src/drvizer/gtf_parser.py` loads transcript models from one or more GTF files.
   - `src/drvizer/bed_parser.py` loads annotation tracks from one or more BED files.
   - `src/drvizer/bam_parser.py` computes coverage arrays from BAM files.
2. Normalize into plotting-ready dictionaries.
3. Compose transcript data and track data in the API layer.
4. Render matplotlib figures in `src/drvizer/visualizer.py`.

### GTF parsing model
`src/drvizer/gtf_parser.py` is the backbone of the repository.

Key behaviors:
- `GTFParser` accepts one or many GTF files.
- Parsing is chunked line-by-line into temporary pandas DataFrames for memory efficiency.
- Only `exon` and `CDS` features are retained.
- Parsed data is cached in memory on the parser instance.
- The parser maintains:
  - `gene_transcripts`
  - `gene_info`
  - `gene_name_to_id`
  - `transcript_to_gene`
- `get_transcript_data()` supports lookup by gene ID, gene name, or transcript ID.
- Multi-file parsing suppresses transcripts already seen in earlier files.

When changing parsing behavior, preserve the output shape expected by the visualizer: transcript dicts with `transcript_id`, `exons`, `cds`, plus top-level gene metadata like `gene_id`, `seqname`, `strand`, `identifier_type`, and `original_identifier`.

### BED parsing model
`src/drvizer/bed_parser.py` handles auxiliary annotation tracks.

Important behaviors:
- Accepts one or many BED files.
- Stores parsed records grouped by chromosome in `anno_data`.
- Supports genomic BED coordinates and transcript-coordinate BED inputs converted back to genomic space using `GTFParser`.
- Supports parser types `distribution` and `score`.
- Region extraction is done with `get_anno_in_region()` or `get_grouped_anno_in_region()`.

### BAM coverage model
`src/drvizer/bam_parser.py` adds coverage-style tracks.

Important behaviors:
- Uses `pysam` to read one or multiple BAM files.
- Aggregates per-base coverage over a requested genomic window.
- Supports `sum` or `mean` aggregation semantics.
- Dynamically bins long regions to keep plotting manageable.
- Supports transcript-coordinate BAM projection through `GTFParser`.

### Plotting model
`src/drvizer/visualizer.py` contains rendering logic.

Key concepts:
- Transcript structures render as exon/CDS tracks.
- Introns are drawn as strand-aware arrows.
- Prepared track data from the API layer is rendered as distribution, score, or coverage subplots.

## Repository-specific notes

- The package import path is `drvizer`.
- Packaging metadata lives in `setup.py`.
- The installed console entry point is `drvizer=drvizer.cli:main`.
- Worktrees for software changes should be created from this `source/` repository.
- For non-trivial implementation work, prefer `source/` worktrees plus subagents to keep software changes isolated and fast to iterate.
- Software worktrees are disposable by default and should be kept through implementation, review, and validation while still active, then removed after acceptance unless the user explicitly asks to retain them.
- Local-only plans, specs, and migration docs belong in `/datf/hanxi/software/drVizer/workflow`, not in this repository.
