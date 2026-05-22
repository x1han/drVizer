# drVizer Codebase Summary

## Repository role

This repository is the formal Python package source for drVizer. It contains package code, packaging metadata, tests, and public documentation. Local workflow material belongs outside this repository in the surrounding workspace.

## Package layout

```text
src/drvizer/
├── __init__.py
├── api.py
├── gtf_parser.py
├── bed_parser.py
├── bam_parser.py
├── _parallel.py
├── _track_build.py
├── visualizer.py
├── utils.py
├── _cython_projection.pyx
├── _cython_gtf.pyx
└── _cython_bed.pyx
```

## Main modules

### `src/drvizer/__init__.py`

Defines package version and exports `DrViz` as the public API surface.

### `src/drvizer/api.py`

High-level builder and reusable plotting API.

Key classes and functions:

- `DrViz`: chainable public builder.
- `PreparedDataSource`: internal combiner for GTF data and prepared tracks.
- `ReusableParser`: reusable object returned by `DrViz.build()`.
- `_validate_split_by_transcript(...)`: validates split transcript options.
- `_build_right_label_groups(...)`: computes split-track transcript label groups.

Typical flow:

1. `DrViz.load_gtf(...)` builds and parses a `GTFParser`.
2. `DrViz.add_bed_track(...)` and `DrViz.add_bam_track(...)` register deferred track specs.
3. `DrViz.build()` prepares tracks through `_track_build.prepare_tracks_parallel(...)`.
4. `ReusableParser.plot(...)` requests normalized gene data and calls `visualize_gene_transcripts(...)`.

### `src/drvizer/gtf_parser.py`

GTF model parser and transcript-coordinate projection logic.

Key behavior:

- Accepts one or more GTF files.
- Reads gzip or plain text files.
- Processes lines in chunks.
- Keeps exon and CDS features.
- Maintains gene, transcript, and name lookup indexes.
- Supports lookup by gene ID, gene name, or transcript ID.
- Converts transcript-coordinate intervals into genomic segments.
- Uses Cython helpers when available and Python fallbacks otherwise.

Important output shape:

- top-level gene metadata: `gene_id`, `seqname`, `strand`, `identifier_type`, `original_identifier`
- transcript entries: `transcript_id`, `strand`, `exons`, `cds`

### `src/drvizer/bed_parser.py`

BED parser for distribution and score annotation tracks.

Key behavior:

- Accepts one or more BED files.
- Supports BED3 and wider BED records.
- Supplies defaults for missing optional BED fields.
- Groups records by chromosome in `anno_data`.
- Uses half-open overlap filtering.
- Supports transcript-coordinate projection via `GTFParser`.
- Supports split transcript-coordinate BED payloads.
- Uses `_cython_bed.parse_bed_records` when available for genomic BED parsing.

### `src/drvizer/bam_parser.py`

BAM coverage parser for genomic and transcript-coordinate coverage tracks.

Key behavior:

- Requires `pysam`.
- Accepts one or more BAM paths.
- Ensures BAM indexes exist by calling `pysam.index(...)` when missing.
- Computes coverage arrays through `compute_region_coverage(...)`.
- Uses parallel aggregation for multi-BAM genomic coverage.
- Supports `sum` and `mean` aggregation.
- Uses int64 coverage accumulators.
- Supports transcript-coordinate BAM projection and split transcript tracks.
- Bins long coverage regions down to target bins for plotting.

### `src/drvizer/_parallel.py`

Low-level helpers for genomic BAM coverage.

Key behavior:

- `compute_region_coverage(...)` reads one BAM region and returns per-base coverage.
- `aggregate_region_coverages_parallel(...)` aggregates multiple BAMs.
- Worker count is capped by BAM count, CPU count, and 32.
- Parallel failures are normalized as `ParallelCoverageError` so `BAMParser` can fall back to serial aggregation.

### `src/drvizer/_track_build.py`

Build-stage track preparation.

Key behavior:

- Converts deferred track specs from `DrViz` into prepared parser objects.
- Prepares process-safe genomic BED tracks in a process pool.
- Prepares transcript-coordinate and BAM tracks serially because they depend on live parser state or pysam behavior.
- Preserves registration order when returning prepared tracks.

### `src/drvizer/visualizer.py`

Matplotlib rendering layer.

Key behavior:

- Builds transcript and track subplot layout.
- Renders exons, CDS blocks, introns, distribution tracks, score tracks, and coverage tracks.
- Applies shared y-axis scaling for numeric tracks.
- Renders right-side labels for split transcript tracks.
- Returns the generated matplotlib figure.

## Tests

Tests live in `tests/` and cover parser acceleration, GTF parsing, projection, BAM parallel behavior, split transcript APIs, split transcript visualization, and build-stage track preparation.

Run full validation with:

```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest -q
```

## Packaging

- `setup.py` defines package metadata, dependencies, and optional Cython extensions.
- `pyproject.toml` declares setuptools, wheel, and Cython build requirements.
- `requirements.txt` lists pandas, matplotlib, and numpy.
- BAM support requires `pysam`, which is imported by `bam_parser.py` but not listed in `requirements.txt`.

## Public documentation map

- [Project overview / PDR](project-overview-pdr.md)
- [System architecture](system-architecture.md)
- [API reference](api-reference.md)
- [Testing guide](testing-guide.md)
- [Code standards](code-standards.md)
- [Changelog](changelog.md)
