# drVizer Product Requirements Document

## Product summary

drVizer is an API-first Python library for building transcript-structure figures from GTF gene models, BED annotation intervals, and BAM coverage tracks. The primary output is a reusable matplotlib figure workflow for direct RNA sequencing and transcriptomics analysis.

## Goals

- Let users build publication-ready transcript structure plots from scripts or notebooks.
- Support one or more GTF files as the transcript model source.
- Overlay BED annotations as interval distributions or numeric score tracks.
- Overlay BAM coverage as sum or mean coverage tracks.
- Support transcript-coordinate BED/BAM inputs by projecting them back to genomic coordinates through the loaded GTF model.
- Allow users to prepare inputs once and plot many genes without repeating parser setup.

## Non-goals

- drVizer is not maintained as a CLI workflow in this repository.
- drVizer does not own read alignment, transcript discovery, or differential analysis.
- drVizer does not manage external workflow orchestration or deployment.

## Users

- Bioinformatics researchers producing transcript-centric figures.
- Analysts comparing annotation or coverage tracks across samples.
- Pipeline authors who need scriptable plotting in Python workflows.

## Main workflow

1. Import `DrViz` from `drvizer`.
2. Load GTF annotations with `load_gtf(...)`.
3. Add optional BED tracks with `add_bed_track(...)`.
4. Add optional BAM tracks with `add_bam_track(...)`.
5. Call `build()` for a reusable parser or `plot(...)` for one-shot plotting.
6. Save or reuse the returned matplotlib figure.

## Functional requirements

### GTF model loading

- Accept a path or list of GTF paths.
- Parse exon and CDS features.
- Support gzip-compressed GTF files.
- Support gene lookup by gene ID, gene name, or transcript ID.
- Cache parsed transcript models on the parser instance.

### BED tracks

- Accept a path or list of BED paths.
- Support BED3 and wider BED formats.
- Use half-open interval overlap semantics for region filtering.
- Support `parser_type="distribution"` for interval blocks.
- Support `parser_type="score"` for numeric bar tracks.
- Support transcript-coordinate BED projection when `transcript_coord=True`.
- Support split transcript-coordinate tracks with `split_by_transcript`.

### BAM tracks

- Accept a path or list of BAM paths.
- Require `pysam` at runtime.
- Read or create BAM index files through `pysam.index` when needed.
- Compute per-base coverage over genomic or projected transcript-coordinate regions.
- Aggregate multiple BAM files with `aggregate_method="sum"` or `aggregate_method="mean"`.
- Use int64 accumulators for coverage sums.
- Cap parallel workers for multi-BAM genomic coverage.

### Rendering

- Render transcript structures as exon/CDS tracks with strand-aware intron arrows.
- Render distribution, score, and coverage subplots under the transcript model.
- Support automatic layout sizing from transcript count and track count.
- Support shared y-axis scaling through `y_axis_group` for numeric tracks.
- Return a `matplotlib.figure.Figure` and optionally save it.

## Constraints and validations

- `load_gtf(...)` must run before `build()` or `plot(...)`.
- `split_by_transcript` must be `None`, `"nc"`, or `"cn"`.
- `split_by_transcript` requires `transcript_coord=True`.
- Split transcript tracks do not support multi-gene plotting.
- `y_axis_group` is valid only for numeric tracks.
- BED/BAM color and alpha lists must match file counts.
- BAM aggregation mode must be `"sum"` or `"mean"`.
- Multi-gene plotting requires all genes to be on the same chromosome.

## Success metrics

- Users can generate expected figures through the documented `DrViz` workflow.
- Public examples run with normal Python imports and no CLI dependency.
- Formal pytest validation passes in the DRS environment.
- Parser behavior remains stable for GTF, BED, BAM, transcript-coordinate, and split-track workflows.
