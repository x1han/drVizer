# drVizer

`drVizer` is a Python library for building transcript-structure figures from GTF models, BED annotations, and BAM coverage tracks.

It is designed for direct RNA sequencing and transcriptomics workflows that need reusable, scriptable, publication-ready matplotlib figures.

## Features

- Build transcript structure figures from one or more GTF files.
- Overlay BED annotation tracks as interval blocks or numeric score bars.
- Add BAM-derived coverage tracks with `sum` or `mean` aggregation.
- Reuse parsed state across many plotting calls through a builder-style API.
- Project transcript-coordinate BED and BAM inputs back into genomic space.
- Split transcript-coordinate tracks into transcript-specific subtracks.
- Share automatic y-axis scaling across matched numeric tracks with `y_axis_group`.
- Export matplotlib figures for publication workflows.

## Installation

Editable install for local development:

```bash
git clone https://github.com/x1han/drVizer.git
cd drVizer
pip install -e .
```

Regular local install:

```bash
pip install .
```

BAM coverage tracks require `pysam`:

```bash
pip install pysam
```

## Quick start

The public workflow is centered on `DrViz`:

```python
from drvizer import DrViz

parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed", label="TE")
    .add_bam_track("reads.bam", label="Coverage")
    .build()
)

fig = parser.plot("TP53", show=False)
```

`build()` freezes the configured GTF and tracks into a reusable parser, so plotting many genes does not repeat setup work.

### Reusable parser with auto-cleanup (Phase 2.2)

`DrViz` is also a context manager; using it as such hands back the
prepared `ReusableParser` and tears down the cache + lazy
`ProcessPool` on exit:

```python
from drvizer import DrViz

with DrViz().load_gtf("genes.gtf").add_bed_track("repeats.bed", label="TE") as parser:
    for gene in ("TP53", "BRCA1", "MYC"):
        fig = parser.plot(gene, show=False)
# cache cleared, pool shut down with wait=True
```

Cache and pool configuration:

- `DrViz(cache_maxsize=N)` — LRU cache capacity (default 128).
  Caches are per-parser-instance; each new `build()` constructs a
  fresh cache.
- `DrViz(adaptive_threshold=N)` — minimum total BED record count
  that opens a `ProcessPool` during build (default 20_000). Below
  the threshold, build is purely sequential.

## Core API

### `load_gtf(...)`

Load one or more GTF files. GTF parsing is required before tracks can be built.

```python
.load_gtf("genes.gtf")
.load_gtf(["reference.gtf", "novel_transcripts.gtf"])
```

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `gtf_files` | `str` or `list[str]` | required | Path to one GTF file or ordered list of GTF files. Exon/CDS features are parsed for transcript rendering. Loading a new GTF resets previously added tracks because tracks depend on the active GTF model. |

drVizer uses exon/CDS features and supports gene ID, gene name, or transcript ID lookup during plotting.

### `add_bed_track(...)`

Add BED-backed annotation tracks.

```python
.add_bed_track("repeats.bed", label="TE", color="tomato")
```

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `bed_files` | `str` or `list[str]` | required | BED file path or multiple BED files grouped into one logical track. BED3 and wider BED records are supported. |
| `label` | `str` or `None` | auto `Track_N` | Track label shown on the left side of the figure. Duplicate labels are made unique in registration order. |
| `color` | `str` or `list[str]` | `"orange"` | Matplotlib color for the track. When multiple BED files are passed, a list gives per-file colors. |
| `alpha` | `float` or `list[float]` | `0.8` | Track transparency from `0` to `1`. When multiple BED files are passed, a list gives per-file alpha values. |
| `parser_type` | `"distribution"` or `"score"` | `"distribution"` | Rendering mode. `distribution` draws interval blocks; `score` draws numeric BED scores as bars. |
| `y_axis_range` | `float` or `None` | `None` | Fixed y-axis maximum for `score` tracks. Takes precedence over automatic scaling and `y_axis_group`. |
| `y_axis_group` | `str` or `None` | `None` | Shared automatic y-axis scaling group for numeric BED `score` tracks. Invalid for `distribution` tracks. |
| `transcript_coord` | `bool` | `False` | Treat BED `chrom` field as transcript ID and project transcript coordinates back to genomic coordinates through the loaded GTF. |
| `layer_order` | `None`, `"ascending"`, or `"descending"` | `"ascending"` | Controls drawing order for layered BED elements. |
| `split_by_transcript` | `None`, `"nc"`, or `"cn"` | `None` | Split transcript-coordinate BED data into transcript-specific subtracks. Requires `transcript_coord=True`. |

Score tracks can share automatic y-axis scaling:

```python
.add_bed_track("control_m6a.bed", label="Control m6A", parser_type="score", y_axis_group="m6A")
.add_bed_track("treated_m6a.bed", label="Treated m6A", parser_type="score", y_axis_group="m6A")
```

### `add_bam_track(...)`

Add BAM-backed coverage tracks.

```python
.add_bam_track("reads.bam", label="Coverage", color="steelblue")
```

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `bam_files` | `str` or `list[str]` | required | BAM file path or multiple BAM files grouped into one logical coverage track. BAM support requires `pysam`. |
| `label` | `str` | `"Coverage"` | Track label shown on the left side of the figure. Duplicate labels are made unique in registration order. |
| `color` | `str` or `list[str]` | `"steelblue"` | Matplotlib color for coverage. When multiple BAM files are rendered as per-file series, a list gives per-file colors. |
| `alpha` | `float` or `list[float]` | `0.6` | Coverage transparency from `0` to `1`. When multiple BAM files are rendered as per-file series, a list gives per-file alpha values. |
| `aggregate_method` | `"sum"` or `"mean"` | `"sum"` | Combines multiple BAM files by summed coverage or average coverage. |
| `y_axis_range` | `float` or `None` | `None` | Fixed coverage y-axis maximum. Takes precedence over automatic scaling and `y_axis_group`. |
| `y_axis_group` | `str` or `None` | `None` | Shared automatic y-axis scaling group for numeric coverage tracks. |
| `transcript_coord` | `bool` | `False` | Treat BAM reference names as transcript IDs and project coverage back to genomic coordinates through the loaded GTF. |
| `layer_order` | `None`, `"ascending"`, or `"descending"` | `"ascending"` | Controls drawing order for per-file coverage series where individual series are rendered. |
| `split_by_transcript` | `None`, `"nc"`, or `"cn"` | `None` | Split transcript-coordinate BAM coverage into transcript-specific subtracks. Requires `transcript_coord=True`. |

Multiple BAM files can be combined:

```python
.add_bam_track(
    ["sample_a.bam", "sample_b.bam"],
    label="Reads",
    aggregate_method="mean",
)
```

Supported aggregation modes:

- `aggregate_method="sum"`: sum coverage across BAM files.
- `aggregate_method="mean"`: average coverage across BAM files.

### `build()` and `plot(...)`

Use `build()` when plotting multiple genes from the same inputs:

```python
parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed", label="TE")
    .add_bam_track("reads.bam", label="Coverage")
    .build()
)

fig1 = parser.plot("TP53", show=False)
fig2 = parser.plot("MYC", show=False)
```

Use one-shot `plot(...)` for quick figures:

```python
fig = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed", label="TE")
    .plot("TP53", show=False)
)
```

`ReusableParser.plot(...)` accepts these parameters:

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `gene` | `str` or `list[str]` | required | Gene ID, gene name, transcript ID, or same-chromosome list of identifiers to plot. |
| `transcript_to_show` | `str`, `list[str]`, or `None` | `None` | Restrict output to one transcript or selected transcripts from the requested gene. |
| `output` | `str` or `None` | `None` | Optional output path. When set, figure is saved with matplotlib using tight bounding box and 300 DPI. |
| `figsize` | `tuple` or `None` | `None` | Explicit final figure size in inches. Overrides automatically computed size. |
| `figfact` | `tuple` or `None` | `None` | Multiplicative width/height factor applied to automatically computed figure size. Ignored when `figsize` is set. |
| `show` | `bool` | `True` | Display figure through matplotlib. If `False`, figure is closed after creation but still returned. |
| `close` | `bool` | `False` | Close figure after showing it. Applies when `show=True`. |
| `**kwargs` | any | | Forwarded to the visualizer, including transcript sorting and layout options. |

## Transcript-coordinate workflows

Set `transcript_coord=True` when BED or BAM records use transcript IDs instead of genomic chromosome names.

```python
.add_bed_track("mods.transcript.bed", transcript_coord=True)
.add_bam_track("reads.transcript.bam", transcript_coord=True)
```

drVizer projects transcript-coordinate intervals and coverage back into genomic plotting space through the loaded GTF model.

## Split transcript tracks

Transcript-coordinate BED and BAM tracks can be split by transcript with `split_by_transcript`:

```python
split_by_transcript="nc"  # transcript-major order
split_by_transcript="cn"  # track-major order
```

Supported modes:

- `None`: keep combined track behavior.
- `"nc"`: transcript-major ordering; each transcript groups its split tracks together.
- `"cn"`: track-major ordering; each track groups its transcript-specific subtracks together.

Split transcript tracks require `transcript_coord=True` and do not support multi-gene plotting.

## Numeric y-axis control

Numeric tracks support two y-axis controls:

- `y_axis_range`: manually fix the y-axis maximum.
- `y_axis_group`: share automatic y-axis scaling across numeric tracks with the same group name.

`y_axis_range` takes precedence over `y_axis_group`.

```python
parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bam_track("control.bam", label="Control", y_axis_group="reads")
    .add_bam_track("treated.bam", label="Treated", y_axis_group="reads")
    .add_bed_track("m6a.bed", label="m6A", parser_type="score", y_axis_group="mod_score")
    .build()
)
```

## Output behavior

drVizer renders matplotlib figures and returns the generated `Figure` object from `plot(...)`.

```python
fig = parser.plot("TP53", show=False, output="tp53.pdf")
fig.set_size_inches((10, 6))
fig.savefig("tp53.png", dpi=300, bbox_inches="tight")
```

## Documentation

Detailed docs live in `docs/`:

- [Project overview / PDR](docs/project-overview-pdr.md)
- [Codebase summary](docs/codebase-summary.md)
- [System architecture](docs/system-architecture.md)
- [API reference](docs/api-reference.md)
- [Testing guide](docs/testing-guide.md)
- [Code standards](docs/code-standards.md)
- [Changelog](docs/changelog.md)

## Testing

```bash
python -m pytest -q
```

## Public entry point

```python
from drvizer import DrViz
```
