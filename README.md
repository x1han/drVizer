# drVizer

**drVizer** is a Python library for visualizing transcript structures from GTF annotations with optional BED and BAM-derived tracks.

It is designed for API-first use in analysis scripts and reusable plotting workflows.

## Installation

```bash
git clone https://github.com/x1han/drVizer.git
cd drVizer
pip install -e .
```

For a regular local install:

```bash
pip install .
```

If you want to use BAM coverage tracks, you may also need:

```bash
pip install pysam
```

## Core workflow

The main API is centered on `DrViz`:

- `load_gtf(...)`
- `add_bed_track(...)`
- `add_bam_track(...)`
- `build()`
- `plot(...)`

## Basic example

```python
from drvizer import DrViz

parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed", label="TE")
    .build()
)

fig = parser.plot("TP53")
```

## Add multiple BED tracks

```python
from drvizer import DrViz

parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("te.bed", label="TE", color="tomato")
    .add_bed_track("mod_sites.bed", label="m6A", color="royalblue")
    .build()
)

fig = parser.plot("TP53")
```

## Add BAM coverage

```python
from drvizer import DrViz

parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bam_track("reads.bam", label="Coverage", color="steelblue")
    .build()
)

fig = parser.plot("TP53")
```

## One-shot plotting

```python
from drvizer import DrViz

fig = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed", label="TE")
    .plot("TP53")
)
```

## Transcript-coordinate tracks

You can treat BED or BAM inputs as transcript-coordinate data by setting:

```python
transcript_coord=True
```

This is useful when records are aligned to transcript IDs rather than genomic chromosome coordinates.

## Split transcript tracks

Transcript-coordinate tracks can be split by transcript with:

```python
split_by_transcript="nc"  # or "cn"
```

Supported modes:

- `None`: keep the combined legacy track behavior
- `"nc"`: transcript-major ordering
- `"cn"`: track-major ordering

When split mode is enabled, drVizer adds transcript labels on the right side of split subtracks.

## Example: split transcript-coordinate BAM and BED tracks

```python
from drvizer import DrViz

parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed", label="TE")
    .add_bed_track(
        ["mod_control.bed", "mod_treated.bed"],
        label="m6A",
        color=["gray", "royalblue"],
        alpha=[0.15, 0.3],
        transcript_coord=True,
        split_by_transcript="nc",
        parser_type="score",
    )
    .add_bam_track(
        ["sample_a.bam", "sample_b.bam"],
        label="Reads",
        color=["#4a98c9", "#f14432"],
        alpha=[0.25, 0.25],
        transcript_coord=True,
        split_by_transcript="nc",
    )
    .build()
)

fig = parser.plot("TP53", show=False)
```

## Common parameters

Common options for `add_bed_track(...)` and `add_bam_track(...)` include:

- `label`: display name for the track
- `color`: track color
- `alpha`: transparency
- `transcript_coord`: treat the input as transcript-coordinate data
- `split_by_transcript`: split transcript-coordinate tracks by transcript (`None`, `"nc"`, `"cn"`)
- `y_axis_range`: fix the y-axis maximum instead of auto-scaling

BAM-specific options:

- `aggregate_method`: combine multiple BAM files using `"sum"` or `"mean"`

## Capabilities

- Parse one or multiple GTF files
- Access genes by gene ID, gene name, or transcript ID
- Add BED annotation tracks with `add_bed_track(...)`
- Add BAM coverage tracks with `add_bam_track(...)`
- Export figures as PNG, PDF, or SVG
- Reuse parsed data efficiently across many genes with `build()`

## Public API

```python
from drvizer import DrViz
```
