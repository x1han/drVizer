# drVizer

**drVizer** is an API-first Python library for visualizing transcript structures from GTF annotations with optional BED annotations and BAM-derived coverage tracks.

It is designed for reusable, scriptable figure generation in transcriptomics and related genomic analysis workflows.

## Key features

- Build transcript structure figures from one or multiple GTF files
- Overlay BED annotation tracks and BAM coverage tracks
- Support transcript-coordinate BED and BAM inputs
- Split transcript-coordinate tracks into transcript-specific subtracks
- Reuse parsed state across many plotting calls through a builder-style API
- Export publication-ready figures through matplotlib

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

The public workflow is centered on `DrViz`:

1. `load_gtf(...)`
2. `add_bed_track(...)`
3. `add_bam_track(...)`
4. `build()`
5. `plot(...)`

### Minimal example

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

### BED and BAM example

```python
from drvizer import DrViz

parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed", label="TE", color="tomato")
    .add_bam_track("reads.bam", label="Coverage", color="steelblue")
    .build()
)

fig = parser.plot("TP53")
```

### One-shot plotting

```python
from drvizer import DrViz

fig = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed", label="TE")
    .plot("TP53")
)
```

## Transcript-coordinate workflows

Set `transcript_coord=True` when BED or BAM records are aligned to transcript IDs rather than genomic chromosome coordinates.

```python
transcript_coord=True
```

This enables drVizer to project transcript-coordinate inputs back into genomic plotting space.

## Split transcript tracks

Transcript-coordinate tracks can be split by transcript with:

```python
split_by_transcript="nc"  # or "cn"
```

Supported modes:

- `None`: keep the combined legacy track behavior
- `"nc"`: transcript-major ordering
- `"cn"`: track-major ordering

When split mode is enabled, drVizer adds transcript labels on the right side of split subtracks and keeps left-side track labels aligned with the rendered subplot order.

### Example: split transcript-coordinate BED and BAM tracks

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

## Key API concepts

### Track configuration

Common options for `add_bed_track(...)` and `add_bam_track(...)` include:

- `label`: display name for the track
- `color`: track color
- `alpha`: transparency
- `transcript_coord`: treat the input as transcript-coordinate data
- `split_by_transcript`: split transcript-coordinate tracks by transcript (`None`, `"nc"`, `"cn"`)
- `y_axis_range`: fix the y-axis maximum instead of auto-scaling

BAM-specific options:

- `aggregate_method`: combine multiple BAM files using `"sum"` or `"mean"`

BED-specific options:

- `parser_type`: render BED data as `"distribution"` or `"score"`

### Reusable plotting

`build()` freezes the configured tracks into a reusable parser object. This is the preferred workflow when plotting multiple genes from the same GTF and track configuration.

```python
parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed", label="TE")
    .build()
)

fig1 = parser.plot("TP53")
fig2 = parser.plot("MYC")
```

## Output behavior

drVizer renders matplotlib figures and returns the resulting `Figure` object from `plot(...)`.

Figures can be:
- displayed directly
- resized after creation
- saved through matplotlib using PNG, PDF, SVG, or other supported formats

## Public entry point

```python
from drvizer import DrViz
```
