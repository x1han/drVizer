# drVizer 🧬

**drVizer** is a Python library for building transcript-structure figures from GTF models, BED annotations, and BAM coverage tracks.

It is designed for direct RNA sequencing and transcriptomics workflows where you need reusable, scriptable, publication-ready matplotlib figures.

## ✨ Key features

- 🧱 Build transcript structure figures from one or multiple GTF files
- 🛏️ Overlay BED annotation tracks as interval blocks or numeric score bars
- 📈 Add BAM-derived coverage tracks with sum or mean aggregation
- 🔁 Reuse parsed state across many plotting calls through a builder-style API
- 🧭 Project transcript-coordinate BED and BAM inputs back into genomic space
- 🧩 Split transcript-coordinate tracks into transcript-specific subtracks
- 📏 Share automatic y-axis scaling across matched numeric tracks with `y_axis_group`
- 🎨 Export publication-ready matplotlib figures

## 🚀 Installation

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

## 🧪 Quick start

The public workflow is centered on `DrViz`:

1. `load_gtf(...)`
2. `add_bed_track(...)`
3. `add_bam_track(...)`
4. `build()`
5. `plot(...)`

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

## 🧭 Core API

### `load_gtf(...)`

Load one or more GTF files. GTF parsing is the required first step.

```python
.load_gtf("genes.gtf")
.load_gtf(["reference.gtf", "novel_transcripts.gtf"])
```

drVizer uses exon/CDS features and supports gene ID, gene name, or transcript ID lookup during plotting.

### `add_bed_track(...)`

Add BED-backed annotation tracks.

```python
.add_bed_track("repeats.bed", label="TE", color="tomato")
```

BED tracks can render in two modes:

- `parser_type="distribution"` — draw intervals as annotation blocks
- `parser_type="score"` — draw interval scores as numeric bars

Score tracks can participate in shared numeric y-axis scaling:

```python
.add_bed_track("control_m6a.bed", label="Control m6A", parser_type="score", y_axis_group="m6A")
.add_bed_track("treated_m6a.bed", label="Treated m6A", parser_type="score", y_axis_group="m6A")
```

`y_axis_group` is rejected for BED distribution tracks because distribution tracks do not use numeric y-axis scaling.

### `add_bam_track(...)`

Add BAM-backed coverage tracks.

```python
.add_bam_track("reads.bam", label="Coverage", color="steelblue")
```

Multiple BAM files can be combined:

```python
.add_bam_track(
    ["sample_a.bam", "sample_b.bam"],
    label="Reads",
    aggregate_method="mean",
)
```

Supported aggregation modes:

- `aggregate_method="sum"` — sum coverage across BAM files
- `aggregate_method="mean"` — average coverage across BAM files

Coverage tracks can share automatic y-axis scaling:

```python
.add_bam_track(control_bams, label="Control Reads", aggregate_method="mean", y_axis_group="reads")
.add_bam_track(treated_bams, label="Treated Reads", aggregate_method="mean", y_axis_group="reads")
```

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

## 📏 Numeric y-axis control

Numeric tracks support two y-axis controls:

- `y_axis_range`: manually fix the y-axis maximum
- `y_axis_group`: share automatic y-axis scaling across numeric tracks with the same group name

`y_axis_range` takes precedence. If one track has `y_axis_range=10`, that track uses `0..10` even if it also has `y_axis_group`.

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

In split transcript tracks, grouping is applied within each transcript. Matched subtracks for the same transcript share a y-axis maximum, but different transcripts keep independent scales.

## 🧬 Transcript-coordinate workflows

Set `transcript_coord=True` when BED or BAM records are aligned to transcript IDs rather than genomic chromosome coordinates.

```python
.add_bed_track("mods.transcript.bed", transcript_coord=True)
.add_bam_track("reads.transcript.bam", transcript_coord=True)
```

drVizer projects transcript-coordinate intervals and coverage back into genomic plotting space using the loaded GTF model.

## 🧩 Split transcript tracks

Transcript-coordinate BED and BAM tracks can be split by transcript with `split_by_transcript`.

```python
split_by_transcript="nc"  # transcript-major order
split_by_transcript="cn"  # track-major order
```

Supported modes:

- `None`: keep combined track behavior
- `"nc"`: transcript-major ordering; each transcript groups its split tracks together
- `"cn"`: track-major ordering; each track groups its transcript-specific subtracks together

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
        y_axis_group="m6A",
    )
    .add_bam_track(
        ["sample_a.bam", "sample_b.bam"],
        label="Reads",
        color=["#4a98c9", "#f14432"],
        alpha=[0.25, 0.25],
        aggregate_method="mean",
        transcript_coord=True,
        split_by_transcript="nc",
        y_axis_group="reads",
    )
    .build()
)

fig = parser.plot("TP53", show=False)
```

## 🛠️ Helper / API reference

This section lists the supported builder helpers and their parameters.

### `DrViz().load_gtf(...)`

```python
.load_gtf(gtf_files)
```

| Parameter | Type | Default | Function |
| --- | --- | --- | --- |
| `gtf_files` | `str` or `list[str]` | required | Load one or more GTF files. Resets previously added tracks because tracks depend on the loaded GTF model. |

Notes:

- GTF parsing keeps exon/CDS features for transcript rendering.
- Lookup during plotting supports gene ID, gene name, and transcript ID.

### `DrViz().add_bed_track(...)`

```python
.add_bed_track(
    bed_files,
    label=None,
    color="orange",
    alpha=0.8,
    parser_type="distribution",
    y_axis_range=None,
    y_axis_group=None,
    transcript_coord=False,
    split_by_transcript=None,
)
```

| Parameter | Type | Default | Function |
| --- | --- | --- | --- |
| `bed_files` | `str` or `list[str]` | required | BED file path or multiple BED files for one logical track. |
| `label` | `str` | auto `Track_N` | Track label shown on the left side of the figure. Duplicate labels are made unique in registration order. |
| `color` | `str` or `list[str]` | `"orange"` | Track color. When multiple BED files are passed, a list gives per-file colors. |
| `alpha` | `float` or `list[float]` | `0.8` | Track transparency. When multiple BED files are passed, a list gives per-file alpha values. |
| `parser_type` | `"distribution"` or `"score"` | `"distribution"` | Rendering mode. `distribution` draws interval blocks; `score` draws numeric score bars. |
| `y_axis_range` | `float` or `None` | `None` | Fixed numeric y-axis maximum for score tracks. Overrides automatic scaling and `y_axis_group`. |
| `y_axis_group` | `str` or `None` | `None` | Share automatic y-axis scaling with other score tracks in the same group. Valid only when `parser_type="score"`. |
| `transcript_coord` | `bool` | `False` | Treat BED `chrom` field as transcript ID and project transcript coordinates back to genomic coordinates. |
| `split_by_transcript` | `None`, `"nc"`, or `"cn"` | `None` | Split transcript-coordinate BED data into transcript-specific subtracks. Requires `transcript_coord=True`. |

BED examples:

```python
# Interval blocks
.add_bed_track("repeats.bed", label="TE", parser_type="distribution")

# Numeric score bars
.add_bed_track("m6a.bed", label="m6A", parser_type="score")

# Multi-file score track with per-file styles
.add_bed_track(
    ["control_m6a.bed", "treated_m6a.bed"],
    label="m6A",
    parser_type="score",
    color=["gray", "royalblue"],
    alpha=[0.2, 0.35],
)

# Transcript-coordinate split score track
.add_bed_track(
    "m6a.transcript.bed",
    label="m6A",
    parser_type="score",
    transcript_coord=True,
    split_by_transcript="nc",
    y_axis_group="m6A",
)
```

### `DrViz().add_bam_track(...)`

```python
.add_bam_track(
    bam_files,
    label="Coverage",
    color="steelblue",
    alpha=0.6,
    aggregate_method="sum",
    y_axis_range=None,
    y_axis_group=None,
    transcript_coord=False,
    split_by_transcript=None,
)
```

| Parameter | Type | Default | Function |
| --- | --- | --- | --- |
| `bam_files` | `str` or `list[str]` | required | BAM file path or multiple BAM files for one logical coverage track. |
| `label` | `str` | `"Coverage"` | Track label shown on the left side of the figure. Duplicate labels are made unique in registration order. |
| `color` | `str` or `list[str]` | `"steelblue"` | Coverage color. When multiple BAM files are passed, a list gives per-file colors where per-file series are rendered. |
| `alpha` | `float` or `list[float]` | `0.6` | Coverage transparency. When multiple BAM files are passed, a list gives per-file alpha values where per-file series are rendered. |
| `aggregate_method` | `"sum"` or `"mean"` | `"sum"` | Combine multiple BAM files by summed coverage or average coverage. |
| `y_axis_range` | `float` or `None` | `None` | Fixed coverage y-axis maximum. Overrides automatic scaling and `y_axis_group`. |
| `y_axis_group` | `str` or `None` | `None` | Share automatic y-axis scaling with other coverage tracks in the same group. |
| `transcript_coord` | `bool` | `False` | Treat BAM references as transcript IDs and project coverage back to genomic coordinates. |
| `split_by_transcript` | `None`, `"nc"`, or `"cn"` | `None` | Split transcript-coordinate BAM coverage into transcript-specific subtracks. Requires `transcript_coord=True`. |

BAM examples:

```python
# One BAM coverage track
.add_bam_track("reads.bam", label="Reads")

# Mean coverage across replicates
.add_bam_track(["rep1.bam", "rep2.bam"], label="Mean Reads", aggregate_method="mean")

# Matched tracks with shared y-axis scaling
.add_bam_track(control_bams, label="Control", aggregate_method="mean", y_axis_group="reads")
.add_bam_track(treated_bams, label="Treated", aggregate_method="mean", y_axis_group="reads")

# Transcript-coordinate split coverage
.add_bam_track(
    "reads.transcript.bam",
    label="Reads",
    transcript_coord=True,
    split_by_transcript="nc",
    y_axis_group="reads",
)
```

### `DrViz().build()`

```python
parser = DrViz().load_gtf("genes.gtf").add_bed_track("repeats.bed").build()
```

| Function | Details |
| --- | --- |
| Freezes builder state | Converts deferred track specs into prepared parser objects. |
| Reusable plotting | Returns `ReusableParser`, which can plot many genes from the same configured inputs. |
| Preparation behavior | Genomic BED tracks may prepare in a process pool; transcript-coordinate tracks and BAM tracks keep GTF-aware serial preparation. |

### `DrViz().plot(...)`

```python
fig = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed")
    .plot("TP53", show=False)
)
```

This is a one-shot helper: it calls `build()` internally, then plots one gene.

### `ReusableParser.plot(...)`

```python
fig = parser.plot(
    gene,
    transcript_to_show=None,
    output=None,
    figsize=None,
    figfact=None,
    show=True,
    close=False,
)
```

| Parameter | Type | Default | Function |
| --- | --- | --- | --- |
| `gene` | `str` or `list[str]` | required | Gene ID/name/transcript ID to plot. Lists are supported when genes are on the same chromosome. |
| `transcript_to_show` | `str`, `list[str]`, or `None` | `None` | Restrict plot to one transcript or selected transcripts. |
| `output` | `str` or `None` | `None` | Save the figure to this path using matplotlib. |
| `figsize` | `tuple` or `None` | `None` | Set final figure size directly. |
| `figfact` | `tuple` or `None` | `None` | Scale generated figure width/height by this factor. |
| `show` | `bool` | `True` | Display the figure through matplotlib. If `False`, the figure is closed after creation but still returned. |
| `close` | `bool` | `False` | Close the figure after showing it. Applies when `show=True`. |


## 🖼️ Output behavior

drVizer renders matplotlib figures and returns the generated `Figure` object from `plot(...)`.

Figures can be:

- displayed directly with `show=True`
- returned without display using `show=False`
- resized with `figsize` or `figfact`
- saved through `output="figure.pdf"`
- saved manually through matplotlib as PNG, PDF, SVG, or other supported formats

```python
fig = parser.plot("TP53", show=False, output="tp53.pdf")
fig.set_size_inches((10, 6))
fig.savefig("tp53.png", dpi=300, bbox_inches="tight")
```

## ✅ Testing

Use the project DRS environment for validation:

```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest -q
```

## 📦 Public entry point

```python
from drvizer import DrViz
```
