# drVizer API Reference

## Import

```python
from drvizer import DrViz
```

`DrViz` is the supported public API. Lower-level parser classes are implementation details unless a maintainer is extending the package.

## `DrViz`

Chainable builder for transcript visualization workflows.

### `DrViz().load_gtf(gtf_files)`

Load one or more GTF files and reset any previously added tracks.

```python
viz = DrViz().load_gtf("genes.gtf")
```

Parameters:

| Name | Type | Required | Description |
| --- | --- | --- | --- |
| `gtf_files` | `str` or `list[str]` | yes | GTF file path or list of GTF file paths. |

Behavior:

- Instantiates `GTFParser`.
- Parses the GTF data immediately.
- Clears existing track specs and configs.
- Returns the same `DrViz` instance.

### `DrViz().add_bed_track(...)`

Register a BED annotation track.

```python
viz = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed", label="Repeats", color="tomato")
)
```

Parameters:

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `bed_files` | `str` or `list[str]` | required | BED file path or paths. |
| `label` | `str` | generated | Display label. Duplicate labels get numeric suffixes. |
| `color` | `str` or `list[str]` | `"orange"` | Track color or one color per BED file. |
| `alpha` | `float` or `list[float]` | `0.8` | Track alpha or one alpha per BED file. |
| `parser_type` | `str` | `"distribution"` | `"distribution"` for interval blocks or `"score"` for score bars. |
| `y_axis_range` | `float` or `None` | `None` | Fixed y-axis maximum for numeric tracks. |
| `y_axis_group` | `str` or `None` | `None` | Shared automatic y-axis group for numeric tracks. |
| `transcript_coord` | `bool` | `False` | Treat BED reference names as transcript IDs. |
| `layer_order` | `str` or `None` | `"ascending"` | Draw order for layered data; valid values are `None`, `"ascending"`, `"descending"`. |
| `split_by_transcript` | `str` or `None` | `None` | Optional keyword-only value from `None`, `"nc"`, `"cn"`. |

Validation:

- `split_by_transcript` requires `transcript_coord=True`.
- `split_by_transcript` must be `None`, `"nc"`, or `"cn"`.
- `y_axis_group` requires `parser_type="score"`.
- color and alpha lists must match the BED file count.

### `DrViz().add_bam_track(...)`

Register a BAM coverage track.

```python
viz = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bam_track("reads.bam", label="Coverage", aggregate_method="sum")
)
```

Parameters:

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `bam_files` | `str` or `list[str]` | required | BAM file path or paths. |
| `label` | `str` | `"Coverage"` | Display label. Duplicate labels get numeric suffixes. |
| `color` | `str` or `list[str]` | `"steelblue"` | Track color or one color per BAM file. |
| `alpha` | `float` or `list[float]` | `0.6` | Track alpha or one alpha per BAM file. |
| `aggregate_method` | `str` | `"sum"` | `"sum"` or `"mean"`. |
| `y_axis_range` | `float` or `None` | `None` | Fixed y-axis maximum. |
| `y_axis_group` | `str` or `None` | `None` | Shared automatic y-axis group. |
| `transcript_coord` | `bool` | `False` | Treat BAM references as transcript IDs. |
| `layer_order` | `str` or `None` | `"ascending"` | Draw order for layered data; valid values are `None`, `"ascending"`, `"descending"`. |
| `split_by_transcript` | `str` or `None` | `None` | Optional keyword-only value from `None`, `"nc"`, `"cn"`. |

Validation:

- Requires `pysam`; otherwise raises `ImportError`.
- `aggregate_method` must be `"sum"` or `"mean"`.
- `split_by_transcript` requires `transcript_coord=True`.
- color and alpha lists must match BAM file count.

### `DrViz().build()`

Prepare all registered tracks and return a reusable parser.

```python
parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed")
    .add_bam_track("reads.bam")
    .build()
)
```

Returns: `ReusableParser`.

Behavior:

- Requires a loaded GTF parser.
- Prepares deferred track specs with `prepare_tracks_parallel(...)`.
- Preserves track registration order.
- Propagates underlying preparation errors when possible.

### `DrViz().get_transcript_data(gene, transcript_to_show=None)`

Return normalized plotting data without rendering.

```python
data = parser_builder.get_transcript_data("TP53")
```

Parameters:

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `gene` | `str` or `list[str]` | required | Gene ID, gene name, transcript ID, or same-chromosome list of identifiers. |
| `transcript_to_show` | `str` or `list[str]` or `None` | `None` | Optional transcript filter. |

Returns: dictionary with transcript model data and prepared tracks.

### `DrViz().plot(...)`

One-shot plotting helper.

```python
fig = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed")
    .plot("TP53", show=False)
)
```

Parameters:

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `gene` | `str` or `list[str]` | required | Gene ID, gene name, transcript ID, or same-chromosome list of identifiers. Multi-gene plotting is not supported when split transcript tracks are enabled. |
| `output` | `str` or `None` | `None` | Optional output image path. |
| `figsize` | `tuple` or `None` | `None` | Explicit figure size. |
| `figfact` | `tuple` or `None` | `None` | Multiplicative size factor. |
| `show` | `bool` | `True` | Whether to display figure. |
| `close` | `bool` | `False` | Close figure after showing. |
| `**kwargs` | any | | Forwarded to `visualize_gene_transcripts(...)`. |

Returns: `matplotlib.figure.Figure`.

## `ReusableParser`

Returned by `DrViz.build()`.

### `ReusableParser.plot(...)`

Plot one gene or multiple same-chromosome genes using already prepared parser state.

```python
parser = DrViz().load_gtf("genes.gtf").add_bed_track("repeats.bed").build()
fig1 = parser.plot("TP53", show=False)
fig2 = parser.plot("MYC", show=False)
```

Parameters match `DrViz.plot(...)`, plus `transcript_to_show` for filtering.

## Split transcript modes

`split_by_transcript` controls ordering for transcript-coordinate tracks.

| Value | Meaning |
| --- | --- |
| `None` | Combined track behavior. |
| `"nc"` | Transcript-major order: each transcript groups its split tracks together. |
| `"cn"` | Track-major order: each track groups its transcript-specific subtracks together. |

Split transcript tracks require `transcript_coord=True` and do not support multi-gene plotting.

## Numeric y-axis control

Numeric track types are BAM coverage tracks and BED score tracks.

- `y_axis_range`: fixed maximum for one track.
- `y_axis_group`: shared automatic maximum across tracks in the same group.

`y_axis_range` takes precedence over `y_axis_group`.

## Lower-level parser classes

Maintainers may work directly with these classes when extending internals:

- `GTFParser` in `src/drvizer/gtf_parser.py`
- `BEDParser` in `src/drvizer/bed_parser.py`
- `BAMParser` in `src/drvizer/bam_parser.py`

The stable public workflow should continue to use `DrViz`.
