# drVizer Public API Contract (v0.1.0)

This document is the authoritative public API contract for the `drvizer` package
as of v0.1.0. It documents the symbols exported from `drvizer.__init__`,
their signatures, the public methods that downstream callers may rely on,
and the error / version contract.

## Version

```python
import drvizer
assert drvizer.__version__ == "0.1.0"
```

`drvizer.__version__` matches `importlib.metadata.version("drvizer")` so
both introspection paths return the same string.

## Re-exports

```python
from drvizer import DrViz, ParallelCoverageError
```

- `DrViz`: chainable builder for transcript visualization.
- `ParallelCoverageError`: `RuntimeError` subclass raised when a BAM
  coverage worker raises inside the parallel aggregator.

Anything under `drvizer.*` that is NOT in this list is implementation
detail and may change between minor versions without notice.

## `DrViz`

```python
DrViz(*, adaptive_threshold: int = 20_000, cache_maxsize: int = 128)
```

| Parameter | Type | Default | Valid range | Description |
| --- | --- | --- | --- | --- |
| `adaptive_threshold` | `int` (keyword-only) | `20_000` | `[1_000, 1_000_000]` | Minimum total BED record count that opens a `ProcessPool` during build. Below the threshold the build path is pure sequential and no pool is allocated. A `warnings.warn` is emitted if the value is outside the recommended range. |
| `cache_maxsize` | `int` (keyword-only) | `128` | `[16, 1024]` | LRU cache capacity on the prepared data source. A `warnings.warn` is emitted if the value is outside the recommended range. |

Both parameters are keyword-only — positional arguments are not accepted.

### `DrViz.load_gtf(gtf_files)`

Load one or more GTF files and reset previously added tracks. Returns `self`
for chaining. `gtf_files` is a `str` or `list[str]`.

### `DrViz.add_bed_track(...)`

Add one BED-backed track. See `docs/api-reference.md` for the full parameter
table; the public contract is the keyword arguments listed there
(`bed_files`, `label`, `color`, `alpha`, `parser_type`, `y_axis_range`,
`y_axis_group`, `transcript_coord`, `layer_order`, `split_by_transcript`).

### `DrViz.add_bam_track(...)`

Add one BAM-backed coverage track. Requires `pysam` to be installed (the
`[bam]` extra).

### `DrViz.build() -> ReusableParser`

Freeze the current builder state into a reusable plotting object. Successive
calls with the same builder state return the same `ReusableParser` instance
so the LRU cache survives across `build()` cycles.

### `DrViz.get_transcript_data(gene, transcript_to_show=None) -> dict`

Return normalized plotting data for one gene (or same-chromosome gene list)
without rendering. `gene` is `str` or `list[str]`; `transcript_to_show` is
`str`, `list[str]`, or `None`.

### `DrViz.plot(gene, output=None, figsize=None, figfact=None, show=True, close=False, **kwargs) -> matplotlib.figure.Figure`

Plot one gene from the builder without explicitly calling `build()`. Returns
the generated `matplotlib.figure.Figure`.

### `DrViz.__enter__() -> ReusableParser`

`with DrViz() as parser:` hands back the prepared `ReusableParser` (not
the builder) so the LRU cache and lazy `ProcessPool` are scoped to the
block.

### `DrViz.__exit__(exc_type, exc_val, exc_tb) -> bool`

Returns `False` (exceptions propagate). The cached `ReusableParser` is
closed on context exit.

## `ReusableParser`

Returned by `DrViz.build()`.

### `ReusableParser.__enter__() -> ReusableParser`

`with parser:` is self-binding; use `parser.plot(...)` inside.

### `ReusableParser.__exit__(...) -> bool`

Runs `close()` on context exit.

### `ReusableParser.close() -> None`

Idempotent. Runs in the mandated order: `clear_cache()` first, then
`pool.shutdown(wait=True)` so worker subprocesses drain, then flip
`_is_closed` so subsequent calls short-circuit. Failures during any of
the three steps surface as `ResourceWarning` instead of being swallowed.

### `ReusableParser.clear_cache() -> None`

Public passthrough that drops the prepared data source cache.

### `ReusableParser.plot(gene, transcript_to_show=None, output=None, figsize=None, figfact=None, show=True, close=False, **kwargs) -> matplotlib.figure.Figure`

Plot one gene (or same-chromosome gene list) using already prepared parser
state. `**kwargs` forwards to the visualizer.

## `ParallelCoverageError`

```python
class ParallelCoverageError(RuntimeError):
    ...
```

Inheritance chain: `RuntimeError` -> `ParallelCoverageError`.

`ParallelCoverageError` is raised when a worker inside
`drvizer._parallel.aggregate_region_coverages_parallel` raises. The
exception is constructed as
`ParallelCoverageError(f"Failed to compute coverage for {bam_path}: {exc}")`
and chained from the worker exception via `raise ... from exc`, so:

- `str(exc)` includes the failing `bam_path`.
- `exc.__cause__` is the original worker exception (e.g. `OSError` for
  a corrupt BAM, `KeyError` for a missing `.bai` reference, etc.).
- `exc.__cause__.__traceback__` preserves the original traceback.

`ParallelCoverageError` is also re-exported from the top-level `drvizer`
package so callers can write:

```python
from drvizer import ParallelCoverageError

try:
    coverage = drvizer._parallel.aggregate_region_coverages_parallel(...)
except ParallelCoverageError as exc:
    print(f"worker failed for {exc.__cause__!r}: {exc}")
```

## Stability guarantee

Within the v0.1.x series, the public symbols listed above will NOT change
their signatures or remove existing keyword arguments without a major
version bump. New keyword arguments may be added; new methods may be added
to `DrViz` and `ReusableParser`. Internal modules
(`drvizer._parallel`, `drvizer._track_build`, `drvizer._cython_*`) are NOT
part of the public contract.