# drVizer Changelog

This changelog summarizes documented repository behavior. It is not a replacement for git history.

## Phase 1 behavior fix (current)

- **Projection coordinates are now 0-based half-open.** The Cython and
  Python projection paths emit the BED-standard 0-based half-open
  genomic interval. Widths are invariant under this change, so
  matplotlib renders are unchanged; any user silently piping projection
  output into BEDTools / pybedtools / GenomicRanges will see a +1
  correction on the right edge.
- **BED and projection coordinates are int64.** Cython BED parser and
  Cython projection now use `int64_t`, so records with `start >= 2**31`
  (very large contigs) no longer silently wrap or overflow.
- **Empty placeholder track contract.** `_expand_split_tracks` pads
  missing (track, transcript) combinations with an `empty: True`
  placeholder. The visualizer short-circuits on `track.get('empty')` to
  render a label-only axis (`No signal / unexpressed`) instead of
  crashing on `np.max([])` and `bar([], [])`.
- **Python and Cython BED / projection paths are parity-locked** via the
  new `test_int64_coordinate_round_trip` and
  `test_64bit_and_zero_based_parity` tests.

## Current documented state

- drVizer is maintained as an API-first Python library.
- The public workflow is centered on `from drvizer import DrViz`.
- Formal validation uses pytest in the DRS conda environment.
- Public docs live in `docs/` and README.

## Recent confirmed behavior

- BED3 records are accepted with default optional BED fields.
- BED region filtering uses half-open overlap semantics.
- Python and Cython BED parsing paths are expected to stay behaviorally aligned.
- BAM coverage accumulation uses int64 arrays to avoid int32 overflow.
- Multi-BAM genomic coverage caps worker count by BAM count, CPU count, and 32.
- Split transcript tracks reject multi-gene plotting with a clear `ValueError`.

## Known product decision pending

Transcript-coordinate BAM `mean` aggregation currently averages over observed series when split by transcript. A product decision is still needed for missing transcript references: average over configured BAMs with missing transcript as zero, or keep observed-series mean and document that behavior.

## Compatibility notes

- BAM support requires `pysam` at runtime via `pip install drvizer[bam]`.
- Cython acceleration is optional at import time. If Cython is missing at
  build, the build emits a `warnings.warn(ImportWarning, ...)` instead
  of silently skipping the extensions. Python fallbacks remain important.
- Python versions 3.7 and 3.8 are no longer classified; minimum is 3.9.
- CLI workflows are not part of the supported current public interface.
