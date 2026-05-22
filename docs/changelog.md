# drVizer Changelog

This changelog summarizes documented repository behavior. It is not a replacement for git history.

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

- BAM support requires `pysam` at runtime.
- Cython acceleration is optional at import time; Python fallbacks remain important.
- CLI workflows are not part of the supported current public interface.
