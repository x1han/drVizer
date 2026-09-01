# Changelog

All notable changes to drVizer are documented here. Versions follow
[Semantic Versioning](https://semver.org/). The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added
- **LRU cache in `PreparedDataSource`** (Phase 2.2): an
  `OrderedDict`-backed LRU cache keyed on
  `(track_id, target_id, chrom, start, end)` so repeat
  `get_transcript_data` / `plot` calls skip BED parse and BAM
  coverage work. `cache_maxsize` defaults to 128; configure via
  `DrViz(cache_maxsize=N)`.
- **`ReusableParser` lifecycle lock** (Phase 2.2): `__enter__`,
  `__exit__`, `close()`, and `clear_cache()` methods. `close()` is
  idempotent and runs the mandated order: `clear_cache` first, then
  `pool.shutdown(wait=True)`, then flip `_is_closed`. The data
  source rejects post-close inserts with `RuntimeError`.
- **`DrViz` context manager**: `with DrViz() as parser:` returns the
  prepared `ReusableParser` (not the builder) so the cache and pool
  are scoped to the block.
- **Adaptive `ProcessPool` threshold** (Phase 2.2): the build-stage
  process pool is only opened when the genomic-BED record estimate
  exceeds `adaptive_threshold` (default 20_000). Under the threshold
  the path is pure sequential and no pool is allocated. Configure via
  `DrViz(adaptive_threshold=N)`.
- **BAM dtype compression**: cached coverage payloads are stored as
  `np.uint32` for `aggregate_method='sum'` and `np.float32` for
  `aggregate_method='mean'`.
- `benchmarks/bench_cache_vs_no_cache.py`: stdlib-only benchmark
  measuring the cache hit/miss speedup across `small/medium/large`
  BED workloads.

### Changed
- `DrViz.build()` is now cached: successive calls with the same
  builder state return the same `ReusableParser` instance so the
  LRU cache survives. Any mutating method (`load_gtf`,
  `add_bed_track`, `add_bam_track`) calls `_mark_dirty()` and forces
  a recompile on the next `build()`.
- `DrViz.plot()` and `DrViz.get_transcript_data()` route through
  `self.build()` instead of constructing a fresh `ReusableParser`
  per call.
- `drvizer._track_build.prepare_tracks_parallel` accepts
  `adaptive_threshold` and `data_source` keyword arguments; the
  ProcessPool executor is owned by the data source so its lifecycle
  is locked to the `ReusableParser.close()` path.
- The build-stage outer `ThreadPoolExecutor` wrapper around
  process + serial paths has been removed; the two paths now run
  sequentially within `prepare_tracks_parallel` (process first, then
  serial), and the pool is gated by the threshold.

### Phase 3 — Cython/Python parity & hardening

#### Added
- **`ParallelCoverageError` re-exported from `drvizer` package**:
  callers can now `from drvizer import ParallelCoverageError` for
  typed error handling around parallel coverage failures.
- **`tests/fixtures/crlf.bed`**: Windows-style CRLF-terminated BED
  fixture (PARSER-011) exercising comment + blank + multi-chrom
  records.
- **`GTFParser._chunk_parse_degradation` counter**: incremented each
  time the Cython chunk parser falls back to per-row Python parsing
  (P1-3).
- **Five new test files** (`tests/test_phase_3_*.py`) covering CRLF
  BED parity, GTF chunk degradation, BAM projection kernel
  consistency, `cpu_count()` None guard, and `ParallelCoverageError`
  exception chain. Two additional tests (`test_close_failure_emits_*`,
  `test_get_transcript_data_none_coord_safe`) extend the Phase 2.2
  LRU-cache suite.

#### Changed
- **P1-8** (`_parallel.py`): the `Pool.imap_unordered` worker now
  catches any `Exception`, re-raises as `ParallelCoverageError` chained
  from the worker exception (`raise ... from exc`), and includes the
  failing BAM path in the message. The original traceback is preserved
  on `__cause__` so corrupt BAM / missing `.bai` failures stay
  diagnosable.
- **P1-3** (`gtf_parser.py`, `_cython_gtf.pyx`): pattern (a) of the
  Cython chunk try/except decision is implemented. The `_process_chunk`
  wrapper calls `_parse_gtf_chunk_impl` inside `try/except`; on
  chunk-level failure it falls back to per-row Python parsing with
  per-row `try/except (ValueError, TypeError)`, increments
  `self._chunk_parse_degradation`, and emits a single
  `RuntimeWarning` at `parse_gtf()` return time. The Cython kernel
  remains a pure best-effort transform with no C-level try/except.
- **P1-6** (`bam_parser.py`): extracted the
  `_project_aligned_blocks_to_transcript(transcript_id, blocks,
  gtf_parser, target_chrom, region_start, region_len, coverage)`
  pure-numeric kernel; both the single-transcript and the batch
  transcript-coverage paths now route through this single helper so
  spliced-read / int64 / `np.add.at` fixes live in one place.
- **P1-5** (`_parallel.py`, `_track_build.py`): `cpu_count()` is now
  guarded with `try/except NotImplementedError` plus a
  `isinstance(int)` fallback (and `os.cpu_count() or 1` in
  `_track_build.py`). The original audit's `min(N, None, 32)` →
  `TypeError` failure mode is impossible.
- **PARSER-011** (`_cython_bed.pyx`, `bed_parser.py`): BED open mode
  unified. Cython path uses `raw_line.rstrip(b'\r\n')`; Python
  fallback opens in `'rb'` and decodes UTF-8 with `errors='replace'`
  before `rstrip('\r\n')` and tab-splitting. CRLF and LF BEDs now
  produce byte-identical records across both paths.
- **Cleanup B** (`api.py`): the three `try/except Exception: pass`
  blocks in `ReusableParser.close()` and the `__del__` finalizer now
  surface failures as `warnings.warn(..., ResourceWarning,
  stacklevel=2)` instead of silently swallowing them.
- **`PreparedDataSource._fetch_coverage_in_region` /
  `_fetch_grouped_anno_in_region`** (`api.py`): the cache key now
  preserves `None` for missing coordinates instead of raising
  `TypeError` from `int(None)`; repeated calls with `start=None,
  end=None` are now cache-hits.

## [1.0.0] - prior

Pre-Phase 2.2 release.
