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

## [1.0.0] - prior

Pre-Phase 2.2 release.
