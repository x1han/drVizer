# Changelog

All notable changes to drVizer are documented here. Versions follow
[Semantic Versioning](https://semver.org/). The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [0.1.1] - 2026-09-02

### Added
- **CI matrix workflow** (`.github/workflows/ci.yml`): ubuntu-latest
  on push to `main` and on `pull_request`, with a Python 3.8 - 3.12
  matrix. Steps include a dedicated Cython build-prerequisites install
  (mirroring `.github/workflows/benchmarks.yml`) so the editable
  install compiles its Cython extensions on every runner, then
  `pytest -q` for the test suite. macOS is intentionally excluded
  from this matrix because Cython compilation on macOS is unvalidated
  in the maintainer's Linux/DSR environment and risks a workflow that
  fails on first push; it can be added in a follow-up phase once the
  Linux matrix is proven stable across several push cycles.
- **Release workflow** (`.github/workflows/release.yml`): triggered by
  push tags matching `v*`. Three sequential jobs: (1) `build` runs
  `python -m build` and `twine check dist/*`; (2) `test-pypi` always
  uploads the built artifact to TestPyPI as a dry-run / canary
  validation step using the `TESTPYPI_API_TOKEN` secret; (3)
  `publish-pypi` uploads to production PyPI using the
  `PYPI_API_TOKEN` secret, gated on a `pypi` GitHub environment that
  requires manual approval via environment protection rules. **No
  auto-publish on first run.** The upload-artifact / download-artifact
  pair is wired so the `dist/` directory flows from `build` to both
  publish jobs.
- **Editable install gotcha section** in `CLAUDE.md`: documents that
  switching git branches (especially after manual ref manipulation
  such as `git update-ref refs/heads/main <sha>` or `git symbolic-ref
  HEAD refs/heads/main`) leaves the editable install pointing at the
  old `src/`. Mandated command:
  `/datf/hanxi/software/miniconda3/envs/DRS/bin/pip install
  --force-reinstall --no-deps -e .` after every merge.

### Changed
- **Tuple-return worker contract** (P1-8, Phase 7 PARTIALLY-FIXED ->
  FIXED; `src/drvizer/_parallel.py`): the
  `_compute_region_coverage_with_path` worker now returns
  `Tuple[str, np.ndarray]` -- `(bam_path, coverage_array)` -- so
  per-bam-path attribution is encoded in the data shape instead of in
  the closure binding that Phase 7 relied on. The master loop in
  `aggregate_region_coverages_parallel` unpacks the tuple:
  `for bam_path, coverage in pool.imap_unordered(...)`. This is a
  **refactor for clarity and test maintainability**, not a
  correctness fix: the Phase 7 closure binding was functionally
  correct (the `bam_path` was bound into the worker closure, the
  error message included the path, and master aggregation was
  order-independent). v0.1.0's "Fixed: imap_unordered attribution
  invariant" entry describes the original closure-binding contract;
  v0.1.1's "Changed: tuple-return worker contract" describes the same
  invariant re-expressed as a data shape. The worker still raises
  `ParallelCoverageError(... from exc)` inside the worker so the
  `__cause__` chain survives the `imap_unordered` boundary intact.
- **Test mocks updated** to match the new contract: 4 FakePool mocks
  across `tests/test_bam_parallel.py` and
  `tests/test_phase_3_cpu_guard.py` now return `(bam_path,
  coverage_array)` tuples derived from `args[i][0]`, and the Phase 7
  `tests/test_phase_7_imap_unordered.py` order-independence test
  also returns `(bam_path, coverage_array)` tuples in reversed
  arrival order. Net: the worker contract is exercised end-to-end in
  every mock that previously returned a raw ndarray.

### Notes
- v0.1.0 was tagged locally only (commit `dd33199`) and was **never
  pushed to PyPI**. v0.1.1 is therefore the first PyPI release of
  drVizer without breaking any existing install.
- No throughput numbers or speedup claims appear in this release.
  The GENCODE + ENCODE + pyGenomeTracks comparison is research-grade
  work that belongs to a separate Stage 3 project, not to v0.1.1.

## [0.1.0] - 2026-09-01

### Added
- **`pyproject.toml` PEP 621 `[project]` section** (PKG-002): `name=drvizer`,
  `version=0.1.0`, `requires-python>=3.8`, MIT license text, README,
  Python 3.8 - 3.12 classifiers (incl. `Development Status :: 4 - Beta`),
  `dependencies` from `requirements.txt`, optional `[bam]=[pysam]` and
  `[test]=[pytest]` extras, and `project_urls`. Authoritative single
  source of install / metadata truth.
- **`docs/api.md`**: public API contract document describing `DrViz`
  constructor signature (`adaptive_threshold`, `cache_maxsize`), `plot()`,
  `get_transcript_data()`, `ReusableParser` context manager,
  `ParallelCoverageError` inheritance (`RuntimeError` -> subclass),
  and the `__version__` contract.
- **Three new test files** (`tests/test_phase_7_*.py`):
  - `test_phase_7_imap_unordered.py` verifies that the existing
    `_compute_region_coverage_with_path` worker binds `bam_path` via
    closure, the `ParallelCoverageError` raised on worker failure
    carries the failing `bam_path` in its message, and the master
    consumer's `imap_unordered` loop pairs `bam_path` to its coverage
    correctly even when results arrive out of order.
  - `test_phase_7_degradation_guard.py` asserts no `RuntimeWarning` is
    emitted when zero chunks degraded and exactly one `RuntimeWarning`
    is emitted with the degradation count in the message when chunks
    degraded.
  - `test_phase_7_version_and_metadata.py` asserts `drvizer.__version__`
    matches `importlib.metadata.version("drvizer")` and both equal
    `"0.1.0"`, and the `pyproject.toml` classifiers include Python
    3.8 through 3.12 plus `Development Status :: 4 - Beta`.

### Changed
- **`setup.py` skeleton** (P1-16): `setup.py` is now a thin C-extension
  compile-only wrapper. `extra_compile_args=["-std=c99"]` is preserved on
  every Cython `Extension` block (panel template did not mention this;
  dropping it would silently violate C99 compliance). All metadata
  (name, version, classifiers, install_requires, extras_require,
  project_urls, long_description) was moved to `pyproject.toml [project]`
  so the two sources no longer drift.
- **`__version__` bumped `1.0.0` -> `0.1.0`** (PKG-002): first public
  release tag.
- **GTF chunk-degradation warning gate** (P1-3): the end-of-parse
  `RuntimeWarning` in `GTFParser.parse_gtf` is now Python-gated on
  `if self._chunk_parse_degradation:` so clean parses emit no warning.

### Fixed
- **`imap_unordered` attribution invariant** (P1-8 verification, P0-8
  follow-up): the worker `_compute_region_coverage_with_path` already
  binds `bam_path = args[0]` via closure; the consumer's
  `imap_unordered` loop relies on the per-`bam_path` `ParallelCoverageError`
  message (which now embeds the path) instead of result order. The
  invariant is locked in `tests/test_phase_7_imap_unordered.py`.

### Hygiene
- **Dead `BEDParser` import removed** from `src/drvizer/api.py:31`
  (OCR-006): `BEDParser` is not referenced anywhere else in `api.py`;
  tests import directly from `drvizer.bed_parser`. Net: zero callers
  depend on the `api.py` re-export.
- **Dead `defaultdict` import removed** from `src/drvizer/utils.py:10`
  (OCR-005): `defaultdict` is not referenced anywhere else in
  `utils.py`.
- **Stale `docs/superpowers/` `.gitignore` rule removed** (DOC-001):
  `docs/` does not contain a `superpowers/` directory; the rule was
  leftover from a prior layout. The corresponding `workspace/`
  superpowers material now lives under `/datf/hanxi/software/drVizer/workspace/`.

## [Unreleased] — Phase 2.2 (LRU cache + lifecycle)

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
