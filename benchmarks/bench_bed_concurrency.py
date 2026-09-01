"""Phase 2.2.0 BED-concurrency data-collection benchmark.

This is a data-collection-only harness used to decide whether the
build-stage ``_prepare_process_batch`` in ``drvizer._track_build`` should
move toward a module-level ``ProcessPoolExecutor`` singleton (audit
P1-7), stick with sequential / thread-pool execution, or migrate to an
altogether different concurrency primitive.

It runs three workloads (small-region batch, whole-track ingestion,
mixed notebook-style submission) through three primitives
(sequential, ``ThreadPoolExecutor``, ``ProcessPoolExecutor``) at four
worker counts (1, 2, 4, 8).  Each (workload, primitive, workers)
configuration runs one uncounted warmup batch then ``--repetitions``
timed batches.  Wall time is reported as min / median / max, peak RSS
via ``resource.getrusage``, and records-per-second as total_records /
median.

The harness intentionally imports the private
``drvizer._track_build._prepare_genomic_bed_track`` so the
``ProcessPoolExecutor`` path pays the same pickle / unpickle cost as
the audited ``_prepare_process_batch`` call site.  Zero files under
``src/drvizer/`` are edited by this benchmark.

Memory caveat: ``ru_maxrss`` is a per-process high-water mark that
never decreases, so for sequential / thread_pool the value reported
for a late configuration is contaminated by the peak of every earlier
configuration in the same interpreter.  This limitation is recorded in
the JSON ``meta.peak_rss_note`` block.

Cython availability is reported in ``meta.cython_bed_available`` but
not swept by this script -- toggling Cython at runtime requires moving
``src/drvizer/_cython_*.so`` aside and re-running under a clean
interpreter (see Phase 2.2.0 plan runbook).
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import random
import statistics
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
from resource import RUSAGE_SELF, RUSAGE_CHILDREN, getrusage

# Make ``drvizer`` importable when running from a checkout that has not
# been pip-installed in editable mode yet.  This matches the conftest.py
# bootstrap used by the formal tests.
_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC_ROOT = _REPO_ROOT / "src"
if _SRC_ROOT.is_dir() and str(_SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(_SRC_ROOT))

import drvizer._track_build as _track_build  # noqa: E402  (sys.path edit above)


# ---------------------------------------------------------------------------
# Fixture generation
# ---------------------------------------------------------------------------

# BED6 columns: chrom, start, end, name, score, strand.  Tab-separated.
_BED_COLUMNS = 6


def _write_bed(path: Path, n_records: int, span_bp: int, rng: random.Random) -> None:
    """Write ``n_records`` BED6 rows into ``path`` packed into ``span_bp``.

    Rows are deterministic given the seed in ``rng``.  Each row occupies
    a 100 bp slot with no overlap so chromosome-position bookkeeping
    stays trivially correct regardless of record count.
    """
    stride = 100
    with path.open("wt") as handle:
        for index in range(n_records):
            start = index * stride
            end = start + stride
            name = f"feat_{index:06d}"
            score = rng.random()
            strand = "+" if index % 2 == 0 else "-"
            handle.write(
                f"chr1\t{start}\t{end}\t{name}\t{score:.4f}\t{strand}\n"
            )


def _build_small_region_workload(root: Path) -> tuple[list[dict], int]:
    """10 BED files x 100 records, all inside a 10 KB window on chr1."""
    rng = random.Random(0)
    specs: list[dict] = []
    total = 0
    for file_index in range(10):
        path = root / f"small_region_{file_index}.bed"
        _write_bed(path, n_records=100, span_bp=10_000, rng=rng)
        specs.append(_build_spec(path, label=f"small_region_{file_index}"))
        total += 100
    return specs, total


def _build_whole_track_workload(root: Path) -> tuple[list[dict], int]:
    """5 BED files x 100_000 records, tiled across chr1 at 100 bp stride."""
    rng = random.Random(1)
    specs: list[dict] = []
    total = 0
    for file_index in range(5):
        path = root / f"whole_track_{file_index}.bed"
        _write_bed(path, n_records=100_000, span_bp=100, rng=rng)
        specs.append(_build_spec(path, label=f"whole_track_{file_index}"))
        total += 100_000
    return specs, total


def _build_mixed_workload(
    small_root: Path, whole_root: Path
) -> tuple[list[dict], int]:
    """5 small-region files + 5 whole-track files submitted as one batch."""
    small_specs, small_total = _build_small_region_workload(small_root)
    whole_specs, whole_total = _build_whole_track_workload(whole_root)
    return small_specs + whole_specs, small_total + whole_total


def _build_spec(path: Path, label: str) -> dict:
    """Return a spec dict shaped exactly like ``_prepare_genomic_bed_track`` reads."""
    return {
        "files": [str(path)],
        "label": label,
        "parser_type": "distribution",
        "color": "orange",
        "alpha": 0.8,
        "file_colors": ["orange"],
        "file_alphas": [0.8],
        "y_axis_range": None,
    }


WORKLOAD_BUILDERS = {
    "small_region": lambda small, _whole: _build_small_region_workload(small),
    "whole_track": lambda _small, whole: _build_whole_track_workload(whole),
    "mixed": lambda small, whole: _build_mixed_workload(small, whole),
}


# ---------------------------------------------------------------------------
# Worker (module-level so ProcessPoolExecutor can pickle it)
# ---------------------------------------------------------------------------

def _load(spec: dict):
    """Call the audited build helper and return the prepared BEDParser.

    Returning the parser keeps the ``ProcessPoolExecutor`` path
    faithful to ``_prepare_process_batch``: each worker sends back a
    several-MB BEDParser over the pipe, which is the suspected
    pickle / unpickle cost we want to measure.
    """
    return _track_build._prepare_genomic_bed_track(spec)


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------

def run_sequential(specs: list[dict]) -> list:
    return [_load(spec) for spec in specs]


def run_threads(specs: list[dict], workers: int) -> list:
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_load, spec) for spec in specs]
        return [future.result() for future in futures]


def run_processes(specs: list[dict], workers: int) -> list:
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_load, spec) for spec in specs]
        return [future.result() for future in futures]


PRIMITIVE_RUNNERS = {
    "sequential": lambda specs, _workers: run_sequential(specs),
    "thread_pool": run_threads,
    "process_pool": run_processes,
}


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------

def _peak_rss_mb() -> float:
    """Return max(self.ru_maxrss, children.ru_maxrss) in MiB.

    On Linux both numbers are already in KiB; on macOS they are in
    bytes.  The drvizer benchmark host is Linux so the divisor is
    1024 to convert KiB to MiB.  See ``meta.platform`` in the output
    JSON for the runtime platform.
    """
    self_kib = getrusage(RUSAGE_SELF).ru_maxrss
    children_kib = getrusage(RUSAGE_CHILDREN).ru_maxrss
    if platform.system() == "Darwin":
        # macOS reports bytes.
        return max(self_kib, children_kib) / (1024 * 1024)
    # Linux reports KiB.
    return max(self_kib, children_kib) / 1024


def _time_once(fn, *args, **kwargs) -> float:
    start = time.perf_counter()
    fn(*args, **kwargs)
    return time.perf_counter() - start


def _measure(
    primitive_name: str, specs: list[dict], workers: int, repetitions: int
) -> dict:
    """Run one configuration: warmup + N timed batches, return metrics."""
    runner = PRIMITIVE_RUNNERS[primitive_name]

    def _invoke():
        return runner(specs, workers)

    # Uncounted warmup: drives pool startup, Cython import caching, etc.
    _invoke()

    wall_times: list[float] = []
    for _ in range(repetitions):
        wall_times.append(_time_once(_invoke))

    return {
        "wall_time_s": {
            "min": min(wall_times),
            "median": statistics.median(wall_times),
            "max": max(wall_times),
        },
        "peak_rss_mb": _peak_rss_mb(),
        # records_per_s is filled in by the caller once it knows the
        # workload's total record count.
    }


def _records_per_s(median_seconds: float, total_records: int) -> float:
    if median_seconds <= 0:
        return float("inf")
    return total_records / median_seconds


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def _format_table(results: dict) -> str:
    header = (
        f"{'workload':<13} {'primitive':<13} {'workers':>7} "
        f"{'median s':>10} {'peak MB':>9} {'records/s':>11}"
    )
    rule = "-" * len(header)
    lines = [header, rule]
    for workload, primitives in results.items():
        if workload == "meta":
            continue
        for primitive, workers_map in primitives.items():
            for workers, metrics in workers_map.items():
                lines.append(
                    f"{workload:<13} {primitive:<13} {workers:>7} "
                    f"{metrics['wall_time_s']['median']:>10.3f} "
                    f"{metrics['peak_rss_mb']:>9.1f} "
                    f"{metrics['records_per_s']:>11.0f}"
                )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _build_meta() -> dict:
    try:
        from drvizer.bed_parser import _CYTHON_BED_AVAILABLE  # type: ignore
        cython_available = bool(_CYTHON_BED_AVAILABLE)
    except Exception:
        cython_available = False
    try:
        import multiprocessing as mp
        start_method = mp.get_start_method()
    except Exception:
        start_method = "unknown"
    return {
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "mp_start_method": start_method,
        "cython_bed_available": cython_available,
        "peak_rss_note": (
            "ru_maxrss is a per-process high-water mark that never "
            "decreases, so values for sequential / thread_pool are "
            "monotonically non-decreasing across configurations in "
            "the same interpreter; treat cross-primitive comparisons "
            "as ordinal, not absolute."
        ),
    }


def _resolve_workloads(args: argparse.Namespace) -> list[str]:
    selected = []
    if args.small_region or (not args.small_region and not args.whole_track):
        selected.append("small_region")
    if args.whole_track or (not args.small_region and not args.whole_track):
        selected.append("whole_track")
    if not args.small_region and not args.whole_track:
        # "all three" implies mixed as well.
        selected.append("mixed")
    elif args.small_region and args.whole_track:
        selected.append("mixed")
    return selected


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Phase 2.2.0 BED-concurrency data-collection benchmark. "
            "Drives drvizer._track_build._prepare_genomic_bed_track over "
            "three workloads x three primitives x four worker counts and "
            "writes JSON plus a human-readable table."
        )
    )
    parser.add_argument(
        "--output",
        default="benchmarks/results.json",
        help="Path to write JSON results (default: benchmarks/results.json).",
    )
    parser.add_argument(
        "--repetitions",
        type=int,
        default=3,
        help="Number of timed batches per configuration (default: 3).",
    )
    parser.add_argument(
        "--small-region",
        action="store_true",
        help="Run only the small-region workload.",
    )
    parser.add_argument(
        "--whole-track",
        action="store_true",
        help="Run only the whole-track workload.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv if argv is not None else sys.argv[1:])
    workloads = _resolve_workloads(args)

    results: dict = {"meta": _build_meta()}

    try:
        for workload_name in workloads:
            print(f"[bench] generating fixtures for {workload_name} ...", flush=True)
            small_root = Path(tempfile.mkdtemp(prefix="bench_small_"))
            whole_root = Path(tempfile.mkdtemp(prefix="bench_whole_"))
            try:
                specs, total_records = WORKLOAD_BUILDERS[workload_name](
                    small_root, whole_root
                )
            finally:
                # Fixtures are read-only after generation; tear them down
                # immediately.  Re-create them per repetition? No -- the
                # benchmark reads from the same files every batch.
                pass

            workload_results: dict = {}
            for primitive_name in PRIMITIVE_RUNNERS:
                primitive_results: dict = {}
                worker_counts = [1] if primitive_name == "sequential" else [1, 2, 4, 8]
                for workers in worker_counts:
                    key = str(workers)
                    print(
                        f"[bench]   {workload_name} / {primitive_name} / "
                        f"workers={workers} ... ",
                        end="",
                        flush=True,
                    )
                    metrics = _measure(
                        primitive_name, specs, workers, args.repetitions
                    )
                    metrics["records_per_s"] = _records_per_s(
                        metrics["wall_time_s"]["median"], total_records
                    )
                    primitive_results[key] = metrics
                    print(
                        f"median={metrics['wall_time_s']['median']:.3f}s "
                        f"rss={metrics['peak_rss_mb']:.1f}MB",
                        flush=True,
                    )
                workload_results[primitive_name] = primitive_results
            results[workload_name] = workload_results

            # Tear down fixture directories now that this workload is done.
            for root in (small_root, whole_root):
                for child in root.iterdir():
                    try:
                        child.unlink()
                    except FileNotFoundError:
                        pass
                try:
                    root.rmdir()
                except OSError:
                    pass
    except Exception:
        raise

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wt") as handle:
        json.dump(results, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"[bench] wrote {output_path}", flush=True)

    if sys.stdout.isatty():
        print(_format_table(results))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())