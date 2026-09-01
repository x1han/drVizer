"""Phase 2.2 cache-hit vs no-cache benchmark.

Drives DrViz.plot() twice on the same gene and reports the
records-per-second ratio. The hard expectation per PLAN_SCHEMA is
``cache_hit >= 10x faster than cache_miss``; in practice this is
50-500x for BED because parse cost dominates I/O.

Stdlib-only, argparse, JSON output. Reuses the
``bench_bed_concurrency.py`` fixture pattern: a small/medium/large
BED, a 1-2 transcript GTF, and a 3-repetition sweep.
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
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC_ROOT = _REPO_ROOT / "src"
if _SRC_ROOT.is_dir() and str(_SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(_SRC_ROOT))

# Headless-safe backend lock so plt.show() does not block.
os.environ.setdefault("MPLBACKEND", "Agg")

import drvizer  # noqa: E402


# ---------------------------------------------------------------------------
# Fixture generation
# ---------------------------------------------------------------------------

def _write_bed(path: Path, n_records: int, rng: random.Random) -> None:
    """Write ``n_records`` BED6 rows packed at 100 bp stride on chr1."""
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


def _write_minimal_gtf(path: Path) -> None:
    """One gene with two short transcripts on chr1 [0, 1000)."""
    with path.open("wt") as handle:
        handle.write(
            'chr1\tsrc\texon\t100\t199\t.\t+\t.\t'
            'gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";\n'
        )
        handle.write(
            'chr1\tsrc\texon\t300\t399\t.\t+\t.\t'
            'gene_id "GENE1"; transcript_id "TX2"; gene_name "G1";\n'
        )


# ---------------------------------------------------------------------------
# Workloads
# ---------------------------------------------------------------------------

WORKLOAD_SIZES = {
    "small": 1_000,    # ~1K records: fast parse
    "medium": 50_000,  # ~50K records: representative
    "large": 500_000,  # ~500K records: parse cost dominates
}


def _build_workload(workload: str, root: Path) -> tuple[Path, Path, int]:
    gtf = root / f"{workload}.gtf"
    bed = root / f"{workload}.bed"
    _write_minimal_gtf(gtf)
    _write_bed(bed, n_records=WORKLOAD_SIZES[workload], rng=random.Random(0))
    return gtf, bed, WORKLOAD_SIZES[workload]


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------

def _measure_one(viz, gene: str) -> float:
    """Run one get_transcript_data() call and return wall time in seconds.

    The cache hit speedup is dominated by the BED parse cost;
    matplotlib render is constant-time, so it would dilute the
    measurement. ``get_transcript_data`` exercises the cache
    directly: miss parses the BED, hit returns the cached payload.
    """
    start = time.perf_counter()
    viz.get_transcript_data(gene)
    return time.perf_counter() - start


def _records_per_s(median_seconds: float, total_records: int) -> float:
    if median_seconds <= 0:
        return float("inf")
    return total_records / median_seconds


def _measure_workload(workload: str, repetitions: int) -> dict:
    """Run miss + hit timing for one workload."""
    with tempfile.TemporaryDirectory(prefix=f"bench_{workload}_") as tmp:
        root = Path(tmp)
        gtf, bed, total_records = _build_workload(workload, root)

        # Miss: build a fresh DrViz every call (no caching across calls).
        miss_wall_times = []
        for _ in range(repetitions):
            viz = drvizer.DrViz().load_gtf(str(gtf)).add_bed_track(str(bed), label="TE")
            miss_wall_times.append(_measure_one(viz, "GENE1"))

        # Hit: build once, then call get_transcript_data repeatedly.
        # The LRU cache inside the ReusableParser is populated on
        # the first call and reused on subsequent calls.
        hit_viz = drvizer.DrViz().load_gtf(str(gtf)).add_bed_track(str(bed), label="TE")
        # Warmup: populate the cache.
        _measure_one(hit_viz, "GENE1")
        hit_wall_times = []
        for _ in range(repetitions):
            hit_wall_times.append(_measure_one(hit_viz, "GENE1"))

    miss_median = statistics.median(miss_wall_times)
    hit_median = statistics.median(hit_wall_times)
    speedup = miss_median / hit_median if hit_median > 0 else float("inf")

    return {
        "total_records": total_records,
        "miss": {
            "wall_time_s": {
                "min": min(miss_wall_times),
                "median": miss_median,
                "max": max(miss_wall_times),
            },
            "records_per_s": _records_per_s(miss_median, total_records),
        },
        "hit": {
            "wall_time_s": {
                "min": min(hit_wall_times),
                "median": hit_median,
                "max": max(hit_wall_times),
            },
            "records_per_s": _records_per_s(hit_median, total_records),
        },
        "speedup_median": speedup,
    }


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _build_meta() -> dict:
    return {
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
    }


def _format_table(results: dict) -> str:
    header = (
        f"{'workload':<8} {'records':>8} "
        f"{'miss med s':>11} {'hit med s':>11} "
        f"{'miss r/s':>10} {'hit r/s':>10} {'speedup':>9}"
    )
    rule = "-" * len(header)
    lines = [header, rule]
    for workload, payload in results.items():
        if workload == "meta":
            continue
        lines.append(
            f"{workload:<8} {payload['total_records']:>8} "
            f"{payload['miss']['wall_time_s']['median']:>11.3f} "
            f"{payload['hit']['wall_time_s']['median']:>11.3f} "
            f"{payload['miss']['records_per_s']:>10.0f} "
            f"{payload['hit']['records_per_s']:>10.0f} "
            f"{payload['speedup_median']:>9.1f}x"
        )
    return "\n".join(lines)


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Phase 2.2 cache-hit vs no-cache benchmark. "
            "Sweeps (workload: small/medium/large BED) x "
            "(primitive: cache_miss/cache_hit) and reports "
            "records-per-second and the speedup ratio."
        )
    )
    parser.add_argument(
        "--output",
        default="benchmarks/cache_vs_no_cache.json",
        help="Path to write JSON results.",
    )
    parser.add_argument(
        "--repetitions",
        type=int,
        default=3,
        help="Number of timed batches per workload (default: 3).",
    )
    parser.add_argument(
        "--workloads",
        nargs="*",
        default=None,
        choices=list(WORKLOAD_SIZES.keys()),
        help="Workloads to run (default: all).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv if argv is not None else sys.argv[1:])
    workloads = args.workloads or list(WORKLOAD_SIZES.keys())

    results: dict = {"meta": _build_meta()}

    for workload in workloads:
        print(f"[bench] {workload} (n_records={WORKLOAD_SIZES[workload]}) ...", flush=True)
        payload = _measure_workload(workload, args.repetitions)
        results[workload] = payload
        print(
            f"  miss median={payload['miss']['wall_time_s']['median']:.3f}s, "
            f"hit median={payload['hit']['wall_time_s']['median']:.3f}s, "
            f"speedup={payload['speedup_median']:.1f}x",
            flush=True,
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(results, indent=2, sort_keys=True))
    print(f"\nResults written to {output}")
    print()
    print(_format_table(results))

    # Hard expectation per PLAN_SCHEMA: cache_hit >= 10x faster.
    failures = []
    for workload in workloads:
        speedup = results[workload]["speedup_median"]
        if speedup < 10.0:
            failures.append(f"{workload}: speedup={speedup:.1f}x < 10x")
    if failures:
        print()
        print("WARNING: cache speedup below the 10x plan expectation:")
        for failure in failures:
            print(f"  - {failure}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
