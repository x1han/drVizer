from pathlib import Path
from time import perf_counter

from drvizer.bed_parser import BEDParser
from drvizer.gtf_parser import GTFParser


BENCHMARK_ROOT = Path(__file__).resolve().parent.parent


def timed(label, fn):
    started = perf_counter()
    result = fn()
    elapsed = perf_counter() - started
    print(f"{label}: {elapsed:.6f}s")
    return result, elapsed


def main():
    print(f"Benchmark root: {BENCHMARK_ROOT}")

    gtf_candidates = [
        BENCHMARK_ROOT / "tests" / "data" / "sample.gtf",
        BENCHMARK_ROOT / "tests" / "data" / "sample.gtf.gz",
    ]
    bed_candidates = [
        BENCHMARK_ROOT / "tests" / "data" / "sample.bed",
        BENCHMARK_ROOT / "tests" / "data" / "sample.bed.gz",
    ]

    gtf_path = next((path for path in gtf_candidates if path.exists()), None)
    bed_path = next((path for path in bed_candidates if path.exists()), None)

    if gtf_path is None:
        print("Skipping GTF benchmark: no sample GTF data found under tests/data.")
    else:
        timed("GTF parse", lambda: GTFParser(str(gtf_path)).parse_gtf())

    if bed_path is None:
        print("Skipping BED benchmark: no sample BED data found under tests/data.")
    else:
        timed("BED parse", lambda: BEDParser(str(bed_path)).parse_bed())


if __name__ == "__main__":
    main()
