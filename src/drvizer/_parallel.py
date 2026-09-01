import multiprocessing

import numpy as np
import pysam


# P1-8: A worker that crashes (corrupt BAM, missing .bai, KeyError, BrokenPipe,
# etc.) used to be silently swallowed into a generic ParallelCoverageError,
# hiding the original cause. We now catch the full Exception tree so we can
# chain the worker exception via ``raise ... from exc`` and preserve the
# traceback so the root cause is visible to the caller.


class ParallelCoverageError(RuntimeError):
    pass


def compute_region_coverage(path, chrom, start, end, contained_only=True):
    region_len = end - start
    coverage = np.zeros(region_len, dtype=np.int64)

    with pysam.AlignmentFile(path, "rb") as sam:
        for read in sam.fetch(chrom, start, end):
            if contained_only:
                if read.reference_start < start or read.reference_end > end:
                    continue

            for b_start, b_end in read.get_blocks():
                idx_s = max(0, b_start - start)
                idx_e = min(region_len, b_end - start)
                if idx_s < idx_e:
                    coverage[idx_s:idx_e] += 1

    return coverage


def _compute_region_coverage(args):
    return compute_region_coverage(*args)


# P1-8: Per-bam-path wrapper -- the panel mandate says we must catch the full
# worker exception, chain it, and surface a per-path message. The Pool
# worker returns the coverage ndarray on success and raises the original
# exception on failure; we re-raise here as ParallelCoverageError chained
# from the worker exception so callers see both the BAM path AND the root
# cause. The pool's imap_unordered iteration then surfaces a single
# ParallelCoverageError with the original traceback attached via __cause__.
def _compute_region_coverage_with_path(args):
    bam_path = args[0]
    try:
        return compute_region_coverage(*args)
    except Exception as exc:
        raise ParallelCoverageError(
            f"Failed to compute coverage for {bam_path}: {exc}"
        ) from exc


def aggregate_region_coverages_parallel(bam_paths, chrom, start, end, contained_only=True):
    if len(bam_paths) <= 1:
        return compute_region_coverage(bam_paths[0], chrom, start, end, contained_only)

    # P1-5: ``cpu_count`` may return ``None`` on BSDs / restricted containers;
    # the original audit flagged this as a TypeError that bypassed the
    # except clause. Use multiprocessing.cpu_count() (which is what
    # existing tests monkeypatch) plus the ``or 1`` / isinstance fallback
    # so the arithmetic stays well-defined.
    try:
        cpu_count = multiprocessing.cpu_count()
    except NotImplementedError:
        cpu_count = 1
    if not isinstance(cpu_count, int):
        cpu_count = 1
    worker_count = min(len(bam_paths), cpu_count, 32)
    total_coverage = np.zeros(end - start, dtype=np.int64)
    with multiprocessing.Pool(processes=worker_count) as pool:
        for coverage in pool.imap_unordered(
            _compute_region_coverage_with_path,
            [
                (path, chrom, start, end, contained_only)
                for path in bam_paths
            ],
        ):
            total_coverage += coverage.astype(np.int64, copy=False)

    return total_coverage
