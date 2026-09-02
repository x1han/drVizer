"""Parallel coverage aggregation for drVizer BAM tracks.

Note: ``multiprocessing.Pool.imap_unordered`` is used by
``aggregate_region_coverages_parallel`` to distribute per-BAM coverage
computes across worker processes. Callers MUST NOT assume any specific
ordering of the yielded ``(bam_path, coverage_array)`` tuples -- the
unordered variant deliberately allows out-of-order completion so slow
workers do not block the fast ones. The current aggregation step only
reads each worker's ``coverage_array`` (the ``bam_path`` is preserved
for downstream attribution but unused here), which makes the reduction
order-independent; if a future change ever needs ordered attribution it
must switch to ``imap`` (or sort explicitly) and document the contract
change.
"""

from typing import Tuple

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


# Phase 8 tuple-return contract: the worker now returns
# ``Tuple[str, np.ndarray]`` -- (bam_path, coverage_array) -- so per-bam
# attribution is encoded in the data shape rather than in a closure
# binding. The exception path still raises ``ParallelCoverageError`` chained
# from the worker exception INSIDE the worker so the ``__cause__`` chain
# survives the imap_unordered boundary intact.
def _compute_region_coverage_with_path(args) -> Tuple[str, np.ndarray]:
    bam_path = args[0]
    try:
        coverage = compute_region_coverage(*args)
    except Exception as exc:
        raise ParallelCoverageError(
            f"Failed to compute coverage for {bam_path}: {exc}"
        ) from exc
    return bam_path, coverage


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
        # Phase 8 tuple contract: each imap_unordered result is
        # (bam_path, coverage_array). Attribution is now a property of
        # the returned data; aggregation is still order-independent
        # because we only read coverage (bam_path is unused here, but
        # preserved in the return value for downstream attribution
        # needs).
        for bam_path, coverage in pool.imap_unordered(
            _compute_region_coverage_with_path,
            [
                (path, chrom, start, end, contained_only)
                for path in bam_paths
            ],
        ):
            total_coverage += coverage.astype(np.int64, copy=False)

    return total_coverage