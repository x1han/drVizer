import multiprocessing

import numpy as np
import pysam


class ParallelCoverageError(RuntimeError):
    pass


def compute_region_coverage(path, chrom, start, end, contained_only=True):
    region_len = end - start
    coverage = np.zeros(region_len, dtype=np.int32)

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


def aggregate_region_coverages_parallel(bam_paths, chrom, start, end, contained_only=True):
    if len(bam_paths) <= 1:
        return compute_region_coverage(bam_paths[0], chrom, start, end, contained_only)

    try:
        with multiprocessing.Pool(processes=len(bam_paths)) as pool:
            coverages = pool.map(
                _compute_region_coverage,
                [
                    (path, chrom, start, end, contained_only)
                    for path in bam_paths
                ],
            )
    except (OSError, ValueError) as exc:
        raise ParallelCoverageError(str(exc)) from exc

    total_coverage = np.zeros(end - start, dtype=np.int32)
    for coverage in coverages:
        total_coverage += coverage
    return total_coverage
