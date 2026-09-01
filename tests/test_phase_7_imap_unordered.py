"""Phase 7 -- imap_unordered attribution invariant (P1-8 / P0-8 verification).

The panel mandate for Phase 7 is that the worker inside
``aggregate_region_coverages_parallel`` must keep the BAM path bound to
its result so that the master loop can always pair coverage to the
correct BAM even though ``imap_unordered`` iterates out of order.

The current implementation in ``drvizer._parallel`` satisfies the
mandate via the closure binding ``bam_path = args[0]`` inside
``_compute_region_coverage_with_path``: when a worker raises, the
re-raised ``ParallelCoverageError`` carries the failing ``bam_path``
in its message; when it returns, the master loop sums the coverage
into the int64 accumulator (order-independent). The tests in this file
lock both behaviors down so a future refactor cannot silently regress
either path.

Tests:

* ``test_worker_binds_bam_path_in_closure``
    The ``_compute_region_coverage_with_path`` wrapper reads
    ``args[0]`` and embeds it in the ``ParallelCoverageError`` message
    on worker failure.

* ``test_imap_unordered_aggregation_is_order_independent``
    A monkeypatched ``Pool.imap_unordered`` returns coverage arrays
    in reversed order from the input BAM paths; the aggregated
    int64 result equals the per-path sum regardless of order.

* ``test_parallel_coverage_error_raises_with_bam_path``
    When the worker raises, the resulting ``ParallelCoverageError``
    carries the failing BAM path in its message and the original
    worker exception is preserved on ``__cause__``.
"""

import numpy as np
import pytest

from drvizer import _parallel
from drvizer._parallel import ParallelCoverageError


def test_worker_binds_bam_path_in_closure(monkeypatch):
    """``_compute_region_coverage_with_path`` must embed ``bam_path`` in
    the ``ParallelCoverageError`` message via the closure binding."""

    def fake_compute_region_coverage(path, chrom, start, end, contained_only=True):
        raise OSError("corrupt BAM index")

    monkeypatch.setattr(_parallel, "compute_region_coverage", fake_compute_region_coverage)

    args = ("/data/sample_corrupt.bam", "chr1", 100, 200, True)
    with pytest.raises(ParallelCoverageError) as excinfo:
        _parallel._compute_region_coverage_with_path(args)

    message = str(excinfo.value)
    assert "/data/sample_corrupt.bam" in message, (
        f"ParallelCoverageError message must include the failing bam_path; got {message!r}"
    )


def test_imap_unordered_aggregation_is_order_independent(monkeypatch):
    """``imap_unordered`` returns results out of order; the master loop
    must still produce the per-path sum. Because the accumulator is a
    single int64 zero array and ``+=`` is order-independent, this works
    regardless of arrival order. The test pins the invariant against a
    shuffled arrival sequence."""

    captured_workers = []

    class FakePool:
        def __init__(self, processes):
            captured_workers.append(processes)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def imap_unordered(self, function, args):
            # Phase 8 tuple-return contract: each task returns
            # (bam_path, coverage_array). We deliberately return the
            # (bam_path, array) tuples in REVERSE order from the input
            # args so the consumer sees results out of order. The
            # master loop's tuple-unpacking handles the reversal
            # transparently.
            ordered_results = [
                (args[i][0], np.full(args[i][3] - args[i][2], i + 1, dtype=np.int64))
                for i in range(len(args))
            ]
            return list(reversed(ordered_results))

    monkeypatch.setattr(_parallel.multiprocessing, "Pool", FakePool)
    monkeypatch.setattr(_parallel.multiprocessing, "cpu_count", lambda: 4)

    paths = [f"/tmp/sample_{i}.bam" for i in range(4)]
    region_len = 50

    # Per-path: path i contributes (i+1) coverage across the whole region.
    # Sum is order-independent.
    expected = np.zeros(region_len, dtype=np.int64)
    for i, _path in enumerate(paths):
        expected += np.full(region_len, i + 1, dtype=np.int64)

    coverage = _parallel.aggregate_region_coverages_parallel(paths, "chr1", 0, region_len)

    assert captured_workers == [4], "Pool must spin up 4 workers"
    assert np.array_equal(coverage, expected), (
        f"aggregated coverage must equal the per-path sum regardless of "
        f"imap_unordered arrival order; got {coverage} vs expected {expected}"
    )


def test_parallel_coverage_error_raises_with_bam_path(monkeypatch):
    """When a worker raises, the master loop must raise
    ``ParallelCoverageError`` with the failing BAM path in its
    message and the original worker exception preserved on
    ``__cause__``."""

    failing_path = "/data/very/specific/missing_index.bam"

    def fake_compute_region_coverage(path, chrom, start, end, contained_only=True):
        if path == failing_path:
            raise KeyError("BAM index not found")
        return np.zeros(end - start, dtype=np.int64)

    class FakePool:
        def __init__(self, processes):
            pass

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def imap_unordered(self, function, args):
            # Each entry is (path, chrom, start, end, contained_only).
            # The worker raises ParallelCoverageError on the failing path;
            # other paths return zero coverage. We invoke the real
            # worker so the raise propagates the way it does in
            # production.
            return [function(a) for a in args]

    monkeypatch.setattr(_parallel.multiprocessing, "Pool", FakePool)
    monkeypatch.setattr(_parallel.multiprocessing, "cpu_count", lambda: 2)
    monkeypatch.setattr(_parallel, "compute_region_coverage", fake_compute_region_coverage)

    paths = ["/data/ok_1.bam", failing_path, "/data/ok_2.bam"]

    with pytest.raises(ParallelCoverageError) as excinfo:
        _parallel.aggregate_region_coverages_parallel(paths, "chr1", 0, 10)

    message = str(excinfo.value)
    assert failing_path in message, (
        f"ParallelCoverageError message must include failing bam_path; got {message!r}"
    )
    assert isinstance(excinfo.value.__cause__, KeyError), (
        "original worker exception must be preserved on __cause__"
    )