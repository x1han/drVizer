"""Phase 8 -- tuple-return worker contract (P1-8 / Phase 7 PARTIALLY-FIXED -> FIXED).

The Phase 7 implementation satisfied the imap_unordered attribution
invariant via a closure binding ``bam_path = args[0]`` inside the
worker. Phase 8 makes the same invariant explicit in the DATA shape:
the worker returns ``Tuple[str, np.ndarray]`` -- (bam_path,
coverage_array) -- so attribution is encoded in the return value
instead of a closure binding.

Tests:

* ``test_worker_returns_bam_path_and_coverage_tuple``
    The ``_compute_region_coverage_with_path`` wrapper returns a
    (bam_path, ndarray) tuple with ``bam_path`` equal to ``args[0]``
    and a coverage array of the correct length on the success path.

* ``test_master_loop_unpacks_tuples_out_of_order``
    A monkeypatched ``Pool.imap_unordered`` returns
    (bam_path, reversed_array) tuples in reverse arrival order; the
    aggregated int64 result equals the per-path sum regardless of
    arrival order.

* ``test_parallel_coverage_error_preserves_original_cause``
    When the worker raises, the ``ParallelCoverageError`` raised by
    the worker (and re-surfaced by the master loop) preserves the
    original worker exception on ``__cause__`` across the
    ``imap_unordered`` boundary.
"""

import numpy as np
import pytest

from drvizer import _parallel
from drvizer._parallel import ParallelCoverageError


def test_worker_returns_bam_path_and_coverage_tuple(monkeypatch):
    """``_compute_region_coverage_with_path`` returns ``(bam_path, ndarray)``
    on the success path with ``bam_path`` equal to ``args[0]``."""

    expected_length = 17

    def fake_compute_region_coverage(path, chrom, start, end, contained_only=True):
        assert path == "/data/sample_A.bam"
        assert chrom == "chr1"
        assert start == 100
        assert end == 100 + expected_length
        return np.arange(expected_length, dtype=np.int64)

    monkeypatch.setattr(_parallel, "compute_region_coverage", fake_compute_region_coverage)

    args = ("/data/sample_A.bam", "chr1", 100, 100 + expected_length, True)
    result = _parallel._compute_region_coverage_with_path(args)

    assert isinstance(result, tuple), (
        f"worker must return a tuple, got {type(result).__name__}"
    )
    assert len(result) == 2, (
        f"worker tuple must have length 2 (bam_path, coverage), got {len(result)}"
    )
    bam_path, coverage = result
    assert bam_path == "/data/sample_A.bam", (
        f"first element of worker tuple must be the bam_path from args[0]; got {bam_path!r}"
    )
    assert isinstance(coverage, np.ndarray), (
        f"second element of worker tuple must be a numpy ndarray; got {type(coverage).__name__}"
    )
    assert coverage.shape == (expected_length,), (
        f"coverage array must have the expected region length; got shape {coverage.shape}"
    )


def test_master_loop_unpacks_tuples_out_of_order(monkeypatch):
    """``aggregate_region_coverages_parallel`` must unpack each
    ``(bam_path, coverage_array)`` tuple from ``imap_unordered`` even
    when results arrive in reverse order, and produce the per-path sum
    (int64 accumulator, order-independent)."""

    captured_workers = []

    class FakePool:
        def __init__(self, processes):
            captured_workers.append(processes)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def imap_unordered(self, function, args):
            # Phase 8 tuple-return contract: each result is
            # (bam_path, coverage_array). Return tuples in REVERSE
            # order so the consumer's tuple-unpacking is exercised
            # with shuffled arrival. The coverage payload carries a
            # distinctive fill value (i+1) so an order bug would
            # produce a visibly wrong array; here it does not, because
            # the int64 accumulator is order-independent.
            ordered_results = [
                (args[i][0], np.full(args[i][3] - args[i][2], i + 1, dtype=np.int64))
                for i in range(len(args))
            ]
            return list(reversed(ordered_results))

    monkeypatch.setattr(_parallel.multiprocessing, "Pool", FakePool)
    monkeypatch.setattr(_parallel.multiprocessing, "cpu_count", lambda: 4)

    paths = [f"/tmp/sample_{i}.bam" for i in range(3)]
    region_len = 12
    expected = np.zeros(region_len, dtype=np.int64)
    for i, _path in enumerate(paths):
        expected += np.full(region_len, i + 1, dtype=np.int64)

    coverage = _parallel.aggregate_region_coverages_parallel(paths, "chr1", 0, region_len)

    assert captured_workers == [3], (
        f"Pool must spin up 3 workers (capped by len(paths)=3 vs cpu_count=4); got {captured_workers}"
    )
    assert np.array_equal(coverage, expected), (
        f"aggregated coverage must equal the per-path sum regardless of "
        f"imap_unordered arrival order; got {coverage} vs expected {expected}"
    )


def test_parallel_coverage_error_preserves_original_cause(monkeypatch):
    """When the worker raises, the resulting ``ParallelCoverageError``
    must chain the original worker exception via ``raise ... from exc``
    so ``__cause__`` is preserved across the ``imap_unordered``
    boundary. Regression protection for P1-8 / Phase 7 invariant."""

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
            # Invoke the real worker so the raise propagates through
            # the tuple-return code path the way it does in production.
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
        f"original worker exception must be preserved on __cause__ across the "
        f"imap_unordered boundary; got {type(excinfo.value.__cause__).__name__}"
    )