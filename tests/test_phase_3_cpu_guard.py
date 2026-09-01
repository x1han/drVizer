"""Phase 3 -- P1-5: cpu_count guard tests.

Two tests:

* test_cpu_count_none_falls_back_to_one
    monkeypatches ``multiprocessing.cpu_count`` to return None and
    asserts ``worker_count`` falls back to 1 with no TypeError.
* test_cpu_guard_present_at_all_entry_points
    Static grep over the four files flagged by P1-5 confirms the
    ``or 1`` (or equivalent isinstance(int)) fallback is present at the
    only two real concurrency entry points: ``_parallel.py`` and
    ``_track_build.py``. ``BedParser`` and ``BAMParser`` have no
    concurrency entry points and are confirmed to NOT contain
    ``multiprocessing.Pool`` / ``ProcessPoolExecutor`` / ``cpu_count``.
"""

import inspect
import multiprocessing

import numpy as np
import pytest

from drvizer import _parallel, _track_build


def test_cpu_count_none_falls_back_to_one(monkeypatch):
    """When ``multiprocessing.cpu_count`` returns None, worker_count must be 1."""

    pool_ref = {}

    class _CapturePool:
        def __init__(self, processes):
            pool_ref["processes"] = processes

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def imap_unordered(self, function, args):
            # Return a coverage ndarray per task so the caller can sum.
            return [
                np.zeros(args[0][3] - args[0][2], dtype=np.int64) for _ in args
            ]

    monkeypatch.setattr(_parallel.multiprocessing, "cpu_count", lambda: None)
    monkeypatch.setattr(_parallel.multiprocessing, "Pool", _CapturePool)

    # Two BAM paths, cpu_count returns None -> worker_count must be 1.
    _parallel.aggregate_region_coverages_parallel(["a.bam", "b.bam"], "chr1", 0, 1)
    assert pool_ref["processes"] == 1, (
        f"cpu_count None must fall back to 1 worker, got {pool_ref['processes']}"
    )


def test_cpu_guard_present_at_all_entry_points():
    """Static check: the ``or 1`` (or isinstance fallback) lives at every
    concurrency entry point. The two real concurrency entry points are
    ``_parallel.py`` and ``_track_build.py``. BedParser / BAMParser do
    not have concurrency code paths, so they are checked for the
    ABSENCE of any cpu_count reference.
    """
    parallel_src = inspect.getsource(_parallel)
    track_build_src = inspect.getsource(_track_build)

    assert "cpu_count" in parallel_src, "_parallel.py must reference cpu_count"
    assert "isinstance" in parallel_src or "or 1" in parallel_src, (
        "_parallel.py must guard cpu_count() for None"
    )

    assert "cpu_count" in track_build_src, "_track_build.py must reference cpu_count"
    assert "or 1" in track_build_src, (
        "_track_build.py must use the `os.cpu_count() or 1` fallback"
    )

    # Confirm BedParser/BAMParser do NOT have any concurrency code paths
    # (they were the pre-review 'over-claim' concern in P1-5).
    from drvizer import bed_parser, bam_parser

    bed_src = inspect.getsource(bed_parser)
    bam_src = inspect.getsource(bam_parser)

    for token in ("multiprocessing.Pool", "ProcessPoolExecutor", "cpu_count"):
        assert token not in bed_src, (
            f"BedParser has no concurrency code path; found {token}"
        )
        assert token not in bam_src, (
            f"BAMParser has no concurrency code path; found {token}"
        )
