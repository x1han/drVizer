"""Phase 3 -- P1-8: ParallelCoverageError chain tests.

One test:

* test_error_chain_preserves_cause_and_message
    Raise ParallelCoverageError from a synthetic OSError inside the
    per-bam-path worker, catch the wrapper, and assert:
    - ``exc.__cause__`` is the original OSError (chained)
    - the path and 'corrupt' substring both appear in ``str(exc.__cause__)``
"""

import pytest

from drvizer._parallel import ParallelCoverageError


def test_error_chain_preserves_cause_and_message(monkeypatch):
    """The worker wrapper must chain the original exception and include
    the bam_path context in the raised ParallelCoverageError.

    We monkeypatch ``compute_region_coverage`` to raise a synthetic
    OSError with 'corrupt' in its message, then call the wrapper and
    inspect the resulting chain.
    """
    from drvizer import _parallel

    def fake_compute_region_coverage(path, chrom, start, end, contained_only=True):
        raise OSError("BAM file is corrupt")

    monkeypatch.setattr(_parallel, "compute_region_coverage", fake_compute_region_coverage)

    args = ("/tmp/missing_or_corrupt.bam", "chr1", 0, 100, True)
    with pytest.raises(ParallelCoverageError) as excinfo:
        _parallel._compute_region_coverage_with_path(args)

    exc = excinfo.value
    # The wrapper chained the original OSError.
    assert exc.__cause__ is not None
    assert isinstance(exc.__cause__, OSError)
    assert "corrupt" in str(exc.__cause__)
    # The bam_path appears in the wrapper's message.
    assert "/tmp/missing_or_corrupt.bam" in str(exc)
