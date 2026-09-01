"""Real BAM/BAI integration test for PARSER-015.

This test exercises ``BAMParser.get_coverage_in_region`` against a real
``pysam.AlignmentFile`` written via ``pysam.AlignmentFile(..., mode='wb')``
and indexed via ``pysam.index``. We deliberately avoid fake/mock
``AlignmentFile`` objects so the test catches genuine integration
regressions in the Cython build, header parsing, and index lookup paths.

If pysam cannot be imported in the DRS env, the test module is skipped at
import time via ``pytest.importorskip('pysam')``. This matches the audit
hygiene contract: the file may be skipped in CI but should still catch real
BAM regressions when run locally.
"""

import pytest

pytest.importorskip("pysam")

import numpy as np

from drvizer.bam_parser import BAMParser


def test_bam_parser_get_coverage_in_region_real_bam(synthetic_bam):
    """Real pysam BAM/BAI integration for get_coverage_in_region.

    The fixture ``synthetic_bam`` writes 2 records to ``chr1``:
      - read1: blocks [100, 103) -> 3-base coverage at 100, 101, 102.
      - read2: blocks [200, 202) -> 2-base coverage at 200, 201.

    Query the full chromosome extent [0, 500) and assert:
      - The x array spans [0, 500).
      - coverage at positions 100, 101, 102 is 1 (read1).
      - coverage at positions 200, 201 is 1 (read2).
      - coverage elsewhere is 0.
    """
    bam_path = str(synthetic_bam)

    parser = BAMParser(bam_path)
    x, y = parser.get_coverage_in_region("chr1", 0, 500)

    # x is np.arange(0, 500) and y is the per-base coverage.
    assert len(x) == 500
    assert len(y) == 500

    # Read1 contributes 3-base coverage at positions 100, 101, 102.
    assert int(y[100]) == 1
    assert int(y[101]) == 1
    assert int(y[102]) == 1

    # Read2 contributes 2-base coverage at positions 200, 201.
    assert int(y[200]) == 1
    assert int(y[201]) == 1

    # Positions outside either record should be zero coverage.
    assert int(y[0]) == 0
    assert int(y[150]) == 0
    assert int(y[300]) == 0

    # Query just the read1 window: [99, 104) -> positions 99, 100, 101, 102, 103.
    x_sub, y_sub = parser.get_coverage_in_region("chr1", 99, 104)
    assert len(x_sub) == 5
    assert x_sub.tolist() == [99, 100, 101, 102, 103]
    assert y_sub.tolist() == [0, 1, 1, 1, 0]

    # Query just the read2 window: [199, 203) -> positions 199, 200, 201, 202.
    x_sub2, y_sub2 = parser.get_coverage_in_region("chr1", 199, 203)
    assert x_sub2.tolist() == [199, 200, 201, 202]
    assert y_sub2.tolist() == [0, 1, 1, 0]


def test_bam_parser_empty_region_returns_empty_arrays(synthetic_bam):
    """A zero-width region should return empty arrays without error."""
    bam_path = str(synthetic_bam)
    parser = BAMParser(bam_path)

    x, y = parser.get_coverage_in_region("chr1", 300, 300)
    assert len(x) == 0
    assert len(y) == 0