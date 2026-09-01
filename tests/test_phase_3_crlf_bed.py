"""Phase 3 -- PARSER-011: CRLF BED parity tests.

Two tests, each anchored to the PARSER-011 mandate:

* test_crlf_bed_parity
    Byte-identical records across the Cython and Python BED paths when
    the fixture is a Windows-style CRLF BED.
* test_crlf_fixture_exists
    The fixture is present and uses ``\\r\\n`` line endings (we check
    raw bytes so a text-mode rewrite never silently strips the CRLF).
"""

from pathlib import Path

import pytest

from drvizer import bed_parser
from drvizer.bed_parser import parse_bed_records_python

FIXTURE_PATH = Path(__file__).resolve().parent / "fixtures" / "crlf.bed"


def test_crlf_fixture_exists():
    """tests/fixtures/crlf.bed must exist and use CRLF line endings."""
    assert FIXTURE_PATH.exists(), f"missing fixture: {FIXTURE_PATH}"
    raw = FIXTURE_PATH.read_bytes()
    # A real CRLF BED must contain at least one \r\n sequence.
    assert b"\r\n" in raw, "fixture is not CRLF-terminated"
    # And it must not contain a bare LF (every \n must be preceded by \r).
    bare_lf = sum(1 for i, byte in enumerate(raw) if byte == 0x0A and (i == 0 or raw[i - 1] != 0x0D))
    assert bare_lf == 0, f"fixture contains {bare_lf} bare-LF characters"


def test_crlf_bed_parity():
    """Python and Cython paths must produce byte-identical records on a CRLF BED.

    The fixture carries records on ``chr1`` and ``chr2`` plus a comment and a
    blank line so the comment/blank skipper is exercised in both paths.
    """
    python_records = parse_bed_records_python([str(FIXTURE_PATH)])
    python_records = dict(python_records)

    cython_records = None
    if bed_parser._CYTHON_BED_AVAILABLE and bed_parser.parse_bed_records is not None:
        cython_records = bed_parser.parse_bed_records([str(FIXTURE_PATH)])
        cython_records = dict(cython_records)

    assert "chr1" in python_records
    assert len(python_records["chr1"]) == 2
    assert python_records["chr1"][0]["start"] == 105
    assert python_records["chr1"][0]["end"] == 120
    assert python_records["chr1"][0]["name"] == "peak1"
    assert python_records["chr2"][0]["start"] == 500
    assert python_records["chr2"][0]["end"] == 550

    if cython_records is not None:
        # Same record count per chromosome
        assert set(python_records.keys()) == set(cython_records.keys())
        for chrom in python_records:
            assert len(python_records[chrom]) == len(cython_records[chrom])
            for p_rec, c_rec in zip(python_records[chrom], cython_records[chrom]):
                assert p_rec == c_rec, (
                    f"parity mismatch on {chrom}: python={p_rec} cython={c_rec}"
                )
