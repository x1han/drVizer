"""Phase 3 -- P1-3: GTF chunk degradation tests.

Two tests, each anchored to the P1-3 mandate:

* test_malformed_row_preserves_partial_data
    A chunk with 9 valid + 1 non-integer coordinate line must surface
    ``>=9`` valid rows and increment ``self._chunk_parse_degradation``.
* test_malformed_only_emits_warning
    A chunk of all-malformed rows still parses 0 rows but increments the
    counter without raising -- the wrapper must NEVER abort a 10k-line
    chunk on a single bad row.
"""

import textwrap
import warnings

import pytest

from drvizer.gtf_parser import GTFParser


# A valid GTF line template (the only thing we vary is column 3/4).
def _valid_line(idx):
    return (
        f'chr1\tsrc\texon\t{100 + idx * 10}\t{120 + idx * 10}\t.\t+\t.'
        f'\tgene_id "GENE1"; transcript_id "TX{idx}"; gene_name "G1";\n'
    )


def _malformed_line():
    return 'chr1\tsrc\texon\tNOT_AN_INT\t200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX_BAD"; gene_name "G1";\n'


def _write_gtf(tmp_path, lines):
    path = tmp_path / "test.gtf"
    path.write_text("".join(lines))
    return path


def test_malformed_row_preserves_partial_data(tmp_path):
    """9 valid + 1 malformed row in a single chunk: at least 9 rows survive."""
    lines = [_valid_line(i) for i in range(9)] + [_malformed_line()]
    gtf_path = _write_gtf(tmp_path, lines)

    parser = GTFParser(str(gtf_path))
    # The wrapper degrades the chunk on the Cython path; the
    # per-row Python fallback keeps the 9 valid rows. With the Cython
    # kernel active and raising on the bad row, the wrapper must fall
    # back so the chunk is NOT silently dropped.
    parser.parse_gtf("GENE1")

    assert "GENE1" in parser.gene_transcripts, "chunk was silently dropped"
    # At least 9 valid transcripts survive (the bad row must not abort the chunk).
    assert len(parser.gene_transcripts["GENE1"]) >= 9, (
        f"expected >=9 valid transcripts, got {len(parser.gene_transcripts['GENE1'])}"
    )
    # The degradation counter increments for the chunk that fell back.
    assert parser._chunk_parse_degradation >= 1, (
        f"expected _chunk_parse_degradation >= 1, got {parser._chunk_parse_degradation}"
    )


def test_malformed_only_emits_warning(tmp_path):
    """All-malformed chunk: 0 rows survive but counter increments; no raise."""
    lines = [_malformed_line() for _ in range(5)]
    gtf_path = _write_gtf(tmp_path, lines)

    parser = GTFParser(str(gtf_path))
    # Must not raise. The wrapper falls back per-row and the bad rows
    # are skipped, so the chunk degrades but never aborts.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        parser.parse_gtf("GENE1")

    # No valid rows survived (every row was malformed).
    assert "GENE1" not in parser.gene_transcripts
    # Counter still incremented because the chunk-level fallback fired.
    assert parser._chunk_parse_degradation >= 1
