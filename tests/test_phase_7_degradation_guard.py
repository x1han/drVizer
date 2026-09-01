"""Phase 7 -- GTF chunk-degradation warning guard (P1-3).

The Cython chunk parser in ``drvizer._cython_gtf`` may raise on a malformed
coordinate. ``GTFParser._parse_chunk_rows`` catches the raise, increments
``self._chunk_parse_degradation``, and falls back to per-row Python parsing.

The end-of-parse ``RuntimeWarning`` that summarizes the degradation count
must only fire when ``self._chunk_parse_degradation > 0``. Clean parses
must not emit the warning.

Tests:

* ``test_no_warning_when_zero_chunks_degraded``
    Parsing a valid GTF with the Cython kernel present does not increment
    ``_chunk_parse_degradation`` and emits no ``RuntimeWarning``.

* ``test_warning_when_chunks_degraded``
    When the Cython kernel raises on every chunk, ``_chunk_parse_degradation``
    is non-zero and exactly one ``RuntimeWarning`` is emitted whose message
    contains the degradation count.
"""

import textwrap
import warnings

import pytest

from drvizer.gtf_parser import GTFParser


def test_no_warning_when_zero_chunks_degraded(tmp_path):
    """Clean GTF: ``_chunk_parse_degradation`` stays at 0, no warning."""
    gtf_path = tmp_path / "clean.gtf"
    gtf_path.write_text(textwrap.dedent(
        """\
        chr1	src	exon	100	149	.	+	.	gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1	src	CDS	110	140	.	+	0	gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1	src	exon	200	249	.	+	.	gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        """
    ))

    parser = GTFParser(str(gtf_path))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        parser.parse_gtf()

    runtime_warnings = [w for w in caught if issubclass(w.category, RuntimeWarning)]
    assert parser._chunk_parse_degradation == 0
    assert runtime_warnings == [], (
        f"clean parse must emit no RuntimeWarning; got {[str(w.message) for w in runtime_warnings]}"
    )


def test_warning_when_chunks_degraded(tmp_path, monkeypatch):
    """When the Cython kernel raises, exactly one ``RuntimeWarning`` is
    emitted and its message includes the degradation count."""
    gtf_path = tmp_path / "degrades.gtf"
    # Write enough rows that the parser produces at least one chunk
    # (chunk_size is 10_000 lines, so we need >=10_000 lines or rely on
    # the final-chunk flush; we use a small-but-valid file and force
    # degradation by monkeypatching the Cython kernel).
    gtf_path.write_text(textwrap.dedent(
        """\
        chr1	src	exon	100	149	.	+	.	gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1	src	CDS	110	140	.	+	0	gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1	src	exon	200	249	.	+	.	gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        """
    ))

    # Force the Cython kernel to raise so the parser falls back to
    # per-row Python parsing and increments _chunk_parse_degradation.
    def fake_raise(chunk_lines):
        raise ValueError("simulated Cython chunk failure")

    monkeypatch.setattr(
        "drvizer.gtf_parser._parse_gtf_chunk_impl",
        fake_raise,
        raising=False,
    )

    parser = GTFParser(str(gtf_path))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        parser.parse_gtf()

    runtime_warnings = [w for w in caught if issubclass(w.category, RuntimeWarning)]
    assert parser._chunk_parse_degradation >= 1, (
        f"_chunk_parse_degradation must be >=1 when Cython kernel raises; "
        f"got {parser._chunk_parse_degradation}"
    )
    assert len(runtime_warnings) == 1, (
        f"exactly one RuntimeWarning must fire when chunks degraded; got "
        f"{len(runtime_warnings)}: {[str(w.message) for w in runtime_warnings]}"
    )
    message = str(runtime_warnings[0].message)
    assert str(parser._chunk_parse_degradation) in message, (
        f"degradation count must appear in warning message; got {message!r}"
    )