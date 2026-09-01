from drvizer import gtf_parser as gtf_parser_module
from drvizer.gtf_parser import GTFParser


# Convention note: the projection path emits 0-based half-open genomic
# intervals (BED standard). The fixture tmp_gtf declares exons
# [100, 149] and [200, 249] in 1-based inclusive GTF coordinates, so the
# expected 0-based half-open halves are [99, 149] and [199, 249]. Width
# is preserved across the convention switch.

def test_segment_projection_splits_across_exons(tmp_gtf):
    parser = GTFParser(str(tmp_gtf))
    parser.parse_gtf()

    chrom, strand, segments = parser.convert_transcript_to_genomic_segments("TX1", 25, 75)

    assert chrom == "chr1"
    assert strand == "+"
    assert segments == [(124, 149), (199, 224)]


def test_legacy_projection_wrapper_returns_outer_bounds(tmp_gtf):
    parser = GTFParser(str(tmp_gtf))
    parser.parse_gtf()

    result = parser.convert_transcript_to_genomic("TX1", 25, 75)

    assert result == ("chr1", 124, 224)


def test_find_transcript_returns_cached_mapping(parsed_gtf):
    gene_id, transcript = parsed_gtf.find_transcript("TX1")

    assert gene_id == "GENE1"
    assert transcript["seqname"] == "chr1"


def test_projection_exons_cache_is_populated(parsed_gtf):
    parsed_gtf.convert_transcript_to_genomic_segments("TX1", 25, 75)

    transcript = parsed_gtf.gene_transcripts["GENE1"]["TX1"]
    assert "_projection_exons" in transcript
    assert len(transcript["_projection_exons"]) == 2


def test_64bit_and_zero_based_parity(tmp_path, monkeypatch):
    """Lock int64 promotion + 0-based half-open + Python/Cython parity.

    Fixture GTF exons:
      exon1: GTF 1-based inclusive [2300000000, 2300000049] -> 0-based half-open [2300000000, 2300000050)
      exon2: GTF 1-based inclusive [2300000100, 2300000149] -> 0-based half-open [2300000100, 2300000150)

    Transcript coordinates:
      exon1 spans transcript [0, 50)
      exon2 spans transcript [50, 100)

    Case A (single exon): convert_transcript_to_genomic_segments('TX1', 10, 30)
      only overlaps transcript [10, 30) within exon1.
      Expected (0-based half-open): [(2299999999 + 10, 2299999999 + 30)]
                                  = [(2300000009, 2300000039)]

    Case B (junctions): convert_transcript_to_genomic_segments('TX1', 40, 60)
      hits the exon1/exon2 seam, transcript [40, 50) in exon1 and [50, 60) in exon2.
      Expected:
        exon1: (2299999999 + 40, 2299999999 + 50) = (2300000039, 2300000049)
        exon2: (2299999999 + 50, 2299999999 + 60) = (2300000049, 2300000059)
      Sorted: [(2300000039, 2300000049), (2300000049, 2300000059)]

    The test runs both arms (Cython fast path + Python fallback) via
    monkeypatch and asserts identical results.
    """
    import textwrap

    gtf = tmp_path / "int64.gtf"
    gtf.write_text(textwrap.dedent(
        """\
        chr1\tsrc\texon\t2300000000\t2300000049\t.\t+\t.\tgene_id "GENE_INT64"; transcript_id "TX1"; gene_name "GINT64";
        chr1\tsrc\texon\t2300000100\t2300000149\t.\t+\t.\tgene_id "GENE_INT64"; transcript_id "TX1"; gene_name "GINT64";
        """
    ))

    def run_case(transcript_start, transcript_end):
        # First arm: whatever is currently loaded (Cython if built, else Python).
        parser = GTFParser(str(gtf))
        parser.parse_gtf()
        fast_result = parser.convert_transcript_to_genomic_segments("TX1", transcript_start, transcript_end)

        # Second arm: force the Python fallback by nulling the module binding.
        fallback = GTFParser(str(gtf))
        fallback.parse_gtf()
        monkeypatch.setattr(gtf_parser_module, "_project_segments_fast", None)
        python_result = fallback.convert_transcript_to_genomic_segments("TX1", transcript_start, transcript_end)

        return fast_result, python_result

    fast_a, python_a = run_case(10, 30)
    assert fast_a == python_a
    assert fast_a == ("chr1", "+", [(2300000009, 2300000039)])

    fast_b, python_b = run_case(40, 60)
    assert fast_b == python_b
    assert fast_b == ("chr1", "+", [(2300000039, 2300000049), (2300000049, 2300000059)])
