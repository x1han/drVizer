from drvizer.gtf_parser import GTFParser


def test_segment_projection_splits_across_exons(tmp_gtf):
    parser = GTFParser(str(tmp_gtf))
    parser.parse_gtf()

    chrom, strand, segments = parser.convert_transcript_to_genomic_segments("TX1", 25, 75)

    assert chrom == "chr1"
    assert strand == "+"
    assert segments == [(125, 150), (200, 225)]


def test_legacy_projection_wrapper_returns_outer_bounds(tmp_gtf):
    parser = GTFParser(str(tmp_gtf))
    parser.parse_gtf()

    result = parser.convert_transcript_to_genomic("TX1", 25, 75)

    assert result == ("chr1", 125, 225)


def test_find_transcript_returns_cached_mapping(parsed_gtf):
    gene_id, transcript = parsed_gtf.find_transcript("TX1")

    assert gene_id == "GENE1"
    assert transcript["seqname"] == "chr1"


def test_projection_exons_cache_is_populated(parsed_gtf):
    parsed_gtf.convert_transcript_to_genomic_segments("TX1", 25, 75)

    transcript = parsed_gtf.gene_transcripts["GENE1"]["TX1"]
    assert "_projection_exons" in transcript
    assert len(transcript["_projection_exons"]) == 2
