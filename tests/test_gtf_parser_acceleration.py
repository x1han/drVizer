from drvizer.gtf_parser import GTFParser


def test_parse_gtf_builds_gene_and_transcript_indexes(tmp_gtf):
    parser = GTFParser(str(tmp_gtf))
    parser.parse_gtf()

    assert parser.get_gene_by_name("G1") == "GENE1"
    assert parser.transcript_to_gene["TX1"] == "GENE1"
    assert parser.identify_gene_identifier_type("TX1") == ("transcript_id", "GENE1")


def test_multi_gtf_deduplicates_transcripts_from_earlier_files(tmp_gtf, tmp_gtf_second):
    parser = GTFParser([str(tmp_gtf), str(tmp_gtf_second)])
    parser.parse_gtf()

    transcript_ids = set(parser.gene_transcripts["GENE1"].keys())
    assert transcript_ids == {"TX1", "TX2", "TX3"}


def test_get_transcript_data_keeps_expected_output_shape(tmp_gtf):
    parser = GTFParser(str(tmp_gtf))
    parser.parse_gtf()

    data = parser.get_transcript_data("GENE1")

    assert data["gene_id"] == "GENE1"
    assert data["seqname"] == "chr1"
    assert [tx["transcript_id"] for tx in data["transcripts"]] == ["TX1", "TX2"]
