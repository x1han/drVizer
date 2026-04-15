from drvizer import gtf_parser as gtf_parser_module
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


def test_parse_gtf_uses_cython_chunk_rows_without_changing_indexes(tmp_gtf, monkeypatch):
    parser = GTFParser(str(tmp_gtf))
    observed_chunks = []

    def fake_parse_gtf_chunk(chunk_lines):
        observed_chunks.append(list(chunk_lines))
        return [
            ("chr1", "exon", 100, 149, "+", 'gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";'),
            ("chr1", "CDS", 110, 140, "+", 'gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";'),
            ("chr1", "exon", 200, 249, "+", 'gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";'),
            ("chr1", "exon", 300, 349, "+", 'gene_id "GENE1"; transcript_id "TX2"; gene_name "G1";'),
        ]

    monkeypatch.setattr(gtf_parser_module, "_parse_gtf_chunk_impl", fake_parse_gtf_chunk)

    parser.parse_gtf()

    assert len(observed_chunks) == 1
    assert parser.gene_info["GENE1"]["transcripts"] == ["TX1", "TX2"]
    assert parser.transcript_to_gene == {"TX1": "GENE1", "TX2": "GENE1"}
    assert parser.gene_transcripts["GENE1"]["TX1"]["exons"] == [
        {"start": 100, "end": 149, "number": 1},
        {"start": 200, "end": 249, "number": 2},
    ]
    assert parser.gene_transcripts["GENE1"]["TX1"]["cds"] == [{"start": 110, "end": 140}]
