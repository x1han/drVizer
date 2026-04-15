import builtins
from pathlib import Path

import pytest

from drvizer.bed_parser import BEDParser


def test_parse_bed_uses_cython_fast_path_when_available(tmp_bed, monkeypatch):
    calls = []

    def fake_parse(bed_file_paths, region):
        calls.append((tuple(bed_file_paths), region))
        return {
            "chr1": [
                {
                    "chrom": "chr1",
                    "start": 105,
                    "end": 120,
                    "name": "peak1",
                    "score": 10.0,
                    "strand": "+",
                    "thickStart": 105,
                    "thickEnd": 120,
                    "itemRgb": "0",
                }
            ]
        }

    monkeypatch.setattr("drvizer.bed_parser._CYTHON_BED_AVAILABLE", True, raising=False)
    monkeypatch.setattr("drvizer.bed_parser.parse_bed_records", fake_parse, raising=False)

    parser = BEDParser(str(tmp_bed))
    result = parser.parse_bed()

    assert calls == [((str(tmp_bed),), None)]
    assert result["chr1"][0]["name"] == "peak1"
    assert parser._parsed is True


def test_parse_bed_falls_back_when_cython_extension_unavailable(tmp_bed, monkeypatch):
    calls = []

    def fake_python_parse(bed_file_paths, region=None):
        calls.append((tuple(bed_file_paths), region))
        return {
            "chr1": [
                {
                    "chrom": "chr1",
                    "start": 105,
                    "end": 120,
                    "name": "peak1",
                    "score": 10.0,
                    "strand": "+",
                    "thickStart": 105,
                    "thickEnd": 120,
                    "itemRgb": "0",
                }
            ]
        }

    monkeypatch.setattr("drvizer.bed_parser._CYTHON_BED_AVAILABLE", False)
    monkeypatch.setattr("drvizer.bed_parser.parse_bed_records_python", fake_python_parse)

    parser = BEDParser(str(tmp_bed))
    result = parser.parse_bed()

    assert calls == [((str(tmp_bed),), None)]
    assert result["chr1"][0]["name"] == "peak1"
    assert result["chr1"][0]["score"] == 10.0


def test_parse_bed_cython_fast_path_matches_python_region_filter(tmp_bed, tmp_bed_second, monkeypatch):
    parser = BEDParser([str(tmp_bed), str(tmp_bed_second)])
    expected = parser.parse_bed(chrom="chr1", start=200, end=330)

    monkeypatch.setattr("drvizer.bed_parser._CYTHON_BED_AVAILABLE", True, raising=False)
    parser = BEDParser([str(tmp_bed), str(tmp_bed_second)])
    actual = parser.parse_bed(chrom="chr1", start=200, end=330)

    assert actual == expected


def test_transcript_coordinate_bed_uses_python_parser_even_when_cython_available(parsed_gtf, tmp_path, monkeypatch):
    bed_path = tmp_path / "transcript_coords.bed"
    bed_path.write_text("TX1\t10\t70\tpeak_tx\t7\t+\n")

    def unexpected_fast_path(*args, **kwargs):
        raise AssertionError("Cython genomic parser should not be used for transcript-coordinate BED")

    monkeypatch.setattr("drvizer.bed_parser._CYTHON_BED_AVAILABLE", True, raising=False)
    monkeypatch.setattr("drvizer.bed_parser.parse_bed_records", unexpected_fast_path, raising=False)

    parser = BEDParser(str(bed_path), transcript_coord=True, gtf_parser=parsed_gtf)
    result = parser.parse_bed()

    assert result == {
        "chr1": [
            {
                "chrom": "chr1",
                "start": 110,
                "end": 150,
                "name": "peak_tx",
                "score": 7.0,
                "strand": "+",
                "thickStart": 10,
                "thickEnd": 70,
                "itemRgb": "0",
            },
            {
                "chrom": "chr1",
                "start": 200,
                "end": 220,
                "name": "peak_tx",
                "score": 7.0,
                "strand": "+",
                "thickStart": 10,
                "thickEnd": 70,
                "itemRgb": "0",
            },
        ]
    }


def test_cython_and_python_paths_report_same_missing_file_error_message(monkeypatch):
    missing_path = "/tmp/definitely-missing-parser-acceleration.bed"

    monkeypatch.setattr("drvizer.bed_parser._CYTHON_BED_AVAILABLE", False)
    parser = BEDParser(missing_path)
    with pytest.raises(ValueError, match=r"^Error reading BED file /tmp/definitely-missing-parser-acceleration\.bed: .*No such file or directory"):
        parser.parse_bed()

    monkeypatch.setattr("drvizer.bed_parser._CYTHON_BED_AVAILABLE", True)
    parser = BEDParser(missing_path)
    with pytest.raises(ValueError, match=r"^Error reading BED file /tmp/definitely-missing-parser-acceleration\.bed: .*No such file or directory"):
        parser.parse_bed()


def test_cython_and_python_paths_report_same_gzip_decode_error_message(tmp_path, monkeypatch):
    bad_gzip_path = tmp_path / "bad.bed.gz"
    bad_gzip_path.write_bytes(b"not-a-gzip-stream")

    monkeypatch.setattr("drvizer.bed_parser._CYTHON_BED_AVAILABLE", False)
    parser = BEDParser(str(bad_gzip_path))
    with pytest.raises(ValueError, match=rf"^Error reading BED file {Path(str(bad_gzip_path))}: .*Not a gzipped file"):
        parser.parse_bed()

    monkeypatch.setattr("drvizer.bed_parser._CYTHON_BED_AVAILABLE", True)
    parser = BEDParser(str(bad_gzip_path))
    with pytest.raises(ValueError, match=rf"^Error reading BED file {Path(str(bad_gzip_path))}: .*Not a gzipped file"):
        parser.parse_bed()




def test_parse_bed_region_filter_accepts_zero_start(tmp_path, monkeypatch):
    bed_path = tmp_path / "zero_start.bed"
    bed_path.write_text("chr1\t0\t10\tzero\t1\t+\nchr1\t20\t30\tlater\t2\t+\n")

    monkeypatch.setattr("drvizer.bed_parser._CYTHON_BED_AVAILABLE", False)

    parser = BEDParser(str(bed_path))
    result = parser.parse_bed(chrom="chr1", start=0, end=15)

    assert result == {
        "chr1": [
            {
                "chrom": "chr1",
                "start": 0,
                "end": 10,
                "name": "zero",
                "score": 1.0,
                "strand": "+",
                "thickStart": 0,
                "thickEnd": 10,
                "itemRgb": "0",
            }
        ]
    }
