from pathlib import Path
import sys
import textwrap

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from drvizer.gtf_parser import GTFParser


@pytest.fixture
def tmp_gtf(tmp_path):
    path = tmp_path / "test.gtf"
    path.write_text(textwrap.dedent(
        """\
        chr1	src	exon	100	149	.	+	.	gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1	src	CDS	110	140	.	+	0	gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1	src	exon	200	249	.	+	.	gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1	src	exon	300	349	.	+	.	gene_id "GENE1"; transcript_id "TX2"; gene_name "G1";
        """
    ))
    return path


@pytest.fixture
def tmp_gtf_second(tmp_path):
    path = tmp_path / "test_second.gtf"
    path.write_text(textwrap.dedent(
        """\
        chr1	src	exon	100	149	.	+	.	gene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1	src	exon	400	449	.	+	.	gene_id "GENE1"; transcript_id "TX3"; gene_name "G1";
        """
    ))
    return path


@pytest.fixture
def tmp_bed(tmp_path):
    path = tmp_path / "test.bed"
    path.write_text(textwrap.dedent(
        """\
        chr1	105	120	peak1	10	+
        chr1	205	220	peak2	5	+
        """
    ))
    return path


@pytest.fixture
def tmp_bed_second(tmp_path):
    path = tmp_path / "test_second.bed"
    path.write_text(textwrap.dedent(
        """\
        chr1	305	320	peak3	8	+
        chr1	405	420	peak4	4	+
        """
    ))
    return path


@pytest.fixture
def parsed_gtf(tmp_gtf):
    parser = GTFParser(str(tmp_gtf))
    parser.parse_gtf()
    return parser
