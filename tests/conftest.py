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
def transcript_split_gtf(tmp_path):
    path = tmp_path / "transcript_split.gtf"
    path.write_text(textwrap.dedent(
        """\
        chr1	test	exon	100	149	.	+	.	gene_id "gene1"; transcript_id "ENST00000111111"; gene_name "GENE1";
        chr1	test	CDS	110	140	.	+	0	gene_id "gene1"; transcript_id "ENST00000111111"; gene_name "GENE1";
        chr1	test	exon	200	249	.	+	.	gene_id "gene1"; transcript_id "ENST00000111111"; gene_name "GENE1";
        chr1	test	CDS	205	235	.	+	0	gene_id "gene1"; transcript_id "ENST00000111111"; gene_name "GENE1";
        chr1	test	exon	300	339	.	+	.	gene_id "gene1"; transcript_id "ENST00000999999"; gene_name "GENE1";
        chr1	test	CDS	305	330	.	+	0	gene_id "gene1"; transcript_id "ENST00000999999"; gene_name "GENE1";
        chr1	test	exon	360	399	.	+	.	gene_id "gene1"; transcript_id "ENST00000999999"; gene_name "GENE1";
        chr1	test	CDS	365	390	.	+	0	gene_id "gene1"; transcript_id "ENST00000999999"; gene_name "GENE1";
        """
    ))
    return path


@pytest.fixture
def transcript_split_bed_a(tmp_path):
    path = tmp_path / "track_a.bed"
    path.write_text(textwrap.dedent(
        """\
        ENST00000111111	5	20	peakA1	10	+
        ENST00000999999	8	18	peakA2	12	+
        """
    ))
    return path


@pytest.fixture
def transcript_split_bed_b(tmp_path):
    path = tmp_path / "track_b.bed"
    path.write_text(textwrap.dedent(
        """\
        ENST00000111111	30	40	peakB1	5	+
        ENST00000999999	25	35	peakB2	7	+
        """
    ))
    return path
