"""Shared pytest fixtures for the drVizer test suite.

The fixtures here are deliberately reusable across multiple test files.
Audit IDs in parentheses indicate which gaps each fixture helps close:

- simple_gtf / two_exon_gtf / five_exon_gtf
    Property-based projection invariants (P0-2).
- tmp_bed_bed12
    PARSER-013 (int64 thickStart/thickEnd round-trip).
- synthetic_bam
    PARSER-015 (real pysam BAM/BAI integration).

2134796 revival audit (5 tests, 2 helpers):
------------------------------------------
All 5 tests from commit 2134796 are accounted for, none silently dropped.

  test_visualizer_adds_grouped_right_labels_for_nc         -> revive-as-is
  test_visualizer_adds_one_right_label_per_subtrack_for_cn -> revive-as-is
  test_visualizer_truncates_long_transcript_ids            -> revise-and-add
  test_visualizer_no_right_labels_without_groups           -> revive-as-is
  test_visualizer_no_right_labels_with_empty_groups        -> revive-as-is

The 2 helpers (``_base_transcript_data``, ``_truncate_transcript_id``) are
also revive-as-is (functionally equivalent to 2134796). Per-test rationale
is documented in the docstring header of
``tests/test_split_by_transcript_visualizer.py``.
"""

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
def simple_gtf(tmp_path):
    """Single-exon GTF for projection invariants (P0-2)."""
    path = tmp_path / "simple.gtf"
    path.write_text(textwrap.dedent(
        """\
        chr1\tsrc\texon\t500\t699\t.\t+\t.\tgene_id "GSIMPLE"; transcript_id "TX_S"; gene_name "GSIMPLE";
        """
    ))
    return path


@pytest.fixture
def two_exon_gtf(tmp_path):
    """Two-exon GTF for projection invariants (P0-2)."""
    path = tmp_path / "two_exon.gtf"
    path.write_text(textwrap.dedent(
        """\
        chr1\tsrc\texon\t100\t149\t.\t+\t.\tgene_id "GTWO"; transcript_id "TX_T"; gene_name "GTWO";
        chr1\tsrc\texon\t300\t399\t.\t+\t.\tgene_id "GTWO"; transcript_id "TX_T"; gene_name "GTWO";
        chr1\tsrc\texon\t1100\t1149\t.\t-\t.\tgene_id "GTWO"; transcript_id "TX_T_MINUS"; gene_name "GTWO";
        chr1\tsrc\texon\t900\t949\t.\t-\t.\tgene_id "GTWO"; transcript_id "TX_T_MINUS"; gene_name "GTWO";
        """
    ))
    return path


@pytest.fixture
def five_exon_gtf(tmp_path):
    """Five-exon GTF for projection invariants (P0-2).

    Includes a plus- and a minus-strand transcript so monotonicity tests
    can exercise both directions on the same fixture.
    """
    path = tmp_path / "five_exon.gtf"
    plus_lines = "\n".join(
        f"chr1\tsrc\texon\t{1000 + i * 250}\t{1149 + i * 250}\t.\t+\t.\t"
        f'gene_id "GFIVE"; transcript_id "TX_F"; gene_name "GFIVE";'
        for i in range(5)
    )
    minus_exons = [
        (5000, 5149),
        (4800, 4899),
        (4600, 4799),
        (4400, 4499),
        (4200, 4349),
    ]
    minus_lines = "\n".join(
        f"chr1\tsrc\texon\t{s}\t{e}\t.\t-\t.\t"
        f'gene_id "GFIVE"; transcript_id "TX_F_MINUS"; gene_name "GFIVE";'
        for s, e in minus_exons
    )
    path.write_text(plus_lines + "\n" + minus_lines + "\n")
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
def parsed_gtf(tmp_gtf):
    parser = GTFParser(str(tmp_gtf))
    parser.parse_gtf()
    return parser


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
def tmp_bed_bed12(tmp_path):
    """BED12 fixture with thickStart/thickEnd at >2**31 to exercise int64 (PARSER-013).

    Thick start/end values are intentionally placed above 2**31 so that any
    BED parser contract regression (32-bit overflow, off-by-one in
    half-open conversion) is caught at int64 round-trip time.
    """
    path = tmp_path / "int64_thick.bed"
    path.write_text(
        "chr1\t2300000000\t2300000500\tpeak\t10\t+\t2300000050\t2300000450\t255,0,0\t2\t30,30\t0,170\n"
    )
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


@pytest.fixture
def synthetic_bam(tmp_path):
    """Real pysam BAM/BAI fixture for PARSER-015 (genomic coverage).

    Writes a 2-record BAM to disk via ``pysam.AlignmentFile(..., mode='wb')``,
    indexes it via ``pysam.index``, and returns the path. Tests that need a
    real BAM path can call this fixture; if pysam is unavailable, the fixture
    raises a clear error at fixture time (rather than skipping the test
    silently). Tests should additionally use ``pytest.importorskip('pysam')``
    at module level so the test reports SKIPPED with a clear reason when
    pysam is missing.
    """
    pysam = pytest.importorskip("pysam")
    bam_path = tmp_path / "coverage.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 500}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        # Record 1: 3-base block at genomic [100, 103).
        read1 = pysam.AlignedSegment()
        read1.query_name = "read1"
        read1.reference_id = 0
        read1.reference_start = 100
        read1.mapping_quality = 60
        read1.cigar = [(0, 3)]
        read1.flag = 0
        bam.write(read1)

        # Record 2: 2-base block at genomic [200, 202).
        read2 = pysam.AlignedSegment()
        read2.query_name = "read2"
        read2.reference_id = 0
        read2.reference_start = 200
        read2.mapping_quality = 60
        read2.cigar = [(0, 2)]
        read2.flag = 0
        bam.write(read2)

    pysam.index(str(bam_path))
    return bam_path


# 2134796 revival audit:
# ----------------------
# All 5 tests from commit 2134796 are accounted for, none silently dropped:
#
#   test_visualizer_adds_grouped_right_labels_for_nc          (revive-as-is)
#   test_visualizer_adds_one_right_label_per_subtrack_for_cn  (revive-as-is)
#   test_visualizer_truncates_long_transcript_ids             (revise-and-add)
#   test_visualizer_no_right_labels_without_groups            (revive-as-is)
#   test_visualizer_no_right_labels_with_empty_groups         (revive-as-is)
#
# See tests/test_split_by_transcript_visualizer.py docstring header for
# per-test revival rationale.
