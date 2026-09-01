"""Phase 3 -- P1-6: BAM projection kernel consistency tests.

Two tests, each anchored to the P1-6 mandate:

* test_single_transcript_equals_batch
    Identical transcript-coverage arrays between the single-transcript
    path (``_get_coverage_for_transcript``) and the batch path
    (``get_coverage_for_transcripts``) on the same synthetic BAM-shaped
    input.
* test_kernel_pure_numeric
    ``_project_aligned_blocks_to_transcript`` accumulates into the
    coverage array with no I/O when given transcript blocks, the GTF
    parser, target_chrom, and the genomic region window.
"""

import numpy as np
import pytest

from drvizer import bam_parser


class _FakeRead:
    def __init__(self, blocks):
        self._blocks = list(blocks)

    def get_blocks(self):
        return list(self._blocks)


class _FakeAlignmentFile:
    """Minimal pysam.AlignmentFile stub for transcript-coord coverage tests.

    The header advertises a single transcript ``TX1`` of length 200; reads
    returned by ``fetch`` are supplied by the test.
    """

    header = {"SQ": [{"SN": "TX1", "LN": 200}]}

    def __init__(self, path, mode):
        self._path = path
        self._reads = _FakeAlignmentFile._READS_BY_PATH.get(path, [])

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def fetch(self, transcript_id, start, end):
        return list(self._reads)


class _FakeGTFParser:
    """Single-exon transcript: chr1 [99, 110] half-open."""

    def __init__(self, exon_start=100, exon_end=110):
        # 0-based half-open: exon GTF start 100 maps to genomic 99.
        self._start = exon_start - 1
        self._end = exon_end
        self._exon_start = exon_start  # GTF 1-based inclusive
        self._exon_end = exon_end      # GTF 1-based inclusive

    def convert_transcript_to_genomic_segments(self, transcript_id, t_start, t_end):
        s = max(self._start, self._start + t_start)
        e = min(self._end, self._start + t_end)
        if s >= e:
            return "chr1", "+", []
        return "chr1", "+", [(s, e)]

    def get_transcript_data(self, gene_identifier):
        return {
            "seqname": "chr1",
            "transcripts": [
                {
                    "transcript_id": gene_identifier,
                    "exons": [{"start": self._exon_start, "end": self._exon_end}],
                    "cds": [],
                    "strand": "+",
                }
            ],
        }


class _ReadPool(list):
    pass


_FakeAlignmentFile._READS_BY_PATH = {}


def test_kernel_pure_numeric():
    """``_project_aligned_blocks_to_transcript`` accumulates into a
    pre-allocated coverage array with no I/O."""
    gtf_parser = _FakeGTFParser(exon_start=100, exon_end=110)
    # Region window [99, 110) -> region_len = 11.
    region_start = 99
    region_len = 11
    coverage = np.zeros(region_len, dtype=np.int64)

    # Transcript block (0, 11) maps to genomic (99, 110) -> indices (0, 11).
    blocks = [(0, 11)]
    bam_parser._project_aligned_blocks_to_transcript(
        transcript_id="TX1",
        blocks=blocks,
        gtf_parser=gtf_parser,
        target_chrom="chr1",
        region_start=region_start,
        region_len=region_len,
        coverage=coverage,
    )
    # Each base in [99, 110) gets one read.
    assert int(coverage.sum()) == 11
    assert int(coverage[0]) == 1
    assert int(coverage[-1]) == 1

    # Mismatched chrom filter: target_chrom != projection chrom skips.
    coverage[:] = 0
    bam_parser._project_aligned_blocks_to_transcript(
        transcript_id="TX1",
        blocks=blocks,
        gtf_parser=gtf_parser,
        target_chrom="chrX",
        region_start=region_start,
        region_len=region_len,
        coverage=coverage,
    )
    assert int(coverage.sum()) == 0, "non-target chrom should accumulate nothing"


def test_single_transcript_equals_batch(monkeypatch, tmp_path):
    """Single-transcript path and batch path must produce identical arrays.

    We point the two paths at the same BAM-shaped fake file. Reads cover
    transcript positions (0, 5) and (5, 10), so the cumulative genomic
    coverage [99, 109] should be exactly [1, 1, 1, 1, 1, 1, 1, 1, 1, 1].
    """
    bam_path = tmp_path / "fake.bam"
    bam_path.write_text("")
    # Empty .bai file so BAMParser.__init__ doesn't try to index.
    (tmp_path / "fake.bam.bai").write_text("")
    _FakeAlignmentFile._READS_BY_PATH = {
        str(bam_path): [_FakeRead([(0, 5)]), _FakeRead([(5, 10)])]
    }

    monkeypatch.setattr(
        "drvizer.bam_parser.pysam.AlignmentFile",
        _FakeAlignmentFile,
    )

    from drvizer.bam_parser import BAMParser

    parser = BAMParser(
        str(bam_path),
        transcript_coord=True,
        gtf_parser=_FakeGTFParser(exon_start=100, exon_end=110),
    )

    transcript_info = {
        "transcript_id": "TX1",
        "exons": [{"start": 100, "end": 110}],
    }
    target_chrom = "chr1"

    # Single-transcript path
    x_single, y_single, _series = parser._get_coverage_for_transcript(transcript_info, target_chrom)
    # Batch path
    x_batch, y_batch = parser.get_coverage_for_transcripts("TX1")

    assert np.array_equal(x_single, x_batch), (
        f"x arrays differ: single={x_single} batch={x_batch}"
    )
    assert np.array_equal(y_single, y_batch), (
        f"y arrays differ: single={y_single} batch={y_batch}"
    )
    assert int(y_single.sum()) == 10
