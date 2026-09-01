import numpy as np
import pytest

from drvizer.bam_parser import BAMParser
from drvizer import _parallel
from drvizer._parallel import ParallelCoverageError


class FakeRead:
    def __init__(self, start, end):
        self.reference_start = start
        self.reference_end = end
        self._blocks = [(start, end)]

    def get_blocks(self):
        return list(self._blocks)


class FakeAlignmentFile:
    def __init__(self, reads_by_path, path, mode):
        self._reads = reads_by_path[path]
        self.path = path
        self.mode = mode

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def fetch(self, chrom, start, end):
        return [
            read for read in self._reads
            if read.reference_start < end and read.reference_end > start
        ]


@pytest.fixture
def fake_bam_paths(tmp_path):
    paths = []
    for name in ("a.bam", "b.bam"):
        path = tmp_path / name
        path.write_text("")
        (tmp_path / f"{name}.bai").write_text("")
        paths.append(str(path))
    return paths


def test_get_coverage_in_region_serially_aggregates_multiple_bams(monkeypatch, fake_bam_paths):
    reads_by_path = {
        fake_bam_paths[0]: [FakeRead(10, 13)],
        fake_bam_paths[1]: [FakeRead(11, 14)],
    }

    monkeypatch.setattr(
        "drvizer.bam_parser.pysam.AlignmentFile",
        lambda path, mode: FakeAlignmentFile(reads_by_path, path, mode),
    )

    parser = BAMParser(fake_bam_paths)

    x, coverage = parser.get_coverage_in_region("chr1", 10, 14)

    assert np.array_equal(x, np.array([10, 11, 12, 13]))
    assert np.array_equal(coverage, np.array([1, 2, 2, 1]))


def test_get_coverage_in_region_returns_mean_for_multiple_bams(monkeypatch, fake_bam_paths):
    reads_by_path = {
        fake_bam_paths[0]: [FakeRead(10, 13)],
        fake_bam_paths[1]: [FakeRead(11, 14)],
    }

    monkeypatch.setattr(
        "drvizer.bam_parser.pysam.AlignmentFile",
        lambda path, mode: FakeAlignmentFile(reads_by_path, path, mode),
    )

    parser = BAMParser(fake_bam_paths, aggregate_method="mean")

    _, coverage = parser.get_coverage_in_region("chr1", 10, 14)

    assert np.array_equal(coverage, np.array([0.5, 1.0, 1.0, 0.5]))


def test_get_coverage_in_region_uses_parallel_only_for_multiple_bams(monkeypatch, fake_bam_paths):
    reads_by_path = {
        fake_bam_paths[0]: [FakeRead(10, 13)],
        fake_bam_paths[1]: [FakeRead(11, 14)],
    }
    parallel_calls = []

    monkeypatch.setattr(
        "drvizer.bam_parser.pysam.AlignmentFile",
        lambda path, mode: FakeAlignmentFile(reads_by_path, path, mode),
    )

    def fake_parallel(paths, chrom, start, end, contained_only):
        parallel_calls.append((tuple(paths), chrom, start, end, contained_only))
        return np.array([1, 2, 2, 1], dtype=np.int32)

    monkeypatch.setattr(
        "drvizer.bam_parser.aggregate_region_coverages_parallel",
        fake_parallel,
    )

    multi_parser = BAMParser(fake_bam_paths)
    _, multi_coverage = multi_parser.get_coverage_in_region("chr1", 10, 14)

    single_parser = BAMParser(fake_bam_paths[0])
    _, single_coverage = single_parser.get_coverage_in_region("chr1", 10, 14)

    assert parallel_calls == [((fake_bam_paths[0], fake_bam_paths[1]), "chr1", 10, 14, True)]
    assert np.array_equal(multi_coverage, np.array([1, 2, 2, 1]))
    assert np.array_equal(single_coverage, np.array([1, 1, 1, 0]))


def test_parallel_coverage_caps_worker_count(monkeypatch):
    worker_counts = []

    class FakePool:
        def __init__(self, processes):
            worker_counts.append(processes)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def imap_unordered(self, function, args):
            # Phase 8 tuple-return contract: each result is
            # (bam_path, coverage_array). Derive bam_path from
            # args[i][0] so attribution is exercised end-to-end.
            return [(args[i][0], np.zeros(args[i][3] - args[i][2], dtype=np.int64)) for i in range(len(args))]

    monkeypatch.setattr(_parallel.multiprocessing, "Pool", FakePool)
    monkeypatch.setattr(_parallel.multiprocessing, "cpu_count", lambda: 4)

    _parallel.aggregate_region_coverages_parallel([f"sample_{i}.bam" for i in range(10)], "chr1", 0, 5)

    assert worker_counts == [4]


def test_parallel_coverage_uses_int64_accumulator(monkeypatch):
    class FakePool:
        def __init__(self, processes):
            pass

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def imap_unordered(self, function, args):
            # Phase 8 tuple-return contract: (bam_path, coverage_array).
            return [
                (args[0][0], np.array([np.iinfo(np.int32).max], dtype=np.int32)),
                (args[1][0], np.array([1], dtype=np.int32)),
            ]

    monkeypatch.setattr(_parallel.multiprocessing, "Pool", FakePool)

    coverage = _parallel.aggregate_region_coverages_parallel(["a.bam", "b.bam"], "chr1", 0, 1)

    assert coverage.dtype == np.int64
    assert coverage[0] == np.iinfo(np.int32).max + 1


def test_transcript_coverage_uses_int64_accumulator(monkeypatch, fake_bam_paths):
    class FakeHeaderAlignmentFile:
        header = {"SQ": [{"SN": "TX1", "LN": 2}]}

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def fetch(self, transcript_id, start, end):
            return [FakeRead(0, 1), FakeRead(1, 2)]

    class FakeGTFParser:
        def convert_transcript_to_genomic_segments(self, transcript_id, start, end):
            # 0-based half-open: subtract 1 from the GTF 'start' to form the
            # origin (matches the audit-fixed projection convention).
            return "chr1", "+", [(99 + start, 99 + end)]

    monkeypatch.setattr(
        "drvizer.bam_parser.pysam.AlignmentFile",
        lambda path, mode: FakeHeaderAlignmentFile(),
    )

    parser = BAMParser(fake_bam_paths[0], transcript_coord=True, gtf_parser=FakeGTFParser())
    _, coverage = parser._get_coverage_for_transcript_from_path(
        fake_bam_paths[0],
        {"transcript_id": "TX1", "exons": [{"start": 100, "end": 101}]},
        "chr1",
    )

    assert coverage.dtype == np.int64
    assert np.array_equal(coverage, np.array([1, 1], dtype=np.int64))


def test_parallel_coverage_defaults_cpu_count_when_unavailable(monkeypatch):
    worker_counts = []

    class FakePool:
        def __init__(self, processes):
            worker_counts.append(processes)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def imap_unordered(self, function, args):
            # Phase 8 tuple-return contract: (bam_path, coverage_array).
            return [(args[i][0], np.zeros(args[i][3] - args[i][2], dtype=np.int64)) for i in range(len(args))]

    monkeypatch.setattr(_parallel.multiprocessing, "Pool", FakePool)
    monkeypatch.setattr(
        _parallel.multiprocessing,
        "cpu_count",
        lambda: (_ for _ in ()).throw(NotImplementedError("unavailable")),
    )

    coverage = _parallel.aggregate_region_coverages_parallel(["a.bam", "b.bam"], "chr1", 0, 2)

    assert worker_counts == [1]
    assert coverage.dtype == np.int64


def test_get_coverage_in_region_serial_uses_int64_accumulator(monkeypatch, fake_bam_paths):
    coverages = iter([
        np.array([np.iinfo(np.int32).max], dtype=np.int32),
        np.array([1], dtype=np.int32),
    ])

    monkeypatch.setattr("drvizer.bam_parser.compute_region_coverage", lambda *args: next(coverages))

    parser = BAMParser(fake_bam_paths)
    coverage = parser._get_coverage_in_region_serial("chr1", 0, 1)

    assert coverage.dtype == np.int64
    assert coverage[0] == np.iinfo(np.int32).max + 1


def test_split_transcript_coverage_uses_int64_sum_accumulator(fake_bam_paths):
    parser = BAMParser(fake_bam_paths, transcript_coord=True, gtf_parser=object())
    coverages = iter([
        np.array([np.iinfo(np.int32).max], dtype=np.int32),
        np.array([1], dtype=np.int32),
    ])
    x = np.array([100])
    parser._get_coverage_for_transcript_from_path = lambda *args: (x, next(coverages))

    _, coverage, series = parser._get_coverage_for_transcript(
        {"transcript_id": "TX1", "exons": [{"start": 100, "end": 100}]},
        "chr1",
    )

    assert coverage.dtype == np.int64
    assert coverage[0] == np.iinfo(np.int32).max + 1
    assert len(series) == 2


def test_get_coverage_in_region_falls_back_to_serial_when_parallel_fails(monkeypatch, fake_bam_paths):
    reads_by_path = {
        fake_bam_paths[0]: [FakeRead(10, 13)],
        fake_bam_paths[1]: [FakeRead(11, 14)],
    }

    monkeypatch.setattr(
        "drvizer.bam_parser.pysam.AlignmentFile",
        lambda path, mode: FakeAlignmentFile(reads_by_path, path, mode),
    )
    monkeypatch.setattr(
        "drvizer.bam_parser.aggregate_region_coverages_parallel",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            ParallelCoverageError("parallel unavailable")
        ),
    )

    parser = BAMParser(fake_bam_paths)

    _, coverage = parser.get_coverage_in_region("chr1", 10, 14)

    assert np.array_equal(coverage, np.array([1, 2, 2, 1]))
