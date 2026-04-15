import numpy as np
import pytest

from drvizer.bam_parser import BAMParser
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
