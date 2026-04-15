import pytest

from drvizer.api import DrViz
from drvizer.bed_parser import BEDParser
from drvizer._track_build import BEDParser as BUILD_BED_PARSER


class FakePreparedTrack:
    def __init__(self, label, fail_on_prepare=False):
        self.track_label = label
        self.fail_on_prepare = fail_on_prepare
        self.prepare_calls = 0
        self.parser_type = "distribution"
        self.color = "orange"
        self.alpha = 0.8
        self.file_colors = ["orange"]
        self.file_alphas = [0.8]

    def prepare_track(self, gtf_parser):
        self.prepare_calls += 1
        if self.fail_on_prepare:
            raise RuntimeError(f"failed to prepare {self.track_label}")

    def get_grouped_anno_in_region(self, chrom, start, end):
        return {self.track_label: []}


def test_add_bed_track_defers_preparation_until_build(tmp_gtf, tmp_bed, monkeypatch):
    viz = DrViz().load_gtf(str(tmp_gtf))
    observed = []

    def fake_parse_bed(self, chrom=None, start=None, end=None):
        observed.append((self.track_label, chrom, start, end))
        return {}

    monkeypatch.setattr(BUILD_BED_PARSER, "parse_bed", fake_parse_bed)

    viz.add_bed_track(str(tmp_bed), label="first")

    assert observed == []


def test_build_invokes_planned_track_preparation_path(tmp_gtf, tmp_bed, tmp_bed_second, monkeypatch):
    viz = DrViz().load_gtf(str(tmp_gtf))
    prepared_tracks = []

    assert hasattr(BEDParser, "prepare_track")

    def fake_parse_bed(self, chrom=None, start=None, end=None):
        return {}

    def fake_prepare_track(self, gtf_parser):
        prepared_tracks.append((self.track_label, gtf_parser))

    monkeypatch.setattr(BUILD_BED_PARSER, "parse_bed", fake_parse_bed)
    monkeypatch.setattr(BEDParser, "prepare_track", fake_prepare_track)
    monkeypatch.setattr("drvizer._track_build._is_process_safe_genomic_bed", lambda spec: False)

    viz.add_bed_track(str(tmp_bed), label="first")
    viz.add_bed_track(str(tmp_bed_second), label="second")

    parser = viz.build()

    assert [label for label, _ in prepared_tracks] == ["first", "second"]
    assert all(gtf_parser is viz.gtf_parser for _, gtf_parser in prepared_tracks)
    assert [track.track_label for track in parser.data_source.tracks] == ["first", "second"]


def test_build_runs_only_genomic_bed_preparation_in_processes(tmp_gtf, tmp_bed, tmp_bed_second, monkeypatch):
    viz = DrViz().load_gtf(str(tmp_gtf))
    observed_calls = []

    def fake_parallel_prepare(specs, gtf_parser):
        observed_calls.extend((spec["label"], spec["kind"], spec.get("transcript_coord", False)) for spec in specs)
        return [FakePreparedTrack(spec["label"]) for spec in specs]

    monkeypatch.setattr("drvizer.api.prepare_tracks_parallel", fake_parallel_prepare)

    viz.add_bed_track(str(tmp_bed), label="genomic-bed")
    viz.add_bed_track(str(tmp_bed_second), label="tx-bed", transcript_coord=True)

    parser = viz.build()

    assert observed_calls == [
        ("genomic-bed", "bed", False),
        ("tx-bed", "bed", True),
    ]
    assert [track.track_label for track in parser.data_source.tracks] == ["genomic-bed", "tx-bed"]


def test_build_strictly_fails_when_any_track_preparation_fails(tmp_gtf, tmp_bed, tmp_bed_second, monkeypatch):
    viz = DrViz().load_gtf(str(tmp_gtf))
    prepared_labels = set()

    assert hasattr(BEDParser, "prepare_track")

    def fake_parse_bed(self, chrom=None, start=None, end=None):
        return {}

    def fake_prepare_track(self, gtf_parser):
        prepared_labels.add(self.track_label)
        if self.track_label == "second":
            raise RuntimeError("failed to prepare second")

    monkeypatch.setattr(BUILD_BED_PARSER, "parse_bed", fake_parse_bed)
    monkeypatch.setattr(BEDParser, "prepare_track", fake_prepare_track)
    monkeypatch.setattr("drvizer._track_build._is_process_safe_genomic_bed", lambda spec: False)

    viz.add_bed_track(str(tmp_bed), label="first")
    viz.add_bed_track(str(tmp_bed_second), label="second")

    with pytest.raises(RuntimeError, match="failed to prepare second"):
        viz.build()

    assert prepared_labels == {"first", "second"}


def test_build_preserves_registration_order_when_genomic_bed_uses_process_pool(tmp_gtf, tmp_bed, tmp_bed_second, monkeypatch):
    viz = DrViz().load_gtf(str(tmp_gtf))

    def fake_parallel_prepare(specs, gtf_parser):
        return [FakePreparedTrack("second"), FakePreparedTrack("first")]

    monkeypatch.setattr("drvizer.api.prepare_tracks_parallel", fake_parallel_prepare)

    viz.add_bed_track(str(tmp_bed), label="first")
    viz.add_bed_track(str(tmp_bed_second), label="second")

    parser = viz.build()

    assert [track.track_label for track in parser.data_source.tracks] == ["first", "second"]


def test_build_preserves_cython_gtf_indexes_for_transcript_coordinate_tracks(tmp_gtf, tmp_bed, monkeypatch):
    viz = DrViz().load_gtf(str(tmp_gtf))
    prepare_calls = []

    def fake_parse_bed(self, chrom=None, start=None, end=None):
        return {}

    def fake_prepare_track(self, gtf_parser):
        prepare_calls.append(self.track_label)
        assert gtf_parser.identify_gene_identifier_type("TX1") == ("transcript_id", "GENE1")
        assert gtf_parser.get_gene_by_name("G1") == "GENE1"
        assert gtf_parser.transcript_to_gene["TX1"] == "GENE1"

    monkeypatch.setattr(BUILD_BED_PARSER, "parse_bed", fake_parse_bed)
    monkeypatch.setattr(BEDParser, "prepare_track", fake_prepare_track)

    viz.add_bed_track(str(tmp_bed), label="tx-bed", transcript_coord=True)

    parser = viz.build()

    assert prepare_calls == ["tx-bed"]
    assert [track.track_label for track in parser.data_source.tracks] == ["tx-bed"]


def test_build_requires_gtf_to_be_loaded_first(tmp_bed, monkeypatch):
    viz = DrViz()

    def fake_parse_bed(self, chrom=None, start=None, end=None):
        return {}

    monkeypatch.setattr(BUILD_BED_PARSER, "parse_bed", fake_parse_bed)

    viz.add_bed_track(str(tmp_bed), label="first")

    with pytest.raises(ValueError, match=r"GTF file must be loaded first using load_gtf\(\)"):
        viz.build()


def test_build_makes_duplicate_track_labels_unique_in_registration_order(tmp_gtf, tmp_bed, tmp_bed_second, monkeypatch):
    viz = DrViz().load_gtf(str(tmp_gtf))

    def fake_parallel_prepare(specs, gtf_parser):
        return [FakePreparedTrack(spec["label"]) for spec in reversed(specs)]

    monkeypatch.setattr("drvizer.api.prepare_tracks_parallel", fake_parallel_prepare)

    viz.add_bed_track(str(tmp_bed), label="repeat")
    viz.add_bed_track(str(tmp_bed_second), label="repeat")
    viz.add_bed_track(str(tmp_bed), label="repeat")

    parser = viz.build()

    assert [spec["label"] for spec in viz.track_specs] == ["repeat", "repeat.1", "repeat.2"]
    assert [config["label"] for config in viz.track_configs] == ["repeat", "repeat.1", "repeat.2"]
    assert [track.track_label for track in parser.data_source.tracks] == ["repeat", "repeat.1", "repeat.2"]
