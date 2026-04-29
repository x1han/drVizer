import numpy as np
import pytest

import drvizer.api as api
import drvizer._track_build as track_build
from drvizer.api import DrViz
from drvizer.bam_parser import BAMParser
from drvizer.bed_parser import BEDParser


@pytest.fixture(autouse=True)
def stub_bam_parser(monkeypatch):
    class DummyBamParser:
        def __init__(self, bam_paths, track_label='BAM Coverage', contained_only=True,
                     color='steelblue', y_axis_range=None, aggregate_method='sum',
                     transcript_coord=False, gtf_parser=None):
            self.bam_paths = [bam_paths] if isinstance(bam_paths, str) else list(bam_paths)
            self.track_label = track_label
            self.contained_only = contained_only
            self.color = color
            self.alpha = 0.6
            self.y_axis_range = y_axis_range
            self.aggregate_method = aggregate_method
            self.parser_type = 'coverage'
            self.transcript_coord = transcript_coord
            self.gtf_parser = gtf_parser

        def get_coverage_in_region(self, chrom, start, end, target_bins=2000):
            return np.array([start, min(start + 1, end)]), np.array([1, 2])

        def get_coverage_for_transcripts(self, gene_identifier, target_bins=2000):
            return np.array([100, 101]), np.array([7, 9])

        def get_coverage_by_transcript(self, gene_identifier):
            return {
                'ENST00000111111': {
                    'x': np.array([100, 101]),
                    'y': np.array([2, 3]),
                },
                'ENST00000999999': {
                    'x': np.array([300, 301]),
                    'y': np.array([5, 1]),
                },
            }

    monkeypatch.setattr(api, 'BAMParser', DummyBamParser)
    monkeypatch.setattr(track_build, 'BAMParser', DummyBamParser)


@pytest.mark.parametrize(
    "bad_value",
    ["", "tx", "NN"],
    ids=["empty", "unknown_short", "unknown_caps"],
)
def test_add_bed_track_rejects_invalid_split_by_transcript(transcript_split_gtf, transcript_split_bed_a, bad_value):
    viz = DrViz().load_gtf(str(transcript_split_gtf))

    with pytest.raises(ValueError, match="must be one of None, 'nc', or 'cn'"):
        viz.add_bed_track(
            str(transcript_split_bed_a),
            label="A",
            transcript_coord=True,
            split_by_transcript=bad_value,
        )


def test_add_bam_track_rejects_invalid_split_by_transcript(transcript_split_gtf):
    viz = DrViz().load_gtf(str(transcript_split_gtf))

    with pytest.raises(ValueError, match="must be one of None, 'nc', or 'cn'"):
        viz.add_bam_track(
            "fake.bam",
            label="Coverage",
            transcript_coord=True,
            split_by_transcript="bad",
        )


def test_add_bam_track_requires_transcript_coord_when_split_is_enabled(transcript_split_gtf):
    viz = DrViz().load_gtf(str(transcript_split_gtf))

    with pytest.raises(ValueError, match="split_by_transcript requires transcript_coord=True"):
        viz.add_bam_track(
            "fake.bam",
            label="Coverage",
            transcript_coord=False,
            split_by_transcript="nc",
        )


def test_add_bam_track_preserves_per_file_color_and_alpha_lists(transcript_split_gtf):
    viz = DrViz().load_gtf(str(transcript_split_gtf))

    viz.add_bam_track(
        ["a.bam", "b.bam"],
        label="Reads",
        color=["#f14432", "#4a98c9"],
        alpha=[0.6, 0.4],
        transcript_coord=True,
        split_by_transcript="nc",
    )

    spec = viz.track_specs[0]
    config = viz.track_configs[0]

    assert spec["color"] == "steelblue"
    assert spec["alpha"] == 0.6
    assert spec["file_colors"] == ["#f14432", "#4a98c9"]
    assert spec["file_alphas"] == [0.6, 0.4]
    assert config["color"] == "steelblue"
    assert config["alpha"] == 0.6


def test_split_by_transcript_none_preserves_legacy_track_behavior(transcript_split_gtf, transcript_split_bed_a):
    viz = DrViz().load_gtf(str(transcript_split_gtf))
    bed_parser = BEDParser(
        str(transcript_split_bed_a),
        transcript_coord=True,
        gtf_parser=viz.gtf_parser,
        track_label="TrackA",
    )

    grouped = bed_parser.get_grouped_anno_in_region("chr1", 90, 410)

    assert sorted(grouped) == ["peakA1", "peakA2"]
    assert grouped["peakA1"][0]["chrom"] == "chr1"
    assert grouped["peakA1"][0]["start"] == 105
    assert grouped["peakA2"][0]["start"] == 308


def test_bed_parser_groups_annotations_by_transcript(transcript_split_gtf, transcript_split_bed_a, transcript_split_bed_b):
    viz = DrViz().load_gtf(str(transcript_split_gtf))
    parser = BEDParser(
        [str(transcript_split_bed_a), str(transcript_split_bed_b)],
        transcript_coord=True,
        gtf_parser=viz.gtf_parser,
    )

    grouped = parser.get_grouped_anno_by_transcript("gene1")

    assert list(grouped) == ["ENST00000111111", "ENST00000999999"]
    assert sorted(grouped["ENST00000111111"]) == ["peakA1", "peakB1"]
    assert grouped["ENST00000111111"]["peakA1"][0]["chrom"] == "chr1"
    assert grouped["ENST00000111111"]["peakA1"][0]["start"] == 105
    assert grouped["ENST00000999999"]["peakB2"][0]["end"] == 335


def test_mixed_split_modes_are_rejected(transcript_split_gtf, transcript_split_bed_a, transcript_split_bed_b):
    parser = (
        DrViz()
        .load_gtf(str(transcript_split_gtf))
        .add_bed_track(str(transcript_split_bed_a), label="TrackA", transcript_coord=True, split_by_transcript="nc")
        .add_bed_track(str(transcript_split_bed_b), label="TrackB", transcript_coord=True, split_by_transcript="cn")
        .build()
    )

    with pytest.raises(ValueError, match="must be consistent across all split tracks"):
        parser.data_source.get_transcript_data("gene1")


def test_prepared_tracks_expand_in_nc_order_for_bed(transcript_split_gtf, transcript_split_bed_a, transcript_split_bed_b):
    parser = (
        DrViz()
        .load_gtf(str(transcript_split_gtf))
        .add_bed_track(str(transcript_split_bed_a), label="TrackA", transcript_coord=True, split_by_transcript="nc")
        .add_bed_track(str(transcript_split_bed_b), label="TrackB", transcript_coord=True, split_by_transcript="nc")
        .build()
    )

    payload = parser.data_source.get_transcript_data("gene1")

    assert [track["label"] for track in payload["prepared_tracks"]] == ["TrackA", "TrackB", "TrackA", "TrackB"]
    assert [track["transcript_id"] for track in payload["prepared_tracks"]] == [
        "ENST00000111111",
        "ENST00000111111",
        "ENST00000999999",
        "ENST00000999999",
    ]
    assert payload["right_label_groups"] == [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 1},
        {"transcript_id": "ENST00000999999", "start_index": 2, "end_index": 3},
    ]


def test_prepared_tracks_expand_in_cn_order_for_bed(transcript_split_gtf, transcript_split_bed_a, transcript_split_bed_b):
    parser = (
        DrViz()
        .load_gtf(str(transcript_split_gtf))
        .add_bed_track(str(transcript_split_bed_a), label="TrackA", transcript_coord=True, split_by_transcript="cn")
        .add_bed_track(str(transcript_split_bed_b), label="TrackB", transcript_coord=True, split_by_transcript="cn")
        .build()
    )

    payload = parser.data_source.get_transcript_data("gene1")

    assert [track["label"] for track in payload["prepared_tracks"]] == ["TrackA", "TrackA", "TrackB", "TrackB"]
    assert [track["transcript_id"] for track in payload["prepared_tracks"]] == [
        "ENST00000111111",
        "ENST00000999999",
        "ENST00000111111",
        "ENST00000999999",
    ]
    assert payload["right_label_groups"] == [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
        {"transcript_id": "ENST00000999999", "start_index": 1, "end_index": 1},
        {"transcript_id": "ENST00000111111", "start_index": 2, "end_index": 2},
        {"transcript_id": "ENST00000999999", "start_index": 3, "end_index": 3},
    ]


def test_bam_parser_returns_coverage_by_transcript(monkeypatch, transcript_split_gtf):
    class DummyRead:
        def __init__(self, blocks):
            self._blocks = blocks

        def get_blocks(self):
            return self._blocks

    class DummyAlignmentFile:
        def __init__(self, path, mode):
            self.header = {'SQ': [{'SN': 'ENST00000111111', 'LN': 100}, {'SN': 'ENST00000999999', 'LN': 80}]}

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def fetch(self, transcript_id, start, end):
            if transcript_id == 'ENST00000111111':
                return [DummyRead([(5, 8)]), DummyRead([(10, 12)])]
            if transcript_id == 'ENST00000999999':
                return [DummyRead([(2, 5)])]
            return []

    monkeypatch.setattr(api, 'BAMParser', BAMParser)
    monkeypatch.setattr(track_build, 'BAMParser', BAMParser)
    monkeypatch.setattr('drvizer.bam_parser.pysam.AlignmentFile', DummyAlignmentFile)
    monkeypatch.setattr('drvizer.bam_parser.os.path.exists', lambda path: True)

    viz = DrViz().load_gtf(str(transcript_split_gtf))
    parser = BAMParser('fake.bam', transcript_coord=True, gtf_parser=viz.gtf_parser)

    grouped = parser.get_coverage_by_transcript('gene1')

    assert list(grouped) == ['ENST00000111111', 'ENST00000999999']
    assert grouped['ENST00000111111']['x'][5:12].tolist() == [105, 106, 107, 108, 109, 110, 111]
    assert grouped['ENST00000111111']['y'][5:12].tolist() == [1, 1, 1, 0, 0, 1, 1]
    assert grouped['ENST00000999999']['x'][2:5].tolist() == [302, 303, 304]
    assert grouped['ENST00000999999']['y'][2:5].tolist() == [1, 1, 1]


def test_prepared_tracks_expand_in_nc_order_for_bam(monkeypatch, transcript_split_gtf):
    class DummyRead:
        def __init__(self, blocks):
            self._blocks = blocks

        def get_blocks(self):
            return self._blocks

    class DummyAlignmentFile:
        def __init__(self, path, mode):
            self.path = path
            self.header = {'SQ': [{'SN': 'ENST00000111111', 'LN': 100}, {'SN': 'ENST00000999999', 'LN': 80}]}

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def fetch(self, transcript_id, start, end):
            if self.path == 'copd.bam':
                if transcript_id == 'ENST00000111111':
                    return [DummyRead([(5, 8)])]
                if transcript_id == 'ENST00000999999':
                    return [DummyRead([(2, 5)])]
            if self.path == 'control.bam':
                if transcript_id == 'ENST00000111111':
                    return [DummyRead([(10, 12)])]
                if transcript_id == 'ENST00000999999':
                    return [DummyRead([(20, 22)])]
            return []

    monkeypatch.setattr(api, 'BAMParser', BAMParser)
    monkeypatch.setattr(track_build, 'BAMParser', BAMParser)
    monkeypatch.setattr('drvizer.bam_parser.pysam.AlignmentFile', DummyAlignmentFile)
    monkeypatch.setattr('drvizer.bam_parser.os.path.exists', lambda path: True)

    parser = (
        DrViz()
        .load_gtf(str(transcript_split_gtf))
        .add_bam_track(
            ['copd.bam', 'control.bam'],
            label='Coverage',
            color=['#f14432', '#4a98c9'],
            alpha=[0.6, 0.6],
            transcript_coord=True,
            split_by_transcript='nc',
        )
        .build()
    )

    payload = parser.data_source.get_transcript_data('gene1')

    assert [track['transcript_id'] for track in payload['prepared_tracks']] == [
        'ENST00000111111',
        'ENST00000999999',
    ]
    assert [track['label'] for track in payload['prepared_tracks']] == [
        'Coverage',
        'Coverage',
    ]
    assert [track['kind'] for track in payload['prepared_tracks']] == [
        'coverage',
        'coverage',
    ]
    assert [series['source_label'] for series in payload['prepared_tracks'][0]['data']['series']] == ['copd.bam', 'control.bam']
    assert payload['prepared_tracks'][0]['data']['series'][0]['y'][5:8].tolist() == [1, 1, 1]
    assert payload['prepared_tracks'][0]['data']['series'][1]['y'][10:12].tolist() == [1, 1]
    assert payload['prepared_tracks'][1]['data']['series'][0]['y'][2:5].tolist() == [1, 1, 1]
    assert payload['prepared_tracks'][1]['data']['series'][1]['y'][20:22].tolist() == [1, 1]


def test_split_bam_mean_uses_aggregated_track_instead_of_per_sample_series(monkeypatch, transcript_split_gtf):
    class DummyRead:
        def __init__(self, blocks):
            self._blocks = blocks

        def get_blocks(self):
            return self._blocks

    class DummyAlignmentFile:
        def __init__(self, path, mode):
            self.path = path
            self.header = {'SQ': [{'SN': 'ENST00000111111', 'LN': 100}, {'SN': 'ENST00000999999', 'LN': 80}]}

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def fetch(self, transcript_id, start, end):
            if self.path == 'copd.bam' and transcript_id == 'ENST00000111111':
                return [DummyRead([(5, 8)])]
            if self.path == 'control.bam' and transcript_id == 'ENST00000111111':
                return [DummyRead([(10, 12)])]
            return []

    monkeypatch.setattr(api, 'BAMParser', BAMParser)
    monkeypatch.setattr(track_build, 'BAMParser', BAMParser)
    monkeypatch.setattr('drvizer.bam_parser.pysam.AlignmentFile', DummyAlignmentFile)
    monkeypatch.setattr('drvizer.bam_parser.os.path.exists', lambda path: True)

    parser = (
        DrViz()
        .load_gtf(str(transcript_split_gtf))
        .add_bam_track(
            ['copd.bam', 'control.bam'],
            label='Coverage',
            color=['#f14432', '#4a98c9'],
            alpha=[0.6, 0.6],
            aggregate_method='mean',
            transcript_coord=True,
            split_by_transcript='nc',
        )
        .build()
    )

    payload = parser.data_source.get_transcript_data('gene1')
    track = payload['prepared_tracks'][0]

    assert 'series' not in track['data']
    assert track['data']['y'][5:8].tolist() == [0.5, 0.5, 0.5]
    assert track['data']['y'][10:12].tolist() == [0.5, 0.5]


def test_mixed_genomic_and_split_transcript_tracks_preserve_genomic_coordinates(monkeypatch, transcript_split_gtf, tmp_bed_second):
    class DummyRead:
        def __init__(self, blocks):
            self._blocks = blocks

        def get_blocks(self):
            return self._blocks

    class DummyAlignmentFile:
        def __init__(self, path, mode):
            self.header = {'SQ': [{'SN': 'ENST00000111111', 'LN': 100}, {'SN': 'ENST00000999999', 'LN': 80}]}

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def fetch(self, transcript_id, start, end):
            if transcript_id == 'ENST00000111111':
                return [DummyRead([(5, 8)]), DummyRead([(10, 12)])]
            if transcript_id == 'ENST00000999999':
                return [DummyRead([(2, 5)])]
            return []

    monkeypatch.setattr(api, 'BAMParser', BAMParser)
    monkeypatch.setattr(track_build, 'BAMParser', BAMParser)
    monkeypatch.setattr('drvizer.bam_parser.pysam.AlignmentFile', DummyAlignmentFile)
    monkeypatch.setattr('drvizer.bam_parser.os.path.exists', lambda path: True)

    parser = (
        DrViz()
        .load_gtf(str(transcript_split_gtf))
        .add_bed_track(str(tmp_bed_second), label='TE')
        .add_bam_track(
            'fake.bam',
            label='Coverage',
            transcript_coord=True,
            split_by_transcript='nc',
        )
        .build()
    )

    payload = parser.data_source.get_transcript_data('gene1')

    te_track = payload['prepared_tracks'][0]
    first_split_track = payload['prepared_tracks'][1]
    second_split_track = payload['prepared_tracks'][2]

    assert te_track['label'] == 'TE'
    assert te_track['data']['peak3'][0]['chrom'] == 'chr1'
    assert te_track['data']['peak3'][0]['start'] == 305
    assert te_track['data']['peak4'][0]['start'] == 405

    assert first_split_track['transcript_id'] == 'ENST00000111111'
    assert first_split_track['data']['x'][5:12].tolist() == [105, 106, 107, 108, 109, 110, 111]
    assert first_split_track['data']['y'][5:12].tolist() == [1, 1, 1, 0, 0, 1, 1]

    assert second_split_track['transcript_id'] == 'ENST00000999999'
    assert second_split_track['data']['x'][2:5].tolist() == [302, 303, 304]
    assert second_split_track['data']['y'][2:5].tolist() == [1, 1, 1]


def test_prepared_tracks_keep_score_kind_for_split_bed(transcript_split_gtf, transcript_split_bed_a, transcript_split_bed_b):
    parser = (
        DrViz()
        .load_gtf(str(transcript_split_gtf))
        .add_bed_track(
            [str(transcript_split_bed_a), str(transcript_split_bed_b)],
            label='ScoreTrack',
            transcript_coord=True,
            split_by_transcript='nc',
            parser_type='score',
            color=['gray', 'blue'],
            alpha=[0.1, 0.25],
        )
        .build()
    )

    payload = parser.data_source.get_transcript_data('gene1')

    assert [track['kind'] for track in payload['prepared_tracks']] == ['score', 'score']
    assert payload['prepared_tracks'][0]['file_colors'] == ['gray', 'blue']
    assert payload['prepared_tracks'][0]['file_alphas'] == [0.1, 0.25]


def test_reusable_parser_skips_labels_for_empty_split_tracks(transcript_split_gtf, tmp_bed_second):
    parser = (
        DrViz()
        .load_gtf(str(transcript_split_gtf))
        .add_bed_track(str(tmp_bed_second), label='TE')
        .add_bed_track(
            str(tmp_bed_second),
            label='m6A',
            transcript_coord=True,
            split_by_transcript='nc',
            parser_type='score',
        )
        .build()
    )

    gene_data = parser.data_source.get_transcript_data('gene1')

    assert [track['label'] for track in gene_data['prepared_tracks']] == ['TE']
    assert [config['label'] for config in parser.track_configs] == ['TE', 'm6A']


def test_split_bed_tracks_are_projected_to_genomic_coordinates(transcript_split_gtf, transcript_split_bed_a):
    parser = (
        DrViz()
        .load_gtf(str(transcript_split_gtf))
        .add_bed_track(
            str(transcript_split_bed_a),
            label='m6A',
            transcript_coord=True,
            split_by_transcript='nc',
            parser_type='score',
        )
        .build()
    )

    payload = parser.data_source.get_transcript_data('gene1')

    first_track_elements = next(iter(payload['prepared_tracks'][0]['data'].values()))
    assert first_track_elements[0]['chrom'] == 'chr1'
    assert first_track_elements[0]['start'] == 105
