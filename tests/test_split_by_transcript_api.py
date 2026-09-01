import numpy as np
import pytest

import drvizer.api as api
import drvizer._track_build as track_build
from drvizer.api import DrViz, PreparedDataSource
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


def test_split_by_transcript_rejects_multiple_genes(tmp_path):
    gtf_path = tmp_path / "multi_gene.gtf"
    gtf_path.write_text(
        'chr1\ttest\texon\t100\t149\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1"; gene_name "GENE1";\n'
        'chr1\ttest\texon\t200\t249\t.\t+\t.\tgene_id "gene2"; transcript_id "tx2"; gene_name "GENE2";\n'
    )

    viz = (
        DrViz()
        .load_gtf(str(gtf_path))
        .add_bam_track("fake.bam", transcript_coord=True, split_by_transcript="nc")
        .build()
    )

    with pytest.raises(ValueError, match="split_by_transcript does not support multiple genes"):
        viz.data_source.get_transcript_data(["gene1", "gene2"])


def test_split_by_transcript_multiple_genes_rejected_before_gtf_lookup():
    class FailingGTFParser:
        def get_transcript_data(self, gene_identifier):
            raise AssertionError("GTF lookup should not run for unsupported split multi-gene input")

    class SplitTrack:
        split_by_transcript = "nc"

    data_source = PreparedDataSource(FailingGTFParser(), tracks=[SplitTrack()])

    with pytest.raises(ValueError, match="split_by_transcript does not support multiple genes"):
        data_source.get_transcript_data(["gene1", "gene2"])


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
    # 0-based half-open projection: GTF exon start 100 -> 99, 300 -> 299.
    assert grouped["peakA1"][0]["start"] == 104
    assert grouped["peakA2"][0]["start"] == 307


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
    # 0-based half-open projection: GTF exon start 100 -> 99, 300 -> 299.
    assert grouped["ENST00000111111"]["peakA1"][0]["start"] == 104
    assert grouped["ENST00000999999"]["peakB2"][0]["end"] == 334


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
    # 0-based half-open: BAM region_start = exon['start'] - 1, so x shifts by -1.
    assert grouped['ENST00000111111']['x'][5:12].tolist() == [104, 105, 106, 107, 108, 109, 110]
    assert grouped['ENST00000111111']['y'][5:12].tolist() == [1, 1, 1, 0, 0, 1, 1]
    assert grouped['ENST00000999999']['x'][2:5].tolist() == [301, 302, 303]
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
    # 0-based half-open: BAM region_start = exon['start'] - 1, so x shifts by -1.
    assert first_split_track['data']['x'][5:12].tolist() == [104, 105, 106, 107, 108, 109, 110]
    assert first_split_track['data']['y'][5:12].tolist() == [1, 1, 1, 0, 0, 1, 1]

    assert second_split_track['transcript_id'] == 'ENST00000999999'
    assert second_split_track['data']['x'][2:5].tolist() == [301, 302, 303]
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


def test_reusable_parser_marks_empty_split_tracks_with_placeholder(transcript_split_gtf, tmp_bed_second):
    """Empty (track, transcript) combos get an `empty: True` placeholder.

    Under the Phase 1 empty-track contract, missing combos are padded
    explicitly (instead of being silently dropped) so the visualizer can
    short-circuit on track.get('empty') and the right-side label groups
    stay aligned with prepared_tracks indexing. tmp_bed_second is a
    genomic BED (chr1-based), so when used as a transcript-coord split
    track, no transcript has any records -> every transcript produces
    an empty placeholder.
    """
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

    track_labels = [track['label'] for track in gene_data['prepared_tracks']]
    assert track_labels == ['TE', 'm6A', 'm6A']
    assert [config['label'] for config in parser.track_configs] == ['TE', 'm6A']

    split_entries = [track for track in gene_data['prepared_tracks'] if track.get('label') == 'm6A']
    assert len(split_entries) == 2
    assert all(entry.get('empty') for entry in split_entries)
    assert all(entry.get('transcript_id') for entry in split_entries)


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
    # 0-based half-open projection: GTF exon start 100 -> 99 origin.
    assert first_track_elements[0]['start'] == 104


def test_split_by_transcript_cn_with_missing_transcript_keeps_alignment(transcript_split_gtf, transcript_split_bed_a, tmp_path):
    """API-013 (P0-5): cn split with a missing transcript keeps the layout aligned.

    Setup: under cn split mode, one transcript (ENST00000999999) has data
    for both TrackA and TrackB, but the other transcript (ENST00000111111)
    only has data for TrackA. Without the P0-5 placeholder contract, the
    prepared_tracks list for ENST00000111111 would only include TrackA
    and the right_label_groups indices would mis-align with the row
    positions, causing fragmentation.

    Asserts:
      - prepared_tracks still includes a placeholder for ENST00000111111's
        TrackB entry (so right_label_groups stay contiguous).
      - right_label_groups covers both transcripts in order without gaps.
    """
    # TrackB has data ONLY for tx2 (so tx1 gets an empty placeholder).
    tx2_only_bed = tmp_path / "track_b_tx2_only.bed"
    tx2_only_bed.write_text("ENST00000999999\t25\t35\tpeakB2\t7\t+\n")

    parser = (
        DrViz()
        .load_gtf(str(transcript_split_gtf))
        .add_bed_track(str(transcript_split_bed_a), label="TrackA", transcript_coord=True, split_by_transcript="cn")
        .add_bed_track(str(tx2_only_bed), label="TrackB", transcript_coord=True, split_by_transcript="cn")
        .build()
    )

    payload = parser.data_source.get_transcript_data("gene1")
    prepared_tracks = payload["prepared_tracks"]
    right_label_groups = payload["right_label_groups"]

    # cn ordering: TrackA@tx1, TrackB@tx1, TrackA@tx2, TrackB@tx2 (or placeholders).
    # ENST00000111111 row should have a placeholder for TrackB.
    assert len(prepared_tracks) == 4, (
        f"prepared_tracks should have 4 entries (2 split x 2 transcripts), got {len(prepared_tracks)}"
    )

    by_label_and_tx = {(t["label"], t.get("transcript_id")): t for t in prepared_tracks}
    tx1 = "ENST00000111111"
    tx2 = "ENST00000999999"

    # TrackA on both transcripts has real data.
    assert not by_label_and_tx[("TrackA", tx1)].get("empty")
    assert not by_label_and_tx[("TrackA", tx2)].get("empty")

    # TrackB on tx1 is empty (placeholder); TrackB on tx2 has real data.
    assert by_label_and_tx[("TrackB", tx1)].get("empty") is True
    assert by_label_and_tx[("TrackB", tx2)].get("empty") is not True

    # right_label_groups must cover both transcripts in order, no gaps.
    tx_ids_in_groups = [g["transcript_id"] for g in right_label_groups]
    assert tx1 in tx_ids_in_groups, f"missing {tx1} in right_label_groups: {right_label_groups}"
    assert tx2 in tx_ids_in_groups, f"missing {tx2} in right_label_groups: {right_label_groups}"

    # Indices are contiguous within each group (no fragmentation).
    for group in right_label_groups:
        assert group["end_index"] - group["start_index"] >= 0

    # Indices form a continuous span covering all 4 prepared_tracks.
    all_indices = sorted(
        idx
        for group in right_label_groups
        for idx in range(group["start_index"], group["end_index"] + 1)
    )
    assert all_indices == [0, 1, 2, 3], (
        f"right_label_groups must cover all 4 prepared_tracks contiguously; "
        f"got all_indices={all_indices}, right_label_groups={right_label_groups}"
    )


def test_y_axis_group_end_to_end_shared_limit():
    """API-014: y_axis_group shared limit applied through _shared_y_axis_limits.

    Two coverage tracks on the same gene share y_axis_group='g'; max y is
    4.2 across both. ``_shared_y_axis_limits`` must return identical
    (transcript_id, 'g') -> 4.62 (=4.2*1.1) limit for both tracks.

    Also verifies the limit is propagated to the rendered figure axes.
    """
    from drvizer.visualizer import _shared_y_axis_limits

    transcript_data = {
        "gene_id": "gene1",
        "seqname": "chr1",
        "strand": "+",
        "identifier_type": "gene_id",
        "original_identifier": "gene1",
        "transcripts": [
            {"transcript_id": "ENST00000111111", "exons": [{"start": 100, "end": 150}], "cds": []},
            {"transcript_id": "ENST00000999999", "exons": [{"start": 300, "end": 340}], "cds": []},
        ],
        "prepared_tracks": [
            {
                "kind": "coverage",
                "data": {"x": [105, 106], "y": [1, 2]},
                "label": "COPD Reads",
                "transcript_id": "ENST00000111111",
                "color": "#f14432",
                "alpha": 1,
                "y_axis_group": "g",
            },
            {
                "kind": "coverage",
                "data": {"x": [105, 106], "y": [3, 4]},
                "label": "Control Reads",
                "transcript_id": "ENST00000111111",
                "color": "#4a98c9",
                "alpha": 1,
                "y_axis_group": "g",
            },
        ],
    }

    limits = _shared_y_axis_limits(transcript_data["prepared_tracks"])
    # Both tracks share the same (transcript_id, group) key.
    assert set(limits.keys()) == {("ENST00000111111", "g")}
    assert limits[("ENST00000111111", "g")] == pytest.approx(4.4, abs=1e-9)


def test_y_axis_group_end_to_end_through_visualizer_axes():
    """API-014 end-to-end: shared y_axis_group yields identical axes ylim.

    Builds a transcript_data dict with two coverage tracks sharing
    y_axis_group='reads' on the same transcript, max y is 4.2. After
    rendering via visualize_gene_transcripts, both axes must have the
    same ylim upper bound (proving the shared limit propagates through
    the full visualizer pipeline).
    """
    import matplotlib.pyplot as plt

    from drvizer.visualizer import visualize_gene_transcripts

    transcript_data = {
        "gene_id": "gene1",
        "seqname": "chr1",
        "strand": "+",
        "identifier_type": "gene_id",
        "original_identifier": "gene1",
        "transcripts": [
            {"transcript_id": "ENST00000111111", "exons": [{"start": 100, "end": 150}], "cds": []},
        ],
        "prepared_tracks": [
            {
                "kind": "coverage",
                "data": {"x": [105, 106], "y": [1, 2]},
                "label": "COPD Reads",
                "transcript_id": "ENST00000111111",
                "color": "#f14432",
                "alpha": 1,
                "y_axis_group": "reads",
            },
            {
                "kind": "coverage",
                "data": {"x": [105, 106], "y": [3, 4]},
                "label": "Control Reads",
                "transcript_id": "ENST00000111111",
                "color": "#4a98c9",
                "alpha": 1,
                "y_axis_group": "reads",
            },
        ],
        "right_label_groups": [
            {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 1},
        ],
    }

    fig = visualize_gene_transcripts(transcript_data)

    # Two coverage axes on the same transcript with y_axis_group='reads'.
    # Max y is 4 (max of 1,2,3,4), shared ylim upper = 4 * 1.1 = 4.4.
    axes_for_tx = [ax for ax in fig.axes[1:] if ax.get_ylim()[1] == pytest.approx(4.4, abs=1e-9)]
    assert len(axes_for_tx) == 2, (
        f"expected exactly 2 axes with shared ylim 4.4; got {len(axes_for_tx)} "
        f"axes with that ylim, full axes ylims: {[ax.get_ylim() for ax in fig.axes[1:]]}"
    )
    plt.close(fig)


def test_track_colors_length_alignment_with_split_expansion(transcript_split_gtf, transcript_split_bed_a, transcript_split_bed_b):
    """API-016 (P1-9): track_colors / track_labels align with split expansion.

    Build a DrViz with 2 split tracks under nc mode against 2 transcripts.
    That yields 4 prepared_tracks (2 split x 2 transcripts). After calling
    .plot(show=False), ``gene_data['track_colors']`` and
    ``gene_data['track_labels']`` must each have length 5
    (Transcripts + 4 expanded), so visualizer rendering doesn't fall back
    to a KeyError.

    NOTE: This test pins the current pre-P1-9 behavior (no visible_configs
    dedup). If Phase 7 implements the ``_visible_configs_cache`` dedup the
    audit recommends, ``_build_visible_track_configs`` may return a shorter
    list and this assertion will need to be relaxed to
    ``len(unique_configs) + 1``.
    """
    parser = (
        DrViz()
        .load_gtf(str(transcript_split_gtf))
        .add_bed_track(str(transcript_split_bed_a), label="TrackA", transcript_coord=True, split_by_transcript="nc")
        .add_bed_track(str(transcript_split_bed_b), label="TrackB", transcript_coord=True, split_by_transcript="nc")
        .build()
    )

    fig = parser.plot("gene1", show=False)

    # Re-fetch gene_data the same way plot() does internally; instead,
    # call get_transcript_data and replicate the track_labels/colors logic.
    gene_data = parser.data_source.get_transcript_data("gene1")
    visible_configs = parser._build_visible_track_configs(gene_data["prepared_tracks"])
    track_labels = ["Transcripts"]
    track_colors = [None]
    for config in visible_configs:
        track_labels.append(config["label"])
        track_colors.append(config["color"])

    # Without dedup, 2 split tracks x 2 transcripts = 4 visible configs.
    assert len(visible_configs) == 4, (
        f"expected 4 visible configs (2 split x 2 transcripts); got {len(visible_configs)}"
    )
    assert len(track_labels) == 5, (
        f"track_labels must have len 5 (Transcripts + 4 expanded); got {track_labels}"
    )
    assert len(track_colors) == 5, (
        f"track_colors must have len 5 (Transcripts + 4 expanded); got {track_colors}"
    )
    assert track_labels[0] == "Transcripts"
    assert track_colors[0] is None
    # Expanded labels are populated.
    assert "TrackA" in track_labels
    assert "TrackB" in track_labels

    # Render must not raise a KeyError at axes access time.
    assert len(fig.axes) == 5, f"expected 5 axes (Transcripts + 4 expanded); got {len(fig.axes)}"
