import matplotlib.pyplot as plt

from drvizer.api import ReusableParser
from drvizer.visualizer import visualize_gene_transcripts


def _base_transcript_data():
    return {
        "gene_id": "gene1",
        "seqname": "chr1",
        "strand": "+",
        "identifier_type": "gene_id",
        "original_identifier": "gene1",
        "transcripts": [
            {"transcript_id": "ENST00000111111", "exons": [{"start": 100, "end": 150}], "cds": []},
            {"transcript_id": "ENST00000999999", "exons": [{"start": 300, "end": 340}], "cds": []},
        ],
    }


def test_visualizer_adds_grouped_right_labels_for_nc():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {"kind": "distribution", "data": {"peak": [{"start": 105, "end": 120}]}, "label": "TrackA", "transcript_id": "ENST00000111111"},
        {"kind": "distribution", "data": {"peak": [{"start": 125, "end": 140}]}, "label": "TrackB", "transcript_id": "ENST00000111111"},
        {"kind": "distribution", "data": {"peak": [{"start": 305, "end": 320}]}, "label": "TrackA", "transcript_id": "ENST00000999999"},
        {"kind": "distribution", "data": {"peak": [{"start": 325, "end": 338}]}, "label": "TrackB", "transcript_id": "ENST00000999999"},
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 1},
        {"transcript_id": "ENST00000999999", "start_index": 2, "end_index": 3},
    ]

    fig = visualize_gene_transcripts(transcript_data)

    right_texts = [text.get_text() for text in fig.texts if text.get_text().startswith("ENST")]
    assert right_texts == ["ENST00000111111", "ENST00000999999"]
    plt.close(fig)


def test_visualizer_adds_one_right_label_per_subtrack_for_cn():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {"kind": "distribution", "data": {"peak": [{"start": 105, "end": 120}]}, "label": "TrackA", "transcript_id": "ENST00000111111"},
        {"kind": "distribution", "data": {"peak": [{"start": 305, "end": 320}]}, "label": "TrackA", "transcript_id": "ENST00000999999"},
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
        {"transcript_id": "ENST00000999999", "start_index": 1, "end_index": 1},
    ]

    fig = visualize_gene_transcripts(transcript_data)

    right_texts = [text.get_text() for text in fig.texts if text.get_text().startswith("ENST")]
    assert right_texts == ["ENST00000111111", "ENST00000999999"]
    plt.close(fig)


def test_visualizer_truncates_long_transcript_ids():
    """Test that long transcript IDs are shortened in right-side labels."""
    long_id = "Novel_transcript_with_a_very_long_identifier_suffix_000123456789"
    transcript_data = {
        "gene_id": "gene1",
        "seqname": "chr1",
        "strand": "+",
        "identifier_type": "gene_id",
        "original_identifier": "gene1",
        "transcripts": [
            {"transcript_id": long_id, "exons": [{"start": 100, "end": 150}], "cds": []},
        ],
    }
    transcript_data["prepared_tracks"] = [
        {"kind": "distribution", "data": {"peak": [{"start": 105, "end": 120}]}, "label": "TrackA", "transcript_id": long_id},
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": long_id, "start_index": 0, "end_index": 0},
    ]

    fig = visualize_gene_transcripts(transcript_data)

    right_texts = [text.get_text() for text in fig.texts]
    assert "000123456789" in right_texts
    assert long_id not in right_texts
    plt.close(fig)


def test_visualizer_no_right_labels_without_groups():
    """Test that no right-side labels are added when right_label_groups is absent."""
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {"kind": "distribution", "data": {"peak": [{"start": 105, "end": 120}]}, "label": "TrackA"},
    ]

    fig = visualize_gene_transcripts(transcript_data)

    right_texts = [text.get_text() for text in fig.texts if text.get_text().startswith("ENST")]
    assert right_texts == []
    plt.close(fig)


def test_visualizer_no_right_labels_with_empty_groups():
    """Test that no right-side labels are added when right_label_groups is empty."""
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {"kind": "distribution", "data": {"peak": [{"start": 105, "end": 120}]}, "label": "TrackA"},
    ]
    transcript_data["right_label_groups"] = []

    fig = visualize_gene_transcripts(transcript_data)

    right_texts = [text.get_text() for text in fig.texts if text.get_text().startswith("ENST")]
    assert right_texts == []
    plt.close(fig)


def test_reusable_parser_closes_figure_when_show_is_false():
    class DummyDataSource:
        def get_transcript_data(self, gene, transcript_to_show=None):
            return {
                "gene_id": gene,
                "seqname": "chr1",
                "strand": "+",
                "identifier_type": "gene_id",
                "original_identifier": gene,
                "transcripts": [
                    {"transcript_id": "ENST00000111111", "exons": [{"start": 100, "end": 150}], "cds": []},
                ],
                "prepared_tracks": [],
            }

    parser = ReusableParser(DummyDataSource(), [])
    fig = parser.plot("gene1", show=False)

    assert fig.number not in plt.get_fignums()


def test_visualizer_positions_nc_right_labels_after_non_split_track():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {"kind": "distribution", "data": {"te": [{"start": 90, "end": 130}]}, "label": "TE"},
        {
            "kind": "coverage",
            "data": {"x": [105, 106], "y": [1, 2]},
            "label": "COPD",
            "transcript_id": "ENST00000111111",
            "file_colors": ["#f14432", "#4a98c9"],
            "file_alphas": [0.6, 0.6],
        },
        {
            "kind": "coverage",
            "data": {"x": [105, 106], "y": [2, 1]},
            "label": "Control",
            "transcript_id": "ENST00000111111",
            "file_colors": ["#f14432", "#4a98c9"],
            "file_alphas": [0.6, 0.6],
        },
        {
            "kind": "coverage",
            "data": {"x": [305, 306], "y": [1, 2]},
            "label": "COPD",
            "transcript_id": "ENST00000999999",
            "file_colors": ["#f14432", "#4a98c9"],
            "file_alphas": [0.6, 0.6],
        },
        {
            "kind": "coverage",
            "data": {"x": [305, 306], "y": [2, 1]},
            "label": "Control",
            "transcript_id": "ENST00000999999",
            "file_colors": ["#f14432", "#4a98c9"],
            "file_alphas": [0.6, 0.6],
        },
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 1, "end_index": 2},
        {"transcript_id": "ENST00000999999", "start_index": 3, "end_index": 4},
    ]

    fig = visualize_gene_transcripts(transcript_data)

    axes = fig.axes
    label_map = {text.get_text(): text for text in fig.texts if text.get_text().startswith("ENST")}

    first_label_y = label_map["ENST00000111111"].get_position()[1]
    te_center_y = (axes[1].get_position().y0 + axes[1].get_position().y1) / 2
    first_group_center_y = (axes[2].get_position().y1 + axes[3].get_position().y0) / 2

    assert abs(first_label_y - first_group_center_y) < 1e-9
    assert abs(first_label_y - te_center_y) > 1e-3
    plt.close(fig)


def test_visualizer_renders_score_tracks_with_multiple_file_colors():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "score",
            "data": {
                "peakA": [
                    {"start": 105, "end": 115, "score": 0.1},
                    {"start": 120, "end": 130, "score": 0.2},
                    {"start": 135, "end": 145, "score": 0.3},
                    {"start": 150, "end": 160, "score": 0.4},
                ]
            },
            "label": "m6A",
            "transcript_id": "ENST00000111111",
            "file_colors": ["gray", "gray", "blue", "red"],
            "file_alphas": [0.1, 0.1, 0.25, 0.25],
        }
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
    ]

    fig = visualize_gene_transcripts(transcript_data)
    ax = fig.axes[1]

    patch_facecolors = [
        tuple(round(value, 3) for value in patch.get_facecolor()[:4])
        for patch in ax.patches
    ]

    assert patch_facecolors == [
        (0.502, 0.502, 0.502, 0.1),
        (0.502, 0.502, 0.502, 0.1),
        (0.0, 0.0, 1.0, 0.25),
        (1.0, 0.0, 0.0, 0.25),
    ]
    plt.close(fig)


def test_visualizer_renders_distribution_tracks_with_track_color():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "distribution",
            "data": {"peakA": [{"start": 105, "end": 115}]},
            "label": "TE",
            "transcript_id": "ENST00000111111",
            "color": "purple",
            "alpha": 0.3,
        }
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
    ]

    fig = visualize_gene_transcripts(transcript_data)
    ax = fig.axes[1]

    patch_facecolors = [
        tuple(round(value, 3) for value in patch.get_facecolor()[:4])
        for patch in ax.patches
    ]

    assert patch_facecolors == [
        (0.502, 0.0, 0.502, 0.3),
    ]
    plt.close(fig)


class DummyDataSourceWithEmptySplitTracks:
    def get_transcript_data(self, gene, transcript_to_show=None):
        return {
            "gene_id": gene,
            "seqname": "chr1",
            "strand": "+",
            "identifier_type": "gene_id",
            "original_identifier": gene,
            "transcripts": [
                {"transcript_id": "ENST00000111111", "exons": [{"start": 100, "end": 150}], "cds": []},
            ],
            "prepared_tracks": [
                {"kind": "distribution", "data": {"peak": [{"start": 105, "end": 120}]}, "label": "TE", "color": "orange", "alpha": 0.8},
                {"kind": "coverage", "data": {"x": [100, 101], "y": [1, 2]}, "label": "Reads", "transcript_id": "ENST00000111111", "color": "steelblue", "alpha": 0.6},
            ],
        }


def test_reusable_parser_skips_empty_track_labels_in_plot():
    parser = ReusableParser(
        DummyDataSourceWithEmptySplitTracks(),
        [
            {"label": "TE", "color": "orange"},
            {"label": "m6A", "color": "orange"},
            {"label": "m5C", "color": "orange"},
            {"label": "Reads", "color": "steelblue"},
        ],
    )

    fig = parser.plot("gene1", show=False)

    axis_labels = [ax.get_ylabel() for ax in fig.axes]
    assert axis_labels == ["Transcripts", "TE", "Reads"]


class DummyDataSourceWithRepeatedVisibleLabels:
    def get_transcript_data(self, gene, transcript_to_show=None):
        return {
            "gene_id": gene,
            "seqname": "chr1",
            "strand": "+",
            "identifier_type": "gene_id",
            "original_identifier": gene,
            "transcripts": [
                {"transcript_id": "tx1", "exons": [{"start": 100, "end": 150}], "cds": []},
                {"transcript_id": "tx2", "exons": [{"start": 300, "end": 350}], "cds": []},
            ],
            "prepared_tracks": [
                {"kind": "distribution", "data": {"peak": [{"start": 90, "end": 120}]}, "label": "TE", "color": "orange", "alpha": 0.8},
                {"kind": "coverage", "data": {"x": [105, 106], "y": [1, 2]}, "label": "Reads", "transcript_id": "tx1", "color": "steelblue", "alpha": 0.6},
                {"kind": "score", "data": {"peak": [{"start": 110, "end": 115, "score": 0.3}]}, "label": "m6A", "transcript_id": "tx1", "color": "orange", "alpha": 0.8},
                {"kind": "coverage", "data": {"x": [305, 306], "y": [2, 1]}, "label": "Reads", "transcript_id": "tx2", "color": "steelblue", "alpha": 0.6},
            ],
            "right_label_groups": [
                {"transcript_id": "tx1", "start_index": 1, "end_index": 2},
                {"transcript_id": "tx2", "start_index": 3, "end_index": 3},
            ],
        }




class DummyDataSourceWithSparseTranscriptTrackCoverage:
    def get_transcript_data(self, gene, transcript_to_show=None):
        return {
            "gene_id": gene,
            "seqname": "chr1",
            "strand": "+",
            "identifier_type": "gene_id",
            "original_identifier": gene,
            "transcripts": [
                {"transcript_id": "tx1", "exons": [{"start": 100, "end": 150}], "cds": []},
                {"transcript_id": "tx2", "exons": [{"start": 300, "end": 350}], "cds": []},
            ],
            "prepared_tracks": [
                {"kind": "distribution", "data": {"peak": [{"start": 90, "end": 120}]}, "label": "TE", "color": "orange", "alpha": 0.8},
                {"kind": "coverage", "data": {"x": [105, 106], "y": [1, 2]}, "label": "Reads", "transcript_id": "tx1", "color": "steelblue", "alpha": 0.6},
                {"kind": "coverage", "data": {"x": [305, 306], "y": [2, 1]}, "label": "Reads", "transcript_id": "tx2", "color": "steelblue", "alpha": 0.6},
                {"kind": "score", "data": {"peak": [{"start": 310, "end": 315, "score": 0.3}]}, "label": "m6A", "transcript_id": "tx2", "color": "orange", "alpha": 0.8},
                {"kind": "coverage", "data": {"x": [307, 308], "y": [3, 2]}, "label": "Reads", "transcript_id": "tx2", "color": "steelblue", "alpha": 0.6},
            ],
            "right_label_groups": [
                {"transcript_id": "tx1", "start_index": 1, "end_index": 1},
                {"transcript_id": "tx2", "start_index": 2, "end_index": 4},
            ],
        }







def test_reusable_parser_repeats_track_configs_for_each_prepared_track_label():
    parser = ReusableParser(
        DummyDataSourceWithSparseTranscriptTrackCoverage(),
        [
            {"label": "TE", "color": "orange"},
            {"label": "m6A", "color": "orange"},
            {"label": "m5C", "color": "orange"},
            {"label": "pseudoU", "color": "orange"},
            {"label": "A-to-I", "color": "orange"},
            {"label": "Reads", "color": "steelblue"},
        ],
    )

    visible_configs = parser._build_visible_track_configs(
        DummyDataSourceWithSparseTranscriptTrackCoverage().get_transcript_data("gene1")["prepared_tracks"]
    )

    assert [config["label"] for config in visible_configs] == ["TE", "Reads", "Reads", "m6A", "Reads"]

