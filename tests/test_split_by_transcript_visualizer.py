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


def test_visualizer_uses_multiple_colors_for_coverage_subtracks():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "coverage",
            "data": {
                "series": [
                    {"x": [105, 106], "y": [1, 2], "color": "#f14432", "alpha": 0.6, "source_label": "COPD"},
                    {"x": [105, 106], "y": [2, 1], "color": "#4a98c9", "alpha": 0.4, "source_label": "Control"},
                ]
            },
            "label": "Reads",
            "transcript_id": "ENST00000111111",
        },
        {
            "kind": "coverage",
            "data": {
                "series": [
                    {"x": [305, 306], "y": [2, 1], "color": "#f14432", "alpha": 0.6, "source_label": "COPD"},
                    {"x": [305, 306], "y": [1, 2], "color": "#4a98c9", "alpha": 0.4, "source_label": "Control"},
                ]
            },
            "label": "Reads",
            "transcript_id": "ENST00000999999",
        },
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
        {"transcript_id": "ENST00000999999", "start_index": 1, "end_index": 1},
    ]

    fig = visualize_gene_transcripts(transcript_data)

    coverage_facecolors = [
        tuple(round(value, 3) for value in collection.get_facecolor()[0][:4])
        for axis in fig.axes[1:]
        for collection in axis.collections
        if len(collection.get_facecolor()) > 0
    ]

    assert coverage_facecolors == [
        (0.945, 0.267, 0.196, 0.6),
        (0.29, 0.596, 0.788, 0.4),
        (0.945, 0.267, 0.196, 0.6),
        (0.29, 0.596, 0.788, 0.4),
    ]
    plt.close(fig)
