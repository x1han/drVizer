"""Visualizer tests, including revived tests from commit 2134796.

2134796 revival decisions (per-test rationale):

  test_visualizer_adds_grouped_right_labels_for_nc
    Decision: revive-as-is (functionally equivalent to 2134796).
    No conflict with the P0-5 placeholder contract or with the
    0-based half-open projection (this test never projects; distribution
    tracks use transcript-coord peaks directly).

  test_visualizer_adds_one_right_label_per_subtrack_for_cn
    Decision: revive-as-is (functionally equivalent to 2134796).
    Cn split mode test with start_index == end_index per transcript.
    Renders two right-side labels. No projection; no empty-track
    placeholder interaction (cn + populated tracks).

  test_visualizer_truncates_long_transcript_ids
    Decision: revise-and-add (revised assertion to match the production
    rename ``_truncate_transcript_id`` -> ``_shorten_transcript_label``).
    The original 2134796 assertion expected the first-10 + '...' +
    last-6 rule; the production code now uses a Novel-split + 40-char
    prefix truncation. The test exercises the current rule by passing
    a long "Novel_*" id and asserting the suffix "000123456789"
    appears (the Novel token keeps the suffix as the displayed label)
    while the full long_id does NOT appear.

  test_visualizer_no_right_labels_without_groups
    Decision: revive-as-is (functionally equivalent to 2134796).
    No projection; no empty-track interaction. Survives the helper
    rename unchanged.

  test_visualizer_no_right_labels_with_empty_groups
    Decision: revive-as-is (functionally equivalent to 2134796).
    No projection; no empty-track interaction. Survives the helper
    rename unchanged.

The 2134796 commit also introduced two helper functions:

  _base_transcript_data  -> revive-as-is (functionally equivalent).
  _truncate_transcript_id -> revive-as-is (now an orphan helper; no
    current test calls it but we keep it for the audit trail).

Phase 6 additions:
  - Layer-order 'descending' tests for score and coverage (API-012 gap).
  - VIZ-013 visual-style regression tests (exon edgecolor, CDS facecolor,
    intron arrow direction, grid alpha, xlabel seqname).
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.collections as mcoll

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


def _score_track_facecolors(ax):
    """Collect per-bar facecolors for a score track after the batched-bar refactor.

    Each score track renders one BarContainer inside one Collection of bar
    rectangles. We iterate the container's patches so the assertion target
    matches what the visualizer produces today (one BarContainer holding
    all bars from a single ax.bar() call), not the per-element ax.bar that
    the previous implementation used.
    """
    facecolors = []
    for container in ax.containers:
        for patch in container.patches:
            facecolors.append(tuple(round(value, 3) for value in patch.get_facecolor()[:4]))
    return facecolors


def _distribution_track_facecolors(ax):
    """Collect facecolors for distribution tracks after the PatchCollection refactor.

    The visualizer now groups distribution rectangles into one PatchCollection
    per bed_name group. We iterate the collections on ax and concatenate the
    per-patch facecolors so the test reads element-by-element regardless of
    how many collections were emitted.
    """
    facecolors = []
    for collection in ax.collections:
        for patch in collection.get_paths():
            rgba = patch.get_facecolor() if hasattr(patch, 'get_facecolor') else None
            if rgba is not None and len(rgba) >= 4:
                facecolors.append(tuple(round(value, 3) for value in rgba[:4]))
    # collections may not surface per-patch colors via get_paths for PatchCollection
    # with match_original=True (path geometries are shared). Fall back to
    # collection-level facecolors concatenated with their original count.
    if not facecolors:
        for collection in ax.collections:
            colors = collection.get_facecolor()
            for row in colors:
                facecolors.append(tuple(round(value, 3) for value in row[:4]))
    return facecolors


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
            "layer_order": None,
        }
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
    ]

    fig = visualize_gene_transcripts(transcript_data)
    ax = fig.axes[1]

    facecolors = _score_track_facecolors(ax)

    assert facecolors == [
        (0.502, 0.502, 0.502, 0.1),
        (0.502, 0.502, 0.502, 0.1),
        (0.0, 0.0, 1.0, 0.25),
        (1.0, 0.0, 0.0, 0.25),
    ]
    plt.close(fig)


def test_visualizer_sorts_score_track_layers_with_small_values_on_top_by_default():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "score",
            "data": {
                "peakA": [
                    {"start": 105, "end": 115, "score": 0.1},
                    {"start": 105, "end": 115, "score": 0.9},
                    {"start": 105, "end": 115, "score": 0.4},
                ]
            },
            "label": "m6A",
            "transcript_id": "ENST00000111111",
            "file_colors": ["blue", "red", "green"],
            "file_alphas": [1, 1, 1],
        }
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
    ]

    fig = visualize_gene_transcripts(transcript_data)
    ax = fig.axes[1]

    facecolors = _score_track_facecolors(ax)

    assert facecolors == [
        (1.0, 0.0, 0.0, 1.0),
        (0.0, 0.502, 0.0, 1.0),
        (0.0, 0.0, 1.0, 1.0),
    ]
    plt.close(fig)


def test_visualizer_preserves_score_track_file_order_when_layer_order_is_none():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "score",
            "data": {
                "peakA": [
                    {"start": 105, "end": 115, "score": 0.1},
                    {"start": 105, "end": 115, "score": 0.9},
                    {"start": 105, "end": 115, "score": 0.4},
                ]
            },
            "label": "m6A",
            "transcript_id": "ENST00000111111",
            "file_colors": ["blue", "red", "green"],
            "file_alphas": [1, 1, 1],
            "layer_order": None,
        }
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
    ]

    fig = visualize_gene_transcripts(transcript_data)
    ax = fig.axes[1]

    facecolors = _score_track_facecolors(ax)

    assert facecolors == [
        (0.0, 0.0, 1.0, 1.0),
        (1.0, 0.0, 0.0, 1.0),
        (0.0, 0.502, 0.0, 1.0),
    ]
    plt.close(fig)


def test_visualizer_sorts_coverage_series_layers_with_small_values_on_top_by_default():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "coverage",
            "data": {
                "series": [
                    {"x": [105, 106], "y": [1, 1], "color": "blue", "alpha": 1},
                    {"x": [105, 106], "y": [9, 9], "color": "red", "alpha": 1},
                    {"x": [105, 106], "y": [4, 4], "color": "green", "alpha": 1},
                ]
            },
            "label": "Reads",
            "transcript_id": "ENST00000111111",
        }
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
    ]

    fig = visualize_gene_transcripts(transcript_data)
    ax = fig.axes[1]

    line_colors = [line.get_color() for line in ax.lines]

    assert line_colors == ["red", "green", "blue"]
    plt.close(fig)


def test_visualizer_preserves_coverage_series_order_when_layer_order_is_none():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "coverage",
            "data": {
                "series": [
                    {"x": [105, 106], "y": [1, 1], "color": "blue", "alpha": 1},
                    {"x": [105, 106], "y": [9, 9], "color": "red", "alpha": 1},
                    {"x": [105, 106], "y": [4, 4], "color": "green", "alpha": 1},
                ]
            },
            "label": "Reads",
            "transcript_id": "ENST00000111111",
            "layer_order": None,
        }
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
    ]

    fig = visualize_gene_transcripts(transcript_data)
    ax = fig.axes[1]

    line_colors = [line.get_color() for line in ax.lines]

    assert line_colors == ["blue", "red", "green"]
    plt.close(fig)


def test_visualizer_sorts_score_track_layers_descending_when_layer_order_descending():
    """API-012 (descending case): score bars put the LARGEST value on top.

    Existing tests cover 'ascending' (default) and ``None`` for score.
    This test fills the genuinely missing 'descending' case so the
    _ordered_by_layer_value branch in visualizer.py is locked for all
    three modes.

    Naming note: 'descending' here means "largest value drawn LAST (on
    top)". With layer_order='descending', ``_ordered_by_layer_value``
    sorts without reversing, so the natural ascending output
    [smallest, ..., largest] is drawn in that order, putting largest
    on top.
    """
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "score",
            "data": {
                "peakA": [
                    {"start": 105, "end": 115, "score": 0.1},
                    {"start": 105, "end": 115, "score": 0.9},
                    {"start": 105, "end": 115, "score": 0.4},
                ]
            },
            "label": "m6A",
            "transcript_id": "ENST00000111111",
            "file_colors": ["blue", "red", "green"],
            "file_alphas": [1, 1, 1],
            "layer_order": "descending",
        }
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
    ]

    fig = visualize_gene_transcripts(transcript_data)
    ax = fig.axes[1]

    facecolors = _score_track_facecolors(ax)

    # 'descending' -> smallest drawn first, largest drawn last (on top).
    assert facecolors == [
        (0.0, 0.0, 1.0, 1.0),    # score=0.1 (drawn first, on bottom)
        (0.0, 0.502, 0.0, 1.0),  # score=0.4
        (1.0, 0.0, 0.0, 1.0),    # score=0.9 (drawn last, on top)
    ]
    plt.close(fig)


def test_visualizer_sorts_coverage_series_descending_when_layer_order_descending():
    """API-012 (descending case): coverage series put the LARGEST y on top.

    Mirrors the score track descending test for coverage series.
    """
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "coverage",
            "data": {
                "series": [
                    {"x": [105, 106], "y": [1, 1], "color": "blue", "alpha": 1},
                    {"x": [105, 106], "y": [9, 9], "color": "red", "alpha": 1},
                    {"x": [105, 106], "y": [4, 4], "color": "green", "alpha": 1},
                ]
            },
            "label": "Reads",
            "transcript_id": "ENST00000111111",
            "layer_order": "descending",
        }
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
    ]

    fig = visualize_gene_transcripts(transcript_data)
    ax = fig.axes[1]

    line_colors = [line.get_color() for line in ax.lines]

    # 'descending' -> smallest y drawn first, largest drawn last (on top).
    assert line_colors == ["blue", "green", "red"]
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

    facecolors = _distribution_track_facecolors(ax)

    assert facecolors == [
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


def _capture_intron_arrows(transcript_data):
    """Render ``transcript_data`` while capturing every FancyArrowPatch created.

    The visualizer wraps intron arrows in a ``PatchCollection`` with
    ``match_original=True``, which absorbs them into the collection and removes
    them from the axes' child tree. To verify the arrow direction we monkey-patch
    ``patches.FancyArrowPatch`` to record each ``(posA, posB)`` tuple before
    handing the original class off to the visualizer.

    Returns (figure, captured_positions) where captured_positions is a list of
    ``((posA_x, posA_y), (posB_x, posB_y))`` tuples in construction order.
    """
    captured = []
    original_fancy_arrow_patch = patches.FancyArrowPatch

    class RecordingFancyArrowPatch(original_fancy_arrow_patch):
        def __init__(self, posA, posB, *args, **kwargs):
            super().__init__(posA, posB, *args, **kwargs)
            captured.append((posA, posB))

    patches.FancyArrowPatch = RecordingFancyArrowPatch
    try:
        fig = visualize_gene_transcripts(transcript_data)
    finally:
        patches.FancyArrowPatch = original_fancy_arrow_patch
    return fig, captured


def test_visualizer_exon_edgecolor_is_black():
    """VIZ-013: exon rectangles must have edgecolor='black'."""
    transcript_data = _base_transcript_data()
    fig = visualize_gene_transcripts(transcript_data)
    ax_gtf = fig.axes[0]

    # Find the exon PatchCollection (edgecolor=black, facecolor='none').
    found_black_edge = False
    for collection in ax_gtf.collections:
        # match_original=True preserves per-patch edgecolor.
        edges = collection.get_edgecolors()
        if edges.size == 0:
            continue
        # Any black edgecolor (alpha=1) on a collection means at least one
        # exon rectangle is black-edged.
        if any(c[0] == 0.0 and c[1] == 0.0 and c[2] == 0.0 for c in edges):
            found_black_edge = True
            break
    assert found_black_edge, (
        f"no exon collection has a black edge; collections: "
        f"{[(type(c).__name__, c.get_edgecolors().tolist()) for c in ax_gtf.collections]}"
    )
    plt.close(fig)


def test_visualizer_cds_facecolor_is_lightblue():
    """VIZ-013: CDS rectangles must have facecolor='lightblue'."""
    transcript_data = _base_transcript_data()
    # Add CDS to make CDS collections present.
    transcript_data["transcripts"][0]["cds"] = [{"start": 110, "end": 140}]

    fig = visualize_gene_transcripts(transcript_data)
    ax_gtf = fig.axes[0]

    found_lightblue = False
    for collection in ax_gtf.collections:
        faces = collection.get_facecolors()
        if faces.size == 0:
            continue
        # 'lightblue' is matplotlib's CSS color, which resolves to
        # rgb(0.6784313725490196, 0.8470588235294118, 0.9019607843137255).
        # Compare with a small tolerance.
        for face in faces:
            if (
                abs(face[0] - 0.6784) < 0.01
                and abs(face[1] - 0.8471) < 0.01
                and abs(face[2] - 0.9020) < 0.01
            ):
                found_lightblue = True
                break
        if found_lightblue:
            break
    assert found_lightblue, (
        f"no CDS collection has a lightblue facecolor; collections: "
        f"{[(type(c).__name__, c.get_facecolors().tolist()) for c in ax_gtf.collections]}"
    )
    plt.close(fig)


def test_visualizer_intron_arrows_point_5_to_3_plus_strand():
    """VIZ-013: plus-strand intron arrow head points toward higher genomic x."""
    transcript_data = {
        "gene_id": "gene1",
        "seqname": "chr1",
        "strand": "+",
        "identifier_type": "gene_id",
        "original_identifier": "gene1",
        "transcripts": [
            {
                "transcript_id": "ENST00000111111",
                "exons": [
                    {"start": 100, "end": 150},
                    {"start": 200, "end": 250},
                    {"start": 300, "end": 350},
                ],
                "cds": [],
            },
        ],
    }
    fig, captured = _capture_intron_arrows(transcript_data)
    assert len(captured) == 2, (
        f"expected 2 intron arrows (3 exons -> 2 introns); got {len(captured)}"
    )
    for pos_a, pos_b in captured:
        # Plus strand: arrow constructed (intron_start, intron_end) with
        # arrowstyle='->' so arrow head at posB (higher x).
        assert pos_a[0] < pos_b[0], (
            f"plus-strand arrow head should point toward higher x; "
            f"got posA={pos_a}, posB={pos_b}"
        )
    plt.close(fig)


def test_visualizer_intron_arrows_point_5_to_3_minus_strand():
    """VIZ-013: minus-strand intron arrow head points toward LOWER genomic x."""
    transcript_data = {
        "gene_id": "gene1",
        "seqname": "chr1",
        "strand": "-",
        "identifier_type": "gene_id",
        "original_identifier": "gene1",
        "transcripts": [
            {
                "transcript_id": "ENST00000111111",
                "exons": [
                    {"start": 100, "end": 150},
                    {"start": 200, "end": 250},
                    {"start": 300, "end": 350},
                ],
                "cds": [],
            },
        ],
    }
    fig, captured = _capture_intron_arrows(transcript_data)
    assert len(captured) == 2, (
        f"expected 2 intron arrows (3 exons -> 2 introns); got {len(captured)}"
    )
    for pos_a, pos_b in captured:
        # Minus strand: arrow constructed (intron_end, intron_start) so
        # posA > posB; arrow head at posB (lower x).
        assert pos_a[0] > pos_b[0], (
            f"minus-strand arrow head should point toward lower x; "
            f"got posA={pos_a}, posB={pos_b}"
        )
    plt.close(fig)


def test_visualizer_grid_alpha_is_025():
    """VIZ-013: x-axis gridlines have alpha=0.25."""
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "distribution",
            "data": {"peakA": [{"start": 105, "end": 115}]},
            "label": "TE",
            "transcript_id": "ENST00000111111",
        },
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
    ]
    fig = visualize_gene_transcripts(transcript_data)

    # Inspect x-axis gridlines on the track axes (axes[1:] for prepared_tracks).
    alphas = set()
    for ax in fig.axes[1:]:
        for line in ax.xaxis.get_gridlines():
            alpha = line.get_alpha()
            if alpha is not None:
                alphas.add(round(alpha, 4))

    assert 0.25 in alphas, (
        f"expected at least one x-axis gridline with alpha=0.25; got alphas={alphas}"
    )
    plt.close(fig)


def test_visualizer_xlabel_has_seqname():
    """VIZ-013: bottom axes' xlabel includes the seqname."""
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "distribution",
            "data": {"peakA": [{"start": 105, "end": 115}]},
            "label": "TE",
            "transcript_id": "ENST00000111111",
        },
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
    ]
    fig = visualize_gene_transcripts(transcript_data)

    last_xlabel = fig.axes[-1].get_xlabel()
    assert "chr1" in last_xlabel, (
        f"bottom axes xlabel must contain 'chr1' (seqname); got {last_xlabel!r}"
    )
    plt.close(fig)


def test_visualizer_shares_split_coverage_y_axis_within_each_transcript_group():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
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
        {
            "kind": "coverage",
            "data": {"x": [305, 306], "y": [10, 20]},
            "label": "COPD Reads",
            "transcript_id": "ENST00000999999",
            "color": "#f14432",
            "alpha": 1,
            "y_axis_group": "reads",
        },
        {
            "kind": "coverage",
            "data": {"x": [305, 306], "y": [30, 40]},
            "label": "Control Reads",
            "transcript_id": "ENST00000999999",
            "color": "#4a98c9",
            "alpha": 1,
            "y_axis_group": "reads",
        },
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 1},
        {"transcript_id": "ENST00000999999", "start_index": 2, "end_index": 3},
    ]

    fig = visualize_gene_transcripts(transcript_data)

    assert fig.axes[1].get_ylim() == (0.0, 4.4)
    assert fig.axes[2].get_ylim() == (0.0, 4.4)
    assert fig.axes[3].get_ylim() == (0.0, 44.0)
    assert fig.axes[4].get_ylim() == (0.0, 44.0)
    plt.close(fig)


def test_visualizer_keeps_different_y_axis_groups_independent():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "coverage",
            "data": {"x": [105, 106], "y": [1, 2]},
            "label": "Reads",
            "transcript_id": "ENST00000111111",
            "color": "steelblue",
            "alpha": 1,
            "y_axis_group": "reads",
        },
        {
            "kind": "coverage",
            "data": {"x": [105, 106], "y": [10, 20]},
            "label": "IP",
            "transcript_id": "ENST00000111111",
            "color": "purple",
            "alpha": 1,
            "y_axis_group": "ip",
        },
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 1},
    ]

    fig = visualize_gene_transcripts(transcript_data)

    assert fig.axes[1].get_ylim() == (0.0, 2.2)
    assert fig.axes[2].get_ylim() == (0.0, 22.0)
    plt.close(fig)


def test_visualizer_shares_score_track_y_axis_group():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "score",
            "data": {"peakA": [{"start": 105, "end": 115, "score": 0.2}]},
            "label": "Control m6A",
            "transcript_id": "ENST00000111111",
            "color": "blue",
            "alpha": 0.5,
            "y_axis_group": "m6A",
        },
        {
            "kind": "score",
            "data": {"peakB": [{"start": 120, "end": 130, "score": 0.8}]},
            "label": "COPD m6A",
            "transcript_id": "ENST00000111111",
            "color": "red",
            "alpha": 0.5,
            "y_axis_group": "m6A",
        },
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 1},
    ]

    fig = visualize_gene_transcripts(transcript_data)

    assert fig.axes[1].get_ylim() == (0.0, 0.8800000000000001)
    assert fig.axes[2].get_ylim() == (0.0, 0.8800000000000001)
    plt.close(fig)


def test_visualizer_y_axis_range_overrides_y_axis_group():
    transcript_data = _base_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "coverage",
            "data": {"x": [105, 106], "y": [1, 2]},
            "label": "Manual",
            "transcript_id": "ENST00000111111",
            "color": "steelblue",
            "alpha": 1,
            "y_axis_group": "reads",
            "y_axis_range": 10,
        },
        {
            "kind": "coverage",
            "data": {"x": [105, 106], "y": [30, 40]},
            "label": "Automatic",
            "transcript_id": "ENST00000111111",
            "color": "purple",
            "alpha": 1,
            "y_axis_group": "reads",
        },
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 1},
    ]

    fig = visualize_gene_transcripts(transcript_data)

    assert fig.axes[1].get_ylim() == (0.0, 10.0)
    assert fig.axes[2].get_ylim() == (0.0, 44.0)
    plt.close(fig)

