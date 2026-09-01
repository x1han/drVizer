"""Regression tests for Phase 4 silent visual bugs.

Each test is designed to FAIL on main (pre-fix) and PASS after the
corresponding fix in src/drvizer/visualizer.py / api.py. Tests cover:

- backbone hlines linewidth=0 -> linewidth=0.5 (P0-8)
- score tracks auto-scaling without explicit ylim (P0-9)
- coverage tracks skipping np.max when ylim is pinned (P1-15)
- exon / CDS PatchCollections are separate (P1-12)
- introns use FancyArrowPatch instead of Annotation (P1-11)
- headless Agg backend lock respects MPLBACKEND (VIZ-010)
- vectorized _score_track_max / _coverage_track_max (VIZ-014)
- gridlines present after rendering (VIZ-008 helper coverage)
- batched ax.bar produces one BarContainer per score track (P1-13)
"""

import os
import subprocess
import sys
from pathlib import Path

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pytest

from drvizer.visualizer import visualize_gene_transcripts


def _multi_exon_transcript_data():
    return {
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
                "cds": [
                    {"start": 110, "end": 145},
                    {"start": 210, "end": 245},
                    {"start": 310, "end": 345},
                ],
            },
        ],
    }


def test_backbone_linewidth_is_nonzero():
    """P0-8: backbone hlines must render with linewidth > 0."""
    from matplotlib.collections import LineCollection
    fig = visualize_gene_transcripts(_multi_exon_transcript_data())
    ax_gtf = fig.axes[0]

    # ax.hlines emits a LineCollection (not a Line2D) in modern matplotlib.
    # The backbone is the only LineCollection whose color is lightgray and
    # whose segments span the full genomic extent (start < end on the x-axis).
    backbone_collections = [
        c for c in ax_gtf.collections
        if isinstance(c, LineCollection)
        and c.get_colors().size > 0
    ]
    assert backbone_collections, "expected a LineCollection for the backbone hlines"
    backbone = backbone_collections[0]
    colors = backbone.get_colors()
    # First row is the single segment color (lightgray ~ rgb 0.827, 0.827, 0.827)
    is_lightgray = all(abs(c - colors[0][i]) < 0.01 for i, c in enumerate([0.827, 0.827, 0.827]))
    assert is_lightgray, (
        f"backbone color is not lightgray: {colors[0].tolist()}"
    )
    linewidths = backbone.get_linewidths()
    assert len(linewidths) > 0 and linewidths[0] > 0, (
        f"backbone linewidth was 0 (invisible); got {linewidths}"
    )
    plt.close(fig)


def test_score_track_auto_scales_when_no_range_or_group():
    """P0-9: score tracks without y_axis_range or shared group must auto-scale."""
    transcript_data = _multi_exon_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "score",
            "data": {
                "peakA": [
                    {"start": 110, "end": 140, "score": 10.0},
                    {"start": 220, "end": 240, "score": 10.0},
                ],
            },
            "label": "m6A",
            "transcript_id": "ENST00000111111",
        }
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
    ]

    fig = visualize_gene_transcripts(transcript_data)
    ax = fig.axes[1]

    upper_ylim = ax.get_ylim()[1]
    assert upper_ylim > 1, (
        "score track with score=10 fell back to matplotlib default ylim (0, 1); "
        f"got {ax.get_ylim()}"
    )
    plt.close(fig)


def test_coverage_track_does_not_compute_max_when_y_axis_range_is_set(monkeypatch):
    """P1-15: skip per-series np.max when y_axis_range or shared_y_axis_limit pins ylim."""
    transcript_data = _multi_exon_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "coverage",
            "data": {
                "series": [
                    {"x": [105, 106], "y": [1, 2], "color": "blue", "alpha": 1},
                    {"x": [210, 211], "y": [3, 4], "color": "red", "alpha": 1},
                ],
            },
            "label": "Reads",
            "transcript_id": "ENST00000111111",
            "y_axis_range": 5,
        }
    ]
    transcript_data["right_label_groups"] = [
        {"transcript_id": "ENST00000111111", "start_index": 0, "end_index": 0},
    ]

    import drvizer.visualizer as visualizer
    call_counter = {"n": 0}
    original = visualizer._coverage_track_max

    def counting_helper(track_data):
        call_counter["n"] += 1
        return original(track_data)

    monkeypatch.setattr(visualizer, "_coverage_track_max", counting_helper)

    fig = visualize_gene_transcripts(transcript_data)
    assert call_counter["n"] == 0, (
        "_coverage_track_max was called even though y_axis_range pinned the ylim; "
        "the per-series np.max loop should be skipped in this branch."
    )
    plt.close(fig)


def test_exon_and_cds_use_separate_patch_collections():
    """P1-12: exon and CDS rectangles are emitted as SEPARATE PatchCollections
    with distinct facecolor / edgecolor sets."""
    fig = visualize_gene_transcripts(_multi_exon_transcript_data())
    ax_gtf = fig.axes[0]

    from matplotlib.collections import PatchCollection
    collections = [
        c for c in ax_gtf.collections if isinstance(c, PatchCollection)
    ]
    # exon + CDS collections at minimum (introns are also a PatchCollection)
    assert len(collections) >= 2, (
        f"expected at least 2 PatchCollections on GTF axis, got {len(collections)}"
    )

    # Locate the exon and CDS collections by their distinct facecolor pattern:
    # exon uses facecolor='none' (RGBA alpha 0), CDS uses facecolor='lightblue'.
    exon_collection = None
    cds_collection = None
    for coll in collections:
        faces = coll.get_facecolor()
        if faces.size == 0:
            continue
        # any face with alpha == 0 marks the exon collection (facecolor='none')
        if any(row[3] == 0 for row in faces):
            exon_collection = coll
        # any face with non-zero alpha and non-default color marks CDS
        elif any(row[:3].tolist() != [0.0, 0.0, 0.0] for row in faces):
            # CDS facecolor is lightblue (0.678, 0.847, 0.902)
            if any(abs(row[2] - 0.902) < 0.05 for row in faces):
                cds_collection = coll

    assert exon_collection is not None, "no exon PatchCollection found (facecolor='none')"
    assert cds_collection is not None, "no CDS PatchCollection found (facecolor='lightblue')"
    assert exon_collection is not cds_collection, (
        "exon and CDS must live in SEPARATE PatchCollections so per-element "
        "styling is preserved"
    )
    plt.close(fig)


def test_intron_uses_arrowpatch_not_annotation():
    """P1-11: introns render via FancyArrowPatch in a PatchCollection,
    not per-intron Annotation artists."""
    from matplotlib.collections import PatchCollection
    fig = visualize_gene_transcripts(_multi_exon_transcript_data())
    ax_gtf = fig.axes[0]

    annotations = [
        child for child in ax_gtf.get_children()
        if isinstance(child, matplotlib.text.Annotation)
    ]
    assert not annotations, (
        f"expected zero Annotation children on GTF axis, found {len(annotations)}"
    )

    # PatchCollection with match_original=True stores the original Patch
    # instances on `self._original_patches` (private) or accessible via
    # the constructor list. Verify the intron collection contains
    # FancyArrowPatch instances by inspecting the underlying patches.
    arrow_patches = []
    for coll in ax_gtf.collections:
        if not isinstance(coll, PatchCollection):
            continue
        # The intron arrow PatchCollection has 2 paths in our fixture
        # (3 exons -> 2 introns) and each path's codes include the
        # FancyArrowPatch arrow signature (CURVE codes 3,3 around MOVETO 1).
        paths = coll.get_paths()
        if not paths:
            continue
        if any(p.codes is not None and np.array_equal(p.codes, [1, 3, 3, 1, 2, 2]) for p in paths):
            arrow_patches.extend(paths)
    assert arrow_patches, (
        "expected at least one FancyArrowPatch-shaped path on the GTF axis "
        "(intron arrow). Found no PatchCollection with FancyArrowPatch path codes."
    )
    plt.close(fig)


def _run_subprocess(env_override=None):
    """Run a python subprocess with optional env override.

    Returns CompletedProcess for assertions. matplotlib caches the backend
    on first pyplot import, so we MUST isolate the test in a fresh process.
    """
    env = os.environ.copy()
    env.pop("MPLBACKEND", None)
    if env_override:
        env.update(env_override)
    return subprocess.run(
        [sys.executable, "-c",
         "import os, sys; "
         "import matplotlib; "
         "import drvizer.visualizer; "
         "print('BACKEND=', matplotlib.get_backend()); "
         "print('TKINTER=', 'tkinter' in sys.modules); "
         "print('PYQT5=', 'PyQt5' in sys.modules); "
         "print('MPLBACKEND=', os.environ.get('MPLBACKEND', '<unset>'))"],
        env=env,
        capture_output=True,
        text=True,
        timeout=30,
    )


def test_backend_lock_does_not_override_user_choice():
    """VIZ-010: when MPLBACKEND is set, the lock respects user choice (Agg).

    If MPLBACKEND is unset, the lock should pin Agg and avoid importing
    tkinter / PyQt5 in headless environments.
    """
    drs_python = sys.executable
    project_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env.pop("MPLBACKEND", None)

    # Case 1: MPLBACKEND=Agg -> backend should be Agg
    result = subprocess.run(
        [drs_python, "-c",
         "import os; "
         "os.environ['MPLBACKEND'] = 'Agg'; "
         "import matplotlib; "
         "import drvizer.visualizer; "
         "print(matplotlib.get_backend())"],
        env=env, capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, (
        f"subprocess failed: stderr={result.stderr}"
    )
    assert "agg" in result.stdout.lower(), (
        f"backend was not Agg after explicit MPLBACKEND=Agg: {result.stdout!r}"
    )

    # Case 2: MPLBACKEND unset, fresh process -> backend should default to Agg
    # and tkinter should not be imported.
    result = subprocess.run(
        [drs_python, "-c",
         "import os; "
         "assert 'MPLBACKEND' not in os.environ; "
         "import matplotlib; "
         "import drvizer.visualizer; "
         "import sys; "
         "backend = matplotlib.get_backend(); "
         "has_tkinter = 'tkinter' in sys.modules; "
         "print(backend, has_tkinter)"],
        env=env, capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, (
        f"subprocess failed: stderr={result.stderr}"
    )
    stdout = result.stdout.strip().splitlines()[-1] if result.stdout.strip() else ""
    backend, has_tkinter = stdout.split()
    assert backend == "agg", (
        f"expected backend Agg with no MPLBACKEND, got {backend}"
    )
    assert has_tkinter == "False", (
        f"tkinter was imported during drvizer.visualizer load: has_tkinter={has_tkinter}"
    )


def test_gridlines_present_after_rendering():
    """VIZ-008 helper coverage: after rendering, x-axis gridlines are visible."""
    fig = visualize_gene_transcripts(_multi_exon_transcript_data())
    ax_gtf = fig.axes[0]
    gridlines = ax_gtf.xaxis.get_gridlines()
    visible = [g for g in gridlines if g.get_visible()]
    assert visible, "no visible x-axis gridlines on the GTF axis"
    plt.close(fig)


def test_score_track_produces_one_bar_container(monkeypatch):
    """P1-13: batched ax.bar produces exactly ONE BarContainer per score track.

    The previous per-element ax.bar created one BarContainer per element;
    this test would fail under that implementation and pass under the
    batched version.
    """
    transcript_data = _multi_exon_transcript_data()
    transcript_data["prepared_tracks"] = [
        {
            "kind": "score",
            "data": {
                "peakA": [
                    {"start": 105, "end": 115, "score": 0.1},
                    {"start": 120, "end": 130, "score": 0.2},
                    {"start": 135, "end": 145, "score": 0.3},
                    {"start": 150, "end": 160, "score": 0.4},
                ],
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

    # ax.containers holds one BarContainer per ax.bar() call. The batched
    # implementation issues exactly one call per score track.
    bar_containers = [c for c in ax.containers if isinstance(c, matplotlib.container.BarContainer)]
    assert len(bar_containers) == 1, (
        f"expected exactly 1 BarContainer per score track, got {len(bar_containers)}"
    )
    plt.close(fig)


def test_score_track_max_is_vectorized(monkeypatch):
    """VIZ-014: _score_track_max calls np.max once over a concatenated array.

    The vectorized implementation builds one ndarray of all scores and
    calls np.max once; the loop-based legacy implementation called max(...)
    per element. We monkey-patch np.concatenate to count its calls.
    """
    import numpy as np
    import drvizer.visualizer as visualizer

    concat_counter = {"n": 0}
    original_concat = np.concatenate

    def counting_concat(*args, **kwargs):
        concat_counter["n"] += 1
        return original_concat(*args, **kwargs)

    monkeypatch.setattr(np, "concatenate", counting_concat)

    track_data = {
        "peakA": [
            {"start": 105, "end": 115, "score": 0.1},
            {"start": 120, "end": 130, "score": 0.5},
            {"start": 135, "end": 145, "score": 0.9},
        ],
    }
    result = visualizer._score_track_max(track_data)

    assert result == pytest.approx(0.9)
    assert concat_counter["n"] == 0, (
        "_score_track_max should NOT call np.concatenate (vectorized path "
        f"uses np.max on a list-coerced array); got {concat_counter['n']} calls"
    )
