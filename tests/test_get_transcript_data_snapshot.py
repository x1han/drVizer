"""Golden-file snapshot harness for PreparedDataSource.get_transcript_data.

Pre-refactor baseline (Phase 2.1 / commit 1 of phase2-1-prepare-data-source).

This test loads 15 scenarios that together cover every branch of
PreparedDataSource.get_transcript_data, captures the returned dict via
a deterministic JSON canonicalizer, and compares it against a checked-in
artifact at tests/_artifacts/get_transcript_data_snapshots.json.

Generation mode: set DRVIZER_SNAPSHOT_UPDATE=1 to overwrite the artifact
from the current src/drvizer/api.py. The artifact is generated against
the UNMODIFIED api.py on commit 085c353 (verified by the plan's
`git diff --exit-code 085c353 -- src/drvizer/api.py` gate), which makes
its pre-refactor provenance verifiable from git history alone: the same
artifact must compare byte-identical after the SRP refactor lands.

Why this exists:
- The window (target_chrom, gene_start, gene_end) and the raw
  gene_identifier pass-through are not otherwise observable in the
  returned dict. Two echo fakes (EchoRegionTrack, EchoIdentifierTrack)
  close that hole.
- numpy arrays/scalars are canonicalized to tagged dicts carrying
  tolist()/item() plus dtype and shape.
- sort_keys=False keeps dict insertion order significant, so the fact
  that right_label_groups is inserted before prepared_tracks is itself
  pinned.
- Exception type and message are captured under {error_type, error}.
"""

import json
import os
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = REPO_ROOT / "tests" / "_artifacts" / "get_transcript_data_snapshots.json"

# conftest.py inserts the worktree's src/ at sys.path[0]. This guard makes
# that fact -- which is the entire reason the harness is meaningful --
# observable as a test rather than something we assume. The DRS env has
# an editable install whose .pth points at the MAIN repo src/, so without
# conftest.py's sys.path.insert(0, parents[1]/'src') the test would
# silently import a different api.py and the snapshot would still pass
# but prove nothing about the refactor.
def test_snapshot_module_under_test_is_local():
    import drvizer.api
    resolved = str(Path(drvizer.api.__file__).resolve())
    # startswith, not Path.is_relative_to: DRS is Python 3.8.
    assert resolved.startswith(str(REPO_ROOT) + os.sep), (
        f"drvizer.api resolved to {resolved!r}; expected inside this worktree "
        f"({REPO_ROOT!r}). conftest.py's sys.path manipulation is not winning."
    )


# ---------------------------------------------------------------------------
# Canonicalizer + capture helpers
# ---------------------------------------------------------------------------

def _canonicalize(obj):
    if isinstance(obj, np.ndarray):
        return {
            "__ndarray__": obj.tolist(),
            "dtype": str(obj.dtype),
            "shape": list(obj.shape),
        }
    if isinstance(obj, np.generic):
        return {
            "__npscalar__": obj.item(),
            "dtype": str(obj.dtype),
        }
    if isinstance(obj, dict):
        return {str(k): _canonicalize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_canonicalize(v) for v in obj]
    return obj


def _capture(fn):
    try:
        return {"ok": True, "value": _canonicalize(fn())}
    except Exception as exc:
        return {"ok": False, "error_type": type(exc).__name__, "error": str(exc)}


# ---------------------------------------------------------------------------
# Echo fakes -- make window and identifier pass-through observable
# ---------------------------------------------------------------------------

class EchoRegionTrack:
    """Reflects (chrom, gene_start, gene_end) into the snapshot via
    get_grouped_anno_in_region. Without this, the window is invisible
    to assertions; a refactor could corrupt it and every other test
    would still pass.
    """
    track_label = "EchoRegion"
    parser_type = "distribution"

    def get_grouped_anno_in_region(self, chrom, start, end):
        return {"__region__": [{"chrom": chrom, "start": start, "end": end}]}


class EchoIdentifierTrack:
    """Reflects repr(gene_identifier) and type(gene_identifier).__name__
    into the snapshot via get_coverage_for_transcripts. Locks the raw
    gene_identifier pass-through: today the raw list (not identifiers[0])
    is forwarded for multi-gene input.
    """
    track_label = "EchoIdent"
    transcript_coord = True

    def get_coverage_for_transcripts(self, gene_identifier):
        return [repr(gene_identifier)], [type(gene_identifier).__name__]


class TxCoordCovTrack:
    """Hits the transcript_coord + get_coverage_for_transcripts branch
    with a real numpy payload, so dtype/shape are part of the snapshot.
    """
    track_label = "TxCov"
    color = "purple"
    alpha = 0.4
    transcript_coord = True

    def get_coverage_for_transcripts(self, gene_identifier):
        return np.array([100, 101], dtype=np.int64), np.array([7, 9], dtype=np.int64)


class UnknownTrack:
    """Hits the `else: continue` path (no dispatchable method)."""
    track_label = "Unknown"


class SplitWithoutByTranscript:
    """Hits the split fall-through invariant: split mode is set but the
    track has neither get_coverage_by_transcript nor
    get_grouped_anno_by_transcript, so it is prepared as an ordinary
    non-split track while its mode is still added to split_modes.
    """
    track_label = "SplitFallthrough"
    split_by_transcript = "nc"

    def get_coverage_in_region(self, chrom, start, end):
        return np.array([start, end], dtype=np.int64), np.array([1, 2], dtype=np.int64)


class NoExonGTF:
    """Fake GTF parser returning one transcript with exons=[]. Hits the
    (0, 1000) fallback when all_starts/all_ends are empty.
    """
    def get_transcript_data(self, gene_identifier):
        return {
            "gene_id": "g0",
            "identifier_type": "gene_id",
            "original_identifier": gene_identifier,
            "seqname": "chrZ",
            "strand": "+",
            "transcripts": [
                {"transcript_id": "t0", "strand": "+", "exons": [], "cds": []}
            ],
        }


# ---------------------------------------------------------------------------
# Local GTF fixtures (two-gene same-chrom and two-gene diff-chrom).
# ---------------------------------------------------------------------------

@pytest.fixture
def two_gene_same_chrom_gtf(tmp_path):
    path = tmp_path / "two_gene_same_chrom.gtf"
    path.write_text(textwrap.dedent("""\
        chr1\ttest\texon\t100\t149\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1"; gene_name "GENE1";
        chr1\ttest\texon\t2000\t2049\t.\t+\t.\tgene_id "gene2"; transcript_id "tx2"; gene_name "GENE2";
    """))
    return path


@pytest.fixture
def two_gene_diff_chrom_gtf(tmp_path):
    path = tmp_path / "two_gene_diff_chrom.gtf"
    path.write_text(textwrap.dedent("""\
        chr1\ttest\texon\t100\t149\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1"; gene_name "GENE1";
        chr2\ttest\texon\t100\t149\t.\t+\t.\tgene_id "gene2"; transcript_id "tx2"; gene_name "GENE2";
    """))
    return path


# ---------------------------------------------------------------------------
# 15 scenarios
# ---------------------------------------------------------------------------

def _all_scenarios(
    transcript_split_gtf,
    transcript_split_bed_a,
    transcript_split_bed_b,
    tmp_bed_second,
    two_gene_same_chrom_gtf,
    two_gene_diff_chrom_gtf,
):
    """Return dict of scenario_id -> {"ok": bool, "value"|("error_type","error")}.

    The fixture parametrize above injects tmp_path-derived paths into each
    scenario's lambdas. We deliberately construct the GTF/BED fixtures
    locally (rather than reusing the conftest ones directly in this
    function body) so the fixtures' tmp_path is the worktree's pytest
    tmp_path -- any other path source would leak temp paths into the
    artifact.
    """
    from drvizer.api import DrViz, PreparedDataSource

    # Genomic BED for s02 / s07.
    genomic_bed = tmp_bed_second.parent / "genomic.bed"
    genomic_bed.write_text(textwrap.dedent("""\
        chr1\t305\t320\tpeak3\t8\t+
        chr1\t405\t420\tpeak4\t4\t+
    """))

    # BED with data for only one transcript (for s06 split-cn missing tx).
    bed_b_tx2_only = transcript_split_bed_a.parent / "track_b_tx2_only.bed"
    bed_b_tx2_only.write_text("ENST00000999999\t25\t35\tpeakB2\t7\t+\n")

    scenarios = {}

    # s01: single gene, no tracks -- the minimal happy path.
    def s01():
        viz = DrViz().load_gtf(str(transcript_split_gtf))
        return viz.build().data_source.get_transcript_data("gene1")
    scenarios["s01_single_gene_no_tracks"] = _capture(s01)

    # s02: single gene + genomic BED -- exercises get_grouped_anno_in_region.
    def s02():
        viz = DrViz().load_gtf(str(transcript_split_gtf)).add_bed_track(
            str(genomic_bed), label="TE"
        )
        return viz.build().data_source.get_transcript_data("gene1")
    scenarios["s02_single_gene_genomic_bed"] = _capture(s02)

    # s03: transcript_to_show as a str that matches -- filters to one tx.
    def s03():
        viz = DrViz().load_gtf(str(transcript_split_gtf))
        return viz.build().data_source.get_transcript_data(
            "gene1", transcript_to_show="ENST00000999999"
        )
    scenarios["s03_transcript_to_show_str"] = _capture(s03)

    # s04: transcript_to_show that matches nothing -- filter is silently
    # ignored, window stays on all transcripts.
    def s04():
        viz = DrViz().load_gtf(str(transcript_split_gtf))
        return viz.build().data_source.get_transcript_data(
            "gene1", transcript_to_show="NOPE"
        )
    scenarios["s04_transcript_to_show_unmatched"] = _capture(s04)

    # s05: split nc with two BED tracks -- 4 prepared_tracks + right_label_groups.
    def s05():
        viz = (
            DrViz()
            .load_gtf(str(transcript_split_gtf))
            .add_bed_track(str(transcript_split_bed_a), label="TrackA",
                           transcript_coord=True, split_by_transcript="nc")
            .add_bed_track(str(transcript_split_bed_b), label="TrackB",
                           transcript_coord=True, split_by_transcript="nc")
        )
        return viz.build().data_source.get_transcript_data("gene1")
    scenarios["s05_split_nc_two_bed"] = _capture(s05)

    # s06: split cn where TrackB has data for only one transcript --
    # empty:True placeholder padding exercises P0-5.
    def s06():
        viz = (
            DrViz()
            .load_gtf(str(transcript_split_gtf))
            .add_bed_track(str(transcript_split_bed_a), label="TrackA",
                           transcript_coord=True, split_by_transcript="cn")
            .add_bed_track(str(bed_b_tx2_only), label="TrackB",
                           transcript_coord=True, split_by_transcript="cn")
        )
        return viz.build().data_source.get_transcript_data("gene1")
    scenarios["s06_split_cn_missing_tx"] = _capture(s06)

    # s07: mixed genomic + split score track.
    def s07():
        viz = (
            DrViz()
            .load_gtf(str(transcript_split_gtf))
            .add_bed_track(str(genomic_bed), label="TE")
            .add_bed_track(str(transcript_split_bed_a), label="m6A",
                           transcript_coord=True, split_by_transcript="nc",
                           parser_type="score")
        )
        return viz.build().data_source.get_transcript_data("gene1")
    scenarios["s07_mixed_genomic_and_split"] = _capture(s07)

    # s08: mixed nc/cn split modes -> ValueError.
    def s08():
        viz = (
            DrViz()
            .load_gtf(str(transcript_split_gtf))
            .add_bed_track(str(transcript_split_bed_a), label="TrackA",
                           transcript_coord=True, split_by_transcript="nc")
            .add_bed_track(str(transcript_split_bed_b), label="TrackB",
                           transcript_coord=True, split_by_transcript="cn")
        )
        return viz.build().data_source.get_transcript_data("gene1")
    scenarios["s08_mixed_split_modes_error"] = _capture(s08)

    # s09: two genes on the same chromosome -- merged transcripts,
    # comma-joined original_identifier, widened window.
    def s09():
        viz = DrViz().load_gtf(str(two_gene_same_chrom_gtf)).add_bed_track(
            str(genomic_bed), label="TE"
        )
        return viz.build().data_source.get_transcript_data(["gene1", "gene2"])
    scenarios["s09_multi_gene_same_chrom"] = _capture(s09)

    # s10: two genes on different chromosomes -> ValueError.
    def s10():
        viz = DrViz().load_gtf(str(two_gene_diff_chrom_gtf))
        return viz.build().data_source.get_transcript_data(["gene1", "gene2"])
    scenarios["s10_multi_gene_diff_chrom_error"] = _capture(s10)

    # s11: transcript_coord + get_coverage_for_transcripts branch
    # via PreparedDataSource directly with a numpy payload.
    gtf_parser = DrViz().load_gtf(str(transcript_split_gtf)).gtf_parser
    def s11():
        return PreparedDataSource(gtf_parser, tracks=[TxCoordCovTrack()]
                                  ).get_transcript_data("gene1")
    scenarios["s11_tx_coord_coverage_track"] = _capture(s11)

    # s12: undispatchable track is skipped while a following track is
    # still prepared -- proves the `else: continue` contributes nothing
    # yet does not abort the loop.
    def s12():
        return PreparedDataSource(
            gtf_parser, tracks=[UnknownTrack(), TxCoordCovTrack()]
        ).get_transcript_data("gene1")
    scenarios["s12_unknown_track_skipped"] = _capture(s12)

    # s13: split_by_transcript set but no by_transcript accessor --
    # the fall-through invariant (prepared as ordinary track,
    # split_track_specs empty, no right_label_groups).
    def s13():
        return PreparedDataSource(
            gtf_parser, tracks=[SplitWithoutByTranscript()]
        ).get_transcript_data("gene1")
    scenarios["s13_split_fallthrough_no_by_transcript"] = _capture(s13)

    # s14: transcripts with no exons -> (0, 1000) window fallback.
    def s14():
        return PreparedDataSource(NoExonGTF(), tracks=[TxCoordCovTrack()]
                                  ).get_transcript_data("gZ")
    scenarios["s14_no_exons_window_fallback"] = _capture(s14)

    # s15: split track + multi-gene -> ValueError raised before any
    # GTF lookup.
    def s15():
        return PreparedDataSource(
            gtf_parser, tracks=[SplitWithoutByTranscript()]
        ).get_transcript_data(["gene1", "gene2"])
    scenarios["s15_split_multi_gene_error"] = _capture(s15)

    # Window echoes via EchoRegionTrack -- locks (target_chrom,
    # gene_start, gene_end) values directly.
    def w_all_transcripts():
        return PreparedDataSource(gtf_parser, tracks=[EchoRegionTrack()]
                                  ).get_transcript_data("gene1")
    def w_transcript_to_show():
        return PreparedDataSource(gtf_parser, tracks=[EchoRegionTrack()]
                                  ).get_transcript_data(
            "gene1", transcript_to_show="ENST00000999999"
        )
    def w_multi_gene():
        parser = DrViz().load_gtf(str(two_gene_same_chrom_gtf)).gtf_parser
        return PreparedDataSource(parser, tracks=[EchoRegionTrack()]
                                  ).get_transcript_data(["gene1", "gene2"])
    def w_no_exon_fallback():
        return PreparedDataSource(NoExonGTF(), tracks=[EchoRegionTrack()]
                                  ).get_transcript_data("gZ")

    scenarios["w01_window_all_transcripts"] = _capture(w_all_transcripts)
    scenarios["w02_window_transcript_to_show"] = _capture(w_transcript_to_show)
    scenarios["w03_window_multi_gene"] = _capture(w_multi_gene)
    scenarios["w04_window_no_exon_fallback"] = _capture(w_no_exon_fallback)

    # Identifier echoes via EchoIdentifierTrack -- locks the raw
    # gene_identifier pass-through.
    def ident_single():
        return PreparedDataSource(gtf_parser, tracks=[EchoIdentifierTrack()]
                                  ).get_transcript_data("gene1")
    def ident_multi_gene_raw_list():
        parser = DrViz().load_gtf(str(two_gene_same_chrom_gtf)).gtf_parser
        return PreparedDataSource(parser, tracks=[EchoIdentifierTrack()]
                                  ).get_transcript_data(["gene1", "gene2"])

    scenarios["i01_identifier_single_gene"] = _capture(ident_single)
    scenarios["i02_identifier_multi_gene_raw_list"] = _capture(ident_multi_gene_raw_list)

    return scenarios


SCENARIO_IDS = [
    "s01_single_gene_no_tracks",
    "s02_single_gene_genomic_bed",
    "s03_transcript_to_show_str",
    "s04_transcript_to_show_unmatched",
    "s05_split_nc_two_bed",
    "s06_split_cn_missing_tx",
    "s07_mixed_genomic_and_split",
    "s08_mixed_split_modes_error",
    "s09_multi_gene_same_chrom",
    "s10_multi_gene_diff_chrom_error",
    "s11_tx_coord_coverage_track",
    "s12_unknown_track_skipped",
    "s13_split_fallthrough_no_by_transcript",
    "s14_no_exons_window_fallback",
    "s15_split_multi_gene_error",
    "w01_window_all_transcripts",
    "w02_window_transcript_to_show",
    "w03_window_multi_gene",
    "w04_window_no_exon_fallback",
    "i01_identifier_single_gene",
    "i02_identifier_multi_gene_raw_list",
]


@pytest.fixture
def all_scenarios(
    transcript_split_gtf,
    transcript_split_bed_a,
    transcript_split_bed_b,
    tmp_bed_second,
    two_gene_same_chrom_gtf,
    two_gene_diff_chrom_gtf,
):
    return _all_scenarios(
        transcript_split_gtf,
        transcript_split_bed_a,
        transcript_split_bed_b,
        tmp_bed_second,
        two_gene_same_chrom_gtf,
        two_gene_diff_chrom_gtf,
    )


@pytest.mark.parametrize("scenario_id", SCENARIO_IDS)
def test_scenario_matches_golden(scenario_id, all_scenarios):
    """Each scenario must match the golden artifact byte-for-byte."""
    if os.environ.get("DRVIZER_SNAPSHOT_UPDATE") == "1":
        ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
        # Recompute everything in a single pass so partial regenerations
        # can never happen: parametrize runs each scenario_id in its own
        # test, so writing inside one of them is the natural integration
        # point. We write every scenario at once and skip the comparison
        # tests in update mode.
        ARTIFACT.write_text(
            json.dumps(all_scenarios, indent=2, sort_keys=False)
        )
        pytest.skip(f"snapshot updated; rerun without DRVIZER_SNAPSHOT_UPDATE to compare")

    expected = json.loads(ARTIFACT.read_text())[scenario_id]
    actual = all_scenarios[scenario_id]
    # sort_keys MUST stay False -- dict insertion order is a real
    # observable (right_label_groups is inserted before prepared_tracks).
    assert (
        json.dumps(actual, indent=2, sort_keys=False)
        == json.dumps(expected, indent=2, sort_keys=False)
    ), f"snapshot drift in {scenario_id!r}"


def test_scenario_count_matches_artifact():
    """The artifact must contain every scenario the harness produces."""
    if os.environ.get("DRVIZER_SNAPSHOT_UPDATE") == "1":
        pytest.skip("update mode")
    recorded = json.loads(ARTIFACT.read_text())
    assert set(recorded.keys()) == set(SCENARIO_IDS)
    assert len(recorded) == len(SCENARIO_IDS)


# ---------------------------------------------------------------------------
# Phase 2.2 cache-hit scenario
# ---------------------------------------------------------------------------

def test_cache_hit_returns_same_payload_without_reparse(tmp_gtf, tmp_bed, tmp_bed_second):
    """Same query twice on the same ReusableParser must hit the LRU
    cache: identical prepared_tracks payload, no second parse call.
    """
    # Construct a DrViz with one BED track, build it, and call
    # get_transcript_data twice. The second call must produce the
    # same payload bytes-for-bytes and must NOT trigger a new
    # parse_bed on the BEDParser.
    from drvizer.bed_parser import BEDParser

    parse_calls = {"n": 0}
    original_parse = BEDParser.parse_bed

    def counting_parse(self, *args, **kwargs):
        parse_calls["n"] += 1
        return original_parse(self, *args, **kwargs)

    monkeypatch = pytest.MonkeyPatch()
    try:
        monkeypatch.setattr(BEDParser, "parse_bed", counting_parse)
        from drvizer.api import DrViz
        viz = (
            DrViz()
            .load_gtf(str(tmp_gtf))
            .add_bed_track(str(tmp_bed), label="TE")
        )
        parser = viz.build()
        # Reset the counter: build() is allowed to parse.
        parse_calls["n"] = 0

        first = parser.data_source.get_transcript_data("GENE1")
        first_calls = parse_calls["n"]
        second = parser.data_source.get_transcript_data("GENE1")
        second_calls = parse_calls["n"]

        # No new parse on the second call.
        assert first_calls == 0
        assert second_calls == 0
        # Identical payload (cache hit returns the same dict).
        assert first["prepared_tracks"] == second["prepared_tracks"]
    finally:
        monkeypatch.undo()