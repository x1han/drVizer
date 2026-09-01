"""Property-based invariants for ``GTFParser.convert_transcript_to_genomic_segments``.

The projection path emits 0-based half-open genomic intervals. These tests
lock the three mathematical properties the contract promises:

1. Length conservation: ``sum(|seg|) == T_end - T_start``.
2. Plus-strand monotonicity: ``G_start[i] < G_start[i+1]`` in the order
   returned (which is also transcript 5'→3').
3. Minus-strand monotonicity: in transcript 5'→3' order, ``G_start``
   strictly descends (so the first projected position has the LARGEST
   ``G_start`` on the genome, since minus strand reads right-to-left).

We use a hand-rolled parametrization over GTF shapes (1-exon, 2-exon,
5-exon) with deterministic random intervals (seeded for reproducibility).
This avoids pulling in ``hypothesis`` as a new dependency — the tests run
in the DRS env without extras and still cover enough diversity to catch
real regressions in the Cython fast path or the Python fallback.

Audit IDs: P0-2, ``test_projection_acceleration.py::test_64bit_and_zero_based_parity``.
"""

import random

import pytest

from drvizer.gtf_parser import GTFParser


# Deterministic seed: running the suite produces the same intervals every time.
# If you change this number, re-verify the test still exercises the expected
# cross-exon junctions on the 2-exon and 5-exon fixtures.
_SEED = 20260101


def _random_intervals(rng, transcript_length, count=5):
    """Generate ``count`` valid transcript half-open intervals.

    Each interval ``[a, b)`` satisfies ``0 <= a < b <= transcript_length``.
    Intervals may overlap; that's fine — the projection is well-defined on
    arbitrary intervals.
    """
    intervals = []
    for _ in range(count):
        a = rng.randint(0, max(0, transcript_length - 1))
        b = rng.randint(a + 1, transcript_length)
        intervals.append((a, b))
    return intervals


@pytest.mark.parametrize(
    "fixture_name,transcript_id,transcript_length,intervals_per_run",
    [
        # 1-exon, plus strand: exon [500, 699] (1-based inclusive) -> 200 bp.
        ("simple_gtf", "TX_S", 200, 5),
        # 2-exon, plus strand: exons [100, 149] + [300, 399] -> 150 bp.
        ("two_exon_gtf", "TX_T", 150, 5),
        # 5-exon, plus strand: five 150-bp exons -> 750 bp.
        ("five_exon_gtf", "TX_F", 750, 5),
        # 2-exon, minus strand: exons [1100, 1149] + [900, 949] -> 100 bp.
        ("two_exon_gtf", "TX_T_MINUS", 100, 5),
        # 5-exon, minus strand: 150 + 100 + 200 + 100 + 150 -> 700 bp.
        ("five_exon_gtf", "TX_F_MINUS", 700, 5),
    ],
)
def test_projection_length_conservation(
    request, fixture_name, transcript_id, transcript_length, intervals_per_run
):
    """Sum of projected segment lengths equals ``T_end - T_start``."""
    gtf_path = request.getfixturevalue(fixture_name)
    parser = GTFParser(str(gtf_path))
    parser.parse_gtf()

    rng = random.Random(_SEED)
    intervals = _random_intervals(rng, transcript_length, count=intervals_per_run)

    for t_start, t_end in intervals:
        _, _, segments = parser.convert_transcript_to_genomic_segments(
            transcript_id, t_start, t_end
        )
        total_length = sum(end - start for start, end in segments)
        assert total_length == t_end - t_start, (
            f"length not conserved for transcript={transcript_id} "
            f"[{t_start}, {t_end}): got segments={segments}, total={total_length}"
        )


@pytest.mark.parametrize(
    "fixture_name,transcript_id,transcript_length,interval",
    [
        # Single-exon case (no junction to cross).
        ("simple_gtf", "TX_S", 200, (10, 100)),
        # Two-exon junction: t_start < 50 <= t_end so we MUST split across exons.
        ("two_exon_gtf", "TX_T", 150, (40, 60)),
        # Two-exon inside single exon.
        ("two_exon_gtf", "TX_T", 150, (10, 40)),
        # Two-exon case where t_end sits exactly on the second exon end.
        ("two_exon_gtf", "TX_T", 150, (45, 55)),
        ("five_exon_gtf", "TX_F", 750, (300, 400)),
        # Five-exon inside one exon.
        ("five_exon_gtf", "TX_F", 750, (10, 50)),
    ],
)
def test_projection_plus_strand_monotonicity(
    request, fixture_name, transcript_id, transcript_length, interval
):
    """Plus-strand: returned ``G_start`` values strictly increase.

    The projection sorts segments by ``G_start`` ascending, so any two
    consecutive segments must satisfy ``seg[i].start < seg[i+1].start``.
    """
    del transcript_length  # only used for documentation
    gtf_path = request.getfixturevalue(fixture_name)
    parser = GTFParser(str(gtf_path))
    parser.parse_gtf()

    t_start, t_end = interval
    _, strand, segments = parser.convert_transcript_to_genomic_segments(
        transcript_id, t_start, t_end
    )
    assert strand == "+"
    assert len(segments) >= 1, (
        f"expected at least one segment for transcript={transcript_id} "
        f"[{t_start}, {t_end}); got empty list"
    )
    g_starts = [s[0] for s in segments]
    for i in range(len(g_starts) - 1):
        assert g_starts[i] < g_starts[i + 1], (
            f"plus-strand monotonicity violated: segments={segments}, "
            f"G_start values={g_starts}"
        )


@pytest.mark.parametrize(
    "fixture_name,transcript_id,transcript_length,interval",
    [
        # Two-exon minus-strand junction: 5' end (highest genomic) first.
        ("two_exon_gtf", "TX_T_MINUS", 200, (40, 60)),
        # Five-exon minus-strand, single position.
        ("five_exon_gtf", "TX_F_MINUS", 750, (10, 50)),
        # Five-exon minus-strand, junction-spanning.
        ("five_exon_gtf", "TX_F_MINUS", 750, (300, 400)),
    ],
)
def test_projection_minus_strand_monotonicity(
    request, fixture_name, transcript_id, transcript_length, interval
):
    """Minus-strand: in transcript 5'→3', ``G_start`` strictly descends.

    The projection returns segments sorted by ``G_start`` ascending, which
    for minus strand corresponds to transcript 3'→5' order. Reversing the
    list gives transcript 5'→3' order, and the first segment must have the
    LARGEST ``G_start``.
    """
    del transcript_length
    gtf_path = request.getfixturevalue(fixture_name)
    parser = GTFParser(str(gtf_path))
    parser.parse_gtf()

    t_start, t_end = interval
    _, strand, segments = parser.convert_transcript_to_genomic_segments(
        transcript_id, t_start, t_end
    )
    assert strand == "-"
    assert len(segments) >= 1

    transcript_5_to_3 = list(reversed(segments))
    g_starts = [s[0] for s in transcript_5_to_3]
    for i in range(len(g_starts) - 1):
        assert g_starts[i] > g_starts[i + 1], (
            f"minus-strand monotonicity violated (in transcript 5'->3'): "
            f"segments={segments}, G_start values (5'->3')={g_starts}"
        )

    # The first segment in transcript 5'->3' order has the largest G_start.
    assert g_starts[0] == max(g_starts), (
        f"minus-strand first (5') segment must have largest G_start: "
        f"G_starts={g_starts}, max={max(g_starts)}"
    )