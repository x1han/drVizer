"""Phase 2.2 LRU cache + adaptive pool + lifecycle lock tests.

Eight required tests, each anchored to a panel mandate:

* test_cache_hit_invariant       (P0-3 / P0-7) cache hit skips BED parse
* test_dirty_invalidates_cache   (P0-3)        add_bed_track after build
                                                returns a fresh parser
* test_context_manager_resource_release (P1-7) ``with viz:`` drains
                                                pool + clears cache
* test_eviction_order_lru        (P0-3)        LRU eviction honours
                                                access order
* test_dtype_compression_snapshot (P0-6)       BAM sum -> uint32,
                                                mean -> float32
* test_adaptive_threshold        (panel E)     20K threshold gate
* test_pool_lazy_init            (P1-7)        no subprocess at
                                                DrViz() construction
* test_benchmark_cache_vs_no_cache (Phase 2.2.0) >=10x speedup on
                                                second-call vs first
"""

import time
import textwrap

import numpy as np
import pytest

import drvizer
import drvizer.api as api
import drvizer._track_build as track_build
from drvizer.api import DrViz, PreparedDataSource, ReusableParser
from drvizer.bed_parser import BEDParser


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def stub_bam_parser(monkeypatch):
    """A no-op BAM parser stub so testbed files don't need real BAM/BAI.

    Coverage payloads are deterministic and inspectable; the cache
    compression tests introspect ``payload[1].dtype`` directly.
    """
    class _StubBam:
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
            self._stub_x = np.arange(0, 200, dtype=np.int64)
            self._stub_y = np.arange(0, 200, dtype=np.int64)

        def get_coverage_in_region(self, chrom, start, end, target_bins=2000):
            self._stub_y = (self._stub_y + 1).astype(np.int64)
            return self._stub_x.copy(), self._stub_y.copy()

        def prepare_track(self, gtf_parser):
            return self

    monkeypatch.setattr(api, 'BAMParser', _StubBam)
    monkeypatch.setattr(track_build, 'BAMParser', _StubBam)


@pytest.fixture
def small_gtf(tmp_path):
    path = tmp_path / "small.gtf"
    path.write_text(textwrap.dedent("""\
        chr1\tsrc\texon\t100\t149\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1\tsrc\texon\t200\t249\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX2"; gene_name "G1";
    """))
    return path


@pytest.fixture
def tiny_bed(tmp_path):
    path = tmp_path / "tiny.bed"
    path.write_text("chr1\t105\t120\tpeak1\t10\t+\nchr1\t205\t220\tpeak2\t5\t+\n")
    return path


# ---------------------------------------------------------------------------
# 1. test_cache_hit_invariant
# ---------------------------------------------------------------------------

def test_cache_hit_invariant(small_gtf, tiny_bed, monkeypatch):
    """Second call to get_transcript_data on the same gene must NOT
    re-construct the BEDParser."""
    counter = {"n": 0}
    original_init = BEDParser.__init__

    def counting_init(self, *args, **kwargs):
        counter["n"] += 1
        return original_init(self, *args, **kwargs)

    monkeypatch.setattr(BEDParser, "__init__", counting_init)

    viz = DrViz().load_gtf(str(small_gtf)).add_bed_track(str(tiny_bed), label="TE")
    parser = viz.build()
    counter["n"] = 0  # reset after build() to measure only the queries

    first = parser.data_source.get_transcript_data("GENE1")
    after_first = counter["n"]
    second = parser.data_source.get_transcript_data("GENE1")
    after_second = counter["n"]

    # Parse happens only during build(), not on get_transcript_data.
    assert after_first == 0
    assert after_second == 0
    # Same payload shape from the cache hit.
    assert first["prepared_tracks"][0]["data"] == second["prepared_tracks"][0]["data"]


# ---------------------------------------------------------------------------
# 2. test_dirty_invalidates_cache
# ---------------------------------------------------------------------------

def test_dirty_invalidates_cache(small_gtf, tiny_bed, monkeypatch):
    """Adding a track after build() must produce a fresh parser on the
    next build(), with the BEDParser re-instantiated."""
    counter = {"n": 0}
    original_init = BEDParser.__init__

    def counting_init(self, *args, **kwargs):
        counter["n"] += 1
        return original_init(self, *args, **kwargs)

    monkeypatch.setattr(BEDParser, "__init__", counting_init)

    viz = DrViz().load_gtf(str(small_gtf)).add_bed_track(str(tiny_bed), label="TE")
    first_parser = viz.build()
    after_first_build = counter["n"]
    first_parser_id = id(first_parser)

    # Reuse -- must short-circuit
    second_call = viz.build()
    assert id(second_call) == first_parser_id
    assert counter["n"] == after_first_build

    # Mutate: add a new track
    viz.add_bed_track(str(tiny_bed), label="TE2")
    assert viz._builder_dirty is True
    assert viz._reusable_parser is None

    third_parser = viz.build()
    assert id(third_parser) != first_parser_id
    # Adding a track forces a full rebuild: both the old TE parser
    # and the new TE2 parser are constructed from scratch. We expect
    # +2 (one for each track in the new build), not +1.
    assert counter["n"] == after_first_build + 2


# ---------------------------------------------------------------------------
# 3. test_context_manager_resource_release
# ---------------------------------------------------------------------------

def test_context_manager_resource_release(small_gtf, tiny_bed):
    """``with DrViz() as parser:`` must close the parser, clear the
    cache, and drain the pool on exit. We don't spin up a real pool
    (the synthetic BED is tiny, threshold gate returns sequential),
    so the assertion is that ``_is_closed`` flips and the cache is
    empty after exit."""
    viz = DrViz().load_gtf(str(small_gtf)).add_bed_track(str(tiny_bed), label="TE")
    with viz as parser:
        assert isinstance(parser, ReusableParser)
        assert not parser._is_closed
        # Populate the cache by querying
        parser.data_source.get_transcript_data("GENE1")
        assert len(parser.data_source._cache) > 0

    # After the with block: closed, cache empty
    assert parser._is_closed is True
    assert parser.data_source._cache == {}
    # Post-close insert must raise (no silent staleness).
    with pytest.raises(RuntimeError, match="cache is closed"):
        parser.data_source._cache_insert(("k", "t", "chr1", 0, 1), {"x": [], "y": []})


# ---------------------------------------------------------------------------
# 4. test_eviction_order_lru
# ---------------------------------------------------------------------------

class _FakeTrack:
    """Stand-in for a track with a cacheable get_grouped_anno_in_region."""
    track_label = "TE"

    def __init__(self, idx):
        self.idx = idx
        self.call_count = 0

    def get_grouped_anno_in_region(self, chrom, start, end):
        self.call_count += 1
        return {"name": [{"chrom": chrom, "start": start, "end": end, "score": float(self.idx)}]}


def test_eviction_order_lru():
    """cache_maxsize=4, 9 distinct queries, 5 should be evicted (LRU),
    4 should remain. Re-query the last 4 -> all hits (still cached)."""
    data_source = PreparedDataSource(gtf_parser=None, cache_maxsize=4)
    fake_tracks = [_FakeTrack(i) for i in range(9)]

    # First 9 distinct queries: keys 0..8
    for i in range(9):
        data_source._fetch_grouped_anno_in_region(
            fake_tracks[i], "chr1", i * 100, i * 100 + 50, 0
        )

    # Cache must be capped at 4 entries (LRU: 5,6,7,8 survive)
    assert len(data_source._cache) == 4

    # The first 5 tracks should have been called exactly once each
    # (initial insert only). Re-querying them would miss; instead, we
    # verify the last 4 still hit.
    for i in range(5):
        assert fake_tracks[i].call_count == 1, f"track {i} called >1 time"

    # Re-query the last 4: all should hit (still in cache, not evicted)
    for i in range(5, 9):
        before = fake_tracks[i].call_count
        data_source._fetch_grouped_anno_in_region(
            fake_tracks[i], "chr1", i * 100, i * 100 + 50, 0
        )
        assert fake_tracks[i].call_count == before, f"query {i} should be a hit"


# ---------------------------------------------------------------------------
# 5. test_dtype_compression_snapshot
# ---------------------------------------------------------------------------

def test_dtype_compression_snapshot():
    """BAM coverage payload y must collapse to np.uint32 for
    aggregate_method='sum' and np.float32 for 'mean' (P0-6)."""
    class _BamSum:
        track_label = "BamSum"
        aggregate_method = 'sum'
        parser_type = 'coverage'

        def get_coverage_in_region(self, chrom, start, end, target_bins=2000):
            return np.array([start, start + 1], dtype=np.int64), np.array([1, 2], dtype=np.int64)

    class _BamMean:
        track_label = "BamMean"
        aggregate_method = 'mean'
        parser_type = 'coverage'

        def get_coverage_in_region(self, chrom, start, end, target_bins=2000):
            return np.array([start, start + 1], dtype=np.int64), np.array([1, 2], dtype=np.int64)

    data_source = PreparedDataSource(gtf_parser=None)
    data_source._fetch_coverage_in_region(_BamSum(), "chr1", 0, 200, 0)
    sum_cached = data_source._cache[("BamSum", "*", "chr1", 0, 200)]
    assert sum_cached["y"].dtype == np.uint32
    assert int(sum_cached["y"].max()) < 2**32

    data_source._fetch_coverage_in_region(_BamMean(), "chr1", 0, 200, 0)
    mean_cached = data_source._cache[("BamMean", "*", "chr1", 0, 200)]
    assert mean_cached["y"].dtype == np.float32


# ---------------------------------------------------------------------------
# 6. test_adaptive_threshold
# ---------------------------------------------------------------------------

def test_adaptive_threshold(tmp_path, monkeypatch):
    """Below the threshold: no pool. Above the threshold: pool."""
    gtf = tmp_path / "gtf.gtf"
    gtf.write_text(textwrap.dedent("""\
        chr1\tsrc\texon\t100\t149\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1\tsrc\texon\t200\t249\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX2"; gene_name "G1";
    """))

    # --- Sub-threshold: small BED, no pool ---
    small_bed = tmp_path / "small.bed"
    small_bed.write_text("chr1\t105\t120\tp\t10\t+\n")
    viz = DrViz(adaptive_threshold=20_000)
    viz.load_gtf(str(gtf)).add_bed_track(str(small_bed), label="TE")
    with viz as parser:
        assert parser.data_source._pool is None

    # --- Above threshold: large BED, pool is opened ---
    pool_constructions = {"n": 0}

    class _FakePool:
        def __init__(self, max_workers=None):
            pool_constructions["n"] += 1
            self.max_workers = max_workers

        def submit(self, fn, *args, **kwargs):
            class _Future:
                def result(self_inner):
                    return fn(*args, **kwargs)
            return _Future()

        def shutdown(self, wait=True):
            pass

    import drvizer._track_build as tb
    monkeypatch.setattr(tb, "ProcessPoolExecutor", _FakePool)

    # >1MB BED: ~1_200_000 / 50 = ~24K records, over the 20K threshold.
    bed = tmp_path / "big.bed"
    with bed.open("wt") as fh:
        for i in range(25_000):
            fh.write(f"chr1\t{100 + i*2}\t{150 + i*2}\tp{i}\t10\t+\n")
    bed2 = tmp_path / "big2.bed"
    bed2.write_text(bed.read_text())

    viz2 = DrViz(adaptive_threshold=20_000)
    viz2.load_gtf(str(gtf)).add_bed_track(str(bed), label="A").add_bed_track(str(bed2), label="B")
    with viz2 as parser:
        pass

    assert pool_constructions["n"] >= 1


# ---------------------------------------------------------------------------
# 7. test_pool_lazy_init
# ---------------------------------------------------------------------------

def test_pool_lazy_init(tmp_path):
    """DrViz() construction, load_gtf(), add_bed_track() must NOT
    spawn a subprocess. The pool only materialises inside build()
    when the threshold gate says so.

    We test the Python-level invariant: data_source._pool is None
    until build() opens one.
    """
    gtf = tmp_path / "gtf.gtf"
    gtf.write_text(textwrap.dedent("""\
        chr1\tsrc\texon\t100\t149\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
    """))
    bed = tmp_path / "bed.bed"
    bed.write_text("chr1\t105\t120\tp\t10\t+\n")

    viz = DrViz()
    # No data_source exists yet (it lives on the ReusableParser,
    # which doesn't exist until build()).
    assert getattr(viz, "_reusable_parser", None) is None

    viz.load_gtf(str(gtf))
    assert getattr(viz, "_reusable_parser", None) is None

    viz.add_bed_track(str(bed), label="TE")
    assert getattr(viz, "_reusable_parser", None) is None

    parser = viz.build()
    # Sub-threshold BED => no pool, even after build.
    assert parser.data_source._pool is None


# ---------------------------------------------------------------------------
# 8. test_benchmark_cache_vs_no_cache
# ---------------------------------------------------------------------------

def test_benchmark_cache_vs_no_cache():
    """Second call (cache hit) must be at least 10x faster than the
    first (cache miss). We use a deliberately-slow parse via
    monkey-patched ``get_grouped_anno_in_region`` so the timing
    signal is real, not a timing-resolution artifact."""

    class _SlowTrack:
        track_label = "TE"
        parser_type = "distribution"

        def __init__(self):
            self.calls = 0

        def get_grouped_anno_in_region(self, chrom, start, end):
            self.calls += 1
            # ~50ms of work; well above the cache lookup ceiling
            # (~10us).
            time.sleep(0.05)
            return {"name": [{"chrom": chrom, "start": start, "end": end}]}

    data_source = PreparedDataSource(gtf_parser=None, cache_maxsize=128)
    slow = _SlowTrack()

    t0 = time.perf_counter()
    data_source._fetch_grouped_anno_in_region(slow, "chr1", 100, 200, 0)
    first_elapsed = time.perf_counter() - t0
    assert slow.calls == 1

    t0 = time.perf_counter()
    data_source._fetch_grouped_anno_in_region(slow, "chr1", 100, 200, 0)
    second_elapsed = time.perf_counter() - t0
    # Cache hit: no sleep, no call into the track.
    assert slow.calls == 1

    # Real wall-time ratio: hit should be orders of magnitude faster.
    # We assert >=10x per the plan; in practice this is 100x+.
    assert second_elapsed * 10 < first_elapsed, (
        f"cache hit ({second_elapsed*1e3:.2f}ms) is not >=10x faster "
        f"than miss ({first_elapsed*1e3:.2f}ms)"
    )
