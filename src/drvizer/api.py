"""High-level public API for drVizer.

Lifecycle and pool ownership (Phase 2.2):
    Pool creation is LAZY: it only spins up on the first track-preparation
    path that exceeds the adaptive threshold. The main thread is the
    only ProcessPool consumer (GIL-safe; the pool is not a module
    singleton). Lifecycle is locked via ReusableParser.__enter__/
    __exit__/close(); DrViz.__enter__ returns the ReusableParser so
    ``with DrViz() as parser`` hands back the prepared data, not the
    builder.
"""

import collections
import os
import warnings
import matplotlib as _matplotlib

# VIZ-010: headless-safe backend lock (mirrors visualizer.py).
if "MPLBACKEND" not in os.environ:
    try:
        _matplotlib.use("Agg", force=False)
    except Exception:
        pass

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from typing import Union, List, Dict, Any, Optional, ContextManager

from .gtf_parser import GTFParser
from .bed_parser import BEDParser
from .visualizer import visualize_gene_transcripts
from ._track_build import prepare_tracks_parallel, TrackPreparationError


def _configure_for_illustrator():
    matplotlib.rcParams['font.size'] = 12
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['font.sans-serif'] = "Arial"
    matplotlib.rcParams['font.family'] = 'Arial'


def _make_unique_label(label, existing_labels):
    if label not in existing_labels:
        return label

    suffix = 1
    while True:
        candidate = f"{label}.{suffix}"
        if candidate not in existing_labels:
            return candidate
        suffix += 1


def _register_track_spec(track_specs, track_configs, spec, config):
    track_specs.append(spec)
    track_configs.append(config)


def _validate_split_by_transcript(split_by_transcript, transcript_coord):
    if split_by_transcript not in (None, 'nc', 'cn'):
        raise ValueError("split_by_transcript must be one of None, 'nc', or 'cn'")
    if split_by_transcript is not None and not transcript_coord:
        raise ValueError("split_by_transcript requires transcript_coord=True")


def _validate_layer_order(layer_order):
    if layer_order not in (None, 'ascending', 'descending'):
        raise ValueError("layer_order must be None, 'ascending', or 'descending'")


def _build_right_label_groups(prepared_tracks, start_index=0):
    groups = []
    index = start_index
    end_index = start_index + len(prepared_tracks)
    while index < end_index:
        relative_index = index - start_index
        transcript_id = prepared_tracks[relative_index].get('transcript_id')
        group_end = index
        while (
            group_end + 1 < end_index
            and prepared_tracks[group_end + 1 - start_index].get('transcript_id') == transcript_id
        ):
            group_end += 1
        groups.append({
            'transcript_id': transcript_id,
            'start_index': index,
            'end_index': group_end,
        })
        index = group_end + 1
    return groups


def _cache_get(cache, key):
    """LRU read: return the value for ``key`` and mark it recently used.

    Returns ``None`` on miss. Hits move the entry to the end of the
    OrderedDict so the next popitem(last=False) evicts the true LRU.
    """
    if key in cache:
        cache.move_to_end(key)
        return cache[key]
    return None


def _cache_put(cache, key, value, maxsize):
    """LRU write: insert ``key`` and evict the oldest entry on overflow.

    ``maxsize`` is the configured capacity; cache.put will silently allow
    growth if maxsize is <= 0 to support disabling. LRU eviction uses
    ``popitem(last=False)`` (the FIFO end == oldest entry).
    """
    if key in cache:
        cache.move_to_end(key)
        cache[key] = value
        return
    cache[key] = value
    if maxsize > 0 and len(cache) > maxsize:
        cache.popitem(last=False)


def _clear_cache(cache):
    """Drop every entry; equivalent to ``cache.clear()`` but isolates the
    import path so tests can patch this without touching the type."""
    cache.clear()


def _compress_coverage_payload(payload, aggregate_method):
    """Compress a BAM coverage payload's y-array to the target dtype.

    ``payload`` is ``(x_array, y_array)`` for a single BAM track.
    ``aggregate_method`` is ``'sum'`` or ``'mean'``; sum collapses to
    ``np.uint32`` (covering counts up to 2**32 - 1), mean collapses to
    ``np.float32`` (~7 decimal digits of precision, ample for typical
    coverage means; ultra-deep coverage above 1e7 will lose fractional
    precision but stays in plot-friendly range).
    """
    x, y = payload
    if aggregate_method == 'sum':
        return np.asarray(x), np.asarray(y, dtype=np.uint32)
    if aggregate_method == 'mean':
        return np.asarray(x), np.asarray(y, dtype=np.float32)
    return payload


_configure_for_illustrator()

try:
    from .bam_parser import BAMParser
except ImportError:
    BAMParser = None


class PreparedDataSource:
    """Internal gene+track data source used by the public DrViz workflow.

    The PreparedDataSource owns the LRU cache (P0-3 / P0-6 / P0-7
    remediation) and the lazy ProcessPool lifecycle. Cache key shape is
    a 5-tuple ``(track_id, target_id, chrom, start, end)`` where
    ``track_id`` is the unique spec label and ``target_id`` is the
    gene identifier string. The cache itself is a single
    ``collections.OrderedDict``; LRU eviction runs on insert when the
    cache reaches ``cache_maxsize``.
    """

    def __init__(self, gtf_parser, tracks=None, cache_maxsize=128):
        self.gtf_parser = gtf_parser
        self.tracks = tracks or []
        # LRU cache: OrderedDict so we get O(1) move_to_end / popitem(last=False).
        self._cache: collections.OrderedDict = collections.OrderedDict()
        self._cache_maxsize = int(cache_maxsize)
        # Pool lifecycle: created on demand, scoped to this data source.
        self._pool: Optional[Any] = None
        self._pool_used: bool = False
        # Lazy: probe total records once, on the first BED dispatch.
        self._total_records_estimate: int = 0
        # Lifecycle lock: once closed, no further inserts are allowed.
        self._is_closed: bool = False

    def _build_track_entry(self, track, track_kind, track_data, track_index, transcript_id=None):
        entry = {
            'kind': track_kind,
            'data': track_data,
            'label': getattr(track, 'track_label', f'Track {track_index + 1}'),
            'color': getattr(track, 'color', 'orange'),
            'alpha': getattr(track, 'alpha', 0.8),
            'y_axis_range': getattr(track, 'y_axis_range', None),
            'y_axis_group': getattr(track, 'y_axis_group', None),
            'file_colors': getattr(track, 'file_colors', None),
            'file_alphas': getattr(track, 'file_alphas', None),
            'layer_order': getattr(track, 'layer_order', 'ascending'),
        }
        if transcript_id is not None:
            entry['transcript_id'] = transcript_id
        return entry

    def _build_empty_track_entry(self, track, track_kind, track_index, transcript_id, region):
        """Placeholder entry for (track, transcript) combinations with no data.

        The visualizer short-circuits on entry['empty'] so np.max([]) and
        bar([], []) never run. Region metadata (chrom/start/end) is
        inherited from the parent gene so right-side label groups stay
        aligned with the prepared_tracks index space.
        """
        region_chrom, region_start, region_end = region
        entry = self._build_track_entry(track, track_kind, {}, track_index, transcript_id)
        entry['empty'] = True
        entry['series'] = []
        entry['intervals'] = []
        entry['chrom'] = region_chrom
        entry['start'] = region_start
        entry['end'] = region_end
        return entry

    def _expand_split_tracks(self, split_mode, transcript_ids, split_track_specs, region):
        prepared_tracks = []
        region_chrom, region_start, region_end = region
        if split_mode == 'nc':
            for transcript_id in transcript_ids:
                for track_index, track, all_track_data in split_track_specs:
                    track_kind = (
                        'coverage'
                        if hasattr(track, 'get_coverage_by_transcript')
                        else getattr(track, 'parser_type', 'distribution')
                    )
                    if transcript_id not in all_track_data:
                        prepared_tracks.append(
                            self._build_empty_track_entry(
                                track,
                                track_kind,
                                track_index,
                                transcript_id,
                                region,
                            )
                        )
                        continue
                    track_data = all_track_data[transcript_id]
                    prepared_tracks.append(
                        self._build_track_entry(
                            track,
                            track_kind,
                            track_data,
                            track_index,
                            transcript_id,
                        )
                    )
        else:
            for track_index, track, all_track_data in split_track_specs:
                track_kind = (
                    'coverage'
                    if hasattr(track, 'get_coverage_by_transcript')
                    else getattr(track, 'parser_type', 'distribution')
                )
                for transcript_id in transcript_ids:
                    if transcript_id not in all_track_data:
                        prepared_tracks.append(
                            self._build_empty_track_entry(
                                track,
                                track_kind,
                                track_index,
                                transcript_id,
                                region,
                            )
                        )
                        continue
                    track_data = all_track_data[transcript_id]
                    prepared_tracks.append(
                        self._build_track_entry(
                            track,
                            track_kind,
                            track_data,
                            track_index,
                            transcript_id,
                        )
                    )
        return prepared_tracks

    def _compute_track_window(self, transcripts, target_chrom, flank_pct=0.025):
        """Compute the (chrom, start, end) window covering `transcripts`.

        `transcripts` is the post-filter visible_transcripts list, so the
        window follows transcript_to_show filtering. `target_chrom` is a
        pass-through: it originates from first_gene_data['seqname'] and
        cannot be derived from `transcripts` (transcript dicts carry
        transcript_id/strand/exons/cds only), and get_transcript_data
        already needs it for the same-chromosome validation loop.

        Padding truncates toward zero via int(). The no-exon case falls
        back to (0, 1000) so downstream windows remain well-defined even
        for transcripts with empty exon lists.
        """
        all_starts = [exon['start'] for t in transcripts for exon in t['exons']]
        all_ends = [exon['end'] for t in transcripts for exon in t['exons']]
        if all_starts and all_ends:
            range_min, range_max = min(all_starts), max(all_ends)
            padding = int((range_max - range_min) * flank_pct)
            gene_start = int(range_min - padding)
            gene_end = int(range_max + padding)
        else:
            gene_start, gene_end = 0, 1000
        return target_chrom, gene_start, gene_end

    def _dispatch_track_data(self, track, track_index, gene_identifier, target_chrom, gene_start, gene_end):
        """Dispatch one non-split track and return its prepared entry list.

        Returns a list of length 0 (undispatchable) or 1 (success). Caller
        is expected to use list.extend(...) so the empty case contributes
        nothing and the populated case appends exactly one entry -- this
        preserves the loop's call order across tracks.

        The hasattr chain order matters: transcript_coord +
        get_coverage_for_transcripts shadows get_coverage_in_region for
        transcript-coord BAM tracks, so it is tested first.
        """
        if hasattr(track, 'transcript_coord') and track.transcript_coord and hasattr(track, 'get_coverage_for_transcripts'):
            track_data = self._fetch_coverage_for_transcripts(
                track, gene_identifier, track_index
            )
            track_kind = 'coverage'
        elif hasattr(track, 'get_coverage_in_region'):
            track_data = self._fetch_coverage_in_region(
                track, target_chrom, gene_start, gene_end, track_index
            )
            track_kind = 'coverage'
        elif hasattr(track, 'get_grouped_anno_in_region'):
            track_data = self._fetch_grouped_anno_in_region(
                track, target_chrom, gene_start, gene_end, track_index
            )
            track_kind = getattr(track, 'parser_type', 'distribution')
        else:
            return []
        return [self._build_track_entry(track, track_kind, track_data, track_index)]

    def _track_id(self, track, track_index):
        """Stable per-track cache identifier.

        ``track_label`` is the user-supplied unique label, which is the
        same value that flows through ``DrViz.track_specs`` and survives
        across ``build()`` cycles. We fall back to the index when the
        track has no label attribute (shouldn't happen for prepared
        tracks, but defensive against manually constructed objects).
        """
        return getattr(track, 'track_label', None) or f'Track_{track_index}'

    def _target_id(self, gene_identifier):
        """Coerce the gene identifier into a hashable string.

        The DrViz facade may pass either a string or a list of strings;
        we stringify the list for the cache key so a (gene_a, gene_b)
        multi-gene query produces one key, not N.
        """
        if isinstance(gene_identifier, (list, tuple)):
            return ','.join(str(item) for item in gene_identifier)
        return str(gene_identifier)

    def _fetch_coverage_for_transcripts(self, track, gene_identifier, track_index):
        track_id = self._track_id(track, track_index)
        target_id = self._target_id(gene_identifier)
        key = (track_id, target_id, '*', -1, -1)
        cached = _cache_get(self._cache, key)
        if cached is not None:
            return cached
        payload = track.get_coverage_for_transcripts(gene_identifier)
        track_data = {'x': payload[0], 'y': payload[1]}
        self._cache_insert(key, track_data)
        return track_data

    def _fetch_coverage_in_region(self, track, target_chrom, gene_start, gene_end, track_index):
        track_id = self._track_id(track, track_index)
        key = (
            track_id,
            '*',
            str(target_chrom),
            int(gene_start) if gene_start is not None else None,
            int(gene_end) if gene_end is not None else None,
        )
        cached = _cache_get(self._cache, key)
        if cached is not None:
            return cached
        payload = track.get_coverage_in_region(target_chrom, gene_start, gene_end)
        aggregate_method = getattr(track, 'aggregate_method', None)
        if aggregate_method in ('sum', 'mean'):
            x, y = _compress_coverage_payload(payload, aggregate_method)
        else:
            x, y = payload
        track_data = {'x': x, 'y': y}
        self._cache_insert(key, track_data)
        return track_data

    def _fetch_grouped_anno_in_region(self, track, target_chrom, gene_start, gene_end, track_index):
        track_id = self._track_id(track, track_index)
        key = (
            track_id,
            '*',
            str(target_chrom),
            int(gene_start) if gene_start is not None else None,
            int(gene_end) if gene_end is not None else None,
        )
        cached = _cache_get(self._cache, key)
        if cached is not None:
            return cached
        track_data = track.get_grouped_anno_in_region(target_chrom, gene_start, gene_end)
        self._cache_insert(key, track_data)
        return track_data

    def _fetch_grouped_anno_by_transcript(self, track, gene_identifier, track_index):
        """Cache key for split-by-transcript BED lookups.

        Split results are cacheable by ``(track_id, target_id, '*', -1, -1)``
        so two queries on the same gene skip the BED projection work.
        """
        track_id = self._track_id(track, track_index)
        target_id = self._target_id(gene_identifier)
        key = (track_id, target_id, '*', -1, -1)
        cached = _cache_get(self._cache, key)
        if cached is not None:
            return cached
        result = track.get_grouped_anno_by_transcript(gene_identifier)
        self._cache_insert(key, result)
        return result

    def _fetch_coverage_by_transcript(self, track, gene_identifier, track_index):
        track_id = self._track_id(track, track_index)
        target_id = self._target_id(gene_identifier)
        key = (track_id, target_id, '*', -1, -1)
        cached = _cache_get(self._cache, key)
        if cached is not None:
            return cached
        result = track.get_coverage_by_transcript(gene_identifier)
        self._cache_insert(key, result)
        return result

    def _cache_insert(self, key, value):
        """Insert into the LRU cache; reject post-close inserts.

        Post-close inserts raise ``RuntimeError`` because the parser
        has been torn down and a cache hit would silently serve a
        payload from a stale data source. Better to fail loud.
        """
        if self._is_closed:
            raise RuntimeError(
                "PreparedDataSource cache is closed; cannot insert new entries"
            )
        _cache_put(self._cache, key, value, self._cache_maxsize)

    def clear_cache(self):
        """Drop every cache entry. Idempotent and safe after close."""
        _clear_cache(self._cache)

    @property
    def pool(self):
        """Public accessor for the lazy ProcessPool (None if not yet provisioned)."""
        return self._pool

    def get_transcript_data(self, gene_identifier, transcript_to_show=None):
        identifiers = [gene_identifier] if isinstance(gene_identifier, str) else gene_identifier
        if len(identifiers) > 1 and any(getattr(track, 'split_by_transcript', None) is not None for track in self.tracks):
            raise ValueError("split_by_transcript does not support multiple genes")
        first_gene_data = self.gtf_parser.get_transcript_data(identifiers[0])
        combined_gene_data = first_gene_data.copy()
        all_transcripts = first_gene_data['transcripts'].copy()
        target_chrom = first_gene_data['seqname']
        combined_identifiers = [first_gene_data.get('original_identifier', identifiers[0])]

        for ident in identifiers[1:]:
            current_data = self.gtf_parser.get_transcript_data(ident)
            if target_chrom != current_data['seqname']:
                raise ValueError(
                    f"Error: Genes must be on the same chromosome. Found {target_chrom} and {current_data['seqname']}"
                )
            all_transcripts.extend(current_data['transcripts'])
            combined_identifiers.append(current_data.get('original_identifier', ident))

        combined_gene_data['transcripts'] = all_transcripts
        combined_gene_data['original_identifier'] = ', '.join(combined_identifiers)

        visible_transcripts = combined_gene_data['transcripts']
        if transcript_to_show is not None:
            transcripts_to_use = [transcript_to_show] if isinstance(transcript_to_show, str) else transcript_to_show
            filtered_transcripts = [
                t for t in combined_gene_data['transcripts']
                if t['transcript_id'] in transcripts_to_use
            ]
            if filtered_transcripts:
                combined_gene_data['transcripts'] = filtered_transcripts
                visible_transcripts = filtered_transcripts

        transcript_ids = [transcript['transcript_id'] for transcript in visible_transcripts]
        target_chrom, gene_start, gene_end = self._compute_track_window(visible_transcripts, target_chrom)

        prepared_tracks = []
        split_track_specs = []
        split_modes = set()
        for i, track in enumerate(self.tracks):
            split_mode = getattr(track, 'split_by_transcript', None)
            if split_mode is not None:
                split_modes.add(split_mode)
                if hasattr(track, 'get_coverage_by_transcript'):
                    split_track_specs.append(
                        (i, track, self._fetch_coverage_by_transcript(track, gene_identifier, i))
                    )
                    continue
                if hasattr(track, 'get_grouped_anno_by_transcript'):
                    split_track_specs.append(
                        (i, track, self._fetch_grouped_anno_by_transcript(track, gene_identifier, i))
                    )
                    continue

            prepared_tracks.extend(self._dispatch_track_data(track, i, gene_identifier, target_chrom, gene_start, gene_end))

        if len(split_modes) > 1:
            raise ValueError("split_by_transcript must be consistent across all split tracks")

        if split_modes:
            split_mode = next(iter(split_modes))
            split_tracks = self._expand_split_tracks(
                split_mode,
                transcript_ids,
                split_track_specs,
                (target_chrom, gene_start, gene_end),
            )
            if split_tracks:
                split_track_start = len(prepared_tracks)
                prepared_tracks.extend(split_tracks)
                combined_gene_data['right_label_groups'] = _build_right_label_groups(
                    split_tracks,
                    start_index=split_track_start,
                )

        combined_gene_data['prepared_tracks'] = prepared_tracks
        return combined_gene_data


class DrViz:
    """Chainable public API for building transcript visualizations.

    Phase 2.2: the builder caches its compiled ``ReusableParser`` so
    successive ``build()`` / ``plot()`` / ``get_transcript_data()``
    calls share one ``PreparedDataSource`` (and its LRU cache). Any
    mutating method (``load_gtf``, ``add_bed_track``, ``add_bam_track``)
    invalidates the cache by calling ``_mark_dirty()``.
    """

    def __init__(self, *, adaptive_threshold: int = 20_000, cache_maxsize: int = 128):
        if adaptive_threshold < 1_000 or adaptive_threshold > 1_000_000:
            warnings.warn(
                f"adaptive_threshold={adaptive_threshold} is outside the recommended "
                f"[1_000, 1_000_000] range; BED preparation may be slow.",
                stacklevel=2,
            )
        if cache_maxsize < 16 or cache_maxsize > 1024:
            warnings.warn(
                f"cache_maxsize={cache_maxsize} is outside the recommended "
                f"[16, 1024] range; memory or hit rate may suffer.",
                stacklevel=2,
            )
        self.gtf_parser = None
        self.track_specs = []
        self.track_configs = []
        self._builder_dirty: bool = True
        self._reusable_parser: Optional['ReusableParser'] = None
        self._adaptive_threshold: int = int(adaptive_threshold)
        self._cache_maxsize: int = int(cache_maxsize)

    def _mark_dirty(self) -> None:
        """Invalidate the cached ReusableParser; next build() will recompile.

        Called at the end of every mutating method. We do NOT touch
        ``PreparedDataSource._cache`` here: the cache lives with the
        parser instance, and a fresh ``build()`` constructs a fresh
        ``PreparedDataSource`` (P0-3 invariant: a new build means a new
        data source means a new cache).
        """
        self._builder_dirty = True
        self._reusable_parser = None

    def load_gtf(self, gtf_files: Union[str, List[str]]) -> 'DrViz':
        """Load one or more GTF files and reset any previously added tracks."""
        if isinstance(gtf_files, str):
            gtf_files = [gtf_files]

        self.gtf_parser = GTFParser(gtf_files)
        self.gtf_parser.parse_gtf()

        # Reset tracks and track configs for fresh start
        self.track_specs = []
        self.track_configs = []

        self._mark_dirty()
        return self

    def add_bed_track(self, bed_files: Union[str, List[str]],
                      label: str = None,
                      color: Union[str, List[str]] = 'orange',
                      alpha: Union[float, List[float]] = 0.8,
                      parser_type: str = 'distribution',
                      y_axis_range: float = None,
                      y_axis_group: str = None,
                      transcript_coord: bool = False,
                      layer_order: str = 'ascending',
                      **kwargs) -> 'DrViz':
        """Add one BED-backed track to the current visualization builder."""
        if label is None:
            label = f'Track_{len(self.track_specs) + 1}'
        label = _make_unique_label(label, {spec['label'] for spec in self.track_specs})

        files = [bed_files] if isinstance(bed_files, str) else bed_files
        colors = [color] * len(files) if isinstance(color, str) else color
        alphas = [alpha] * len(files) if isinstance(alpha, (float, int)) else alpha
        split_by_transcript = kwargs.pop('split_by_transcript', None)
        _validate_split_by_transcript(split_by_transcript, transcript_coord)
        _validate_layer_order(layer_order)
        if y_axis_group is not None and parser_type != 'score':
            raise ValueError("y_axis_group requires a numeric y-axis track")

        if len(colors) != len(files) or len(alphas) != len(files):
            raise ValueError("Length of color and alpha lists must match number of BED files")

        resolved_color = colors[0] if len(set(colors)) == 1 else 'orange'
        resolved_alpha = alphas[0] if len(set(alphas)) == 1 else 0.8

        _register_track_spec(
            self.track_specs,
            self.track_configs,
            {
                'kind': 'bed',
                'files': files,
                'label': label,
                'color': resolved_color,
                'alpha': resolved_alpha,
                'file_colors': colors,
                'file_alphas': alphas,
                'parser_type': parser_type,
                'y_axis_range': y_axis_range,
                'y_axis_group': y_axis_group,
                'transcript_coord': transcript_coord,
                'split_by_transcript': split_by_transcript,
                'layer_order': layer_order,
                'parser_kwargs': dict(kwargs),
            },
            {
                'label': label,
                'color': resolved_color,
                'alpha': resolved_alpha,
                'type': parser_type,
                'y_axis_group': y_axis_group,
                'file_colors': colors,
                'file_alphas': alphas,
            },
        )
        self._mark_dirty()
        return self

    def add_bam_track(self, bam_files: Union[str, List[str]],
                      label: str = "Coverage",
                      color: Union[str, List[str]] = 'steelblue',
                      alpha: Union[float, List[float]] = 0.6,
                      aggregate_method: str = 'sum',
                      y_axis_range: float = None,
                      y_axis_group: str = None,
                      transcript_coord: bool = False,
                      layer_order: str = 'ascending',
                      **kwargs) -> 'DrViz':
        """Add one BAM-backed coverage track to the current visualization builder.

        Args:
            bam_files: Path or list of paths to BAM files
            label: Track label for display
            color: Color for the coverage plot
            alpha: Transparency (0-1)
            aggregate_method: 'sum' or 'mean' for multiple BAM files
            y_axis_range: Fixed y-axis maximum (None for auto)
            y_axis_group: Group name for shared automatic y-axis scaling
            transcript_coord: If True, treat BAM as transcript-aligned
        """
        if BAMParser is None:
            raise ImportError("BAM support requires pysam to be installed")

        files = [bam_files] if isinstance(bam_files, str) else bam_files
        colors = [color] * len(files) if isinstance(color, str) else color
        alphas = [alpha] * len(files) if isinstance(alpha, (float, int)) else alpha
        split_by_transcript = kwargs.pop('split_by_transcript', None)
        _validate_split_by_transcript(split_by_transcript, transcript_coord)
        _validate_layer_order(layer_order)
        if len(colors) != len(files) or len(alphas) != len(files):
            raise ValueError("Length of color and alpha lists must match number of BAM files")

        resolved_color = colors[0] if len(set(colors)) == 1 else 'steelblue'
        resolved_alpha = alphas[0] if len(set(alphas)) == 1 else 0.6

        label = _make_unique_label(label, {spec['label'] for spec in self.track_specs})
        _register_track_spec(
            self.track_specs,
            self.track_configs,
            {
                'kind': 'bam',
                'files': files,
                'label': label,
                'color': resolved_color,
                'alpha': resolved_alpha,
                'file_colors': colors,
                'file_alphas': alphas,
                'aggregate_method': aggregate_method,
                'y_axis_range': y_axis_range,
                'y_axis_group': y_axis_group,
                'transcript_coord': transcript_coord,
                'split_by_transcript': split_by_transcript,
                'layer_order': layer_order,
                'parser_kwargs': dict(kwargs),
            },
            {
                'label': label,
                'color': resolved_color,
                'alpha': resolved_alpha,
                'type': 'coverage',
                'y_axis_group': y_axis_group,
                'file_colors': colors,
                'file_alphas': alphas,
            },
        )
        self._mark_dirty()
        return self

    def build(self) -> 'ReusableParser':
        """Freeze the current builder state into a reusable plotting object.

        Cached: when the builder state hasn't changed since the last
        successful ``build()`` and the cached parser hasn't been
        closed, we return the same ``ReusableParser`` instance so the
        LRU cache inside its ``PreparedDataSource`` survives.

        Pool is created LAZILY by ``prepare_tracks_parallel`` when the
        threshold gate is exceeded; under the threshold the
        ``PreparedDataSource._pool`` stays ``None``.
        """
        if self.gtf_parser is None:
            raise ValueError("GTF file must be loaded first using load_gtf()")

        if (
            not self._builder_dirty
            and self._reusable_parser is not None
            and not getattr(self._reusable_parser, '_is_closed', True)
        ):
            return self._reusable_parser

        # Build the PreparedDataSource first so prepare_tracks_parallel
        # can attach the lazy ProcessPool to it. The pool then has a
        # proper lifecycle owner for ReusableParser.close().
        data_source = PreparedDataSource(
            self.gtf_parser, tracks=[], cache_maxsize=self._cache_maxsize
        )
        try:
            prepared_tracks = prepare_tracks_parallel(
                self.track_specs,
                self.gtf_parser,
                adaptive_threshold=self._adaptive_threshold,
                data_source=data_source,
            )
        except TrackPreparationError as exc:
            if exc.__cause__ is not None:
                raise exc.__cause__
            raise

        # Re-home the prepared tracks onto the data source now that
        # the (possibly pool-spawned) build is done. This is the
        # correct ordering because prepare_tracks_parallel needed an
        # already-instantiated data_source to record the pool.
        data_source.tracks = []

        track_by_label = {getattr(track, 'track_label', None): track for track in prepared_tracks}
        ordered_tracks = []
        for spec in self.track_specs:
            label = spec["label"]
            if label not in track_by_label:
                raise RuntimeError(f"Prepared track missing for label {label}")
            track = track_by_label.pop(label)
            if spec.get('split_by_transcript') is not None:
                track.split_by_transcript = spec['split_by_transcript']
            ordered_tracks.append(track)
        if track_by_label:
            raise RuntimeError(
                f"Unexpected prepared tracks returned: {', '.join(sorted(str(label) for label in track_by_label))}"
            )

        data_source.tracks = ordered_tracks

        parser = ReusableParser(data_source, self.track_configs)
        self._reusable_parser = parser
        self._builder_dirty = False
        return parser

    def get_transcript_data(self, gene: Union[str, List[str]], transcript_to_show: Union[str, List[str]] = None) -> Dict[str, Any]:
        """Return normalized plotting data for one gene or a same-chromosome gene list."""
        parser = self.build()
        return parser.data_source.get_transcript_data(gene, transcript_to_show=transcript_to_show)

    def plot(self, gene: str,
             output: str = None,
             figsize: tuple = None,
             figfact: tuple = None,
             show: bool = True,
             close: bool = False,
             **kwargs) -> plt.Figure:
        """Plot one gene directly from the builder without explicitly calling build()."""
        parser = self.build()
        return parser.plot(gene, output=output, figsize=figsize, figfact=figfact, show=show, close=close, **kwargs)

    def __enter__(self) -> 'ReusableParser':
        """Hand the caller the prepared ReusableParser.

        ``with DrViz() as parser:`` -- the parser (not the builder) is
        the right value to render against because it owns the LRU
        cache and the pool lifecycle. This is asymmetric with
        ``ReusableParser.__enter__`` (which returns ``self``); see the
        module docstring's lifecycle section.
        """
        return self.build()

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        """Tear down the cached ReusableParser on context exit.

        Returns ``False`` -- exceptions propagate so the caller's
        ``with`` block can still see them.
        """
        if self._reusable_parser is not None:
            self._reusable_parser.close()
        return False

    def __del__(self):
        """Finalizer: warn if the builder was never closed properly.

        ResourceWarning is the canonical choice; libraries that ship
        ``__del__`` must issue it when state was not torn down.

        We do NOT close the parser from here. The common pattern is
        ``parser = DrViz()...build()`` where the inline DrViz is
        immediately GC'd but the parser is still in use; closing the
        parser from ``__del__`` would invalidate the cache the user
        is about to read from. The ReusableParser owns its own
        lifecycle; users who want auto-cleanup should use the
        ``with DrViz() as parser:`` context manager.
        """
        parser = self._reusable_parser
        if parser is not None and not getattr(parser, '_is_closed', True):
            warnings.warn(
                'DrViz held a ReusableParser that was never closed; '
                'resources may leak. Use "with DrViz() as parser:" or '
                'call parser.close() explicitly.',
                ResourceWarning,
                stacklevel=2,
            )


class ReusableParser:
    """Reusable plotting object backed by one prepared DrViz data source.

    Phase 2.2: the parser owns the lifecycle lock. ``close()`` is
    idempotent and runs in the mandated order -- clear_cache first so
    in-flight cache inserts cannot race with a closed pool, then
    ``pool.shutdown(wait=True)`` so worker subprocesses drain, then
    flip the ``_is_closed`` flag so subsequent ``plot()`` calls
    short-circuit instead of silently using a torn-down data source.
    """

    def __init__(self, data_source, track_configs):
        self.data_source = data_source
        self.track_configs = track_configs
        self._is_closed: bool = False

    def __enter__(self):
        """``with parser:`` is a self-binding; use ``parser.plot()`` inside."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        """Run ``close()`` on context exit; do not suppress exceptions."""
        self.close()
        return False

    def close(self) -> None:
        """Release cache + pool resources.

        Order is mandated by the panel (C): clear_cache FIRST so any
        in-flight insert fails fast against a closed pool, THEN
        ``pool.shutdown(wait=True)`` so worker subprocesses drain
        cleanly, THEN flip ``_is_closed``. ``ProcessPoolExecutor.shutdown``
        is idempotent so the ``is not None`` guard is sufficient even
        if ``close()`` was already called.
        """
        if self._is_closed:
            return
        # 1. clear cache first (panel mandate C)
        try:
            self.data_source.clear_cache()
        except Exception as exc:
            # Cleanup B: surface the failure as a ResourceWarning instead
            # of silently swallowing it; cache-clear failures during
            # shutdown are a leak signal the caller can act on.
            warnings.warn(
                f"ReusableParser.close(): cache clear failed: {exc!r}",
                ResourceWarning,
                stacklevel=2,
            )
        # 2. shutdown pool with wait=True (panel mandate C)
        pool = getattr(self.data_source, '_pool', None)
        if pool is not None:
            try:
                pool.shutdown(wait=True)
            except Exception as exc:
                # Cleanup B: surface pool-shutdown failure instead of
                # swallowing; worker subprocess leaks are a real risk.
                warnings.warn(
                    f"ReusableParser.close(): pool shutdown failed: {exc!r}",
                    ResourceWarning,
                    stacklevel=2,
                )
        # 3. flip the flags (both parser and data_source so cache_insert
        # on either rejects post-close)
        self._is_closed = True
        try:
            self.data_source._is_closed = True
        except Exception as exc:
            # Cleanup B: surface flag-flip failure instead of
            # swallowing; a stuck _is_closed flag would silently serve
            # stale data on the next plot call.
            warnings.warn(
                f"ReusableParser.close(): failed to flip data_source._is_closed: {exc!r}",
                ResourceWarning,
                stacklevel=2,
            )

    def clear_cache(self) -> None:
        """Public passthrough: drop the PreparedDataSource cache."""
        self.data_source.clear_cache()

    def __del__(self):
        if not self._is_closed:
            warnings.warn(
                'ReusableParser was never closed; resources may leak. '
                'Use "with viz.build() as parser:" or call '
                'parser.close() explicitly.',
                ResourceWarning,
                stacklevel=2,
            )
            try:
                self.close()
            except Exception as exc:
                # Cleanup B: surface __del__ close failures instead of
                # swallowing; if finalization fails, the GC cycle ends
                # with leaked cache and pool resources.
                warnings.warn(
                    f"ReusableParser.__del__: close() failed: {exc!r}",
                    ResourceWarning,
                    stacklevel=2,
                )

    def _build_visible_track_configs(self, prepared_tracks):
        if not self.track_configs:
            return []

        config_by_label = {
            config['label']: config
            for config in self.track_configs
            if 'label' in config
        }

        visible_configs = []
        for prepared_track in prepared_tracks:
            label = prepared_track.get('label')
            if label is None or label not in config_by_label:
                continue
            visible_configs.append(config_by_label[label])
        return visible_configs

    def plot(self, gene: Union[str, List[str]],
             transcript_to_show: Union[str, List[str]] = None,
             output: str = None,
             figsize: tuple = None,
             figfact: tuple = None,
             show: bool = True,
             close: bool = False,
             **kwargs) -> plt.Figure:
        """Plot one gene, or multiple genes on the same chromosome, from prepared data."""
        _configure_for_illustrator()
        gene_data = self.data_source.get_transcript_data(gene, transcript_to_show=transcript_to_show)

        visible_track_configs = self._build_visible_track_configs(gene_data.get('prepared_tracks', []))
        track_labels = ['Transcripts']
        track_colors = [None]

        for config in visible_track_configs:
            track_labels.append(config['label'])
            track_colors.append(config['color'])

        gene_data['track_labels'] = track_labels
        gene_data['track_colors'] = track_colors

        fig = visualize_gene_transcripts(gene_data, transcript_to_show=transcript_to_show, **kwargs)

        if figsize is not None:
            fig.set_size_inches(figsize)
        elif figfact is not None:
            current_width, current_height = fig.get_size_inches()
            fig.set_size_inches((current_width * figfact[0], current_height * figfact[1]))

        if output:
            fig.savefig(output, bbox_inches='tight', dpi=300)
            print(f"Plot saved to {output}")

        if show:
            plt.show()
        else:
            plt.close(fig)

        if close and show:
            plt.close(fig)

        return fig
