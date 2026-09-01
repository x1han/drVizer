"""Build-stage helpers for preparing deferred track specifications.

Phase 2.2: the parallel path is gated by an adaptive threshold. Below
the threshold we run a single-process sequential loop (no pool
allocation, no IPC). Above the threshold we open a
``ProcessPoolExecutor`` capped at ``min(len(specs), 8, cpu_count)`` so
the pool size never oversubscribes for small batches. The pool is
owned by the ``PreparedDataSource``; this module returns prepared
parsers to the caller and never stores the executor.
"""

import os
from concurrent.futures import ProcessPoolExecutor

from .bed_parser import BEDParser

try:
    from .bam_parser import BAMParser
except ImportError:
    BAMParser = None


class TrackPreparationError(RuntimeError):
    """Wrap failures that occur while preparing deferred tracks."""


def _is_process_safe_genomic_bed(spec):
    return spec["kind"] == "bed" and not spec.get("transcript_coord", False)


def _estimate_bed_record_count(spec):
    """Probe BED file size on disk to estimate record count for the threshold gate.

    We deliberately do NOT parse the file: the probe is a fast
    byte-size / stride heuristic. BED6 lines average ~50 bytes; a
    100KB BED is ~2K records, a 1MB BED is ~20K records. We use 50
    bytes/line as the canonical BED6 estimate.
    """
    files = spec.get("files", [])
    total_records = 0
    bytes_per_record = 50
    for path in files:
        try:
            total_records += os.path.getsize(path) // bytes_per_record
        except OSError:
            # Missing/unreadable file: fall back to zero so the gate
            # stays on the safe side (the parse will fail later with a
            # clear error).
            pass
    return total_records


def _apply_common_track_settings(parser, spec):
    parser.alpha = spec["alpha"]
    parser.file_colors = spec["file_colors"]
    parser.file_alphas = spec["file_alphas"]
    parser.y_axis_group = spec.get("y_axis_group")
    parser.layer_order = spec.get("layer_order", "ascending")


def _prepare_genomic_bed_track(spec):
    parser_kwargs = dict(spec.get("parser_kwargs", {}))
    parser_kwargs.pop("split_by_transcript", None)
    parser = BEDParser(
        spec["files"],
        track_label=spec["label"],
        parser_type=spec["parser_type"],
        y_axis_range=spec.get("y_axis_range"),
        transcript_coord=False,
        gtf_parser=None,
        **parser_kwargs,
    )
    parser.color = spec["color"]
    _apply_common_track_settings(parser, spec)
    parser.prepare_track(None)
    return parser


def prepare_track(spec, gtf_parser):
    """Instantiate and prepare one deferred track spec."""
    kind = spec["kind"]

    if kind == "bed":
        if _is_process_safe_genomic_bed(spec):
            return _prepare_genomic_bed_track(spec)

        parser_kwargs = dict(spec.get("parser_kwargs", {}))
        parser_kwargs.pop("split_by_transcript", None)
        parser = BEDParser(
            spec["files"],
            track_label=spec["label"],
            parser_type=spec["parser_type"],
            y_axis_range=spec.get("y_axis_range"),
            transcript_coord=spec.get("transcript_coord", False),
            gtf_parser=gtf_parser if spec.get("transcript_coord", False) else None,
            **parser_kwargs,
        )
        parser.color = spec["color"]
        _apply_common_track_settings(parser, spec)
        parser.prepare_track(gtf_parser)
        return parser

    if kind == "bam":
        if BAMParser is None:
            raise ImportError("BAM support requires pysam to be installed")

        parser_kwargs = dict(spec.get("parser_kwargs", {}))
        parser_kwargs.pop("split_by_transcript", None)
        parser = BAMParser(
            spec["files"],
            track_label=spec["label"],
            color=spec["color"],
            aggregate_method=spec.get("aggregate_method", "sum"),
            y_axis_range=spec.get("y_axis_range"),
            transcript_coord=spec.get("transcript_coord", False),
            gtf_parser=gtf_parser if spec.get("transcript_coord", False) else None,
            **parser_kwargs,
        )
        _apply_common_track_settings(parser, spec)
        if hasattr(parser, "prepare_track"):
            parser.prepare_track(gtf_parser)
        return parser

    raise ValueError(f"Unsupported track spec kind: {kind}")


def _prepare_serial(specs, gtf_parser):
    return [prepare_track(spec, gtf_parser) for spec in specs]


def _prepare_process_batch(specs, total_records_estimate, adaptive_threshold, data_source=None):
    """Prepare a batch of process-safe BED tracks.

    Under the threshold: pure sequential loop, no pool.
    At/over the threshold: open a ``ProcessPoolExecutor`` with a
    sensible worker cap. The executor is attached to ``data_source``
    so the lifecycle (shutdown on ReusableParser.close) is owned by
    the data source, not by this helper.
    """
    if not specs:
        return []

    if total_records_estimate < adaptive_threshold:
        # Below threshold: pure sequential. No pool allocation, no IPC.
        return [_prepare_genomic_bed_track(spec) for spec in specs]

    cpu_count = os.cpu_count() or 1
    max_workers = min(min(len(specs), 8), cpu_count)
    if max_workers <= 1:
        return [_prepare_genomic_bed_track(spec) for spec in specs]

    # Above threshold: spin up the pool. The PreparedDataSource is
    # the lifecycle owner; we attach the executor before the first
    # submit so that a data_source reference exists for shutdown().
    executor = ProcessPoolExecutor(max_workers=max_workers)
    if data_source is not None:
        data_source._pool = executor
        data_source._pool_used = True
    try:
        futures = [executor.submit(_prepare_genomic_bed_track, spec) for spec in specs]
        prepared = []
        for spec, future in zip(specs, futures):
            try:
                prepared.append(future.result())
            except Exception as exc:
                raise TrackPreparationError(
                    f"Failed to prepare track {spec.get('label', '<unknown>')}"
                ) from exc
    finally:
        # The pool stays alive; the PreparedDataSource owns shutdown
        # through ReusableParser.close(). We do NOT call
        # executor.shutdown() here -- the executor is supposed to
        # outlive this helper call.
        pass
    return prepared


def prepare_tracks_parallel(specs, gtf_parser, adaptive_threshold=20_000, data_source=None):
    """Prepare all deferred track specs, preserving registration order.

    Threshold-gated concurrency: the process-safe genomic BED batch
    runs in parallel only when its total record estimate exceeds
    ``adaptive_threshold``. The pool is owned by ``data_source`` when
    provided; otherwise it is local to this function and torn down
    on return (caller loses the speedup only on the next parse).
    """
    if not specs:
        return []

    process_specs = [spec for spec in specs if _is_process_safe_genomic_bed(spec)]
    serial_specs = [spec for spec in specs if not _is_process_safe_genomic_bed(spec)]

    # Estimate total records across the process-safe batch. Lazy:
    # we only probe files in the process path, never the serial
    # path, because the serial path cannot use a pool.
    total_records_estimate = sum(
        _estimate_bed_record_count(spec) for spec in process_specs
    )
    if data_source is not None:
        data_source._total_records_estimate = total_records_estimate

    # The threshold gate is per-batch, not per-file. The serial
    # batch (BAM, transcript-coord BED) always runs in the main
    # thread: BAM is excluded from the pool by panel mandate D, and
    # transcript-coord BED requires gtf_parser which cannot be
    # pickled across processes.
    process_prepared = _prepare_process_batch(
        process_specs, total_records_estimate, adaptive_threshold,
        data_source=data_source,
    )
    serial_prepared = _prepare_serial(serial_specs, gtf_parser)

    process_results = {spec["label"]: parser for spec, parser in zip(process_specs, process_prepared)}
    serial_results = {spec["label"]: parser for spec, parser in zip(serial_specs, serial_prepared)}

    prepared = []
    for spec in specs:
        label = spec["label"]
        if label in process_results:
            prepared.append(process_results[label])
        else:
            prepared.append(serial_results[label])

    return prepared
