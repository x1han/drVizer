"""Build-stage helpers for preparing deferred track specifications."""

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

from .bed_parser import BEDParser

try:
    from .bam_parser import BAMParser
except ImportError:
    BAMParser = None


class TrackPreparationError(RuntimeError):
    """Wrap failures that occur while preparing deferred tracks."""


def _is_process_safe_genomic_bed(spec):
    return spec["kind"] == "bed" and not spec.get("transcript_coord", False)


def _prepare_genomic_bed_track(spec):
    parser_kwargs = dict(spec.get("parser_kwargs", {}))
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
    parser.alpha = spec["alpha"]
    parser.file_colors = list(spec["file_colors"])
    parser.file_alphas = list(spec["file_alphas"])
    parser.prepare_track(None)
    return parser


def prepare_track(spec, gtf_parser):
    """Instantiate and prepare one deferred track spec."""
    kind = spec["kind"]

    if kind == "bed":
        if _is_process_safe_genomic_bed(spec):
            return _prepare_genomic_bed_track(spec)

        parser_kwargs = dict(spec.get("parser_kwargs", {}))
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
        parser.alpha = spec["alpha"]
        parser.file_colors = list(spec["file_colors"])
        parser.file_alphas = list(spec["file_alphas"])
        parser.prepare_track(gtf_parser)
        return parser

    if kind == "bam":
        if BAMParser is None:
            raise ImportError("BAM support requires pysam to be installed")

        parser_kwargs = dict(spec.get("parser_kwargs", {}))
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
        parser.alpha = spec["alpha"]
        if hasattr(parser, "prepare_track"):
            parser.prepare_track(gtf_parser)
        return parser

    raise ValueError(f"Unsupported track spec kind: {kind}")


def _prepare_serial(specs, gtf_parser):
    return [prepare_track(spec, gtf_parser) for spec in specs]


def _prepare_process_batch(specs):
    if not specs:
        return []
    max_workers = min(len(specs), 32)
    if max_workers <= 1:
        return [_prepare_genomic_bed_track(spec) for spec in specs]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(_prepare_genomic_bed_track, spec) for spec in specs]
        prepared = []
        for spec, future in zip(specs, futures):
            try:
                prepared.append(future.result())
            except Exception as exc:
                raise TrackPreparationError(
                    f"Failed to prepare track {spec.get('label', '<unknown>')}"
                ) from exc
    return prepared


def prepare_tracks_parallel(specs, gtf_parser):
    """Prepare all deferred track specs, preserving registration order."""
    if not specs:
        return []

    process_specs = [spec for spec in specs if _is_process_safe_genomic_bed(spec)]
    serial_specs = [spec for spec in specs if not _is_process_safe_genomic_bed(spec)]

    process_results = {}
    serial_results = {}

    if process_specs and serial_specs:
        with ThreadPoolExecutor(max_workers=2) as executor:
            process_future = executor.submit(_prepare_process_batch, process_specs)
            serial_future = executor.submit(_prepare_serial, serial_specs, gtf_parser)
            process_prepared = process_future.result()
            serial_prepared = serial_future.result()
    else:
        process_prepared = _prepare_process_batch(process_specs)
        serial_prepared = _prepare_serial(serial_specs, gtf_parser)

    for spec, parser in zip(process_specs, process_prepared):
        process_results[spec["label"]] = parser
    for spec, parser in zip(serial_specs, serial_prepared):
        serial_results[spec["label"]] = parser

    prepared = []
    for spec in specs:
        label = spec["label"]
        if label in process_results:
            prepared.append(process_results[label])
        else:
            prepared.append(serial_results[label])

    return prepared
