"""Build-stage helpers for preparing deferred track specifications."""

from concurrent.futures import ThreadPoolExecutor

from .bed_parser import BEDParser

try:
    from .bam_parser import BAMParser
except ImportError:
    BAMParser = None


class TrackPreparationError(RuntimeError):
    """Wrap failures that occur while preparing deferred tracks."""


def prepare_track(spec, gtf_parser):
    """Instantiate and prepare one deferred track spec."""
    kind = spec["kind"]

    if kind == "bed":
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


def prepare_tracks_parallel(specs, gtf_parser):
    """Prepare all deferred track specs, preserving registration order."""
    if not specs:
        return []

    max_workers = min(len(specs), 32)
    if max_workers <= 1:
        return [prepare_track(spec, gtf_parser) for spec in specs]

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(prepare_track, spec, gtf_parser) for spec in specs]
        prepared = []
        for spec, future in zip(specs, futures):
            try:
                prepared.append(future.result())
            except Exception as exc:
                raise TrackPreparationError(
                    f"Failed to prepare track {spec.get('label', '<unknown>')}"
                ) from exc

    return prepared
