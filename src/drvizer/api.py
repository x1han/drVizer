"""High-level public API for drVizer."""

import matplotlib
import matplotlib.pyplot as plt
from typing import Union, List, Dict, Any

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


_configure_for_illustrator()

try:
    from .bam_parser import BAMParser
except ImportError:
    BAMParser = None


class PreparedDataSource:
    """Internal gene+track data source used by the public DrViz workflow."""

    def __init__(self, gtf_parser, tracks=None):
        self.gtf_parser = gtf_parser
        self.tracks = tracks or []

    def _build_track_entry(self, track, track_kind, track_data, track_index, transcript_id=None):
        entry = {
            'kind': track_kind,
            'data': track_data,
            'label': getattr(track, 'track_label', f'Track {track_index + 1}'),
            'color': getattr(track, 'color', 'orange'),
            'alpha': getattr(track, 'alpha', 0.8),
            'y_axis_range': getattr(track, 'y_axis_range', None),
            'file_colors': getattr(track, 'file_colors', None),
            'file_alphas': getattr(track, 'file_alphas', None),
        }
        if transcript_id is not None:
            entry['transcript_id'] = transcript_id
        return entry

    def _expand_split_tracks(self, split_mode, transcript_ids, split_track_specs):
        prepared_tracks = []
        if split_mode == 'nc':
            for transcript_id in transcript_ids:
                for track_index, track, all_track_data in split_track_specs:
                    if transcript_id not in all_track_data:
                        continue
                    track_data = all_track_data[transcript_id]
                    prepared_tracks.append(
                        self._build_track_entry(
                            track,
                            'coverage' if hasattr(track, 'get_coverage_by_transcript') else getattr(track, 'parser_type', 'distribution'),
                            track_data if hasattr(track, 'get_coverage_by_transcript') else track_data,
                            track_index,
                            transcript_id,
                        )
                    )
        else:
            for track_index, track, all_track_data in split_track_specs:
                for transcript_id in transcript_ids:
                    if transcript_id not in all_track_data:
                        continue
                    track_data = all_track_data[transcript_id]
                    prepared_tracks.append(
                        self._build_track_entry(
                            track,
                            'coverage' if hasattr(track, 'get_coverage_by_transcript') else getattr(track, 'parser_type', 'distribution'),
                            track_data if hasattr(track, 'get_coverage_by_transcript') else track_data,
                            track_index,
                            transcript_id,
                        )
                    )
        return prepared_tracks

    def get_transcript_data(self, gene_identifier, transcript_to_show=None):
        identifiers = [gene_identifier] if isinstance(gene_identifier, str) else gene_identifier
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
        all_starts = [exon['start'] for t in visible_transcripts for exon in t['exons']]
        all_ends = [exon['end'] for t in visible_transcripts for exon in t['exons']]

        if all_starts and all_ends:
            range_min, range_max = min(all_starts), max(all_ends)
            padding = int((range_max - range_min) * 0.025)
            gene_start = int(range_min - padding)
            gene_end = int(range_max + padding)
        else:
            gene_start, gene_end = 0, 1000

        prepared_tracks = []
        split_track_specs = []
        split_modes = set()
        for i, track in enumerate(self.tracks):
            split_mode = getattr(track, 'split_by_transcript', None)
            if split_mode is not None:
                split_modes.add(split_mode)
                if hasattr(track, 'get_coverage_by_transcript'):
                    split_track_specs.append((i, track, track.get_coverage_by_transcript(gene_identifier)))
                    continue
                if hasattr(track, 'get_grouped_anno_by_transcript'):
                    split_track_specs.append((i, track, track.get_grouped_anno_by_transcript(gene_identifier)))
                    continue

            if hasattr(track, 'transcript_coord') and track.transcript_coord and hasattr(track, 'get_coverage_for_transcripts'):
                payload = track.get_coverage_for_transcripts(gene_identifier)
                track_data = {'x': payload[0], 'y': payload[1]}
                track_kind = 'coverage'
            elif hasattr(track, 'get_coverage_in_region'):
                payload = track.get_coverage_in_region(target_chrom, gene_start, gene_end)
                track_data = {'x': payload[0], 'y': payload[1]}
                track_kind = 'coverage'
            elif hasattr(track, 'get_grouped_anno_in_region'):
                track_data = track.get_grouped_anno_in_region(target_chrom, gene_start, gene_end)
                track_kind = getattr(track, 'parser_type', 'distribution')
            else:
                continue

            prepared_tracks.append(self._build_track_entry(track, track_kind, track_data, i))

        if len(split_modes) > 1:
            raise ValueError("split_by_transcript must be consistent across all split tracks")

        if split_modes:
            split_mode = next(iter(split_modes))
            split_tracks = self._expand_split_tracks(split_mode, transcript_ids, split_track_specs)
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
    """Chainable public API for building transcript visualizations."""

    def __init__(self):
        self.gtf_parser = None
        self.track_specs = []
        self.track_configs = []

    def load_gtf(self, gtf_files: Union[str, List[str]]) -> 'DrViz':
        """Load one or more GTF files and reset any previously added tracks."""
        if isinstance(gtf_files, str):
            gtf_files = [gtf_files]

        self.gtf_parser = GTFParser(gtf_files)
        self.gtf_parser.parse_gtf()

        # Reset tracks and track configs for fresh start
        self.track_specs = []
        self.track_configs = []

        return self

    def add_bed_track(self, bed_files: Union[str, List[str]],
                      label: str = None,
                      color: Union[str, List[str]] = 'orange',
                      alpha: Union[float, List[float]] = 0.8,
                      parser_type: str = 'distribution',
                      y_axis_range: float = None,
                      transcript_coord: bool = False,
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
                'transcript_coord': transcript_coord,
                'split_by_transcript': split_by_transcript,
                'parser_kwargs': dict(kwargs),
            },
            {
                'label': label,
                'color': resolved_color,
                'alpha': resolved_alpha,
                'type': parser_type,
                'file_colors': colors,
                'file_alphas': alphas,
            },
        )
        return self

    def add_bam_track(self, bam_files: Union[str, List[str]],
                      label: str = "Coverage",
                      color: Union[str, List[str]] = 'steelblue',
                      alpha: Union[float, List[float]] = 0.6,
                      aggregate_method: str = 'sum',
                      y_axis_range: float = None,
                      transcript_coord: bool = False,
                      **kwargs) -> 'DrViz':
        """Add one BAM-backed coverage track to the current visualization builder.

        Args:
            bam_files: Path or list of paths to BAM files
            label: Track label for display
            color: Color for the coverage plot
            alpha: Transparency (0-1)
            aggregate_method: 'sum' or 'mean' for multiple BAM files
            y_axis_range: Fixed y-axis maximum (None for auto)
            transcript_coord: If True, treat BAM as transcript-aligned
        """
        if BAMParser is None:
            raise ImportError("BAM support requires pysam to be installed")

        files = [bam_files] if isinstance(bam_files, str) else bam_files
        colors = [color] * len(files) if isinstance(color, str) else color
        alphas = [alpha] * len(files) if isinstance(alpha, (float, int)) else alpha
        split_by_transcript = kwargs.pop('split_by_transcript', None)
        _validate_split_by_transcript(split_by_transcript, transcript_coord)
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
                'transcript_coord': transcript_coord,
                'split_by_transcript': split_by_transcript,
                'parser_kwargs': dict(kwargs),
            },
            {
                'label': label,
                'color': resolved_color,
                'alpha': resolved_alpha,
                'type': 'coverage',
                'file_colors': colors,
                'file_alphas': alphas,
            },
        )
        return self

    def build(self) -> 'ReusableParser':
        """Freeze the current builder state into a reusable plotting object."""
        if self.gtf_parser is None:
            raise ValueError("GTF file must be loaded first using load_gtf()")

        try:
            prepared_tracks = prepare_tracks_parallel(self.track_specs, self.gtf_parser)
        except TrackPreparationError as exc:
            if exc.__cause__ is not None:
                raise exc.__cause__
            raise

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

        return ReusableParser(PreparedDataSource(self.gtf_parser, ordered_tracks), self.track_configs)

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
        _configure_for_illustrator()
        parser = self.build()
        return parser.plot(gene, output=output, figsize=figsize, figfact=figfact, show=show, close=close, **kwargs)


class ReusableParser:
    """Reusable plotting object backed by one prepared DrViz data source."""

    def __init__(self, data_source, track_configs):
        self.data_source = data_source
        self.track_configs = track_configs

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

        track_labels = ['Transcripts']
        track_colors = [None]

        if self.track_configs:
            for config in self.track_configs:
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
