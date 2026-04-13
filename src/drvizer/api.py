"""High-level public API for drVizer."""

import matplotlib.pyplot as plt
from typing import Union, List, Dict, Any

from .gtf_parser import GTFParser
from .bed_parser import BEDParser
from .visualizer import visualize_gene_transcripts

try:
    from .bam_parser import BAMParser
except ImportError:
    BAMParser = None


class PreparedDataSource:
    """Internal gene+track data source used by the public DrViz workflow."""

    def __init__(self, gtf_parser, tracks=None):
        self.gtf_parser = gtf_parser
        self.tracks = tracks or []

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
        for i, track in enumerate(self.tracks):
            # Handle transcript-coordinate BAM tracks
            if hasattr(track, 'transcript_coord') and track.transcript_coord:
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

            prepared_tracks.append({
                'kind': track_kind,
                'data': track_data,
                'label': getattr(track, 'track_label', f'Track {i + 1}'),
                'color': getattr(track, 'color', 'orange'),
                'alpha': getattr(track, 'alpha', 0.8),
                'y_axis_range': getattr(track, 'y_axis_range', None),
                'file_colors': getattr(track, 'file_colors', None),
                'file_alphas': getattr(track, 'file_alphas', None),
            })

        combined_gene_data['prepared_tracks'] = prepared_tracks
        return combined_gene_data


class DrViz:
    """Chainable public API for building transcript visualizations."""

    def __init__(self):
        self.gtf_parser = None
        self.parsers = []
        self.track_configs = []

    def load_gtf(self, gtf_files: Union[str, List[str]]) -> 'DrViz':
        """Load one or more GTF files and reset any previously added tracks."""
        if isinstance(gtf_files, str):
            gtf_files = [gtf_files]

        self.gtf_parser = GTFParser(gtf_files)
        self.gtf_parser.parse_gtf()

        # Reset parsers and track configs for fresh start
        self.parsers = []
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
            label = f'Track_{len(self.parsers) + 1}'

        files = [bed_files] if isinstance(bed_files, str) else bed_files
        colors = [color] * len(files) if isinstance(color, str) else color
        alphas = [alpha] * len(files) if isinstance(alpha, (float, int)) else alpha

        if len(colors) != len(files) or len(alphas) != len(files):
            raise ValueError("Length of color and alpha lists must match number of BED files")

        bp = BEDParser(
            files,
            track_label=label,
            parser_type=parser_type,
            y_axis_range=y_axis_range,
            transcript_coord=transcript_coord,
            gtf_parser=self.gtf_parser if transcript_coord else None
        )
        bp.color = colors[0] if len(set(colors)) == 1 else 'orange'
        bp.alpha = alphas[0] if len(set(alphas)) == 1 else 0.8
        bp.file_colors = colors
        bp.file_alphas = alphas
        bp.parse_bed()
        self.parsers.append(bp)
        self.track_configs.append({
            'label': label,
            'color': bp.color,  # Use the actual valid color
            'alpha': bp.alpha,
            'type': parser_type,
            'file_colors': colors,
            'file_alphas': alphas
        })
        return self

    def add_bam_track(self, bam_files: Union[str, List[str]],
                      label: str = "Coverage",
                      color: str = 'steelblue',
                      alpha: float = 0.6,
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

        bmp = BAMParser(
            bam_files,
            track_label=label,
            color=color,
            aggregate_method=aggregate_method,
            y_axis_range=y_axis_range,
            transcript_coord=transcript_coord,
            gtf_parser=self.gtf_parser if transcript_coord else None
        )
        bmp.alpha = alpha
        self.parsers.append(bmp)
        self.track_configs.append({
            'label': label,
            'color': color,
            'alpha': alpha,
            'type': 'coverage'
        })
        return self

    def build(self) -> 'ReusableParser':
        """Freeze the current builder state into a reusable plotting object."""
        if self.gtf_parser is None:
            raise ValueError("GTF file must be loaded first using load_gtf()")

        return ReusableParser(PreparedDataSource(self.gtf_parser, self.parsers), self.track_configs)

    def get_transcript_data(self, gene: Union[str, List[str]], transcript_to_show: Union[str, List[str]] = None) -> Dict[str, Any]:
        """Return normalized plotting data for one gene or a same-chromosome gene list."""
        parser = self.build()
        return parser.data_source.get_transcript_data(gene, transcript_to_show=transcript_to_show)

    def plot(self, gene: str,
             output: str = None,
             figsize: tuple = None,
             figfact: tuple = None,
             show: bool = False,
             close: bool = False,
             **kwargs) -> plt.Figure:
        """Plot one gene directly from the builder without explicitly calling build()."""
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
             show: bool = False,
             close: bool = False,
             **kwargs) -> plt.Figure:
        """Plot one gene, or multiple genes on the same chromosome, from prepared data."""
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

        if close:
            plt.close(fig)

        return fig

    def get_available_genes(self) -> List[str]:
        """Return all available gene IDs from the loaded GTF data."""
        if hasattr(self.data_source, 'gtf_parser'):
            return list(self.data_source.gtf_parser.gene_info.keys())
        return []

    def get_available_gene_names(self) -> List[str]:
        """Return all available gene names from the loaded GTF data."""
        if hasattr(self.data_source, 'gtf_parser'):
            return list(self.data_source.gtf_parser.gene_name_to_id.keys())
        return []

    def batch_plot(self, genes: List[str],
                   output_dir: str,
                   prefix: str = "",
                   **kwargs) -> List[str]:
        """Plot multiple genes and save each figure into one output directory."""
        import os
        os.makedirs(output_dir, exist_ok=True)

        output_files = []
        for gene in genes:
            filename = f"{prefix}{gene}.pdf"
            output_path = os.path.join(output_dir, filename)
            self.plot(gene, output=output_path, show=False, close=True, **kwargs)
            output_files.append(output_path)
            print(f"Plotted {gene} -> {output_path}")

        return output_files



def quick_batch(gtf_file: str,
                genes: List[str],
                output_dir: str,
                bed_files: List[str] = None,
                **kwargs) -> List[str]:
    """Convenience helper for batch plotting from one GTF and optional BED tracks."""
    viz = DrViz().load_gtf(gtf_file)

    if bed_files:
        for bed_file in bed_files:
            viz.add_bed_track(bed_file)

    return viz.batch_plot(genes, output_dir, **kwargs)


def quick_plot(gtf_file: str,
               gene: str,
               bed_files: List[str] = None,
               output: str = None,
               **kwargs) -> plt.Figure:
    """Convenience helper for one-shot plotting from one GTF and optional BED tracks."""
    viz = DrViz().load_gtf(gtf_file)

    if bed_files:
        for bed_file in bed_files:
            viz.add_bed_track(bed_file)

    return viz.plot(gene, output=output, **kwargs)
