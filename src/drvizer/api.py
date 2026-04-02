"""
Simplified API for drVizer - Scanpy-style interface
===================================================

This module provides a simplified, chainable API for easy transcript visualization
in Jupyter notebooks and GUI applications.
"""

import matplotlib.pyplot as plt
from typing import Union, List, Optional, Dict, Any
from pathlib import Path

from .gtf_parser import GTFParser
from .bed_parser import BEDParser
from .visualizer import visualize_gene_transcripts

try:
    from .bam_parser import BAMParser
except ImportError:
    BAMParser = None


class _BuiltParser:
    """
    Internal parser wrapper used by the high-level DrViz API.
    """

    def __init__(self, gtf_parser, other_parsers):
        self.gtf_parser = gtf_parser
        self.other_parsers = other_parsers

    def get_transcript_data(self, gene_identifier, transcript_to_show=None):
        identifiers = [gene_identifier] if isinstance(gene_identifier, str) else gene_identifier

        combined_gene_data = None
        target_chrom = None

        for ident in identifiers:
            current_data = self.gtf_parser.get_transcript_data(ident)

            if target_chrom is None:
                target_chrom = current_data['seqname']
            elif target_chrom != current_data['seqname']:
                raise ValueError(
                    f"Error: Genes must be on the same chromosome. Found {target_chrom} and {current_data['seqname']}"
                )

            if combined_gene_data is None:
                combined_gene_data = current_data
            else:
                combined_gene_data['transcripts'].extend(current_data['transcripts'])
                combined_gene_data['original_identifier'] += f", {current_data['original_identifier']}"

        seqname = target_chrom

        if transcript_to_show is not None:
            transcripts_to_use = [transcript_to_show] if isinstance(transcript_to_show, str) else transcript_to_show
            filtered_transcripts = [
                t for t in combined_gene_data['transcripts']
                if t['transcript_id'] in transcripts_to_use
            ]
            all_starts = [exon['start'] for t in filtered_transcripts for exon in t['exons']]
            all_ends = [exon['end'] for t in filtered_transcripts for exon in t['exons']]
        else:
            all_starts = [exon['start'] for t in combined_gene_data['transcripts'] for exon in t['exons']]
            all_ends = [exon['end'] for t in combined_gene_data['transcripts'] for exon in t['exons']]

        if all_starts and all_ends:
            range_min, range_max = min(all_starts), max(all_ends)
            padding = int((range_max - range_min) * 0.025)
            gene_start = int(range_min - padding)
            gene_end = int(range_max + padding)
        else:
            gene_start, gene_end = 0, 1000

        parsers_data = []
        peak_score_array = []
        for i, parser in enumerate(self.other_parsers):
            if hasattr(parser, 'get_coverage_in_region'):
                x_data, y_data = parser.get_coverage_in_region(seqname, gene_start, gene_end)
                data_payload = {'x': x_data, 'y': y_data}
            elif hasattr(parser, 'get_grouped_anno_in_region'):
                data_payload = parser.get_grouped_anno_in_region(seqname, gene_start, gene_end)
            else:
                continue

            is_score_type = getattr(parser, 'parser_type', 'distribution') == 'score'
            parsers_data.append({
                'data': data_payload,
                'parser_type': getattr(parser, 'parser_type', 'distribution'),
                'use_peak_score': is_score_type,
                'y_axis_range': getattr(parser, 'y_axis_range', None),
                'track_label': getattr(parser, 'track_label', f'Track {i+1}'),
                'color': getattr(parser, 'color', 'orange'),
                'alpha': getattr(parser, 'alpha', 0.8),
                'file_colors': getattr(parser, 'file_colors', None),
                'file_alphas': getattr(parser, 'file_alphas', None)
            })
            peak_score_array.append(is_score_type)

        combined_gene_data['separated_parsers_data'] = parsers_data
        combined_gene_data['peak_score'] = peak_score_array
        return combined_gene_data


class DrViz:
    """
    Main API class for drVizer - Simplified transcript visualization.
    
    Examples:
        >>> # Simple usage - returns merged parser for reuse
        >>> parser = DrViz().load_gtf('ref.gtf').add_bed_track('annotations.bed').build()
        >>> parser.plot('GENE1')  # Reuse the same parser
        >>> parser.plot('GENE2')  # No re-parsing needed
        
        >>> # Direct plotting (one-time use)
        >>> DrViz().load_gtf('ref.gtf').add_bed_track('annotations.bed').plot('GENE1')
        
        >>> # Advanced usage with multiple tracks
        >>> parser = DrViz().load_gtf('ref.gtf') \
        ...    .add_bed_track('te.bed', label='TE', color='red') \
        ...    .add_bed_track('methylation.bed', label='m6A', color='blue') \
        ...    .build()
        >>> parser.plot('TP53', output='tp53_viz.pdf')
        >>> parser.plot('BRCA1', output='brca1_viz.pdf')  # Reuse without re-parsing
    """
    
    def __init__(self):
        self.gtf_parser = None
        self.parsers = []  # List of additional parsers
        self.track_configs = []  # Track configuration for visualization
    
    def load_gtf(self, gtf_files: Union[str, List[str]]) -> 'DrViz':
        """
        Load GTF file(s) for transcript annotation.
        
        Args:
            gtf_files: Path to GTF file or list of paths
            
        Returns:
            self for chaining
        """
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
        """
        Add a BED track for visualization. Multiple BED files will be displayed in a SINGLE track with different colors.
        
        Args:
            bed_files: Path to BED file or list of paths
            label: Track label for legend
            color: Track color (single color or list of colors for multiple files)
            alpha: Track transparency (single alpha or list of alphas for multiple files)
            parser_type: Type of parser ('distribution' or 'score')
            y_axis_range: Y-axis range for score tracks
            transcript_coord: Whether BED files use transcript coordinates (requires GTF to be loaded)
        """
        if label is None:
            label = f'Track_{len(self.parsers) + 1}'
            
        # Handle color and alpha lists for multiple files in SAME TRACK
        files = [bed_files] if isinstance(bed_files, str) else bed_files
        colors = [color] * len(files) if isinstance(color, str) else color
        alphas = [alpha] * len(files) if isinstance(alpha, (float, int)) else alpha
        
        # Validate lengths
        if len(colors) != len(files) or len(alphas) != len(files):
            raise ValueError("Length of color and alpha lists must match number of BED files")
            
        # Create SINGLE parser for ALL files (they go in same track)
        bp = BEDParser(
            files,  # Pass all files to single parser
            track_label=label, 
            parser_type=parser_type, 
            y_axis_range=y_axis_range,
            transcript_coord=transcript_coord,
            gtf_parser=self.gtf_parser if transcript_coord else None
        )
        bp.color = colors[0] if len(set(colors)) == 1 else 'orange'  # Use first color or default orange
        bp.alpha = alphas[0] if len(set(alphas)) == 1 else 0.8  # Use common alpha or default
        bp.file_colors = colors  # Store individual file colors
        bp.file_alphas = alphas  # Store individual file alphas
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
                     **kwargs) -> 'DrViz':
        """
        Add a BAM track for visualization. Multiple BAM files will be aggregated into a single track.
        
        Args:
            bam_files: Path to BAM file or list of paths
            label: Track label for legend
            color: Track color
            alpha: Track transparency
            aggregate_method: How to aggregate multiple BAM files ('sum' or 'mean')
            y_axis_range: Y-axis range for coverage track
        """
        # Create a single BAM parser with all files (they get aggregated)
        if BAMParser is None:
            raise ImportError("BAM support requires pysam to be installed")

        bmp = BAMParser(
            bam_files, 
            track_label=label, 
            color=color,
            aggregate_method=aggregate_method,
            y_axis_range=y_axis_range
        )
        bmp.alpha = alpha  # Set alpha
        self.parsers.append(bmp)
        self.track_configs.append({
            'label': label, 
            'color': color, 
            'alpha': alpha, 
            'type': 'coverage'
        })
        return self
    
    def build(self) -> 'ReusableParser':
        """
        Build and return a reusable merged parser object.
        
        Returns:
            ReusableParser: Parser object that can plot multiple genes without re-parsing
        """
        if self.gtf_parser is None:
            raise ValueError("GTF file must be loaded first using load_gtf()")
        
        if self.parsers:
            built_parser = _BuiltParser(self.gtf_parser, self.parsers)
            return ReusableParser(built_parser, self.track_configs)
        else:
            return ReusableParser(self.gtf_parser, [])
    
    def plot(self, gene: str,
             output: str = None,
             figsize: tuple = None,
             figfact: tuple = None,
             show: bool = False,
             close: bool = False,
             **kwargs) -> plt.Figure:
        """
        Plot transcript visualization for a gene (one-time use).
        
        Args:
            gene: Gene ID or name to visualize
            output: Output file path (optional)
            figsize: Figure size (width, height). If None, use auto sizing.
            figfact: Width/height scaling factor applied to auto size when figsize is not provided.
            show: Whether to display the plot
            close: Whether to close the figure after plotting
            **kwargs: Additional arguments passed to visualize_gene_transcripts
            
        Returns:
            matplotlib Figure object
        """
        parser = self.build()
        return parser.plot(gene, output=output, figsize=figsize, figfact=figfact, show=show, close=close, **kwargs)


class ReusableParser:
    """
    Reusable parser for plotting multiple genes efficiently.
    
    This parser holds pre-parsed data and can plot multiple genes
    without re-parsing the input files.
    """
    
    def __init__(self, merged_parser, track_configs):
        self.merged_parser = merged_parser
        self.track_configs = track_configs
    
    def plot(self, gene: Union[str, List[str]],
             transcript_to_show: Union[str, List[str]] = None,
             output: str = None,
             figsize: tuple = None,
             figfact: tuple = None,
             show: bool = False,
             close: bool = False,
             **kwargs) -> plt.Figure:
        """
        Plot transcript visualization for a gene using pre-parsed data.
        
        Args:
            gene: Gene ID or name to visualize, or list of genes (must be on same chromosome)
            transcript_to_show: Specific transcript ID or list of transcript IDs to show (optional)
            output: Output file path (optional)
            figsize: Figure size (width, height). If None, use auto sizing.
            figfact: Width/height scaling factor applied to auto size when figsize is not provided.
            show: Whether to display the plot (set to False in notebooks to avoid double display)
            close: Whether to close the figure after plotting
            **kwargs: Additional arguments passed to visualize_gene_transcripts
            
        Returns:
            matplotlib Figure object
        """
        # Check if we're using a merged parser (has other parsers) or just GTF
        if hasattr(self.merged_parser, 'other_parsers') and self.merged_parser.other_parsers:
            # MultiMergedParser - supports transcript_to_show parameter
            gene_data = self.merged_parser.get_transcript_data(gene, transcript_to_show=transcript_to_show)
        else:
            # GTFParser only - check if we need to handle multiple genes
            if isinstance(gene, list):
                # For multiple genes with GTFParser, we need to merge the data manually
                # Get data for first gene to establish chromosome reference
                first_gene_data = self.merged_parser.get_transcript_data(gene[0])
                gene_data = first_gene_data.copy()
                
                # Collect all transcripts from all genes
                all_transcripts = gene_data['transcripts'].copy()
                target_chrom = gene_data['seqname']
                
                # Start building the combined identifier
                combined_identifiers = [first_gene_data.get('original_identifier', gene[0])]
                
                # Add transcripts from remaining genes
                for gene_id in gene[1:]:
                    current_data = self.merged_parser.get_transcript_data(gene_id)
                    # Verify chromosome consistency
                    if target_chrom != current_data['seqname']:
                        raise ValueError(f"Error: Genes must be on the same chromosome. Found {target_chrom} and {current_data['seqname']}")
                    all_transcripts.extend(current_data['transcripts'])
                    combined_identifiers.append(current_data.get('original_identifier', gene_id))
                
                # Update the transcripts list
                gene_data['transcripts'] = all_transcripts
                
                # Update the original_identifier to show all genes
                gene_data['original_identifier'] = ', '.join(combined_identifiers)
                # Keep the original identifier_type from the first gene
            else:
                # Single gene case
                gene_data = self.merged_parser.get_transcript_data(gene)
            # transcript_to_show filtering would need to be done manually if needed
            
        # Prepare track labels - ensure GTF track is always labeled
        track_labels = ['Transcripts']  # Default GTF label
        track_colors = [None]  # Default color for GTF
        
        # Add additional track configurations
        if self.track_configs:
            for config in self.track_configs:
                track_labels.append(config['label'])
                track_colors.append(config['color'])
        
        gene_data['track_labels'] = track_labels
        gene_data['track_colors'] = track_colors
        
        # Create visualization with optional transcript_to_show filtering
        fig = visualize_gene_transcripts(gene_data, transcript_to_show=transcript_to_show, **kwargs)

        if figsize is not None:
            fig.set_size_inches(figsize)
        elif figfact is not None:
            current_width, current_height = fig.get_size_inches()
            fig.set_size_inches((current_width * figfact[0], current_height * figfact[1]))

        # Save if output specified
        if output:
            fig.savefig(output, bbox_inches='tight', dpi=300)
            print(f"Plot saved to {output}")

        if show:
            plt.show()

        if close:
            plt.close(fig)

        return fig
    
    def get_available_genes(self) -> List[str]:
        """
        Get list of available genes in the parsed data.
        
        Returns:
            List of gene IDs/names
        """
        # Access gene information through the GTF parser
        if hasattr(self.merged_parser, 'gtf_parser'):
            return list(self.merged_parser.gtf_parser.gene_info.keys())
        elif hasattr(self.merged_parser, 'gene_info'):
            return list(self.merged_parser.gene_info.keys())
        return []

    def get_available_gene_names(self) -> List[str]:
        """
        Get list of available gene names (Symbols) in the parsed data.
        
        Returns:
            List of gene names
        """
        # Access gene name mapping through the GTF parser
        if hasattr(self.merged_parser, 'gtf_parser'):
            return list(self.merged_parser.gtf_parser.gene_name_to_id.keys())
        elif hasattr(self.merged_parser, 'gene_name_to_id'):
            return list(self.merged_parser.gene_name_to_id.keys())
        return []

    def batch_plot(self, genes: List[str], 
                   output_dir: str,
                   prefix: str = "",
                   **kwargs) -> List[str]:
        """
        Batch plot multiple genes.
        
        Args:
            genes: List of gene IDs/names
            output_dir: Directory to save plots
            prefix: Prefix for output filenames
            **kwargs: Arguments passed to plot()
            
        Returns:
            List of output file paths
        """
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
    """
    Quick batch plotting.
    """
    viz = DrViz().load_gtf(gtf_file)
    
    if bed_files:
        for bed_file in bed_files:
            viz.add_bed_track(bed_file)
    
    return viz.batch_plot(genes, output_dir, **kwargs)


# Convenience functions for quick usage
def quick_plot(gtf_file: str, 
               gene: str, 
               bed_files: List[str] = None,
               output: str = None,
               **kwargs) -> plt.Figure:
    """
    Quick one-liner for simple transcript visualization.
    
    Args:
        gtf_file: Path to GTF file
        gene: Gene ID/name
        bed_files: Optional list of BED files
        output: Output file path
        **kwargs: Additional arguments
        
    Returns:
        matplotlib Figure
    """
    viz = DrViz().load_gtf(gtf_file)
    
    if bed_files:
        for bed_file in bed_files:
            viz.add_bed_track(bed_file)
    
    return viz.plot(gene, output=output, **kwargs)