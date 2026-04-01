"""
Genomic Visualization Tools
"""

__version__ = "0.1.0"

# Core classes
from .bed_parser import BEDParser, parse_bed_for_region, parse_bed_all
from .gtf_parser import GTFParser, parse_gtf_for_gene, parse_gtf_all_genes
from .visualizer import visualize_gene_transcripts, save_visualization, merge_parsers, MultiMergedParser

# Utilities
from .utils import (
    convert_to_json, 
    convert_to_csv, 
    get_transcript_stats, 
    filter_transcripts,
    merge_transcript_data
)

# Simplified API (Scanpy-style)
from .api import DrViz, quick_plot, quick_batch

__all__ = [
    # Core classes
    'BEDParser',
    'GTFParser', 
    'parse_bed_for_region',
    'parse_bed_all',
    'parse_gtf_for_gene',
    'parse_gtf_all_genes',
    'visualize_gene_transcripts',
    'save_visualization',
    'merge_parsers',
    'MultiMergedParser',
    
    # Utilities
    'convert_to_json',
    'convert_to_csv',
    'get_transcript_stats',
    'filter_transcripts',
    'merge_transcript_data',
    
    # Simplified API
    'DrViz',
    'quick_plot',
    'quick_batch'
]