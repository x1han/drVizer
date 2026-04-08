#!/usr/bin/env python3
"""
drVizer - Gene Transcript Visualization Tool
============================================

Command-line interface for visualizing gene transcript structures from GTF annotations,
with optional BED tracks.
"""

import argparse
import sys

from .api import DrViz
from .utils import convert_to_json, convert_to_csv, get_transcript_stats, filter_transcripts


def main():
    """Main function to handle command line arguments and execute drVizer."""
    parser = argparse.ArgumentParser(
        description="drVizer - Gene Transcript Visualization Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Visualize a specific gene from a GTF file:
    drvizer --gtf genes.gtf --gene ENSG00000136997 --output gene_plot.png

  Visualize one gene with BED annotations:
    drvizer --gtf genes.gtf --bed repeats.bed --gene TP53 --output merged_plot.png
        """
    )

    parser.add_argument('--gtf', nargs='+', help='Path to GTF file(s)')
    parser.add_argument('--bed', nargs='+', help='Path to BED file(s)')
    parser.add_argument('--gene', help='Gene ID or name to visualize')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--format', choices=['png', 'pdf', 'svg'], default='png', help='Output format (default: png)')
    parser.add_argument('--width', type=float, default=12, help='Figure width in inches (default: 12)')
    parser.add_argument('--height', type=float, default=8, help='Figure height in inches (default: 8)')
    parser.add_argument('--min-exons', type=int, help='Minimum number of exons')
    parser.add_argument('--max-exons', type=int, help='Maximum number of exons')
    parser.add_argument('--min-length', type=int, help='Minimum transcript length')
    parser.add_argument('--max-length', type=int, help='Maximum transcript length')
    parser.add_argument('--json', help='Save data as JSON file')
    parser.add_argument('--csv', help='Save data as CSV file')
    parser.add_argument('--stats', action='store_true', help='Print transcript statistics')

    args = parser.parse_args()

    if not args.gtf:
        parser.error("--gtf is required for the supported CLI workflow")
    if not args.gene:
        parser.error("--gene is required for the supported CLI workflow")

    try:
        print(f"Parsing GTF file(s): {args.gtf}")
        viz = DrViz().load_gtf(args.gtf)
        transcript_data = viz.get_transcript_data(args.gene)

        if any([args.min_exons, args.max_exons, args.min_length, args.max_length]):
            transcript_data = filter_transcripts(
                transcript_data,
                min_exons=args.min_exons,
                max_exons=args.max_exons,
                min_length=args.min_length,
                max_length=args.max_length
            )

        if args.stats:
            stats = get_transcript_stats(transcript_data)
            print("Transcript Statistics:")
            for key, value in stats.items():
                print(f"  {key}: {value}")

        if args.json:
            convert_to_json(transcript_data, args.json)
            print(f"Saved JSON data to {args.json}")

        if args.csv:
            convert_to_csv(transcript_data, args.csv)
            print(f"Saved CSV data to {args.csv}")

        if args.bed:
            print(f"Parsing BED file(s): {args.bed}")
            viz.add_bed_track(args.bed)

        print(f"Creating visualization: {args.output}")
        fig = viz.plot(args.gene, show=False)
        fig.set_size_inches(args.width, args.height)
        fig.savefig(args.output, format=args.format, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {args.output}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
