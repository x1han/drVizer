#!/usr/bin/env python3
"""Compatibility wrapper for running the CLI from a repository checkout."""

import argparse
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from gtf_parser import GTFParser
from bed_parser import BEDParser
from visualizer import visualize_gene_transcripts, save_visualization, merge_parsers
from utils import convert_to_json, convert_to_csv, get_transcript_stats, filter_transcripts


def main():
    """Main function to handle command line arguments and execute drVizer."""
    parser = argparse.ArgumentParser(
        description="drVizer - Gene Transcript Visualization Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Parse and visualize a specific gene from GTF file:
    python drvizer_cli.py --gtf genes.gtf --gene ENSG00000136997 --output gene_plot.png

  Parse BED file for a genomic region:
    python drvizer_cli.py --bed repeats.bed --region chr1:1000000-2000000 --output te_plot.png

  Merge GTF and BED data for comprehensive visualization:
    python drvizer_cli.py --gtf genes.gtf --bed repeats.bed --gene TP53 --output merged_plot.png
        """
    )

    parser.add_argument('--gtf', nargs='+', help='Path to GTF file(s)')
    parser.add_argument('--bed', nargs='+', help='Path to BED file(s)')
    parser.add_argument('--gene', help='Gene ID or name to visualize')
    parser.add_argument('--region', help='Genomic region in format chr:start-end')
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

    if not args.gtf and not args.bed:
        parser.error("Either --gtf or --bed must be specified")
    if args.gene and args.region:
        parser.error("Cannot specify both --gene and --region")
    if args.gene and not args.gtf:
        parser.error("--gene requires --gtf to be specified")
    if args.region and not args.bed:
        parser.error("--region requires --bed to be specified")

    try:
        gtf_parser = None
        transcript_data = None

        if args.gtf:
            print(f"Parsing GTF file(s): {args.gtf}")
            gtf_parser = GTFParser(args.gtf)

            if args.gene:
                gtf_parser.parse_gtf()
                transcript_data = gtf_parser.get_transcript_data(args.gene)

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

        bed_parser = None
        bed_data = None

        if args.bed:
            print(f"Parsing BED file(s): {args.bed}")
            bed_parser = BEDParser(args.bed)
            bed_parser.parse_bed()

            if args.region:
                chrom, pos = args.region.split(':')
                start, end = map(int, pos.split('-'))
                bed_data = bed_parser.get_anno_in_region(chrom, start, end)
                print(f"Found {len(bed_data)} elements in region {args.region}")

        if transcript_data or bed_data:
            print(f"Creating visualization: {args.output}")

            if transcript_data and bed_data and gtf_parser and bed_parser:
                merged_parser = merge_parsers(gtf_parser, bed_parser)
                merged_data = merged_parser.get_transcript_data(args.gene)
                fig = visualize_gene_transcripts(merged_data, transcript_width=args.width / 12)
            elif transcript_data:
                fig = visualize_gene_transcripts(transcript_data, transcript_width=args.width / 12)
            else:
                print("BED-only visualization not yet implemented")
                return

            fig.set_size_inches(args.width, args.height)
            save_visualization(fig, args.output, format=args.format)
            print(f"Visualization saved to {args.output}")
        else:
            print("No data to visualize")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
