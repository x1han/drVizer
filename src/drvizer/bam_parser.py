import pysam
import numpy as np
import os

class BAMParser:
    """
    BAM Parser for drVizer.
    Supports single or multiple BAM files to generate an aggregated coverage profile.
    Supports both genomic-coordinate and transcript-coordinate BAM files.
    """
    def __init__(self, bam_paths, track_label='BAM Coverage', contained_only=True,
                 color='steelblue', y_axis_range=None, aggregate_method='sum',
                 transcript_coord=False, gtf_parser=None):
        if isinstance(bam_paths, str):
            self.bam_paths = [bam_paths]
        elif isinstance(bam_paths, list):
            self.bam_paths = bam_paths
        else:
            raise ValueError("bam_paths must be a string or a list of strings")

        self.track_label = track_label
        self.contained_only = contained_only
        self.color = color
        self.alpha = 0.6
        self.y_axis_range = y_axis_range
        self.aggregate_method = aggregate_method
        self.parser_type = 'coverage'
        self.transcript_coord = transcript_coord
        self.gtf_parser = gtf_parser

        if self.aggregate_method not in ['sum', 'mean']:
            raise ValueError("aggregate_method must be 'sum' or 'mean'")

        for path in self.bam_paths:
            if not os.path.exists(path + ".bai") and not os.path.exists(path[:-4] + ".bai"):
                print(f"Warning: Index file (.bai) not found for {path}. Attempting to index...")
                pysam.index(path)

    def get_coverage_in_region(self, chrom, start, end, target_bins=2000):
        """
        Calculate aggregated coverage array for a specific region across all BAM files.
        """
        region_len = end - start
        if region_len <= 0:
            return np.array([]), np.array([])

        total_coverage = np.zeros(region_len, dtype=np.int32)

        for path in self.bam_paths:
            with pysam.AlignmentFile(path, "rb") as sam:
                for read in sam.fetch(chrom, start, end):
                    if self.contained_only:
                        if read.reference_start < start or read.reference_end > end:
                            continue

                    for b_start, b_end in read.get_blocks():
                        idx_s = max(0, b_start - start)
                        idx_e = min(region_len, b_end - start)
                        if idx_s < idx_e:
                            total_coverage[idx_s:idx_e] += 1

        if self.aggregate_method == 'mean' and len(self.bam_paths) > 1:
            total_coverage = total_coverage.astype(np.float64) / len(self.bam_paths)

        x = np.arange(start, end)

        if region_len > target_bins:
            bin_size = region_len // target_bins
            crop_len = (region_len // bin_size) * bin_size
            coverage_binned = total_coverage[:crop_len].reshape(-1, bin_size).mean(axis=1)
            x_binned = x[:crop_len].reshape(-1, bin_size).mean(axis=1)
            return x_binned, coverage_binned

        return x, total_coverage

    def get_coverage_for_transcripts(self, gene_identifier, target_bins=2000):
        """
        Calculate coverage for transcript-aligned BAM files.
        Converts transcript coordinates to genomic coordinates and merges coverage.

        Args:
            gene_identifier: Gene ID, gene name, or transcript ID to look up
            target_bins: Target number of bins for dynamic binning

        Returns:
            tuple: (genomic_x, coverage_y) arrays in genomic coordinates
        """
        if not self.gtf_parser:
            raise ValueError("gtf_parser is required for transcript coordinate BAM files")

        gene_data = self.gtf_parser.get_transcript_data(gene_identifier)
        transcripts = gene_data['transcripts']
        target_chrom = gene_data['seqname']

        region_start = None
        region_end = None
        for transcript_info in transcripts:
            for exon in transcript_info['exons']:
                exon_start = exon['start']
                exon_end = exon['end'] + 1
                region_start = exon_start if region_start is None else min(region_start, exon_start)
                region_end = exon_end if region_end is None else max(region_end, exon_end)

        if region_start is None or region_end is None or region_end <= region_start:
            return np.array([]), np.array([])

        region_len = region_end - region_start
        coverage = np.zeros(region_len, dtype=np.int32)

        for path in self.bam_paths:
            with pysam.AlignmentFile(path, "rb") as sam:
                transcript_lengths = {ref['SN']: ref['LN'] for ref in sam.header['SQ']}

                for transcript_info in transcripts:
                    transcript_id = transcript_info['transcript_id']
                    transcript_len = transcript_lengths.get(transcript_id)
                    if transcript_len is None:
                        continue

                    for read in sam.fetch(transcript_id, 0, transcript_len):
                        blocks = read.get_blocks()
                        if not blocks:
                            continue

                        for t_start, t_end in blocks:
                            result = self.gtf_parser.convert_transcript_to_genomic_segments(
                                transcript_id, t_start, t_end
                            )
                            if not result:
                                continue

                            chrom, _, segments = result
                            if chrom != target_chrom:
                                continue

                            for g_start, g_end in segments:
                                idx_s = max(0, g_start - region_start)
                                idx_e = min(region_len, g_end - region_start)
                                if idx_s < idx_e:
                                    coverage[idx_s:idx_e] += 1

        if self.aggregate_method == 'mean' and len(self.bam_paths) > 1:
            coverage = coverage.astype(np.float64) / len(self.bam_paths)

        x = np.arange(region_start, region_end)

        if region_len > target_bins:
            bin_size = region_len // target_bins
            crop_len = (region_len // bin_size) * bin_size
            coverage_binned = coverage[:crop_len].reshape(-1, bin_size).mean(axis=1)
            x_binned = x[:crop_len].reshape(-1, bin_size).mean(axis=1)
            return x_binned, coverage_binned

        return x, coverage
