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
        self.alpha = 0.6  # default alpha value
        self.y_axis_range = y_axis_range
        self.aggregate_method = aggregate_method  # 'sum' or 'mean'
        self.parser_type = 'coverage'  # Used by visualizer to trigger area plot
        self.transcript_coord = transcript_coord
        self.gtf_parser = gtf_parser

        # Validate aggregate_method
        if self.aggregate_method not in ['sum', 'mean']:
            raise ValueError("aggregate_method must be 'sum' or 'mean'")

        # Check indices for all BAM files
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

        # Initialize global coverage array
        total_coverage = np.zeros(region_len, dtype=np.int32)

        for path in self.bam_paths:
            with pysam.AlignmentFile(path, "rb") as sam:
                # fetch uses genomic coords
                for read in sam.fetch(chrom, start, end):
                    # Filtering: Contained logic
                    if self.contained_only:
                        if read.reference_start < start or read.reference_end > end:
                            continue

                    # Spliced reads handling
                    for b_start, b_end in read.get_blocks():
                        idx_s = max(0, b_start - start)
                        idx_e = min(region_len, b_end - start)
                        if idx_s < idx_e:
                            total_coverage[idx_s:idx_e] += 1

        # Apply aggregation method
        if self.aggregate_method == 'mean' and len(self.bam_paths) > 1:
            total_coverage = total_coverage.astype(np.float64) / len(self.bam_paths)

        # X-axis generation
        x = np.arange(start, end)

        # Dynamic Binning Optimization
        if region_len > target_bins:
            bin_size = region_len // target_bins
            crop_len = (region_len // bin_size) * bin_size

            # Using mean for smoothing, max for peak preservation (mean is standard for coverage)
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

        # Get gene information and transcript list
        gene_data = self.gtf_parser.get_transcript_data(gene_identifier)
        transcripts = gene_data['transcripts']
        target_chrom = gene_data['seqname']

        # Collect all genomic intervals with coverage counts
        genomic_intervals = []  # list of (chrom, start, end, count)

        for path in self.bam_paths:
            with pysam.AlignmentFile(path, "rb") as sam:
                for transcript_info in transcripts:
                    transcript_id = transcript_info['transcript_id']

                    # Skip if this transcript is not in the BAM file
                    try:
                        # Get transcript length from BAM header
                        transcript_len = None
                        for ref in sam.header['SQ']:
                            if ref['SN'] == transcript_id:
                                transcript_len = ref['LN']
                                break
                        if transcript_len is None:
                            continue

                        # Fetch all reads aligned to this transcript
                        for read in sam.fetch(transcript_id, 0, transcript_len):
                            # Get blocks in transcript coordinates
                            blocks = read.get_blocks()
                            if not blocks:
                                continue

                            # Convert each block to genomic coordinates
                            for t_start, t_end in blocks:
                                result = self._convert_transcript_block_to_genomic(
                                    transcript_id, t_start, t_end
                                )
                                if result:
                                    chrom, g_start, g_end = result
                                    genomic_intervals.append((chrom, g_start, g_end, 1))

                    except ValueError:
                        # Transcript not found in BAM, skip
                        continue

        if not genomic_intervals:
            return np.array([]), np.array([])

        # Find the genomic region to cover
        all_starts = [iv[1] for iv in genomic_intervals]
        all_ends = [iv[2] for iv in genomic_intervals]
        region_start = min(all_starts)
        region_end = max(all_ends)
        region_len = region_end - region_start

        if region_len <= 0:
            return np.array([]), np.array([])

        # Create coverage array
        coverage = np.zeros(region_len, dtype=np.int32)

        # Accumulate coverage
        for chrom, g_start, g_end, count in genomic_intervals:
            if chrom != target_chrom:
                continue
            idx_s = max(0, g_start - region_start)
            idx_e = min(region_len, g_end - region_start)
            if idx_s < idx_e:
                coverage[idx_s:idx_e] += count

        # Apply aggregation method across multiple BAM files
        if self.aggregate_method == 'mean' and len(self.bam_paths) > 1:
            coverage = coverage.astype(np.float64) / len(self.bam_paths)

        # X-axis generation
        x = np.arange(region_start, region_end)

        # Dynamic Binning Optimization
        if region_len > target_bins:
            bin_size = region_len // target_bins
            crop_len = (region_len // bin_size) * bin_size
            coverage_binned = coverage[:crop_len].reshape(-1, bin_size).mean(axis=1)
            x_binned = x[:crop_len].reshape(-1, bin_size).mean(axis=1)
            return x_binned, coverage_binned

        return x, coverage

    def _convert_transcript_block_to_genomic(self, transcript_id, transcript_start, transcript_end):
        """
        Convert a single block from transcript coordinates to genomic coordinates.

        Args:
            transcript_id: The transcript ID
            transcript_start: Start position in transcript coordinates
            transcript_end: End position in transcript coordinates

        Returns:
            tuple: (chrom, genomic_start, genomic_end) or None if not found
        """
        if not self.gtf_parser:
            return None

        # Find transcript structure
        transcript_structure = None
        for gene_id, transcripts in self.gtf_parser.gene_transcripts.items():
            if transcript_id in transcripts:
                transcript_structure = transcripts[transcript_id]
                break

        if not transcript_structure:
            return None

        # Get exons and strand
        exons = sorted(transcript_structure['exons'], key=lambda x: x['start'])
        genomic_strand = transcript_structure['strand']
        seqname = transcript_structure['seqname']

        genomic_start = None
        genomic_end = None

        if genomic_strand == '+':
            # For positive strand
            transcript_pos = 0
            for exon in exons:
                exon_len = exon['end'] - exon['start'] + 1
                exon_end_pos = transcript_pos + exon_len

                # Check if transcript_start falls in this exon
                if transcript_pos <= transcript_start < exon_end_pos:
                    genomic_start = exon['start'] + (transcript_start - transcript_pos)

                # Check if transcript_end falls in this exon
                if transcript_pos <= transcript_end < exon_end_pos:
                    genomic_end = exon['start'] + (transcript_end - transcript_pos)
                    break

                transcript_pos = exon_end_pos

        else:  # genomic_strand == '-'
            # For negative strand, process exons in reverse order
            transcript_pos = 0
            for exon in reversed(exons):
                exon_len = exon['end'] - exon['start'] + 1
                exon_end_pos = transcript_pos + exon_len

                # Check if transcript_start falls in this exon
                if transcript_pos <= transcript_start < exon_end_pos:
                    genomic_start = exon['end'] - (transcript_start - transcript_pos)

                # Check if transcript_end falls in this exon
                if transcript_pos <= transcript_end < exon_end_pos:
                    genomic_end = exon['end'] - (transcript_end - transcript_pos)
                    break

                transcript_pos = exon_end_pos

        if genomic_start is None or genomic_end is None:
            return None

        # Ensure start <= end
        if genomic_start > genomic_end:
            genomic_start, genomic_end = genomic_end, genomic_start

        return seqname, genomic_start, genomic_end
