import pysam
import numpy as np
import os

class BAMParser:
    """
    BAM Parser for drVizer.
    Supports single or multiple BAM files to generate an aggregated coverage profile.
    """
    def __init__(self, bam_paths, track_label='BAM Coverage', contained_only=True, color='steelblue', y_axis_range=None, aggregate_method='sum'):
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