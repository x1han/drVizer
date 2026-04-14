"""
BED Parser for drVizer
======================

This module provides functions to parse BED files and extract annotation elements 
information for visualization along with gene transcripts.
"""

import pandas as pd
from collections import defaultdict
import gzip


class BEDParser:
    """
    A class to parse BED files and extract annotation information for visualization.
    """
    
    def __init__(self, bed_file_path, transcript_coord=False, gtf_parser=None, parser_type='distribution', y_axis_range=None, track_label='Track'):
        """
        Initialize the BEDParser with the path to the BED file or list of BED files.
        
        Args:
            bed_file_path (str or list): Path to the BED file or list of BED file paths
            transcript_coord (bool): Whether the BED file uses transcript coordinates instead of genomic coordinates
            gtf_parser (GTFParser): GTF parser instance to use for coordinate conversion if transcript_coord is True
            parser_type (str): Type of parser ('distribution' or 'score')
            y_axis_range (float, optional): Maximum value for y-axis when plotting as barplot (for score parsers)
            track_label (str): Label for this track in visualization, default is 'Track'
        """
        if isinstance(bed_file_path, str):
            self.bed_file_paths = [bed_file_path]
        elif isinstance(bed_file_path, list):
            self.bed_file_paths = bed_file_path
        else:
            raise ValueError("bed_file_path must be a string or list of strings")
        self.transcript_coord = transcript_coord
        self.gtf_parser = gtf_parser
        self.parser_type = parser_type  # 'distribution' or 'score'
        self.y_axis_range = y_axis_range  # For score parsers, defines max y value
        self.track_label = track_label
        self.alpha = 0.8  # default alpha value
        self.color = 'orange'  # default color (valid matplotlib color)
        self.file_colors = None  # Store colors for different files
        self.file_alphas = None  # Store alphas for different files
        self.anno_data = defaultdict(list)
        self._parsed = False
        
    def convert_transcript_to_genomic_coords(self, transcript_id, transcript_start, transcript_end):
        """
        Convert transcript coordinates to genomic coordinates using GTF parser.
        """
        if not self.gtf_parser:
            raise ValueError("gtf_parser is required for transcript coordinate conversion")

        result = self.gtf_parser.convert_transcript_to_genomic(transcript_id, transcript_start, transcript_end)
        if not result:
            return None, None, None, None

        chrom, genomic_start, genomic_end = result
        _, transcript_structure = self.gtf_parser.find_transcript(transcript_id)
        genomic_strand = transcript_structure['strand'] if transcript_structure else None
        return chrom, genomic_start, genomic_end, genomic_strand

        

    

        
    def parse_bed(self, chrom=None, start=None, end=None):
        """
        Parse the BED file(s) and extract annotation information.
        If region (chrom, start, end) is provided, only parse annotations within that region.
        
        Args:
            chrom (str, optional): Chromosome name to filter
            start (int, optional): Start position to filter
            end (int, optional): End position to filter
            
        Returns:
            dict: Dictionary containing annotation information within the specified region or all annotations
        """
        # Parse if not already parsed, or if a specific region is requested, or if we need to repopulate anno_data
        # We need to re-parse if:
        # 1. Not parsed yet
        # 2. A specific region is requested (chrom, start, end are all provided)
        # 3. anno_data is empty (might happen if previous parse was with a specific region)
        if not self._parsed or (chrom and start and end) or len(self.anno_data) == 0:
            # Clear existing data if we're parsing a specific region or if we need to re-parse everything
            if (chrom and start and end) or not self._parsed:
                self.anno_data.clear()
            self._process_bed_file(chrom, start, end)
            # Only mark as parsed if we're parsing all data (no specific region)
            if not (chrom and start and end):
                self._parsed = True
            
        return dict(self.anno_data)
    
    def _process_bed_file(self, chrom=None, start=None, end=None):
        """
        Process BED file and extract annotation information.
        
        Args:
            chrom (str, optional): Chromosome name to filter
            start (int, optional): Start position to filter
            end (int, optional): End position to filter
        """
        # Process each BED file
        for bed_file_path in self.bed_file_paths:
            # Determine if file is gzipped
            open_func = gzip.open if bed_file_path.endswith('.gz') else open
            
            try:
                with open_func(bed_file_path, 'rt') as f:
                    for line in f:
                        # Skip comment lines
                        if line.startswith('#') or line.strip() == '':
                            continue
                        
                        parts = line.strip().split('\t')
                        # BED files need at least 3 required columns: chrom, start, end
                        # Optional columns include: name (4th), score (5th), strand (6th)
                        # These are at positions 0, 1, 2, 3, 4, 5 (0-based) or 1, 2, 3, 4, 5, 6 (1-based)
                        if len(parts) < 4:  # Require at least chrom, start, end
                            # Skip lines with insufficient columns
                            continue
                            
                        # Parse the BED record with required column handling
                        record = {}
                        record['chrom'] = parts[0]
                        try:
                            record['start'] = int(parts[1])
                            record['end'] = int(parts[2])
                        except (ValueError, IndexError):
                            # Skip records with invalid numeric values
                            continue
                            
                        # Handle name (4th column, 0-based index 3)
                        record['name'] = parts[3] if len(parts) > 3 else '.'
                        
                        # Handle score (5th column, 0-based index 4)
                        if len(parts) >= 5:
                            try:
                                record['score'] = float(parts[4])
                            except ValueError:
                                record['score'] = 0.0
                        else:
                            record['score'] = 0.0
                        
                        # Strand is in the 6th column (0-based index 5)
                        record['strand'] = parts[5] if len(parts) > 5 else '.'
                            
                        # Handle thickStart and thickEnd (7th and 8th columns)
                        if len(parts) >= 7:
                            try:
                                record['thickStart'] = int(parts[6])
                            except ValueError:
                                record['thickStart'] = record['start']
                        else:
                            record['thickStart'] = record['start']
                            
                        if len(parts) >= 8:
                            try:
                                record['thickEnd'] = int(parts[7])
                            except ValueError:
                                record['thickEnd'] = record['end']
                        else:
                            record['thickEnd'] = record['end']
                            
                        # Handle itemRgb (9th column)
                        if len(parts) >= 9:
                            record['itemRgb'] = parts[8]
                        else:
                            record['itemRgb'] = '0'
                        
                        # Filter by region if specified
                        if chrom and start and end:
                            if (record['chrom'] != chrom or 
                                record['end'] < start or 
                                record['start'] > end):
                                continue
                        
                        # If transcript coordinates are used, convert to genomic coordinates
                        if self.transcript_coord:
                            if not self.gtf_parser:
                                raise ValueError("gtf_parser is required when transcript_coord is True")
                            
                            # Extract transcript ID from the first column
                            transcript_id = record['chrom']  # In transcript coordinate BED, chrom field is actually transcript ID
                            transcript_start = record['start']
                            transcript_end = record['end']
                            transcript_strand = record['strand']
                            
                            # Convert transcript coordinates to genomic coordinates
                            chrom, genomic_start, genomic_end, genomic_strand = \
                                self.convert_transcript_to_genomic_coords(transcript_id, transcript_start, transcript_end)
                            
                            if chrom is None:
                                # If transcript not found in GTF, skip this record
                                continue
                            
                            # Update record with genomic coordinates
                            record['chrom'] = chrom
                            record['start'] = genomic_start
                            record['end'] = genomic_end
                            record['strand'] = genomic_strand
                            
                            # Store annotation data grouped by chromosome
                            self.anno_data[record['chrom']].append(record)
                        else:
                            # Store annotation data grouped by chromosome
                            self.anno_data[record['chrom']].append(record)
                    
            except Exception as e:
                raise ValueError(f"Error reading BED file {bed_file_path}: {e}")
    
    def get_anno_in_region(self, chrom, start, end):
        """
        Get annotation elements within a specific genomic region.
        
        Args:
            chrom (str): Chromosome name
            start (int): Start position
            end (int): End position
            
        Returns:
            list: List of annotation elements within the specified region
        """
        # Parse the BED file if not already parsed
        if not self._parsed:
            self.parse_bed()
            
        # Filter annotations within the region
        anno_in_region = []
        if chrom in self.anno_data:
            for anno in self.anno_data[chrom]:
                # Check if annotation overlaps with the region
                if anno['end'] >= start and anno['start'] <= end:
                    anno_in_region.append(anno)
                    
        return anno_in_region
    
    def get_all_chromosomes(self):
        """
        Get all chromosomes in the BED file.
        
        Returns:
            list: List of chromosome names
        """
        if not self._parsed:
            self.parse_bed()
        return list(self.anno_data.keys())
    
    def get_grouped_anno_in_region(self, chrom, start, end):
        """
        Get annotation elements within a specific genomic region, grouped by name.
        
        Args:
            chrom (str): Chromosome name
            start (int): Start position
            end (int): End position
            
        Returns:
            dict: Dictionary with annotation names as keys and list of annotation elements as values
        """
        # Parse the BED file if not already parsed
        if not self._parsed:
            self.parse_bed()
            
        # Group annotations by name within the region
        grouped_anno = {}
        if chrom in self.anno_data:
            for anno in self.anno_data[chrom]:
                # Ensure start <= end for proper region checking
                anno_start, anno_end = anno['start'], anno['end']
                if anno_start > anno_end:
                    anno_start, anno_end = anno_end, anno_start
                
                # Check if annotation overlaps with the region
                if anno_end >= start and anno_start <= end:
                    name = anno['name']
                    if name not in grouped_anno:
                        grouped_anno[name] = []
                    grouped_anno[name].append(anno)
                    
        return grouped_anno


def parse_bed_for_region(bed_file_path, chrom, start, end):
    """
    Parse BED file and extract annotation information for a specific genomic region.
    
    Args:
        bed_file_path (str): Path to the BED file
        chrom (str): Chromosome name
        start (int): Start position
        end (int): End position
        
    Returns:
        list: List of annotation elements within the specified region
    """
    parser = BEDParser(bed_file_path)
    return parser.get_anno_in_region(chrom, start, end)


def parse_bed_all(bed_file_path):
    """
    Parse BED file and extract all annotation information.
    
    Args:
        bed_file_path (str): Path to the BED file
        
    Returns:
        dict: Dictionary of all annotation elements grouped by chromosome
    """
    parser = BEDParser(bed_file_path)
    parser.parse_bed()
    return dict(parser.anno_data)