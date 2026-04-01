"""
GTF Parser for drVizer
======================

This module provides functions to parse GTF files and extract transcript information
for a specific gene, which can then be used for visualization.
"""

import pandas as pd
from collections import defaultdict
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip


class GTFParser:
    """
    A class to parse GTF files and extract transcript information for visualization.
    Optimized for handling multiple genes efficiently by caching parsed data.
    """
    
    def __init__(self, gtf_file_path, track_label='Transcript'):
        """
        Initialize the GTFParser with the path to the GTF file or list of GTF files.
        
        Args:
            gtf_file_path (str or list): Path to the GTF file or list of GTF file paths
            track_label (str): Label for this track in visualization, default is 'Transcript'
        """
        if isinstance(gtf_file_path, str):
            self.gtf_file_paths = [gtf_file_path]
        elif isinstance(gtf_file_path, list):
            self.gtf_file_paths = gtf_file_path
        else:
            raise ValueError("gtf_file_path must be a string or list of strings")
            
        self.track_label = track_label
        self.gene_transcripts = defaultdict(dict)
        self.gene_info = {}
        self._parsed = False
        self.gene_name_to_id = {}  # Map gene names to gene IDs
        self._transcripts_from_previous_files = set()  # Track transcripts from previous files
        
    def parse_gtf(self, gene_id=None):
        """
        Parse the GTF file(s) and extract transcript information.
        If called without gene_id, parses all genes and caches the data.
        If called with gene_id and data is already cached, returns cached data.
        
        Args:
            gene_id (str, optional): Specific gene ID to parse. If None, parse all genes.
            
        Returns:
            dict: Dictionary containing transcript information
        """
        # If we've already parsed everything, just return the requested data
        if self._parsed:
            if gene_id:
                return {gene_id: self.gene_transcripts[gene_id]} if gene_id in self.gene_transcripts else {}
            else:
                return dict(self.gene_transcripts)
        
        # If we're looking for a specific gene and already have it parsed, return it
        if gene_id and gene_id in self.gene_transcripts:
            return {gene_id: self.gene_transcripts[gene_id]}
            
        # Process each GTF file in order
        for i, gtf_file_path in enumerate(self.gtf_file_paths):
            # For files after the first, copy current transcripts to previous file tracking
            if i > 0:
                self._update_transcripts_from_previous_files()
                
            self._process_gtf_file_optimized(gtf_file_path, gene_id)
            
        # Mark as fully parsed if we didn't specify a gene or if we parsed all genes
        if not gene_id:
            self._parsed = True
            
        if gene_id:
            return {gene_id: self.gene_transcripts[gene_id]} if gene_id in self.gene_transcripts else {}
        else:
            return dict(self.gene_transcripts)
    
    def _update_transcripts_from_previous_files(self):
        """
        Update the set of transcripts from previous files.
        This is called before processing each file (except the first) to track what transcripts
        were already present before processing the current file.
        """
        for gene_id, transcripts in self.gene_transcripts.items():
            for transcript_id in transcripts.keys():
                transcript_key = f"{gene_id}:{transcript_id}"
                self._transcripts_from_previous_files.add(transcript_key)
    
    def _process_gtf_file_optimized(self, gtf_file_path, gene_id=None):
        """
        Process GTF file using optimized chunked approach with parallel processing.
        
        Args:
            gtf_file_path (str): Path to the GTF file
            gene_id (str, optional): Specific gene ID to parse
        """
        # Define column names for GTF file
        columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        
        # Determine if file is gzipped
        open_func = gzip.open if gtf_file_path.endswith('.gz') else open
        
        # Process file in chunks for memory efficiency
        chunk_size = 10000
        chunk = []
        
        try:
            with open_func(gtf_file_path, 'rt') as f:
                for line in f:
                    # Skip comment lines
                    if line.startswith('#') or line.strip() == '':
                        continue
                    
                    chunk.append(line)
                    
                    # Process chunk when it reaches the specified size
                    if len(chunk) >= chunk_size:
                        self._process_chunk(chunk, columns, gene_id)
                        chunk = []  # Reset chunk
                
                # Process remaining lines in final chunk
                if chunk:
                    self._process_chunk(chunk, columns, gene_id)
                    
        except Exception as e:
            raise ValueError(f"Error reading GTF file {gtf_file_path}: {e}")
    
    def _process_chunk(self, chunk_lines, columns, gene_id=None):
        """
        Process a chunk of GTF lines.
        
        Args:
            chunk_lines (list): List of GTF file lines
            columns (list): Column names for GTF format
            gene_id (str, optional): Specific gene ID to parse
        """
        # Parse lines into a temporary dataframe for this chunk
        chunk_data = []
        for line in chunk_lines:
            parts = line.strip().split('\t')
            if len(parts) >= 9 and (parts[2] == 'exon' or parts[2] == 'CDS'):  # Process exon and CDS features
                chunk_data.append(parts[:8] + [parts[8]])  # Ensure we have 9 columns
        
        if not chunk_data:
            return
            
        # Create DataFrame for this chunk
        try:
            gtf_data = pd.DataFrame(chunk_data, columns=columns)
            gtf_data['start'] = gtf_data['start'].astype(int)
            gtf_data['end'] = gtf_data['end'].astype(int)
        except Exception as e:
            # If there's an error parsing the chunk, skip it
            return
            
        # Process each row
        for _, row in gtf_data.iterrows():
            # Only process exon and CDS features
            if row['feature'] not in ['exon', 'CDS']:
                continue
                
            # Parse attributes using optimized function
            attributes = self._parse_attributes_fast(row['attribute'])
            
            # Extract gene and transcript IDs
            gene_id_row = attributes.get('gene_id', None)
            transcript_id = attributes.get('transcript_id', None)
            gene_name = attributes.get('gene_name', None)
            
            if not gene_id_row or not transcript_id:
                continue
                
            # If a specific gene is requested, only process that gene
            if gene_id and gene_id_row != gene_id:
                continue
                
            # Check if this transcript was already in previous files
            transcript_key = f"{gene_id_row}:{transcript_id}"
            if transcript_key in self._transcripts_from_previous_files:
                continue  # Skip transcripts from previous files
                
            # Store gene name mapping
            if gene_name and gene_id_row:
                self.gene_name_to_id[gene_name] = gene_id_row
                
            # Store gene information (only if not already stored)
            if gene_id_row not in self.gene_info:
                self.gene_info[gene_id_row] = {
                    'seqname': row['seqname'],
                    'strand': row['strand'],
                    'transcripts': []
                }
                
            # Add transcript to gene's transcript list (avoid duplicates)
            if transcript_id not in self.gene_info[gene_id_row]['transcripts']:
                self.gene_info[gene_id_row]['transcripts'].append(transcript_id)
                
            # Initialize transcript structure if not already present
            if gene_id_row not in self.gene_transcripts:
                self.gene_transcripts[gene_id_row] = {}
                
            # Initialize transcript structure if not already present
            if transcript_id not in self.gene_transcripts[gene_id_row]:
                self.gene_transcripts[gene_id_row][transcript_id] = {
                    'exons': [],
                    'cds': [],  # Add CDS information
                    'strand': row['strand'],
                    'seqname': row['seqname']
                }
                
            # Add feature information based on type
            if row['feature'] == 'exon':
                # Add exon information
                new_exon = {
                    'start': row['start'],
                    'end': row['end'],
                    'number': len(self.gene_transcripts[gene_id_row][transcript_id]['exons']) + 1
                }
                
                # Check if exon already exists to avoid duplicates
                exon_exists = False
                for existing_exon in self.gene_transcripts[gene_id_row][transcript_id]['exons']:
                    if existing_exon['start'] == row['start'] and existing_exon['end'] == row['end']:
                        exon_exists = True
                        break
                        
                if not exon_exists:
                    self.gene_transcripts[gene_id_row][transcript_id]['exons'].append(new_exon)
            elif row['feature'] == 'CDS':
                # Add CDS information
                new_cds = {
                    'start': row['start'],
                    'end': row['end']
                }
                
                # Check if CDS already exists to avoid duplicates
                cds_exists = False
                for existing_cds in self.gene_transcripts[gene_id_row][transcript_id]['cds']:
                    if existing_cds['start'] == row['start'] and existing_cds['end'] == row['end']:
                        cds_exists = True
                        break
                        
                if not cds_exists:
                    self.gene_transcripts[gene_id_row][transcript_id]['cds'].append(new_cds)
                
    def _parse_attributes_fast(self, attribute_string):
        """
        Fast parsing of the attribute column of a GTF file.
        
        Args:
            attribute_string (str): The attribute string from a GTF row
            
        Returns:
            dict: Dictionary of attributes
        """
        attributes = {}
        # Use regex for faster parsing
        # Match key-value pairs like gene_id "ENSG00000123456";
        pattern = r'(\w+) "([^"]*)";?'
        matches = re.findall(pattern, attribute_string)
        for key, value in matches:
            attributes[key] = value
        return attributes
    
    def _parse_attributes(self, attribute_string):
        """
        Parse the attribute column of a GTF file.
        Kept for backward compatibility.
        
        Args:
            attribute_string (str): The attribute string from a GTF row
            
        Returns:
            dict: Dictionary of attributes
        """
        attributes = {}
        # Split by semicolon and process each attribute
        for attr in attribute_string.strip().split(';'):
            attr = attr.strip()
            if not attr:
                continue
            # Split by space and remove quotes
            parts = attr.split(' ')
            if len(parts) >= 2:
                key = parts[0].strip()
                value = ' '.join(parts[1:]).strip().strip('"')
                attributes[key] = value
        return attributes
    
    def get_gene_info(self, gene_id):
        """
        Get information about a specific gene.
        
        Args:
            gene_id (str): The gene ID
            
        Returns:
            dict: Gene information
        """
        return self.gene_info.get(gene_id, {})
    
    def get_transcript_structure(self, gene_id, transcript_id):
        """
        Get the structure of a specific transcript.
        
        Args:
            gene_id (str): The gene ID
            transcript_id (str): The transcript ID
            
        Returns:
            dict: Transcript structure information
        """
        if gene_id in self.gene_transcripts and transcript_id in self.gene_transcripts[gene_id]:
            return self.gene_transcripts[gene_id][transcript_id]
        return None
    
    def get_all_genes(self):
        """
        Get all genes in the GTF file(s).
        
        Returns:
            list: List of gene IDs
        """
        return list(self.gene_info.keys())
    
    def get_gene_by_name(self, gene_name):
        """
        Get gene ID by gene name.
        
        Args:
            gene_name (str): The gene name
            
        Returns:
            str: Gene ID or None if not found
        """
        return self.gene_name_to_id.get(gene_name)
    
    def identify_gene_identifier_type(self, identifier):
        """
        Identify the type of gene identifier (gene_id, gene_name, or transcript_id).
        
        Args:
            identifier (str): The identifier to identify
            
        Returns:
            tuple: (identifier_type, gene_id) where identifier_type is one of 'gene_id', 'gene_name', 'transcript_id'
        """
        # Check if it's a gene ID (exists in gene_transcripts)
        if identifier in self.gene_transcripts:
            return 'gene_id', identifier
            
        # Check if it's a gene name (exists in gene_name_to_id)
        if identifier in self.gene_name_to_id:
            return 'gene_name', self.gene_name_to_id[identifier]
            
        # Check if it's a transcript ID
        for gene_id, transcripts in self.gene_transcripts.items():
            if identifier in transcripts:
                return 'transcript_id', gene_id
                
        # Try to parse if not found
        self.parse_gtf(identifier)
        
        # Check again if it's a gene ID
        if identifier in self.gene_transcripts:
            return 'gene_id', identifier
            
        # Check again if it's a transcript ID
        for gene_id, transcripts in self.gene_transcripts.items():
            if identifier in transcripts:
                return 'transcript_id', gene_id
                
        return None, None
    
    def get_transcript_data(self, gene_identifier):
        """
        Extract all transcripts of a specific gene by gene ID, gene name, or transcript ID for visualization.
        
        Args:
            gene_identifier (str): The gene ID, gene name, or transcript ID to extract
            
        Returns:
            dict: Data for the gene transcripts ready for visualization
        """
        # Identify the type of identifier
        identifier_type, gene_id = self.identify_gene_identifier_type(gene_identifier)
        
        if identifier_type is None or gene_id is None:
            raise ValueError(f"Identifier '{gene_identifier}' not found in GTF file")
            
        # Parse the gene data if not already parsed
        if gene_id not in self.gene_transcripts:
            self.parse_gtf(gene_id)
            
        if gene_id not in self.gene_transcripts:
            raise ValueError(f"Gene {gene_id} not found in GTF file")
            
        # Get gene information
        gene_info = self.get_gene_info(gene_id)
        
        # Prepare visualization data
        visualization_data = {
            'gene_id': gene_id,
            'identifier_type': identifier_type,
            'original_identifier': gene_identifier,
            'seqname': gene_info['seqname'],
            'strand': gene_info['strand'],
            'transcripts': []
        }
        
        # If the identifier is a transcript_id, filter to only that transcript
        if identifier_type == 'transcript_id':
            transcript_data = self.gene_transcripts[gene_id][gene_identifier]
            # Sort exons by start position
            exons = sorted(transcript_data['exons'], key=lambda x: x['start'])
            # Sort CDS by start position
            cds = sorted(transcript_data['cds'], key=lambda x: x['start']) if 'cds' in transcript_data else []
            
            transcript_info = {
                'transcript_id': gene_identifier,
                'strand': transcript_data['strand'],
                'exons': exons,
                'cds': cds  # Include CDS information
            }
            
            visualization_data['transcripts'].append(transcript_info)
        else:
            # Add all transcript data for gene_id or gene_name
            for transcript_id, transcript_data in self.gene_transcripts[gene_id].items():
                # Sort exons by start position
                exons = sorted(transcript_data['exons'], key=lambda x: x['start'])
                # Sort CDS by start position
                cds = sorted(transcript_data['cds'], key=lambda x: x['start']) if 'cds' in transcript_data else []
                
                transcript_info = {
                    'transcript_id': transcript_id,
                    'strand': transcript_data['strand'],
                    'exons': exons,
                    'cds': cds  # Include CDS information
                }
                
                visualization_data['transcripts'].append(transcript_info)
            
        return visualization_data


def parse_gtf_for_gene(gtf_file_path, gene_id):
    """
    Parse GTF file(s) and extract information for a specific gene.
    
    Args:
        gtf_file_path (str or list): Path to the GTF file or list of GTF file paths
        gene_id (str): The gene ID to extract information for
        
    Returns:
        dict: Visualization data for the gene transcripts
    """
    parser = GTFParser(gtf_file_path)
    parser.parse_gtf(gene_id)
    return parser.get_transcript_data(gene_id)


def parse_gtf_all_genes(gtf_file_path):
    """
    Parse GTF file(s) and extract information for all genes.
    
    Args:
        gtf_file_path (str or list): Path to the GTF file or list of GTF file paths
        
    Returns:
        dict: Dictionary of all genes and their transcripts
    """
    parser = GTFParser(gtf_file_path)
    parser.parse_gtf()
    return dict(parser.gene_transcripts)


# Test code for the enhanced functionality
if __name__ == "__main__":
    # This section contains test code that can be run in Jupyter
    # Example usage:
    # 
    # # Test with gene ID
    # parser = GTFParser("path/to/your.gtf")
    # data = parser.get_transcript_data("ENSG00000136997")
    # print(f"Identifier type: {data['identifier_type']}")
    # print(f"Original identifier: {data['original_identifier']}")
    # 
    # # Test with gene name
    # data = parser.get_transcript_data("TP53")
    # print(f"Identifier type: {data['identifier_type']}")
    # print(f"Original identifier: {data['original_identifier']}")
    # 
    # # Test with transcript ID
    # data = parser.get_transcript_data("ENST00000269305")
    # print(f"Identifier type: {data['identifier_type']}")
    # print(f"Original identifier: {data['original_identifier']}")
    #
    # # Visualize the data
    # import matplotlib.pyplot as plt
    # from visualizer import visualize_gene_transcripts
    # fig = visualize_gene_transcripts(data)
    # plt.show()
    pass