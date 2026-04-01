"""
Utility functions for drVizer Python module
===========================================

This module provides utility functions for data processing and format conversion.
"""

import json
import pandas as pd
from collections import defaultdict


def convert_to_json(transcript_data, output_path=None):
    """
    Convert transcript data to JSON format.
    
    Args:
        transcript_data (dict): Transcript data from GTF parser
        output_path (str, optional): Path to save JSON file. If None, returns JSON string.
        
    Returns:
        str: JSON string if output_path is None, otherwise None
    """
    json_data = json.dumps(transcript_data, indent=2)
    
    if output_path:
        with open(output_path, 'w') as f:
            f.write(json_data)
        return None
    else:
        return json_data


def convert_to_csv(transcript_data, output_path):
    """
    Convert transcript data to CSV format.
    
    Args:
        transcript_data (dict): Transcript data from GTF parser
        output_path (str): Path to save CSV file
    """
    # Create a list to hold all exon data
    exon_data = []
    
    gene_id = transcript_data['gene_id']
    strand = transcript_data['strand']
    
    # Extract exon data for each transcript
    for transcript in transcript_data['transcripts']:
        transcript_id = transcript['transcript_id']
        for exon in transcript['exons']:
            exon_data.append({
                'gene_id': gene_id,
                'strand': strand,
                'transcript_id': transcript_id,
                'exon_number': exon.get('number', ''),
                'exon_start': exon['start'],
                'exon_end': exon['end']
            })
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(exon_data)
    df.to_csv(output_path, index=False)


def get_transcript_stats(transcript_data):
    """
    Get statistics for transcript data.
    
    Args:
        transcript_data (dict): Transcript data from GTF parser
        
    Returns:
        dict: Statistics about the transcripts
    """
    stats = {}
    
    gene_id = transcript_data['gene_id']
    transcripts = transcript_data['transcripts']
    
    stats['gene_id'] = gene_id
    stats['num_transcripts'] = len(transcripts)
    
    # Calculate exon statistics
    exon_counts = []
    transcript_lengths = []
    
    for transcript in transcripts:
        exons = transcript['exons']
        exon_counts.append(len(exons))
        
        # Calculate transcript length
        if exons:
            start = min(exon['start'] for exon in exons)
            end = max(exon['end'] for exon in exons)
            transcript_lengths.append(end - start + 1)
    
    stats['exon_counts'] = exon_counts
    stats['transcript_lengths'] = transcript_lengths
    stats['total_exons'] = sum(exon_counts)
    stats['avg_exons_per_transcript'] = sum(exon_counts) / len(exon_counts) if exon_counts else 0
    stats['min_exons'] = min(exon_counts) if exon_counts else 0
    stats['max_exons'] = max(exon_counts) if exon_counts else 0
    stats['avg_transcript_length'] = sum(transcript_lengths) / len(transcript_lengths) if transcript_lengths else 0
    stats['min_transcript_length'] = min(transcript_lengths) if transcript_lengths else 0
    stats['max_transcript_length'] = max(transcript_lengths) if transcript_lengths else 0
    
    return stats


def filter_transcripts(transcript_data, min_exons=None, max_exons=None, 
                      min_length=None, max_length=None):
    """
    Filter transcripts based on various criteria.
    
    Args:
        transcript_data (dict): Transcript data from GTF parser
        min_exons (int, optional): Minimum number of exons
        max_exons (int, optional): Maximum number of exons
        min_length (int, optional): Minimum transcript length
        max_length (int, optional): Maximum transcript length
        
    Returns:
        dict: Filtered transcript data
    """
    filtered_data = transcript_data.copy()
    filtered_transcripts = []
    
    for transcript in transcript_data['transcripts']:
        # Count exons
        exon_count = len(transcript['exons'])
        
        # Calculate transcript length
        transcript_length = 0
        if transcript['exons']:
            start = min(exon['start'] for exon in transcript['exons'])
            end = max(exon['end'] for exon in transcript['exons'])
            transcript_length = end - start + 1
        
        # Apply filters
        if min_exons is not None and exon_count < min_exons:
            continue
        if max_exons is not None and exon_count > max_exons:
            continue
        if min_length is not None and transcript_length < min_length:
            continue
        if max_length is not None and transcript_length > max_length:
            continue
            
        filtered_transcripts.append(transcript)
    
    filtered_data['transcripts'] = filtered_transcripts
    return filtered_data


def merge_transcript_data(data_list):
    """
    Merge multiple transcript data dictionaries.
    
    Args:
        data_list (list): List of transcript data dictionaries
        
    Returns:
        dict: Merged transcript data
    """
    if not data_list:
        return {}
    
    # Use the first entry as the base
    merged_data = data_list[0].copy()
    merged_transcripts = merged_data['transcripts'][:]
    
    # Add transcripts from other entries
    for data in data_list[1:]:
        for transcript in data['transcripts']:
            # Check if transcript already exists
            transcript_id = transcript['transcript_id']
            existing_ids = [t['transcript_id'] for t in merged_transcripts]
            
            if transcript_id not in existing_ids:
                merged_transcripts.append(transcript)
    
    merged_data['transcripts'] = merged_transcripts
    return merged_data