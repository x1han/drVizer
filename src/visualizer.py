"""
Transcript Visualizer for drVizer
=================================

This module provides functions to visualize gene transcript structures using matplotlib.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import numpy as np
from typing import Union, List


def _compute_track_layout(num_transcripts, separated_parsers_data=None, grouped_bed_data=None, bed_data=None,
                          transcript_width=12, transcript_height=0.5,
                          transcript_row_unit=0.5, distribution_row_unit=0.35,
                          score_track_unit=1.6, coverage_track_unit=1.8,
                          min_track_height=1.2, track_gap=0.12,
                          title_gap=0.22, bottom_margin=0.45):
    track_heights = [max(min_track_height, num_transcripts * transcript_row_unit + 0.8)]

    if separated_parsers_data is not None:
        for parser_info in separated_parsers_data:
            parser_type = parser_info.get('parser_type')
            parser_data = parser_info.get('data', {})
            use_peak_score = parser_info.get('use_peak_score', False)

            if parser_type == 'coverage':
                track_heights.append(max(min_track_height, coverage_track_unit))
            elif use_peak_score:
                track_heights.append(max(min_track_height, score_track_unit))
            else:
                row_count = len(parser_data) if isinstance(parser_data, dict) else 1
                track_heights.append(max(min_track_height, row_count * distribution_row_unit + 0.8))
    elif grouped_bed_data is not None:
        for _, bed_elements in grouped_bed_data.items():
            track_heights.append(max(min_track_height, score_track_unit if bed_elements else min_track_height))
    elif bed_data is not None and len(bed_data) > 0:
        track_heights.extend([max(min_track_height, score_track_unit)] * len(bed_data))

    axes_height = sum(track_heights)
    total_height = axes_height + title_gap + bottom_margin + track_gap * max(0, len(track_heights) - 1)

    return {
        'figure_width': transcript_width,
        'figure_height': total_height,
        'track_heights': track_heights,
        'track_gap': track_gap,
        'title_gap': title_gap,
        'bottom_margin': bottom_margin,
        'axes_height': axes_height
    }


def visualize_gene_transcripts(transcript_data, sort_by_exon_order=True, reverse_order=False,
                              transcript_width=12, transcript_height=0.5, y_spacing=0.3,
                              bed_data=None, track_labels=None,
                              transcript_to_show: Union[str, List[str]] = None):
    """
    Visualize transcript structures for a gene with improved formatting.
    
    Args:
        transcript_data (dict): Transcript data from GTF parser
        sort_by_exon_order (bool): Whether to sort transcripts by exon order
        reverse_order (bool): Whether to reverse the sorting order
        transcript_width (float): Width per transcript (in inches)
        transcript_height (float): Height per transcript (in inches)
        y_spacing (float): Vertical spacing between transcripts (in inches)
        bed_data (list, optional): List of BED elements to visualize
        track_labels (list, optional): List of labels for each track, default is ['GTF Transcripts', 'Track 1', 'Track 2', ...]
        transcript_to_show (Union[str, List[str]], optional): If specified, only plot this specific transcript ID or list of transcript IDs
        
    Returns:
        matplotlib.figure.Figure: The generated figure
    """
    # Extract BED data from transcript_data if it exists (for merged parser)
    if 'bed_data' in transcript_data and bed_data is None:
        bed_data = transcript_data['bed_data']
    
    # Extract grouped BED data from transcript_data if it exists (for merged parser)
    grouped_bed_data = None
    if 'grouped_bed_data' in transcript_data and bed_data is None:
        grouped_bed_data = transcript_data['grouped_bed_data']
    
    # Extract separated parsers data from transcript_data if it exists (for updated merged parser)
    separated_parsers_data = None
    if 'separated_parsers_data' in transcript_data and bed_data is None and grouped_bed_data is None:
        separated_parsers_data = transcript_data['separated_parsers_data']
        
    # Extract track_labels from transcript_data if available (priority over passed parameter)
    if 'track_labels' in transcript_data:
        effective_track_labels = transcript_data['track_labels']
    else:
        effective_track_labels = track_labels
        
    # Extract peak_score from transcript_data if it exists
    if 'peak_score' in transcript_data:
        peak_score = transcript_data['peak_score']
    else:
        # Default: all parsers have peak_score=False
        # Initialize num_other_parsers first
        if separated_parsers_data is not None:
            num_other_parsers = len(separated_parsers_data)
        elif grouped_bed_data is not None:
            num_other_parsers = len(grouped_bed_data)
        elif bed_data is not None:
            num_other_parsers = len(bed_data)
        else:
            num_other_parsers = 0
        peak_score = [False] * (num_other_parsers if num_other_parsers > 0 else 1)
        
    # Extract data
    gene_id = transcript_data['gene_id']
    strand = transcript_data['strand']
    transcripts = transcript_data['transcripts']
    
    # New: Filter transcripts by specific transcript_to_show (support both str and list)
    if transcript_to_show:
        # Convert single string to list for uniform processing
        if isinstance(transcript_to_show, str):
            transcript_ids_to_show = [transcript_to_show]
        else:
            transcript_ids_to_show = transcript_to_show
            
        filtered_transcripts = [t for t in transcripts if t['transcript_id'] in transcript_ids_to_show]
        if not filtered_transcripts:
            print(f"Warning: Transcript(s) {transcript_to_show} not found in gene {gene_id}.")
            # Continue with original transcripts if not found
        else:
            transcripts = filtered_transcripts
            # Update gene_id to reflect multiple transcripts mode
            if len(transcript_ids_to_show) == 1:
                gene_id = transcript_ids_to_show[0]
                print(f"Plotting single transcript: {transcript_ids_to_show[0]}")
            else:
                print(f"Plotting {len(transcript_ids_to_show)} transcripts: {', '.join(transcript_ids_to_show)}")
    
    # Sort transcripts by exon order if requested
    if sort_by_exon_order:
        transcripts = sort_transcripts_by_exon_order(transcripts, reverse_order)
    
    # Calculate number of rows needed for transcripts and BED elements
    num_transcripts = len(transcripts)
    
    # For separated parsers data, count the number of parsers
    if separated_parsers_data is not None:
        num_other_parsers = len(separated_parsers_data)
    # For grouped BED data, count the number of groups instead of individual elements
    elif grouped_bed_data is not None:
        num_bed_groups = len(grouped_bed_data)
        num_other_parsers = num_bed_groups
    else:
        num_bed_groups = 0
        num_other_parsers = len(bed_data) if bed_data else 0
    
    # Total number of tracks (GTF + other parsers)
    total_tracks = 1 + num_other_parsers  # 1 for GTF, plus others

    layout = _compute_track_layout(
        num_transcripts,
        separated_parsers_data=separated_parsers_data,
        grouped_bed_data=grouped_bed_data,
        bed_data=bed_data,
        transcript_width=transcript_width,
        transcript_height=transcript_height
    )
    figure_width = layout['figure_width']
    total_height = layout['figure_height']
    track_heights = layout['track_heights']

    # Create figure with subplots - one for each track
    if total_tracks > 1:
        fig, axes = plt.subplots(
            total_tracks,
            1,
            figsize=(figure_width, total_height),
            sharex=True,
            squeeze=False,
            gridspec_kw={'height_ratios': track_heights}
        )
        axes = axes.flatten()  # Ensure axes is a 1D array
    else:
        fig, ax1 = plt.subplots(figsize=(figure_width, total_height))
        axes = [ax1]
    
    # Find global start and end positions
    all_starts = []
    all_ends = []
    for transcript in transcripts:
        for exon in transcript['exons']:
            all_starts.append(exon['start'])
            all_ends.append(exon['end'])
    
    # Include BED data in global range calculation
    if bed_data:
        for bed_element in bed_data:
            all_starts.append(bed_element['start'])
            all_ends.append(bed_element['end'])
    elif grouped_bed_data:
        # For grouped BED data, include all elements from all groups
        for bed_elements in grouped_bed_data.values():
            for bed_element in bed_elements:
                all_starts.append(bed_element['start'])
                all_ends.append(bed_element['end'])
    
    global_start = min(all_starts) if all_starts else 0
    global_end = max(all_ends) if all_ends else 1000
    
    # Add 2.5% padding to each side of the calculated range (based on filtered transcripts)
    if all_starts and all_ends:
        range_span = global_end - global_start
        padding = range_span * 0.025  # 2.5% padding on each side
        global_start -= padding
        global_end += padding
    
    # Plot GTF transcripts on the first track (axes[0])
    if num_transcripts > 0:
        ax_gtf = axes[0]
        transcript_y_positions = []
        for i, transcript in enumerate(transcripts):
            transcript_id = transcript['transcript_id']
            exons = sorted(transcript['exons'], key=lambda x: x['start'])
            cds_list = sorted(transcript.get('cds', []), key=lambda x: x['start'])  # Get CDS information
            # Calculate y position with proper spacing
            y_pos = num_transcripts - i - 0.5  # Start from top
            transcript_y_positions.append(y_pos)
            
            # Plot transcript line
            ax_gtf.hlines(y_pos, global_start, global_end, color='lightgray', linewidth=0, zorder=1)
            
            # Plot introns as arrows to indicate strand direction
            for j in range(len(exons) - 1):
                intron_start = exons[j]['end']
                intron_end = exons[j+1]['start']
                
                # Draw arrow to indicate strand direction
                if strand == '+':
                    # Forward strand: arrow from left to right
                    ax_gtf.annotate('', xy=(intron_end, y_pos), xytext=(intron_start, y_pos),
                               arrowprops=dict(arrowstyle='->', color='gray', lw=1), zorder=2)
                elif strand == '-':
                    # Reverse strand: arrow from right to left
                    ax_gtf.annotate('', xy=(intron_start, y_pos), xytext=(intron_end, y_pos),
                               arrowprops=dict(arrowstyle='->', color='gray', lw=1), zorder=2)
                else:
                    # Unknown strand: draw simple line
                    ax_gtf.hlines(y_pos, intron_start, intron_end, color='gray', linewidth=1, zorder=2)
            
            # Plot exons
            for exon in exons:
                start = exon['start']
                end = exon['end']
                
                # Create rectangle for exon (hollow with black border)
                exon_rect = patches.Rectangle(
                    (start, y_pos - transcript_height/2),  # (x, y) - centered vertically
                    end - start,            # width
                    transcript_height,      # height
                    linewidth=1,
                    edgecolor='black',
                    facecolor='none',  # Hollow by default
                    zorder=5  # Ensure black border is above CDS
                )
                ax_gtf.add_patch(exon_rect)
                
                # Check if this exon has any CDS regions
                exon_cds_regions = []
                for cds in cds_list:
                    # Find overlap between exon and CDS
                    cds_start = max(start, cds['start'])
                    cds_end = min(end, cds['end'])
                    if cds_start < cds_end:  # There is an overlap
                        exon_cds_regions.append((cds_start, cds_end))
                
                # Plot CDS regions within this exon (filled without border)
                for cds_start, cds_end in exon_cds_regions:
                    cds_rect = patches.Rectangle(
                        (cds_start, y_pos - transcript_height/2),  # (x, y) - centered vertically
                        cds_end - cds_start,    # width
                        transcript_height,      # height
                        linewidth=0,            # No border
                        edgecolor='none',       # No border
                        facecolor='lightblue',  # Filled for CDS regions
                        zorder=4
                    )
                    ax_gtf.add_patch(cds_rect)
        
        # Set y-axis labels for transcripts
        transcript_labels = [t['transcript_id'] for t in transcripts]
        
        # Process labels: if contains underscore, use the last element after splitting
        processed_labels = []
        for label in transcript_labels:
            if 'Novel' in label:
                processed_labels.append(label.split('_')[-1])
            else:
                processed_labels.append(label)
        
        # Truncate long transcript IDs for better display
        truncated_labels = []
        for label in processed_labels:
            if len(label) > 40:
                truncated_labels.append(label[:40] + '...')
            else:
                truncated_labels.append(label)
        
        ax_gtf.set_yticks(transcript_y_positions)
        ax_gtf.set_yticklabels(truncated_labels)
        
        # Formatting for GTF track
        ax_gtf.tick_params(axis='y', which='major', labelsize=8, pad=2)
        ax_gtf.grid(True, axis='x', alpha=0.25)
        ax_gtf.ticklabel_format(axis='x', style='plain', useOffset=False)
        
        # Set track label for GTF track
        if effective_track_labels and len(effective_track_labels) > 0:
            ax_gtf.set_ylabel(effective_track_labels[0], fontsize=10)
            ax_gtf.yaxis.set_label_coords(-0.2, 0.5)
        else:
            ax_gtf.set_ylabel('GTF Transcripts', fontsize=10)
            ax_gtf.yaxis.set_label_coords(-0.2, 0.5)
    
    # Plot BED elements from other parsers on subsequent tracks
    if separated_parsers_data is not None and len(separated_parsers_data) > 0:
        # Plot separated parsers data - each parser gets its own track
        for i, parser_info in enumerate(separated_parsers_data):
            ax_bed = axes[i + 1]  # +1 because first track is GTF
            
            # Get parser-specific information
            parser_bed_data = parser_info['data']
            use_peak_score = parser_info.get('use_peak_score', False)
            y_axis_range = parser_info.get('y_axis_range', None)
            
            # Calculate BED element height (1/3 of transcript height)
            bed_height = transcript_height / 2.5
            
            # Added: Draw Coverage Area Plot (BAM)
            if parser_info.get('parser_type') == 'coverage':
                x = parser_info['data']['x']
                y = parser_info['data']['y']
                color = parser_info.get('color', 'steelblue')
                
                if len(x) > 0:
                    # Fill background
                    alpha = parser_info.get('alpha', 0.6)
                    ax_bed.fill_between(x, y, color=color, alpha=alpha, step='mid', zorder=3)
                    # Draw outline
                    ax_bed.plot(x, y, color=color, lw=0.7, alpha=alpha*1.2, zorder=4)
                    
                    # Y-axis auto scaling or using preset
                    if parser_info.get('y_axis_range'):
                        ax_bed.set_ylim(0, parser_info['y_axis_range'])
                    else:
                        ax_bed.set_ylim(0, np.max(y) * 1.1 if len(y) > 0 else 1)
                
                ax_bed.tick_params(axis='y', which='major', labelsize=8, pad=2)
                ax_bed.set_ylabel(parser_info['track_label'], fontsize=10)
                ax_bed.yaxis.set_label_coords(-0.2, 0.5)
                ax_bed.grid(True, axis='x', alpha=0.25)
                continue # Skip subsequent BED logic
            
            # Handle score type BED data (use_peak_score = True)
            if use_peak_score and isinstance(parser_bed_data, dict):
                # For score type, data is grouped by name, each group contains elements with scores
                bed_names = list(parser_bed_data.keys())
                
                # Get file-level colors if available (for mixed colors in same track)
                file_colors = parser_info.get('file_colors', [parser_info.get('color', 'orange')])
                file_alphas = parser_info.get('file_alphas', [parser_info.get('alpha', 0.8)])
                
                # Track which file each element came from (simplified: alternate files)
                file_cycle = 0
                
                # Plot each group as bars with alternating colors
                for name_idx, (name, bed_elements) in enumerate(parser_bed_data.items()):
                    for bed_element in bed_elements:
                        start = bed_element['start']
                        end = bed_element['end']
                        score = bed_element.get('score', 0.0)  # Get score value
                        
                        # Alternate between file colors
                        color_idx = file_cycle % len(file_colors)
                        color = file_colors[color_idx]
                        alpha = file_alphas[color_idx]
                        file_cycle += 1
                        
                        # Plot a bar for the score with file-specific color
                        ax_bed.bar((start + end) / 2, score, width=end-start, bottom=0, 
                               color=color, edgecolor='none', linewidth=0, zorder=3, alpha=alpha)
                
                # Set y-axis range if specified
                if parser_info.get('y_axis_range'):
                    ax_bed.set_ylim(0, parser_info['y_axis_range'])
                
                ax_bed.tick_params(axis='y', which='major', labelsize=8, pad=2)
                ax_bed.set_ylabel(parser_info['track_label'], fontsize=10)
                ax_bed.yaxis.set_label_coords(-0.2, 0.5)
                ax_bed.grid(True, axis='x', alpha=0.25)
                continue # Skip subsequent BED logic
            
            # For distribution type, plot elements with same name on same line
            if not use_peak_score:
                # Calculate y positions for different names
                bed_names = list(parser_bed_data.keys())
                num_bed_names = len(bed_names)
                bed_y_positions = []
                bed_labels = []
                
                # Add margin at the top to prevent the first element from touching the subplot boundary
                top_margin = 0.3
                
                # Calculate y positions for each name
                for name_idx, name in enumerate(bed_names):
                    y_pos = num_bed_names - name_idx - 0.5 + top_margin  # Start from top with minimal spacing and margin
                    bed_y_positions.append(y_pos)
                    # Use track_label for distribution type, actual name for score type
                    bed_labels.append(name)
                
                # Plot each element in the parser's data
                for name_idx, (name, bed_elements) in enumerate(parser_bed_data.items()):
                    y_pos = bed_y_positions[name_idx]
                    for bed_element in bed_elements:
                        start = bed_element['start']
                        end = bed_element['end']
                        
                        # Create rectangle for BED element with border and 1/3 height
                        # For BED entries with zero-width (start==end), adjust slightly to be visible
                        width = end - start
                        if width == 0:
                            width = 1  # Make sure even point-like entries are visible
                            start = start - 0.5
                        
                        bed_rect = patches.Rectangle(
                            (start, y_pos - bed_height/2),  # (x, y) - centered vertically
                            width,                  # width
                            bed_height,             # height (1/3 of transcript height)
                            linewidth=1,            # Add border
                            edgecolor=None,         # None border
                            facecolor=parser_info.get('color', 'orange'),     # Use parser color
                            alpha=parser_info.get('alpha', 0.8),  # Use parser alpha
                            zorder=3                # Below CDS but above transcript lines
                        )
                        ax_bed.add_patch(bed_rect)
                
                # Set y-axis labels for BED elements
                # Process labels: truncate if too long
                processed_bed_labels = []
                for label in bed_labels:
                    if len(label) > 40:
                        processed_bed_labels.append(label[:40] + '...')
                    else:
                        processed_bed_labels.append(label)
                
                ax_bed.set_yticks(bed_y_positions)
                ax_bed.set_yticklabels(processed_bed_labels)
                
                # Set y-axis limits to ensure consistent spacing at top and bottom
                if num_bed_names > 0:
                    # The first element position (topmost) with margin
                    first_element_y = num_bed_names - 0.5 + top_margin
                    # The last element position (bottommost) with margin
                    last_element_y = 1 - 0.5 + top_margin  # num_bed_names - (num_bed_names-1) - 0.5 + top_margin
                    # Set limits with 0.3 margin at both top and bottom
                    ax_bed.set_ylim(last_element_y - 0.3, first_element_y + 0.3)
                
                # Set y-axis label for this BED track (distribution type)
                if effective_track_labels and len(effective_track_labels) > i + 1:
                    ax_bed.set_ylabel(effective_track_labels[i + 1], fontsize=10)
                else:
                    ax_bed.set_ylabel(f'Track {i + 1}', fontsize=10)
                ax_bed.yaxis.set_label_coords(-0.2, 0.5)
                
                ax_bed.tick_params(axis='y', which='major', labelsize=8, pad=2)
                ax_bed.grid(True, axis='x', alpha=0.25)
                ax_bed.ticklabel_format(axis='x', style='plain', useOffset=False)
            else:
                # For score type, plot normally (existing logic)
                # Handle case where all elements might have the same name
                all_names = set()
                for name in parser_bed_data.keys():
                    all_names.add(name)
                
                # If all names are the same (like '.'), treat each element individually
                if len(all_names) == 1 and list(all_names)[0] == '.':
                    # Create individual entries for each element
                    individual_data = {}
                    element_counter = 0
                    for name, bed_elements in parser_bed_data.items():
                        for bed_element in bed_elements:
                            element_key = f"element_{element_counter}"
                            individual_data[element_key] = [bed_element]
                            element_counter += 1
                    parser_bed_data = individual_data
                
                for name, bed_elements in parser_bed_data.items():
                    for bed_element in bed_elements:
                        start = bed_element['start']
                        end = bed_element['end']
                        score = bed_element.get('score', 1.0)
                        
                        # For BED entries, ensure start <= end
                        if start > end:
                            start, end = end, start
                        
                        # Use barplot for peak score visualization
                        score = bed_element.get('score', 1.0)  # Use score value for bar height
                        
                        # If y_axis_range is defined and score exceeds it, clip the score
                        if y_axis_range is not None:
                            score = min(score, y_axis_range)
                        
                        # Calculate width (should be positive now)
                        width = end - start
                        # For BED entries with zero-width, use a minimum visible width
                        if width <= 0:
                            width = 1  # Keep actual width logic as requested
                        
                        # Plot hollow bars with red borders only (no fill)
                        ax_bed.bar((start + end) / 2, score, width=width, bottom=0, 
                               color='none', edgecolor='red', linewidth=0.25, 
                               zorder=10)
                        
                        # Set y-axis limit based on y_axis_range if provided
                        if y_axis_range is not None:
                            # First relim to include all data points
                            ax_bed.relim()
                            ax_bed.autoscale_view(scalex=False, scaley=True)
                            # Then constrain to y_axis_range with padding
                            current_ylim = ax_bed.get_ylim()
                            ax_bed.set_ylim(0, min(current_ylim[1], y_axis_range * 1.1))
                        else:
                            # Auto-scale based on maximum score in the track
                            ax_bed.relim()
                            ax_bed.autoscale_view(scalex=False, scaley=True)
                        
                        # Set y-axis limit based on y_axis_range if provided
                        if y_axis_range is not None:
                            ax_bed.set_ylim(0, y_axis_range * 1.1)  # Add 10% padding
                        else:
                            # Auto-scale based on maximum score in the track
                            ax_bed.relim()
                            ax_bed.autoscale_view(scalex=False, scaley=True)
                
                # Set y-axis label for this BED track
                if effective_track_labels and len(effective_track_labels) > i + 1:
                    ax_bed.set_ylabel(effective_track_labels[i + 1], fontsize=10)
                else:
                    ax_bed.set_ylabel(f'Track {i + 1}', fontsize=10)
                ax_bed.yaxis.set_label_coords(-0.2, 0.5)
                
                # Formatting for BED track
                ax_bed.tick_params(axis='y', which='major', labelsize=8, pad=2)
                ax_bed.grid(True, axis='x', alpha=0.25)
                ax_bed.ticklabel_format(axis='x', style='plain', useOffset=False)
    elif grouped_bed_data is not None and num_bed_groups > 0:
        # Plot grouped BED data - each group gets its own track
        for i, (name, bed_elements) in enumerate(grouped_bed_data.items()):
            ax_bed = axes[i + 1]  # +1 because first track is GTF
            
            # Determine if this group should be plotted as peak scores
            use_peak_score = False
            if i < len(peak_score):
                use_peak_score = peak_score[i]
            
            # Calculate BED element height (1/3 of transcript height)
            bed_height = transcript_height / 2.5
            
            # Plot each element in the group
            for bed_element in bed_elements:
                start = bed_element['start']
                end = bed_element['end']
                
                if use_peak_score:
                    # Use barplot for peak score visualization
                    score = bed_element.get('score', 1.0)  # Use score value for bar height
                    # Plot a bar without border
                    color = parser_info.get('color', 'orange')
                    alpha = parser_info.get('alpha', 0.8)
                    ax_bed.bar((start + end) / 2, score, width=end-start, bottom=0, 
                           color=color, edgecolor='none', linewidth=0, zorder=3, alpha=alpha)
                else:
                    # Create rectangle for BED element with border and 1/3 height
                    color = parser_info.get('color', 'orange')
                    alpha = parser_info.get('alpha', 0.8)
                    bed_rect = patches.Rectangle(
                        (start, 0),  # (x, y) - centered vertically at y=0
                        end - start,            # width
                        bed_height,             # height (1/3 of transcript height)
                        linewidth=1,            # Add border
                        edgecolor=None,         # None border
                        facecolor=color,     # Use parser color
                        alpha=alpha,  # Use parser alpha
                        zorder=3                # Below CDS but above transcript lines
                    )
                    ax_bed.add_patch(bed_rect)
            
            # Set y-axis label for this BED track
            if effective_track_labels and len(effective_track_labels) > i + 1:
                ax_bed.set_ylabel(effective_track_labels[i + 1], fontsize=10)
            else:
                ax_bed.set_ylabel(f'Track {i + 1}', fontsize=10)
            ax_bed.yaxis.set_label_coords(-0.2, 0.5)
            
            # Formatting for BED track
            ax_bed.tick_params(axis='y', which='major', labelsize=8, pad=2)
            ax_bed.grid(True, axis='x', alpha=0.25)
            ax_bed.ticklabel_format(axis='x', style='plain', useOffset=False)
            ax_bed.set_ylim(-bed_height/2, bed_height/2)  # Center the track
    
    elif bed_data and len(bed_data) > 0:
        # Plot ungrouped BED data (backward compatibility) - each element gets its own track
        for i, bed_element in enumerate(bed_data):
            ax_bed = axes[i + 1]  # +1 because first track is GTF
            
            # Determine if this element should be plotted as peak scores
            use_peak_score = False
            if i < len(peak_score):
                use_peak_score = peak_score[i]
            
            # Calculate BED element height (1/3 of transcript height)
            bed_height = transcript_height / 2.5
            
            start = bed_element['start']
            end = bed_element['end']
            
            if use_peak_score:
                # Use barplot for peak score visualization
                score = bed_element.get('score', 1.0)  # Use score value for bar height
                # Plot a bar without border
                ax_bed.bar((start + end) / 2, score, width=end-start, bottom=0, 
                       color='orange', edgecolor='none', linewidth=0, zorder=3, alpha=0.8)
            else:
                # Create rectangle for BED element with border and 1/3 height
                bed_rect = patches.Rectangle(
                    (start, 0),  # (x, y) - centered vertically at y=0
                    end - start,            # width
                    bed_height,             # height (1/3 of transcript height)
                    linewidth=1,            # Add border
                    edgecolor=None,         # None border
                    facecolor='orange',     # Filled with single color
                    alpha=0.8,  # Default alpha
                    zorder=3                # Below CDS but above transcript lines
                )
                ax_bed.add_patch(bed_rect)
            
            # Set y-axis label for this BED track
            if effective_track_labels and len(effective_track_labels) > i + 1:
                ax_bed.set_ylabel(effective_track_labels[i + 1], fontsize=10)
            else:
                ax_bed.set_ylabel(f'Track {i + 1}', fontsize=10)
            ax_bed.yaxis.set_label_coords(-0.2, 0.5)
            
            # Formatting for BED track
            ax_bed.tick_params(axis='y', which='major', labelsize=8, pad=2)
            ax_bed.grid(True, axis='x', alpha=0.25)
            ax_bed.ticklabel_format(axis='x', style='plain', useOffset=False)
            ax_bed.set_ylim(-bed_height/2, bed_height/2)  # Center the track
    
    # Set x-axis label only on the bottom subplot
    axes[-1].set_xlabel(f'Genomic Position ({transcript_data["seqname"]})' if 'seqname' in transcript_data else 'Genomic Position', fontsize=10)
    
    # Improved title with strand direction indicators and identifier type
    # Respect the original identifier type detected by the parser
    identifier_type = transcript_data.get('identifier_type', 'gene_id')
    original_identifier = transcript_data.get('original_identifier', gene_id)
    
    # Check if we have multiple genes (comma-separated identifiers)
    gene_identifiers = original_identifier.split(',') if ',' in original_identifier else [original_identifier]
    
    # Override to gene_id only when showing multiple transcripts and input wasn't transcript_id
    if not transcript_to_show and identifier_type == 'transcript_id':
        # This case: user gave transcript_id but we're showing all transcripts of the gene
        # So display as Gene ID for clarity
        identifier_label = f'Gene ID: {original_identifier}'
        identifier_type = 'gene_id'
    elif transcript_to_show:
        # Handle both single transcript and multiple transcripts display
        if isinstance(transcript_to_show, str):
            # Single transcript mode - show the specified transcript
            identifier_label = f'Transcript ID: {transcript_to_show}'
            identifier_type = 'transcript_id'
        else:
            # Multiple transcripts mode - show count and first few transcript IDs
            if len(transcript_to_show) <= 3:
                transcript_list = ', '.join(transcript_to_show)
            else:
                transcript_list = ', '.join(transcript_to_show[:3]) + f' and {len(transcript_to_show) - 3} more'
            identifier_label = f'{len(transcript_to_show)} Transcripts: {transcript_list}'
            identifier_type = 'transcript_id'
    elif len(gene_identifiers) > 1:
        # Multiple genes case - preserve the original identifier type but show count
        gene_count = len(gene_identifiers)
        if gene_count == 2:
            identifier_label = f'{gene_count} Genes: {original_identifier}'
        else:
            # For more than 2 genes, show first 2 and count the rest
            if len(gene_identifiers) <= 3:
                identifier_label = f'{gene_count} Genes: {original_identifier}'
            else:
                identifier_label = f'{gene_count} Genes: {", ".join(gene_identifiers[:2]).strip()} and {len(gene_identifiers) - 2} more'
    elif identifier_type == 'gene_name':
        # Gene name input - preserve the original gene name in display
        identifier_label = f'Gene Name: {original_identifier}'
    elif identifier_type == 'gene_id':
        # Gene ID input - show as Gene ID
        identifier_label = f'Gene ID: {original_identifier}'
    else:
        # Fallback - show as Gene ID
        identifier_label = f'Gene ID: {original_identifier}'
        identifier_type = 'gene_id'
    
    if strand == '+':
        title = f'Transcript Structure for {identifier_label} (5\' → 3\')'
    elif strand == '-':
        title = f'Transcript Structure for {identifier_label} (3\' ← 5\')'
    else:
        title = f'Transcript Structure for {identifier_label} ({strand} strand)'
    
    fig.suptitle(title, fontsize=14)

    figure_height = layout['figure_height']
    title_gap = layout['title_gap']
    bottom_margin = layout['bottom_margin']
    track_gap = layout['track_gap']
    axes_height = layout['axes_height']

    if len(axes) > 1:
        hspace = track_gap * len(axes) / axes_height if axes_height > 0 else 0
    else:
        hspace = 0

    title_height = 0.28
    top_fraction = 1 - (title_gap + title_height) / figure_height
    bottom_fraction = bottom_margin / figure_height
    title_y = 1 - title_height / (2 * figure_height)

    fig.subplots_adjust(top=top_fraction, bottom=bottom_fraction, hspace=hspace)
    fig._suptitle.set_y(title_y)
    
    # 1. Uniform processing of X-axis visibility control for all tracks
    for i, ax in enumerate(axes):
        # If not the last track
        if i < len(axes) - 1:
            # Hide X-axis tick marks (bottom=False) and tick labels (labelbottom=False)
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
            # Remove X-axis spine (optional, if you want the plot to look more coherent)
            # ax.spines['bottom'].set_visible(False) 
        else:
            # Bottommost track ensures ticks and labels are displayed
            ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True)
            ax.ticklabel_format(axis='x', style='plain', useOffset=False)

    # 2. Remove x-axis titles from all tracks except the bottommost one (prevent overlap)
    for ax in axes[:-1]:
        ax.set_xlabel('')

    # 3. Adjust layout and spacing
    # Use fixed physical gaps converted to subplot fractions so spacing stays visually stable.

    # Debug: Check actual xlim after all adjustments
    
    return fig


def sort_transcripts_by_exon_order(transcripts, reverse=False):
    """
    Sort transcripts by the order of their first exon start positions.
    
    Args:
        transcripts (list): List of transcript dictionaries
        reverse (bool): Whether to reverse the sorting order
        
    Returns:
        list: Sorted list of transcript dictionaries
    """
    def get_first_exon_start(transcript):
        if transcript['exons']:
            return min(exon['start'] for exon in transcript['exons'])
        return 0
    
    # Sort and return the transcripts
    sorted_transcripts = sorted(transcripts, key=get_first_exon_start, reverse=reverse)
    
    # Print sorting information for debugging
    # print("Transcript sorting order:")
    # for i, transcript in enumerate(sorted_transcripts):
    #     first_exon_start = get_first_exon_start(transcript)
    #     print(f"  {i+1}. {transcript['transcript_id']}: first exon start at {first_exon_start}")
    
    return sorted_transcripts


def save_visualization(fig, output_path, format='png', dpi=300):
    """
    Save the visualization to a file.
    
    Args:
        fig (matplotlib.figure.Figure): The figure to save
        output_path (str): Path to save the figure
        format (str): Format to save the figure ('png', 'pdf', 'svg', etc.)
        dpi (int): DPI for the saved figure
    """
    fig.savefig(output_path, format=format, dpi=dpi, bbox_inches='tight')
    plt.close(fig)





def merge_parsers(*parsers):
    """
    Merge multiple parsers (GTF, BED, etc.) to create a combined parser that can visualize 
    gene transcripts and multiple annotation elements.
    
    Args:
        *parsers: Variable number of parser objects (first must be GTFParser, followed by BEDParser(s) or other parsers)
        
    Returns:
        MultiMergedParser: A merged parser object that combines multiple data sources
    """
    if len(parsers) < 2:
        raise ValueError("At least 2 parsers are required for merging")
    
    # First parser must be GTF parser
    gtf_parser = parsers[0]
    
    # Check if the first parser is a GTF parser
    from .gtf_parser import GTFParser
    if not isinstance(gtf_parser, GTFParser):
        raise ValueError("First parser must be a GTFParser instance")
    
    other_parsers = parsers[1:]
    
    # Ensure GTF parser is parsed
    if not gtf_parser._parsed:
        gtf_parser.parse_gtf()
    
    # Ensure other parsers are parsed
    for parser in other_parsers:
        if hasattr(parser, 'parse_bed') and hasattr(parser, '_parsed'):
            if not parser._parsed:
                parser.parse_bed()
        elif hasattr(parser, 'parse_gtf') and hasattr(parser, '_parsed'):
            if not parser._parsed:
                parser.parse_gtf()
    
    return MultiMergedParser(gtf_parser, other_parsers)

class MultiMergedParser:
    """
    A class that combines GTF parser with multiple other parsers (BED, BAM, etc.) for visualization.
    """
    
    def __init__(self, gtf_parser, other_parsers):
        """
        Initialize the multi-merged parser with GTF parser and multiple other parsers.
        
        Args:
            gtf_parser (GTFParser): Initialized and parsed GTF parser
            other_parsers (list): List of other initialized and parsed parsers (BEDParser, etc.)
        """
        self.gtf_parser = gtf_parser
        self.other_parsers = other_parsers
        
    def get_transcript_data(self, gene_identifier, y_axis_ranges=None, transcript_to_show=None):
        """
        Support single gene identifier (str) or multiple gene identifiers (list).
        
        Args:
            gene_identifier (Union[str, List[str]]): The gene ID, gene name, or transcript ID to extract,
                                                    or a list of such identifiers
            y_axis_ranges (list): List of maximum y-axis values for score parsers
            transcript_to_show (Union[str, List[str]], optional): If specified, only include data for these specific transcript(s)
            
        Returns:
            dict: Data for the gene transcripts and overlapping annotations from multiple sources ready for visualization
            
        Raises:
            ValueError: If genes are not on the same chromosome
        """
        # Convert input to list for unified processing
        identifiers = [gene_identifier] if isinstance(gene_identifier, str) else gene_identifier
        
        combined_gene_data = None
        target_chrom = None

        for ident in identifiers:
            current_data = self.gtf_parser.get_transcript_data(ident)
            
            # Chromosome consistency check
            if target_chrom is None:
                target_chrom = current_data['seqname']
            elif target_chrom != current_data['seqname']:
                raise ValueError(f"Error: Genes must be on the same chromosome. Found {target_chrom} and {current_data['seqname']}")

            if combined_gene_data is None:
                combined_gene_data = current_data
            else:
                # Merge transcripts
                combined_gene_data['transcripts'].extend(current_data['transcripts'])
                # Update identifier display
                combined_gene_data['original_identifier'] += f", {current_data['original_identifier']}"

        # Subsequent logic for extracting overlapping data remains the same, but using target_chrom
        seqname = target_chrom
        
        # Calculate range - if transcript_to_show is specified, use only those transcripts
        if transcript_to_show is not None:
            # Convert to list if single string
            transcripts_to_use = [transcript_to_show] if isinstance(transcript_to_show, str) else transcript_to_show
            
            # Filter transcripts to only include those specified
            filtered_transcripts = [
                t for t in combined_gene_data['transcripts'] 
                if t['transcript_id'] in transcripts_to_use
            ]
            
            # Calculate coordinates from filtered transcripts
            all_starts = [exon['start'] for t in filtered_transcripts for exon in t['exons']]
            all_ends = [exon['end'] for t in filtered_transcripts for exon in t['exons']]
        else:
            # Use all transcripts (original behavior)
            all_starts = [exon['start'] for t in combined_gene_data['transcripts'] for exon in t['exons']]
            all_ends = [exon['end'] for t in combined_gene_data['transcripts'] for exon in t['exons']]
        
        # Add 2.5% padding to the range
        if all_starts and all_ends:
            range_min, range_max = min(all_starts), max(all_ends)
            padding = int((range_max - range_min) * 0.025)  # Convert to integer
            gene_start = int(range_min - padding)  # Ensure integer
            gene_end = int(range_max + padding)    # Ensure integer
        else:
            # Fallback if no exons found
            gene_start, gene_end = 0, 1000
        
        # Extract other track data (logic same as original MultiMergedParser, but ensure alpha parameters are passed)
        parsers_data = []
        peak_score_array = []  # Track peak_score for each parser
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
                'file_colors': getattr(parser, 'file_colors', None),  # Add file-level colors
                'file_alphas': getattr(parser, 'file_alphas', None)   # Add file-level alphas
            })
            peak_score_array.append(is_score_type)
        
        combined_gene_data['separated_parsers_data'] = parsers_data
        combined_gene_data['peak_score'] = peak_score_array  # Add peak_score array
        return combined_gene_data