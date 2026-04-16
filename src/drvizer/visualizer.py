"""
Transcript Visualizer for drVizer
=================================

This module provides functions to visualize gene transcript structures using matplotlib.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from typing import Union, List


def _has_right_labels(prepared_tracks):
    """Check if right-side transcript labels should be rendered."""
    if not prepared_tracks:
        return False
    return any(track.get("transcript_id") for track in prepared_tracks)


def _shorten_transcript_label(label):
    if 'Novel' in label:
        label = label.split('_')[-1]
    if len(label) > 40:
        return label[:40] + '...'
    return label


def _add_right_side_transcript_labels(fig, axes, right_label_groups):
    """Render vertical transcript labels on the right side of the figure."""
    if not right_label_groups:
        return

    for group in right_label_groups:
        start_idx = group["start_index"]
        end_idx = group["end_index"]
        transcript_id = group["transcript_id"]

        start_axis = axes[start_idx + 1]
        end_axis = axes[end_idx + 1]

        start_box = start_axis.get_position()
        end_box = end_axis.get_position()

        y_center = (start_box.y1 + end_box.y0) / 2
        display_label = _shorten_transcript_label(transcript_id)

        fig.text(
            0.985,
            y_center,
            display_label,
            rotation=270,
            va="center",
            ha="center",
            fontsize=8,
        )


def _compute_track_layout(num_transcripts, prepared_tracks=None,
                          transcript_width=12, transcript_height=0.5,
                          transcript_row_unit=0.5, distribution_row_unit=0.35,
                          score_track_unit=1.6, coverage_track_unit=1.8,
                          min_track_height=1.2, track_gap=0.12,
                          title_gap=0.22, bottom_margin=0.45):
    track_heights = [max(min_track_height, num_transcripts * transcript_row_unit + 0.8)]

    for track in prepared_tracks or []:
        track_kind = track.get('kind', 'distribution')
        track_data = track.get('data', {})

        if track_kind == 'coverage':
            track_heights.append(max(min_track_height, coverage_track_unit))
        elif track_kind == 'score':
            track_heights.append(max(min_track_height, score_track_unit))
        else:
            row_count = len(track_data) if isinstance(track_data, dict) else 1
            track_heights.append(max(min_track_height, row_count * distribution_row_unit + 0.8))

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
        bed_data (list, optional): Deprecated legacy parameter; normalized track data should be passed in transcript_data
        track_labels (list, optional): List of labels for each track, default is ['GTF Transcripts', 'Track 1', 'Track 2', ...]
        transcript_to_show (Union[str, List[str]], optional): If specified, only plot this specific transcript ID or list of transcript IDs

    Returns:
        matplotlib.figure.Figure: The generated figure
    """
    prepared_tracks = transcript_data.get('prepared_tracks', [])

    if 'track_labels' in transcript_data:
        effective_track_labels = transcript_data['track_labels']
    else:
        effective_track_labels = track_labels

    gene_id = transcript_data['gene_id']
    strand = transcript_data['strand']
    transcripts = transcript_data['transcripts']

    if transcript_to_show:
        transcript_ids_to_show = [transcript_to_show] if isinstance(transcript_to_show, str) else transcript_to_show
        filtered_transcripts = [t for t in transcripts if t['transcript_id'] in transcript_ids_to_show]
        if not filtered_transcripts:
            print(f"Warning: Transcript(s) {transcript_to_show} not found in gene {gene_id}.")
        else:
            transcripts = filtered_transcripts
            if len(transcript_ids_to_show) == 1:
                gene_id = transcript_ids_to_show[0]
                print(f"Plotting single transcript: {transcript_ids_to_show[0]}")
            else:
                print(f"Plotting {len(transcript_ids_to_show)} transcripts: {', '.join(transcript_ids_to_show)}")

    if sort_by_exon_order:
        transcripts = sort_transcripts_by_exon_order(transcripts, reverse_order)

    num_transcripts = len(transcripts)
    total_tracks = 1 + len(prepared_tracks)

    layout = _compute_track_layout(
        num_transcripts,
        prepared_tracks=prepared_tracks,
        transcript_width=transcript_width,
        transcript_height=transcript_height
    )
    figure_width = layout['figure_width']
    total_height = layout['figure_height']
    track_heights = layout['track_heights']

    if total_tracks > 1:
        fig, axes = plt.subplots(
            total_tracks,
            1,
            figsize=(figure_width, total_height),
            sharex=True,
            squeeze=False,
            gridspec_kw={'height_ratios': track_heights}
        )
        axes = axes.flatten()
    else:
        fig, ax1 = plt.subplots(figsize=(figure_width, total_height))
        axes = [ax1]

    all_starts = []
    all_ends = []
    for transcript in transcripts:
        for exon in transcript['exons']:
            all_starts.append(exon['start'])
            all_ends.append(exon['end'])

    for track in prepared_tracks:
        track_data = track.get('data', {})
        track_kind = track.get('kind', 'distribution')

        if track_kind == 'coverage':
            positions = list(track_data.get('x', []))
            all_starts.extend(positions)
            all_ends.extend(positions)
        elif isinstance(track_data, dict):
            for elements in track_data.values():
                for element in elements:
                    all_starts.append(element['start'])
                    all_ends.append(element['end'])

    global_start = min(all_starts) if all_starts else 0
    global_end = max(all_ends) if all_ends else 1000

    if all_starts and all_ends:
        range_span = global_end - global_start
        padding = range_span * 0.025
        global_start -= padding
        global_end += padding

    if num_transcripts > 0:
        ax_gtf = axes[0]
        transcript_y_positions = []
        for i, transcript in enumerate(transcripts):
            exons = sorted(transcript['exons'], key=lambda x: x['start'])
            cds_list = sorted(transcript.get('cds', []), key=lambda x: x['start'])
            y_pos = num_transcripts - i - 0.5
            transcript_y_positions.append(y_pos)

            ax_gtf.hlines(y_pos, global_start, global_end, color='lightgray', linewidth=0, zorder=1)

            for j in range(len(exons) - 1):
                intron_start = exons[j]['end']
                intron_end = exons[j + 1]['start']

                if strand == '+':
                    ax_gtf.annotate('', xy=(intron_end, y_pos), xytext=(intron_start, y_pos),
                                    arrowprops=dict(arrowstyle='->', color='gray', lw=1), zorder=2)
                elif strand == '-':
                    ax_gtf.annotate('', xy=(intron_start, y_pos), xytext=(intron_end, y_pos),
                                    arrowprops=dict(arrowstyle='->', color='gray', lw=1), zorder=2)
                else:
                    ax_gtf.hlines(y_pos, intron_start, intron_end, color='gray', linewidth=1, zorder=2)

            for exon in exons:
                start = exon['start']
                end = exon['end']

                exon_rect = patches.Rectangle(
                    (start, y_pos - transcript_height / 2),
                    end - start,
                    transcript_height,
                    linewidth=1,
                    edgecolor='black',
                    facecolor='none',
                    zorder=5
                )
                ax_gtf.add_patch(exon_rect)

                exon_cds_regions = []
                for cds in cds_list:
                    cds_start = max(start, cds['start'])
                    cds_end = min(end, cds['end'])
                    if cds_start < cds_end:
                        exon_cds_regions.append((cds_start, cds_end))

                for cds_start, cds_end in exon_cds_regions:
                    cds_rect = patches.Rectangle(
                        (cds_start, y_pos - transcript_height / 2),
                        cds_end - cds_start,
                        transcript_height,
                        linewidth=0,
                        edgecolor='none',
                        facecolor='lightblue',
                        zorder=4
                    )
                    ax_gtf.add_patch(cds_rect)

        transcript_labels = [t['transcript_id'] for t in transcripts]
        truncated_labels = [_shorten_transcript_label(label) for label in transcript_labels]

        ax_gtf.set_yticks(transcript_y_positions)
        ax_gtf.set_yticklabels(truncated_labels)
        ax_gtf.tick_params(axis='y', which='major', labelsize=8, pad=2)
        ax_gtf.grid(True, axis='x', alpha=0.25)
        ax_gtf.ticklabel_format(axis='x', style='plain', useOffset=False)

        if effective_track_labels and len(effective_track_labels) > 0:
            ax_gtf.set_ylabel(effective_track_labels[0], fontsize=10)
        else:
            ax_gtf.set_ylabel('GTF Transcripts', fontsize=10)
        ax_gtf.yaxis.set_label_coords(-0.2, 0.5)

    for i, track in enumerate(prepared_tracks):
        ax_track = axes[i + 1]
        track_data = track.get('data', {})
        track_kind = track.get('kind', 'distribution')
        track_label = track.get('label')
        track_color = track.get('color', 'orange')
        track_alpha = track.get('alpha', 0.8)
        y_axis_range = track.get('y_axis_range')
        file_colors = track.get('file_colors') or [track_color]
        file_alphas = track.get('file_alphas') or [track_alpha]
        bed_height = transcript_height / 2.5

        if track_kind == 'coverage':
            series = track_data.get('series')
            if series:
                max_y = 0
                for item in series:
                    x = item.get('x', [])
                    y = item.get('y', [])
                    if len(x) == 0:
                        continue
                    color = item.get('color', track_color)
                    alpha = item.get('alpha', track_alpha)
                    ax_track.fill_between(x, y, color=color, alpha=alpha, step='mid', zorder=3)
                    ax_track.plot(x, y, color=color, lw=0.7, alpha=min(alpha * 1.2, 1.0), zorder=4)
                    if len(y) > 0:
                        max_y = max(max_y, float(np.max(y)))
                if y_axis_range:
                    ax_track.set_ylim(0, y_axis_range)
                else:
                    ax_track.set_ylim(0, max_y * 1.1 if max_y > 0 else 1)
            else:
                x = track_data.get('x', [])
                y = track_data.get('y', [])
                if len(x) > 0:
                    color = file_colors[i % len(file_colors)]
                    alpha = file_alphas[i % len(file_alphas)]
                    ax_track.fill_between(x, y, color=color, alpha=alpha, step='mid', zorder=3)
                    ax_track.plot(x, y, color=color, lw=0.7, alpha=min(alpha * 1.2, 1.0), zorder=4)
                    if y_axis_range:
                        ax_track.set_ylim(0, y_axis_range)
                    else:
                        ax_track.set_ylim(0, np.max(y) * 1.1 if len(y) > 0 else 1)
            ax_track.tick_params(axis='y', which='major', labelsize=8, pad=2)
            ax_track.grid(True, axis='x', alpha=0.25)
        elif track_kind == 'score':
            file_cycle = 0
            for _, bed_elements in track_data.items():
                for bed_element in bed_elements:
                    start = bed_element['start']
                    end = bed_element['end']
                    score = bed_element.get('score', 0.0)
                    color_idx = file_cycle % len(file_colors)
                    color = file_colors[color_idx]
                    alpha = file_alphas[color_idx]
                    file_cycle += 1
                    ax_track.bar((start + end) / 2, score, width=end - start, bottom=0,
                                 color=color, edgecolor='none', linewidth=0, zorder=3, alpha=alpha)
            if y_axis_range:
                ax_track.set_ylim(0, y_axis_range)
            ax_track.tick_params(axis='y', which='major', labelsize=8, pad=2)
            ax_track.grid(True, axis='x', alpha=0.25)
        else:
            bed_names = list(track_data.keys())
            num_bed_names = len(bed_names)
            bed_y_positions = []
            bed_labels = []
            top_margin = 0.3

            for name_idx, name in enumerate(bed_names):
                y_pos = num_bed_names - name_idx - 0.5 + top_margin
                bed_y_positions.append(y_pos)
                bed_labels.append(name)

            for name_idx, (_, bed_elements) in enumerate(track_data.items()):
                y_pos = bed_y_positions[name_idx]
                for bed_element in bed_elements:
                    start = bed_element['start']
                    end = bed_element['end']
                    width = end - start
                    if width == 0:
                        width = 1
                        start = start - 0.5

                    bed_rect = patches.Rectangle(
                        (start, y_pos - bed_height / 2),
                        width,
                        bed_height,
                        linewidth=1,
                        edgecolor=None,
                        facecolor=track_color,
                        alpha=track_alpha,
                        zorder=3
                    )
                    ax_track.add_patch(bed_rect)

            processed_bed_labels = []
            for label in bed_labels:
                if len(label) > 40:
                    processed_bed_labels.append(label[:40] + '...')
                else:
                    processed_bed_labels.append(label)

            ax_track.set_yticks(bed_y_positions)
            ax_track.set_yticklabels(processed_bed_labels)
            if num_bed_names > 0:
                first_element_y = num_bed_names - 0.5 + top_margin
                last_element_y = 1 - 0.5 + top_margin
                ax_track.set_ylim(last_element_y - 0.3, first_element_y + 0.3)
            ax_track.tick_params(axis='y', which='major', labelsize=8, pad=2)
            ax_track.grid(True, axis='x', alpha=0.25)
            ax_track.ticklabel_format(axis='x', style='plain', useOffset=False)

        if effective_track_labels and len(effective_track_labels) > i + 1:
            ax_track.set_ylabel(effective_track_labels[i + 1], fontsize=10)
        elif track_label:
            ax_track.set_ylabel(track_label, fontsize=10)
        else:
            ax_track.set_ylabel(f'Track {i + 1}', fontsize=10)
        ax_track.yaxis.set_label_coords(-0.2, 0.5)

    axes[-1].set_xlabel(f'Genomic Position ({transcript_data["seqname"]})' if 'seqname' in transcript_data else 'Genomic Position', fontsize=10)

    identifier_type = transcript_data.get('identifier_type', 'gene_id')
    original_identifier = transcript_data.get('original_identifier', gene_id)
    gene_identifiers = original_identifier.split(',') if ',' in original_identifier else [original_identifier]

    if not transcript_to_show and identifier_type == 'transcript_id':
        identifier_label = f'Gene ID: {original_identifier}'
        identifier_type = 'gene_id'
    elif transcript_to_show:
        if isinstance(transcript_to_show, str):
            identifier_label = f'Transcript ID: {transcript_to_show}'
            identifier_type = 'transcript_id'
        else:
            if len(transcript_to_show) <= 3:
                transcript_list = ', '.join(transcript_to_show)
            else:
                transcript_list = ', '.join(transcript_to_show[:3]) + f' and {len(transcript_to_show) - 3} more'
            identifier_label = f'{len(transcript_to_show)} Transcripts: {transcript_list}'
            identifier_type = 'transcript_id'
    elif len(gene_identifiers) > 1:
        gene_count = len(gene_identifiers)
        if gene_count == 2:
            identifier_label = f'{gene_count} Genes: {original_identifier}'
        else:
            if len(gene_identifiers) <= 3:
                identifier_label = f'{gene_count} Genes: {original_identifier}'
            else:
                identifier_label = f'{gene_count} Genes: {", ".join(gene_identifiers[:2]).strip()} and {len(gene_identifiers) - 2} more'
    elif identifier_type == 'gene_name':
        identifier_label = f'Gene Name: {original_identifier}'
    elif identifier_type == 'gene_id':
        identifier_label = f'Gene ID: {original_identifier}'
    else:
        identifier_label = f'Gene ID: {original_identifier}'
        identifier_type = 'gene_id'

    if strand == '+':
        title = f"Transcript Structure for {identifier_label} (5' → 3')"
    elif strand == '-':
        title = f"Transcript Structure for {identifier_label} (3' ← 5')"
    else:
        title = f'Transcript Structure for {identifier_label} ({strand} strand)'

    fig.suptitle(title, fontsize=14)

    figure_height = layout['figure_height']
    title_gap = layout['title_gap']
    bottom_margin = layout['bottom_margin']
    track_gap = layout['track_gap']
    axes_height = layout['axes_height']

    right_label_groups = transcript_data.get('right_label_groups', [])
    needs_right_labels = _has_right_labels(prepared_tracks)

    if len(axes) > 1:
        hspace = track_gap * len(axes) / axes_height if axes_height > 0 else 0
    else:
        hspace = 0

    title_height = 0.28
    top_fraction = 1 - (title_gap + title_height) / figure_height
    bottom_fraction = bottom_margin / figure_height
    title_y = 1 - title_height / (2 * figure_height)

    right_margin = 0.9 if needs_right_labels else 0.98

    fig.subplots_adjust(top=top_fraction, bottom=bottom_fraction, right=right_margin, hspace=hspace)
    fig._suptitle.set_y(title_y)

    if needs_right_labels:
        _add_right_side_transcript_labels(fig, axes, right_label_groups)

    for i, ax in enumerate(axes):
        ax.patch.set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.spines['right'].set_visible(True)
        ax.spines['top'].set_visible(True)
        ax.spines['bottom'].set_visible(True)

        if i < len(axes) - 1:
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        else:
            ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True)
            ax.ticklabel_format(axis='x', style='plain', useOffset=False)

    for ax in axes[:-1]:
        ax.set_xlabel('')

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

    sorted_transcripts = sorted(transcripts, key=get_first_exon_start, reverse=reverse)
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
