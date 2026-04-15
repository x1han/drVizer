from collections import defaultdict
import gzip


def _is_comment_or_blank(line):
    return line.startswith('#') or line.strip() == ''


def _parse_bed_fields(parts):
    if len(parts) < 4:
        return None

    record = {}
    record['chrom'] = parts[0]
    try:
        record['start'] = int(parts[1])
        record['end'] = int(parts[2])
    except (ValueError, IndexError):
        return None

    record['name'] = parts[3] if len(parts) > 3 else '.'

    if len(parts) >= 5:
        try:
            record['score'] = float(parts[4])
        except ValueError:
            record['score'] = 0.0
    else:
        record['score'] = 0.0

    record['strand'] = parts[5] if len(parts) > 5 else '.'

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

    record['itemRgb'] = parts[8] if len(parts) >= 9 else '0'
    return record


def _record_in_region(record, region):
    if region is None:
        return True

    chrom, start, end = region
    return not (
        record['chrom'] != chrom or
        record['end'] < start or
        record['start'] > end
    )


def _open_text_file(bed_file_path, mode='rt'):
    open_func = gzip.open if bed_file_path.endswith('.gz') else open
    return open_func(bed_file_path, mode)


def _normalize_bed_read_error(bed_file_path, error):
    return ValueError(f"Error reading BED file {bed_file_path}: {error}")


def parse_bed_records_python(bed_file_paths, region=None):
    grouped = defaultdict(list)

    for bed_file_path in bed_file_paths:
        try:
            with _open_text_file(bed_file_path, 'rt') as handle:
                for line in handle:
                    if _is_comment_or_blank(line):
                        continue

                    record = _parse_bed_fields(line.strip().split('\t'))
                    if record is None or not _record_in_region(record, region):
                        continue

                    grouped[record['chrom']].append(record)
        except Exception as e:
            raise _normalize_bed_read_error(bed_file_path, e) from e

    return grouped


try:
    from ._cython_bed import parse_bed_records
except ImportError:
    parse_bed_records = None
    _CYTHON_BED_AVAILABLE = False
else:
    _CYTHON_BED_AVAILABLE = True


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
        if not self._parsed or (chrom is not None and start is not None and end is not None) or len(self.anno_data) == 0:
            if (chrom is not None and start is not None and end is not None) or not self._parsed:
                self.anno_data.clear()
            self._process_bed_file(chrom, start, end)
            if not (chrom is not None and start is not None and end is not None):
                self._parsed = True

        return dict(self.anno_data)

    def prepare_track(self, gtf_parser):
        """Prepare BED-backed track data after the builder has loaded the GTF."""
        if self.transcript_coord:
            self.gtf_parser = gtf_parser
        self.parse_bed()
        return self

    def _process_bed_file(self, chrom=None, start=None, end=None):
        """
        Process BED file and extract annotation information.

        Args:
            chrom (str, optional): Chromosome name to filter
            start (int, optional): Start position to filter
            end (int, optional): End position to filter
        """
        region = (chrom, start, end) if chrom is not None and start is not None and end is not None else None

        if not self.transcript_coord:
            grouped = self._parse_genomic_records(region)
            for record_chrom, records in grouped.items():
                self.anno_data[record_chrom].extend(records)
            return

        grouped = parse_bed_records_python(self.bed_file_paths, region)
        for _, records in grouped.items():
            for record in records:
                if not self.gtf_parser:
                    raise ValueError("gtf_parser is required when transcript_coord is True")

                transcript_id = record['chrom']
                result = self.gtf_parser.convert_transcript_to_genomic_segments(
                    transcript_id, record['start'], record['end']
                )
                if not result:
                    continue

                chrom, genomic_strand, segments = result
                for genomic_start, genomic_end in segments:
                    segment_record = record.copy()
                    segment_record['chrom'] = chrom
                    segment_record['start'] = genomic_start
                    segment_record['end'] = genomic_end
                    segment_record['strand'] = genomic_strand
                    self.anno_data[segment_record['chrom']].append(segment_record)

    def _parse_genomic_records(self, region=None):
        if _CYTHON_BED_AVAILABLE and parse_bed_records is not None:
            try:
                return parse_bed_records(self.bed_file_paths, region)
            except Exception as e:
                if isinstance(e, ValueError) and str(e).startswith("Error reading BED file "):
                    raise
                return parse_bed_records_python(self.bed_file_paths, region)
        return parse_bed_records_python(self.bed_file_paths, region)

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
        if not self._parsed:
            self.parse_bed()

        anno_in_region = []
        if chrom in self.anno_data:
            for anno in self.anno_data[chrom]:
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
        if not self._parsed:
            self.parse_bed()

        grouped_anno = {}
        if chrom in self.anno_data:
            for anno in self.anno_data[chrom]:
                anno_start, anno_end = anno['start'], anno['end']
                if anno_start > anno_end:
                    anno_start, anno_end = anno_end, anno_start

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
