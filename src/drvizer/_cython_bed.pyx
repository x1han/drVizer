from collections import defaultdict
import gzip


def _parse_optional_int(value, default_value):
    try:
        return int(value)
    except ValueError:
        return default_value


def _parse_optional_float(value, default_value):
    try:
        return float(value)
    except ValueError:
        return default_value


def _parse_record(parts):
    if len(parts) < 4:
        return None

    try:
        start = int(parts[1])
        end = int(parts[2])
    except (ValueError, IndexError):
        return None

    name = parts[3] if len(parts) > 3 else '.'
    score = _parse_optional_float(parts[4], 0.0) if len(parts) >= 5 else 0.0
    strand = parts[5] if len(parts) > 5 else '.'
    thick_start = _parse_optional_int(parts[6], start) if len(parts) >= 7 else start
    thick_end = _parse_optional_int(parts[7], end) if len(parts) >= 8 else end
    item_rgb = parts[8] if len(parts) >= 9 else '0'

    return {
        'chrom': parts[0],
        'start': start,
        'end': end,
        'name': name,
        'score': score,
        'strand': strand,
        'thickStart': thick_start,
        'thickEnd': thick_end,
        'itemRgb': item_rgb,
    }


def _record_in_region(record, region):
    if region is None:
        return True

    chrom, start, end = region
    return not (
        record['chrom'] != chrom or
        record['end'] < start or
        record['start'] > end
    )


def parse_bed_records(bed_file_paths, region=None):
    grouped = defaultdict(list)

    for bed_file_path in bed_file_paths:
        open_func = gzip.open if bed_file_path.endswith('.gz') else open

        try:
            with open_func(bed_file_path, 'rt') as handle:
                for raw_line in handle:
                    if raw_line.startswith('#') or raw_line.strip() == '':
                        continue

                    record = _parse_record(raw_line.strip().split('\t'))
                    if record is None or not _record_in_region(record, region):
                        continue

                    grouped[record['chrom']].append(record)
        except Exception as e:
            raise ValueError(f"Error reading BED file {bed_file_path}: {e}") from e

    return grouped
