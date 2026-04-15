from collections import defaultdict
import gzip

from libc.stdlib cimport strtod, strtol
from cpython.bytes cimport PyBytes_AS_STRING, PyBytes_GET_SIZE


cdef inline int _parse_optional_int_text(const char* text, int text_len, int default_value):
    cdef char* endptr
    cdef long value
    if text_len == 0:
        return default_value
    value = strtol(text, &endptr, 10)
    if endptr == text or endptr[0] != 0:
        return default_value
    return <int>value


cdef inline object _parse_optional_float_text(const char* text, int text_len, object default_value):
    cdef char* endptr
    cdef double value
    if text_len == 0:
        return default_value
    value = strtod(text, &endptr)
    if endptr == text or endptr[0] != 0:
        return default_value
    return value


cdef inline object _decode_field(bytes line_bytes, Py_ssize_t start, Py_ssize_t end):
    if end <= start:
        return ''
    return (<bytes>line_bytes[start:end]).decode('utf-8')


cdef inline object _build_record_from_line(bytes line_bytes):
    cdef const char* buffer = PyBytes_AS_STRING(line_bytes)
    cdef Py_ssize_t length = PyBytes_GET_SIZE(line_bytes)
    cdef Py_ssize_t field_start = 0
    cdef Py_ssize_t idx
    cdef Py_ssize_t field_count = 0
    cdef Py_ssize_t starts[9]
    cdef Py_ssize_t ends[9]
    cdef object chrom
    cdef object name
    cdef object strand
    cdef object item_rgb
    cdef object score
    cdef int start
    cdef int end
    cdef int thick_start
    cdef int thick_end

    cdef char tab = b'\t'

    for idx in range(length):
        if buffer[idx] == tab:
            if field_count < 9:
                starts[field_count] = field_start
                ends[field_count] = idx
            field_count += 1
            field_start = idx + 1

    if field_count < 9:
        starts[field_count] = field_start
        ends[field_count] = length
    field_count += 1

    if field_count < 4:
        return None

    start = _parse_optional_int_text(buffer + starts[1], <int>(ends[1] - starts[1]), -1)
    if start < 0:
        return None
    end = _parse_optional_int_text(buffer + starts[2], <int>(ends[2] - starts[2]), -1)
    if end < 0:
        return None

    chrom = _decode_field(line_bytes, starts[0], ends[0])
    name = _decode_field(line_bytes, starts[3], ends[3]) if field_count > 3 else '.'
    score = _parse_optional_float_text(buffer + starts[4], <int>(ends[4] - starts[4]), 0.0) if field_count >= 5 else 0.0
    strand = _decode_field(line_bytes, starts[5], ends[5]) if field_count > 5 else '.'
    thick_start = _parse_optional_int_text(buffer + starts[6], <int>(ends[6] - starts[6]), start) if field_count >= 7 else start
    thick_end = _parse_optional_int_text(buffer + starts[7], <int>(ends[7] - starts[7]), end) if field_count >= 8 else end
    item_rgb = _decode_field(line_bytes, starts[8], ends[8]) if field_count >= 9 else '0'

    return {
        'chrom': chrom,
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
    cdef dict grouped = defaultdict(list)
    cdef object bed_file_path
    cdef object handle
    cdef object raw_line
    cdef bytes line_bytes
    cdef object record

    for bed_file_path in bed_file_paths:
        try:
            handle = gzip.open(bed_file_path, 'rb') if bed_file_path.endswith('.gz') else open(bed_file_path, 'rb')
            with handle:
                for raw_line in handle:
                    line_bytes = raw_line.strip()
                    if not line_bytes or line_bytes.startswith(b'#'):
                        continue

                    record = _build_record_from_line(line_bytes)
                    if record is None or not _record_in_region(record, region):
                        continue

                    grouped[record['chrom']].append(record)
        except Exception as e:
            raise ValueError(f"Error reading BED file {bed_file_path}: {e}") from e

    return grouped
