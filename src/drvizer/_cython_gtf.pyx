# cython: language_level=3
#
# P1-3 (Phase 3) -- pattern (a): Python-side validation BEFORE calling this
# Cython kernel. Cython cannot reliably catch arbitrary C-level exceptions
# inside a tight loop, so the gtf_parser.py wrapper wraps the call to
# parse_gtf_chunk in try/except and falls back to a per-row Python branch
# with try/except (ValueError, TypeError) per row. The kernel below is a
# pure best-effort transform; malformed rows are filtered upstream so the
# whole 10k-line chunk is never silently aborted on one bad row.
import re

_ATTRIBUTE_PATTERN = re.compile(r'(\w+) "([^"]*)";?')


def parse_attributes_fast(str attribute_string):
    return {key: value for key, value in _ATTRIBUTE_PATTERN.findall(attribute_string)}


def parse_gtf_chunk(list chunk_lines):
    cdef list rows = []
    cdef str line
    cdef list parts
    for line in chunk_lines:
        parts = line.strip().split('\t')
        if len(parts) >= 9 and (parts[2] == 'exon' or parts[2] == 'CDS'):
            rows.append((parts[0], parts[2], int(parts[3]), int(parts[4]), parts[6], parts[8]))
    return rows
