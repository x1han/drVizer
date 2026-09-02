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
    cdef str start_str, end_str
    for line in chunk_lines:
        parts = line.strip().split('\t')
        if len(parts) >= 9 and (parts[2] == 'exon' or parts[2] == 'CDS'):
            start_str = parts[3]
            end_str = parts[4]
            # Explicit integer validation. Cython compiles
            # int(typed_str_var) to C-level strtol on some
            # toolchains, which silently coerces non-integer input
            # to 0 instead of raising ValueError -- the wrapper in
            # gtf_parser.py relies on the ValueError to detect
            # malformed rows and fall back to per-row Python
            # parsing. Validate first to make the failure mode
            # deterministic across Cython / setuptools / Python
            # combinations.
            if not start_str.lstrip('-').isdigit() or not end_str.lstrip('-').isdigit():
                raise ValueError(
                    f"GTF line has non-integer coordinate: {start_str!r}-{end_str!r}"
                )
            rows.append((parts[0], parts[2], int(start_str), int(end_str), parts[6], parts[8]))
    return rows
