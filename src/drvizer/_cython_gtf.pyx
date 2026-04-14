# cython: language_level=3
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
