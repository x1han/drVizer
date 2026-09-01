# cython: language_level=3

from libc.stdint cimport int64_t


def project_segments(list projection_exons, str strand, int transcript_start, int transcript_end):
    cdef int64_t interval_start = transcript_start if transcript_start <= transcript_end else transcript_end
    cdef int64_t interval_end = transcript_end if transcript_end >= transcript_start else transcript_start
    cdef list segments = []
    cdef object exon
    cdef int64_t exon_t_start
    cdef int64_t exon_t_end
    cdef int64_t overlap_start
    cdef int64_t overlap_end
    cdef int64_t genomic_start
    cdef int64_t genomic_end

    if interval_start == interval_end:
        return segments

    for exon, exon_t_start, exon_t_end in projection_exons:
        overlap_start = interval_start if interval_start > exon_t_start else exon_t_start
        overlap_end = interval_end if interval_end < exon_t_end else exon_t_end
        if overlap_start < overlap_end:
            if strand == "+":
                # 0-based half-open projection on the GTF start side.
                # GTF 'start' is 1-based inclusive (e.g. 100 -> 1-based start),
                # so the equivalent 0-based half-open origin is (exon['start'] - 1).
                genomic_start = (exon["start"] - 1) + (overlap_start - exon_t_start)
                genomic_end = (exon["start"] - 1) + (overlap_end - exon_t_start)
            else:
                # Minus strand: 0-based half-open projection that mirrors the
                # plus strand by reading the exon right-edge in 0-based half-
                # open coordinates. Width is preserved across the convention
                # switch; the right edge drops the prior + 1 quirk.
                genomic_start = exon["end"] - (overlap_end - exon_t_start)
                genomic_end = exon["end"] - (overlap_start - exon_t_start)
            if genomic_start > genomic_end:
                genomic_start, genomic_end = genomic_end, genomic_start
            segments.append((genomic_start, genomic_end))

    segments.sort(key=lambda x: x[0])
    return segments
