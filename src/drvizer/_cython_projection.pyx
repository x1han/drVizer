# cython: language_level=3

def project_segments(list projection_exons, str strand, int transcript_start, int transcript_end):
    cdef int interval_start = transcript_start if transcript_start <= transcript_end else transcript_end
    cdef int interval_end = transcript_end if transcript_end >= transcript_start else transcript_start
    cdef list segments = []
    cdef object exon
    cdef int exon_t_start
    cdef int exon_t_end
    cdef int overlap_start
    cdef int overlap_end
    cdef int genomic_start
    cdef int genomic_end

    if interval_start == interval_end:
        return segments

    for exon, exon_t_start, exon_t_end in projection_exons:
        overlap_start = interval_start if interval_start > exon_t_start else exon_t_start
        overlap_end = interval_end if interval_end < exon_t_end else exon_t_end
        if overlap_start < overlap_end:
            if strand == "+":
                genomic_start = exon["start"] + (overlap_start - exon_t_start)
                genomic_end = exon["start"] + (overlap_end - exon_t_start)
            else:
                genomic_start = exon["end"] - (overlap_end - exon_t_start) + 1
                genomic_end = exon["end"] - (overlap_start - exon_t_start) + 1
            if genomic_start > genomic_end:
                genomic_start, genomic_end = genomic_end, genomic_start
            segments.append((genomic_start, genomic_end))

    segments.sort(key=lambda x: x[0])
    return segments
