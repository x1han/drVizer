# Parser Acceleration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Cython-backed and parallel acceleration for GTF, BED, and BAM processing while preserving drVizer's public API, correctness, and pure-Python fallback behavior.

**Architecture:** Keep `gtf_parser.py`, `bed_parser.py`, and `bam_parser.py` as the stable public orchestration layer, and add internal optimized modules for projection, parsing, and BAM parallelism. The implementation lands in small slices: shared projection acceleration first, then packaging and tests, then GTF acceleration, then BAM parallelism, then optional BED parsing acceleration if profiling still shows it matters.

**Tech Stack:** Python, Cython, numpy, pysam, setuptools, pytest

---

## File structure

### Files to create
- `docs/superpowers/plans/2026-04-14-parser-acceleration.md` — this implementation plan
- `pyproject.toml` — build-system requirements for Cython-aware builds and editable installs
- `tests/test_projection_acceleration.py` — regression tests for transcript-to-genomic projection correctness and fallback parity
- `tests/test_gtf_parser_acceleration.py` — parser parity tests for GTF parsing and transcript de-duplication
- `tests/test_bam_parallel.py` — serial vs parallel BAM aggregation tests
- `tests/conftest.py` — shared pytest fixtures and path helpers for parser tests
- `src/drvizer/_cython_projection.pyx` — shared accelerated transcript-to-genomic segment projection
- `src/drvizer/_cython_gtf.pyx` — accelerated GTF hot-loop helpers
- `src/drvizer/_parallel.py` — parallel BAM worker orchestration
- `benchmarks/benchmark_parsers.py` — opt-in benchmark script run with the DRS environment

### Files that will be modified
- `setup.py` — add conditional extension build configuration and package data handling
- `requirements.txt` — keep runtime requirements accurate if needed; do not add build-only requirements here if `pyproject.toml` handles them
- `src/drvizer/gtf_parser.py` — add import fallback, integrate `_cython_gtf`, keep canonical parser state in Python
- `src/drvizer/bed_parser.py` — use `_cython_projection` when available and preserve fallback semantics
- `src/drvizer/bam_parser.py` — use `_cython_projection` when available and add multi-BAM parallel orchestration with serial fallback
- `src/drvizer/__init__.py` — export behavior only if needed to keep package imports stable

### Files to inspect during implementation
- `src/drvizer/api.py` — confirm parser construction remains API-compatible
- `README.md` — only if a brief installation/testing note is needed after implementation
- `test_api_fix.py` and `test_transcript_list.py` — preserve existing standalone validation scripts

---

### Task 1: Add pytest-based regression test scaffolding

**Files:**
- Create: `tests/conftest.py`
- Create: `tests/test_projection_acceleration.py`
- Create: `tests/test_gtf_parser_acceleration.py`
- Modify: `setup.py:1-47`

- [ ] **Step 1: Write the failing test files and fixtures**

```python
# tests/conftest.py
from pathlib import Path
import textwrap

import pytest

from drvizer.gtf_parser import GTFParser


@pytest.fixture
def tmp_gtf(tmp_path):
    path = tmp_path / "test.gtf"
    path.write_text(textwrap.dedent(
        """\
        chr1\tsrc\texon\t100\t149\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1\tsrc\tCDS\t110\t140\t.\t+\t0\tgene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1\tsrc\texon\t200\t249\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1\tsrc\texon\t300\t349\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX2"; gene_name "G1";
        """
    ))
    return path


@pytest.fixture
def tmp_gtf_second(tmp_path):
    path = tmp_path / "test_second.gtf"
    path.write_text(textwrap.dedent(
        """\
        chr1\tsrc\texon\t100\t149\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "G1";
        chr1\tsrc\texon\t400\t449\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX3"; gene_name "G1";
        """
    ))
    return path


@pytest.fixture
def parsed_gtf(tmp_gtf):
    parser = GTFParser(str(tmp_gtf))
    parser.parse_gtf()
    return parser
```

```python
# tests/test_projection_acceleration.py
from drvizer.gtf_parser import GTFParser


def test_segment_projection_splits_across_exons(tmp_gtf):
    parser = GTFParser(str(tmp_gtf))
    parser.parse_gtf()

    chrom, strand, segments = parser.convert_transcript_to_genomic_segments("TX1", 25, 75)

    assert chrom == "chr1"
    assert strand == "+"
    assert segments == [(125, 150), (200, 225)]


def test_legacy_projection_wrapper_returns_outer_bounds(tmp_gtf):
    parser = GTFParser(str(tmp_gtf))
    parser.parse_gtf()

    result = parser.convert_transcript_to_genomic("TX1", 25, 75)

    assert result == ("chr1", 125, 225)
```

```python
# tests/test_gtf_parser_acceleration.py
from drvizer.gtf_parser import GTFParser


def test_parse_gtf_builds_gene_and_transcript_indexes(tmp_gtf):
    parser = GTFParser(str(tmp_gtf))
    parser.parse_gtf()

    assert parser.get_gene_by_name("G1") == "GENE1"
    assert parser.transcript_to_gene["TX1"] == "GENE1"
    assert parser.identify_gene_identifier_type("TX1") == ("transcript_id", "GENE1")


def test_multi_gtf_deduplicates_transcripts_from_earlier_files(tmp_gtf, tmp_gtf_second):
    parser = GTFParser([str(tmp_gtf), str(tmp_gtf_second)])
    parser.parse_gtf()

    transcript_ids = set(parser.gene_transcripts["GENE1"].keys())
    assert transcript_ids == {"TX1", "TX2", "TX3"}
```

- [ ] **Step 2: Run the new tests to verify baseline behavior and environment readiness**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest tests/test_projection_acceleration.py tests/test_gtf_parser_acceleration.py -v
```

Expected: either PASS on current behavior or FAIL only due to missing `pytest` / import-path setup that the next step will fix.

- [ ] **Step 3: Make the minimal test-enablement changes in packaging/import setup**

```python
# setup.py
from pathlib import Path
from setuptools import setup, find_packages

ROOT = Path(__file__).parent
long_description = (ROOT / "README.md").read_text(encoding="utf-8")
requirements = [
    line.strip()
    for line in (ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines()
    if line.strip() and not line.startswith("#")
]

setup(
    name="drvizer",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A Python tool for parsing GTF/BED files and visualizing gene transcript structures",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/drvizer",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "drvizer=drvizer.cli:main",
        ],
    },
    keywords="bioinformatics, genomics, transcriptomics, gtf, bed, visualization",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/drvizer/issues",
        "Source": "https://github.com/yourusername/drvizer",
        "Documentation": "https://github.com/yourusername/drvizer#readme",
    },
)
```

- [ ] **Step 4: Install pytest if needed and rerun the tests until they pass**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pip install pytest
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest tests/test_projection_acceleration.py tests/test_gtf_parser_acceleration.py -v
```

Expected: PASS for both test files.

- [ ] **Step 5: Commit the testing scaffold**

```bash
git add tests/conftest.py tests/test_projection_acceleration.py tests/test_gtf_parser_acceleration.py setup.py
git commit -m "test: add parser acceleration regression scaffolding"
```

### Task 2: Add a shared optional Cython projection module with pure-Python fallback

**Files:**
- Create: `src/drvizer/_cython_projection.pyx`
- Modify: `src/drvizer/gtf_parser.py:323-415`
- Modify: `src/drvizer/bed_parser.py:1-220`
- Modify: `src/drvizer/bam_parser.py:1-155`
- Test: `tests/test_projection_acceleration.py`

- [ ] **Step 1: Extend the projection tests with fallback-parity expectations**

```python
# append to tests/test_projection_acceleration.py

def test_find_transcript_returns_cached_mapping(parsed_gtf):
    gene_id, transcript = parsed_gtf.find_transcript("TX1")

    assert gene_id == "GENE1"
    assert transcript["seqname"] == "chr1"


def test_projection_exons_cache_is_populated(parsed_gtf):
    parsed_gtf.convert_transcript_to_genomic_segments("TX1", 25, 75)

    transcript = parsed_gtf.gene_transcripts["GENE1"]["TX1"]
    assert "_projection_exons" in transcript
    assert len(transcript["_projection_exons"]) == 2
```

- [ ] **Step 2: Run the projection tests and verify they pass before introducing the optimized path**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest tests/test_projection_acceleration.py -v
```

Expected: PASS.

- [ ] **Step 3: Add the optional Cython projection module and Python import fallback**

```cython
# src/drvizer/_cython_projection.pyx
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
```

```python
# gtf_parser.py import section
try:
    from ._cython_projection import project_segments as _project_segments_fast
except ImportError:
    _project_segments_fast = None
```

```python
# in GTFParser.convert_transcript_to_genomic_segments
projection_exons = self._get_projection_exons(transcript_structure)
if _project_segments_fast is not None:
    segments = _project_segments_fast(projection_exons, genomic_strand, interval_start, interval_end)
else:
    segments = []
    for exon, exon_t_start, exon_t_end in projection_exons:
        overlap_start = max(interval_start, exon_t_start)
        overlap_end = min(interval_end, exon_t_end)
        if overlap_start < overlap_end:
            if genomic_strand == '+':
                genomic_start = exon['start'] + (overlap_start - exon_t_start)
                genomic_end = exon['start'] + (overlap_end - exon_t_start)
            else:
                genomic_start = exon['end'] - (overlap_end - exon_t_start) + 1
                genomic_end = exon['end'] - (overlap_start - exon_t_start) + 1
            if genomic_start > genomic_end:
                genomic_start, genomic_end = genomic_end, genomic_start
            segments.append((genomic_start, genomic_end))
    segments.sort(key=lambda x: x[0])
return seqname, genomic_strand, segments
```

- [ ] **Step 4: Update BED and BAM to use the unchanged parser API and verify the projection tests still pass**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest tests/test_projection_acceleration.py tests/test_gtf_parser_acceleration.py -v
```

Expected: PASS, with no change in public behavior.

- [ ] **Step 5: Commit the shared projection acceleration layer**

```bash
git add src/drvizer/_cython_projection.pyx src/drvizer/gtf_parser.py src/drvizer/bed_parser.py src/drvizer/bam_parser.py tests/test_projection_acceleration.py
git commit -m "feat: add optional accelerated projection layer"
```

### Task 3: Add Cython-aware build configuration with pure-Python fallback

**Files:**
- Create: `pyproject.toml`
- Modify: `setup.py:1-47`
- Test: `tests/test_projection_acceleration.py`

- [ ] **Step 1: Write the failing build/install expectation as a smoke test note in the plan execution log**

```python
# No repository file change in this step.
# Expected failure to reproduce before the build config change:
# python -m pip install -e . may not know how to build .pyx files.
```

- [ ] **Step 2: Reproduce the editable-install/build gap in the DRS environment**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pip install -e .
```

Expected: either editable install succeeds without extensions, or fails to build `.pyx` support, which justifies the next packaging change.

- [ ] **Step 3: Add `pyproject.toml` and conditional extension building in `setup.py`**

```toml
# pyproject.toml
[build-system]
requires = ["setuptools>=64", "wheel", "Cython>=3.0"]
build-backend = "setuptools.build_meta"
```

```python
# setup.py
from pathlib import Path
from setuptools import Extension, setup, find_packages

ROOT = Path(__file__).parent
long_description = (ROOT / "README.md").read_text(encoding="utf-8")
requirements = [
    line.strip()
    for line in (ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines()
    if line.strip() and not line.startswith("#")
]

extensions = [
    Extension("drvizer._cython_projection", ["src/drvizer/_cython_projection.pyx"]),
]

try:
    from Cython.Build import cythonize
    ext_modules = cythonize(extensions, compiler_directives={"language_level": "3"})
except ImportError:
    ext_modules = []

setup(
    name="drvizer",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A Python tool for parsing GTF/BED files and visualizing gene transcript structures",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/drvizer",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    ext_modules=ext_modules,
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "drvizer=drvizer.cli:main",
        ],
    },
    keywords="bioinformatics, genomics, transcriptomics, gtf, bed, visualization",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/drvizer/issues",
        "Source": "https://github.com/yourusername/drvizer",
        "Documentation": "https://github.com/yourusername/drvizer#readme",
    },
)
```

- [ ] **Step 4: Reinstall in the DRS environment and rerun the projection tests**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pip install -e .
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest tests/test_projection_acceleration.py -v
```

Expected: editable install succeeds and the projection tests still PASS.

- [ ] **Step 5: Commit the packaging support for optional Cython extensions**

```bash
git add pyproject.toml setup.py
git commit -m "build: add optional cython extension packaging"
```

### Task 4: Add accelerated GTF parsing helpers behind the existing parser API

**Files:**
- Create: `src/drvizer/_cython_gtf.pyx`
- Modify: `src/drvizer/gtf_parser.py:1-250`
- Test: `tests/test_gtf_parser_acceleration.py`

- [ ] **Step 1: Add a failing GTF parity test for cached parser output shape**

```python
# append to tests/test_gtf_parser_acceleration.py

def test_get_transcript_data_keeps_expected_output_shape(tmp_gtf):
    parser = GTFParser(str(tmp_gtf))
    parser.parse_gtf()

    data = parser.get_transcript_data("GENE1")

    assert data["gene_id"] == "GENE1"
    assert data["seqname"] == "chr1"
    assert [tx["transcript_id"] for tx in data["transcripts"]] == ["TX1", "TX2"]
```

- [ ] **Step 2: Run the GTF tests to verify the current baseline**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest tests/test_gtf_parser_acceleration.py -v
```

Expected: PASS.

- [ ] **Step 3: Add `_cython_gtf.pyx` and integrate it as an optional hot-loop helper**

```cython
# src/drvizer/_cython_gtf.pyx
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
```

```python
# gtf_parser.py import section
try:
    from ._cython_gtf import parse_attributes_fast as _parse_attributes_fast_impl
    from ._cython_gtf import parse_gtf_chunk as _parse_gtf_chunk_impl
except ImportError:
    _parse_attributes_fast_impl = None
    _parse_gtf_chunk_impl = None
```

```python
# gtf_parser.py
    def _parse_attributes_fast(self, attribute_string):
        if _parse_attributes_fast_impl is not None:
            return _parse_attributes_fast_impl(attribute_string)
        attributes = {}
        pattern = r'(\w+) "([^"]*)";?'
        matches = re.findall(pattern, attribute_string)
        for key, value in matches:
            attributes[key] = value
        return attributes
```

```python
# inside _process_chunk in gtf_parser.py
        chunk_rows = _parse_gtf_chunk_impl(chunk_lines) if _parse_gtf_chunk_impl is not None else None
        if chunk_rows is not None:
            iterable_rows = chunk_rows
        else:
            iterable_rows = []
            for line in chunk_lines:
                parts = line.strip().split('\t')
                if len(parts) >= 9 and (parts[2] == 'exon' or parts[2] == 'CDS'):
                    iterable_rows.append((parts[0], parts[2], int(parts[3]), int(parts[4]), parts[6], parts[8]))

        for seqname, feature, start_i, end_i, strand, attribute_string in iterable_rows:
            attributes = self._parse_attributes_fast(attribute_string)
            gene_id_row = attributes.get('gene_id', None)
            transcript_id = attributes.get('transcript_id', None)
            gene_name = attributes.get('gene_name', None)
            # keep the existing downstream assembly logic unchanged
```

- [ ] **Step 4: Reinstall editable mode and rerun the GTF tests**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pip install -e .
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest tests/test_gtf_parser_acceleration.py tests/test_projection_acceleration.py -v
```

Expected: PASS, with no behavior change.

- [ ] **Step 5: Commit the optional accelerated GTF parser helpers**

```bash
git add src/drvizer/_cython_gtf.pyx src/drvizer/gtf_parser.py tests/test_gtf_parser_acceleration.py setup.py
git commit -m "feat: add optional accelerated gtf parsing"
```

### Task 5: Add multi-BAM parallel aggregation with serial fallback

**Files:**
- Create: `src/drvizer/_parallel.py`
- Create: `tests/test_bam_parallel.py`
- Modify: `src/drvizer/bam_parser.py:1-155`
- Test: `tests/test_bam_parallel.py`

- [ ] **Step 1: Write the failing BAM aggregation and fallback tests**

```python
# tests/test_bam_parallel.py
import numpy as np

from drvizer.bam_parser import BAMParser


class FakeSam:
    def __init__(self, blocks):
        self._blocks = blocks
        self.header = {'SQ': [{'SN': 'TX1', 'LN': 100}]}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def fetch(self, chrom, start, end):
        return self._blocks


class FakeRead:
    def __init__(self, blocks, reference_start=0, reference_end=100):
        self._blocks = blocks
        self.reference_start = reference_start
        self.reference_end = reference_end

    def get_blocks(self):
        return self._blocks


def test_genomic_serial_and_parallel_paths_match(monkeypatch):
    parser = BAMParser(["a.bam", "b.bam"])

    fake_reads = [FakeRead([(10, 20)]), FakeRead([(12, 18)])]

    def fake_alignment_file(path, mode):
        return FakeSam(fake_reads)

    monkeypatch.setattr("drvizer.bam_parser.pysam.AlignmentFile", fake_alignment_file)

    x, y = parser.get_coverage_in_region("chr1", 0, 30)

    assert len(x) == 30
    assert np.max(y) >= 1
```

- [ ] **Step 2: Run the BAM tests to capture the current baseline**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest tests/test_bam_parallel.py -v
```

Expected: FAIL because the parallel orchestration layer does not exist yet.

- [ ] **Step 3: Add the parallel worker module and integrate serial fallback in `BAMParser`**

```python
# src/drvizer/_parallel.py
from concurrent.futures import ProcessPoolExecutor


def run_parallel(worker, tasks, max_workers=None):
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        return list(executor.map(worker, tasks))
```

```python
# bam_parser.py import section
try:
    from ._parallel import run_parallel
except ImportError:
    run_parallel = None
```

```python
# bam_parser.py helper outline
    def _combine_coverages(self, coverages):
        total = coverages[0].copy()
        for coverage in coverages[1:]:
            total += coverage
        if self.aggregate_method == 'mean' and len(coverages) > 1:
            total = total.astype(np.float64) / len(coverages)
        return total
```

```python
# bam_parser.py usage outline inside get_coverage_in_region / transcript path
if len(self.bam_paths) > 1 and run_parallel is not None:
    try:
        coverages = run_parallel(worker_function, task_payloads)
        coverage = self._combine_coverages(coverages)
    except Exception:
        coverage = self._run_serial_coverage(task_payloads)
else:
    coverage = self._run_serial_coverage(task_payloads)
```

- [ ] **Step 4: Run the BAM tests and parser regression tests**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest tests/test_bam_parallel.py tests/test_projection_acceleration.py tests/test_gtf_parser_acceleration.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit the parallel BAM aggregation support**

```bash
git add src/drvizer/_parallel.py src/drvizer/bam_parser.py tests/test_bam_parallel.py
git commit -m "feat: add parallel bam aggregation"
```

### Task 6: Add optional BED parsing acceleration if profiling still shows a bottleneck

**Files:**
- Create: `src/drvizer/_cython_bed.pyx`
- Modify: `src/drvizer/bed_parser.py:1-220`
- Test: `tests/test_projection_acceleration.py`
- Test: `benchmarks/benchmark_parsers.py`

- [ ] **Step 1: Run the benchmark script without BED acceleration and record the baseline**

```python
# benchmarks/benchmark_parsers.py
from pathlib import Path
from time import perf_counter

from drvizer.gtf_parser import GTFParser
from drvizer.bed_parser import BEDParser


def timed(label, fn):
    start = perf_counter()
    fn()
    elapsed = perf_counter() - start
    print(f"{label}: {elapsed:.4f}s")


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    print(f"benchmark root: {root}")
```
```

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python benchmarks/benchmark_parsers.py
```

Expected: a baseline timing printout you can compare after the change.

- [ ] **Step 2: Add `_cython_bed.pyx` only if profiling shows BED parsing is still a meaningful bottleneck**

```cython
# src/drvizer/_cython_bed.pyx
# cython: language_level=3

def parse_bed_fields(str line):
    parts = line.strip().split('\t')
    if len(parts) < 4:
        return None
    try:
        start = int(parts[1])
        end = int(parts[2])
    except (ValueError, IndexError):
        return None
    return {
        "chrom": parts[0],
        "start": start,
        "end": end,
        "name": parts[3] if len(parts) > 3 else ".",
        "score": float(parts[4]) if len(parts) > 4 else 0.0,
        "strand": parts[5] if len(parts) > 5 else ".",
        "thickStart": int(parts[6]) if len(parts) > 6 else start,
        "thickEnd": int(parts[7]) if len(parts) > 7 else end,
        "itemRgb": parts[8] if len(parts) > 8 else "0",
    }
```

```python
# bed_parser.py import section
try:
    from ._cython_bed import parse_bed_fields as _parse_bed_fields_fast
except ImportError:
    _parse_bed_fields_fast = None
```

- [ ] **Step 3: Wire the optional BED fast path while preserving all existing semantics**

```python
# inside bed_parser.py _process_bed_file
record = _parse_bed_fields_fast(line) if _parse_bed_fields_fast is not None else None
if record is None:
    parts = line.strip().split('\t')
    if len(parts) < 4:
        continue
    try:
        record = {
            'chrom': parts[0],
            'start': int(parts[1]),
            'end': int(parts[2]),
            'name': parts[3] if len(parts) > 3 else '.',
            'score': float(parts[4]) if len(parts) > 4 else 0.0,
            'strand': parts[5] if len(parts) > 5 else '.',
            'thickStart': int(parts[6]) if len(parts) > 6 else int(parts[1]),
            'thickEnd': int(parts[7]) if len(parts) > 7 else int(parts[2]),
            'itemRgb': parts[8] if len(parts) > 8 else '0',
        }
    except ValueError:
        continue
```

- [ ] **Step 4: Rerun tests and benchmarks, and keep the change only if it is measurably useful**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest tests/test_projection_acceleration.py tests/test_gtf_parser_acceleration.py tests/test_bam_parallel.py -v
/datf/hanxi/software/miniconda3/envs/DRS/bin/python benchmarks/benchmark_parsers.py
```

Expected: tests PASS, and benchmark results show BED parsing is not slower than the baseline.

- [ ] **Step 5: Commit the optional BED acceleration only if Step 4 shows a real benefit**

```bash
git add src/drvizer/_cython_bed.pyx src/drvizer/bed_parser.py benchmarks/benchmark_parsers.py setup.py
git commit -m "feat: add optional bed parsing acceleration"
```

### Task 7: Final verification in the DRS environment

**Files:**
- Modify: none required unless verification reveals a bug
- Test: `tests/test_projection_acceleration.py`
- Test: `tests/test_gtf_parser_acceleration.py`
- Test: `tests/test_bam_parallel.py`
- Test: `test_api_fix.py`
- Test: `test_transcript_list.py`

- [ ] **Step 1: Run the full targeted pytest suite**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pytest tests/test_projection_acceleration.py tests/test_gtf_parser_acceleration.py tests/test_bam_parallel.py -v
```

Expected: PASS.

- [ ] **Step 2: Run the existing repository validation scripts**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python test_api_fix.py
/datf/hanxi/software/miniconda3/envs/DRS/bin/python test_transcript_list.py
```

Expected: both scripts complete successfully.

- [ ] **Step 3: Reinstall editable mode in the DRS environment and confirm the package still imports**

Run:
```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -m pip install -e .
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -c "from drvizer import DrViz; print('ok')"
```

Expected: editable install succeeds and the import command prints `ok`.

- [ ] **Step 4: If any verification step fails, make the minimal fix and rerun only the affected verification command before continuing**

```python
# No predetermined repository code for this step.
# Apply the smallest fix required by the failing verification step,
# then rerun the exact command that failed.
```

- [ ] **Step 5: Commit the final verified parser acceleration work**

```bash
git add pyproject.toml setup.py src/drvizer tests benchmarks
git commit -m "feat: accelerate parser hot paths with cython and parallelism"
```

---

## Self-review

### Spec coverage
- shared transcript projection acceleration: covered in Task 2
- GTF parser acceleration: covered in Task 4
- BAM multi-file parallel aggregation: covered in Task 5
- optional BED parser acceleration only if profiling justifies it: covered in Task 6
- build and fallback behavior: covered in Tasks 2 and 3
- pytest in DRS environment: covered in Tasks 1, 3, 4, 5, and 7
- regression protection for transcript-coordinate intron correctness: covered in Tasks 1, 2, and 7

### Placeholder scan
- No `TODO`, `TBD`, or deferred placeholders remain in task steps.
- The only intentionally open conditional is Task 6, which the spec explicitly defines as optional based on profiling.

### Type consistency
- Shared projection helper name remains `project_segments`
- GTF helper names remain `parse_attributes_fast` and `parse_gtf_chunk`
- Parallel entry point remains `run_parallel`
- Parser public methods remain `convert_transcript_to_genomic_segments`, `get_coverage_in_region`, and `get_coverage_for_transcripts`
