# drVizer parser acceleration design

Date: 2026-04-14

## Goal

Speed up GTF, BED, and BAM processing in drVizer without changing the public API or changing plotting semantics.

The accelerated design must preserve:
- `DrViz().load_gtf(...).add_bed_track(...).add_bam_track(...).build()` usage
- existing CLI behavior
- transcript-coordinate BED/BAM projection into genomic space
- current `sum` and `mean` aggregation semantics
- fallback behavior when optimized components are unavailable

## Constraints

- The user allows adding compiled dependencies, including Cython.
- Validation should use `/datf/hanxi/software/miniconda3/envs/DRS/bin/python`.
- It is acceptable to install `pytest` into that environment for testing.
- The package must remain runnable even when compiled extensions are unavailable.
- The existing transcript-coordinate intron fix must remain correct: transcript-coordinate BED/BAM intervals must expand exon-by-exon and must not collapse across introns.

## Recommended approach

Use a hybrid design:
- accelerate GTF and BED hot loops with Cython
- keep BAM reading on top of `pysam`
- add Python-level parallelism for multi-BAM coverage aggregation
- add a shared optimized transcript-to-genomic projection layer used by BED and transcript-coordinate BAM

This provides meaningful speedups while avoiding a risky rewrite of BAM decoding.

## Alternatives considered

### Option A: GTF/BED in Cython, BAM via `pysam` plus Python parallelism
Recommended.

Pros:
- targets the hottest Python parsing loops directly
- preserves current architecture well
- avoids re-implementing BAM decoding
- gives clear reuse via a shared projection layer

Cons:
- adds compiled build requirements
- requires fallback logic and packaging updates

### Option B: parallelize everything, keep pure Python
Not recommended.

Pros:
- lower build complexity
- smaller packaging change

Cons:
- GTF/BED text parsing remains Python-bound
- parallel overhead may dilute gains
- lower upside on large annotation files

### Option C: heavy Cythonization including BAM coverage internals
Not recommended for the first pass.

Pros:
- highest theoretical ceiling

Cons:
- significantly higher implementation and debugging cost
- harder `pysam` integration boundary
- more fragile for current project size

## High-level architecture

Keep the current parser classes as the stable API layer:
- `src/drvizer/gtf_parser.py`
- `src/drvizer/bed_parser.py`
- `src/drvizer/bam_parser.py`

Add internal optimization modules:
- `src/drvizer/_cython_gtf.pyx`
- `src/drvizer/_cython_projection.pyx`
- `src/drvizer/_cython_bed.pyx` (optional in the first implementation wave)
- `src/drvizer/_parallel.py`

### Responsibilities

#### Public parser modules
Responsible for:
- argument validation
- selecting optimized vs pure Python paths
- preserving output shapes and behavior
- organizing cached parser state
- handling fallback behavior

#### `_cython_gtf.pyx`
Responsible for:
- hot-path GTF line parsing
- attribute extraction from the GTF attribute column
- fast exon/CDS row processing support
- producing lightweight parsed records for Python-side assembly

#### `_cython_projection.pyx`
Responsible for:
- transcript-to-genomic segmented projection
- exon-offset traversal for transcript intervals
- shared use by transcript-coordinate BED and transcript-coordinate BAM

#### `_cython_bed.pyx`
Responsible for:
- hot-path BED line parsing
- numeric field conversion
- optional fast transcript-coordinate BED segment expansion

#### `_parallel.py`
Responsible for:
- orchestrating multi-BAM parallel execution
- worker lifecycle and result aggregation
- serial fallback when parallel execution fails

## Detailed design by parser

### GTFParser

#### Objective
Reduce time spent parsing large GTF files while keeping the existing cached data model unchanged.

#### Design
- Preserve `gene_transcripts`, `gene_info`, `gene_name_to_id`, and `transcript_to_gene` as the canonical Python data structures.
- Keep `GTFParser` as the orchestration layer.
- Move the hottest line-processing and attribute-parsing work into `_cython_gtf.pyx`.
- Keep multi-file ordering and duplicate-transcript suppression in Python so behavior remains explicit and easy to verify.

#### Notes
- The current parser already has chunk-oriented structure, which is a good integration point for the accelerated path.
- Existing transcript indexes and projection caches should remain owned by `gtf_parser.py`.

### BEDParser

#### Objective
Speed up both plain genomic BED parsing and transcript-coordinate BED expansion.

#### Design
- Keep `BEDParser` as the public entry point and record organizer.
- Use `_cython_bed.pyx` for hot-path BED row parsing and numeric conversion.
- Use `_cython_projection.pyx` for transcript-coordinate segment expansion.
- Preserve `anno_data` layout and current region-query behavior.

#### Notes
- Transcript-coordinate BED must continue to expand into one genomic segment per overlapping exon section.
- No coordinate-collapsing behavior may be reintroduced.

### BAMParser

#### Objective
Improve coverage computation throughput without replacing `pysam`.

#### Design
- Keep `pysam.AlignmentFile.fetch()` for BAM reading.
- Add a worker model in `_parallel.py` that parallelizes across BAM files, not across transcripts.
- Each worker computes a per-file coverage array.
- The main process aggregates per-file arrays using the existing `sum` or `mean` semantics.
- Use `_cython_projection.pyx` for transcript-coordinate block projection where available.

#### Why parallelize by BAM file
- better task granularity than transcript-level splitting
- avoids excessive scheduling overhead
- avoids repeated BAM/header setup for fine-grained transcript tasks
- maps naturally onto the existing aggregation model

#### Notes
- Single-BAM execution can remain serial by default.
- Parallel execution should activate only when it can help, such as multiple BAM inputs.

## Fallback strategy

### Compiled extensions unavailable
- Import optimized modules opportunistically.
- If import fails, fall back automatically to the pure Python implementation.
- The package must remain functional without requiring users to change code.

### Parallel execution failure
- If worker startup or aggregation fails, retry using serial execution.
- Performance degradation is acceptable; correctness regression is not.

## Packaging and build changes

Update packaging so internal extensions can be built without changing public usage.

Required changes:
- extend `setup.py` to build extension modules with `cythonize`
- likely add `pyproject.toml` to declare build-system requirements
- ensure installation still works in pure Python mode when compiled artifacts are absent

## Testing and validation

### Test environment
Use:
- `/datf/hanxi/software/miniconda3/envs/DRS/bin/python`

It is acceptable to install `pytest` in that environment.

### Functional test coverage
Add tests covering:
- single-file GTF parsing matches current behavior
- multi-file GTF transcript de-duplication matches current behavior
- genomic BED parsing matches current behavior
- transcript-coordinate BED segment expansion matches current behavior
- genomic BAM coverage matches current behavior
- transcript-coordinate BAM coverage matches current behavior
- `sum` and `mean` aggregation remain unchanged
- serial and parallel multi-BAM execution return identical results

### Regression checks
Must explicitly verify:
- transcript-coordinate BAM does not fill introns
- transcript-coordinate BED does not collapse exon-separated intervals into continuous genomic spans

### Performance checks
Add lightweight benchmark scripts or opt-in tests for:
- first-pass GTF parse time
- transcript-coordinate BED parse time
- multi-BAM coverage time

## Implementation order

Recommended order:
1. add `_cython_projection.pyx` and switch transcript-coordinate BED/BAM to use it
2. add `_cython_gtf.pyx` for GTF parsing hot loops
3. add `_parallel.py` for multi-BAM parallel aggregation
4. add `_cython_bed.pyx` for additional BED parsing acceleration if profiling still shows it matters

## Success criteria

The work is successful if all of the following are true:
- public API usage is unchanged
- pure Python fallback still works
- transcript-coordinate correctness is preserved
- multi-BAM results are identical between serial and parallel execution
- large GTF/BED parse paths are measurably faster in the DRS environment
- there is a testable path using `pytest` in `/datf/hanxi/software/miniconda3/envs/DRS/bin/python`

## Scope guardrails

This design does not include:
- rewriting plotting code
- replacing `pysam`
- changing user-facing parser semantics
- broad refactors unrelated to parsing and coverage performance
