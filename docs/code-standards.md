# drVizer Code Standards

## General style

- Keep changes surgical and tied to the requested behavior.
- Prefer clear, direct code over new abstractions.
- Match the existing Python style in nearby code.
- Keep public behavior centered on `DrViz`.
- Avoid adding comments unless they explain a non-obvious constraint.

## API design

- Preserve the chainable builder style in `DrViz`.
- Validate user-facing inputs at public API boundaries.
- Keep parser internals behind the public builder workflow unless exposing a new public API is intentional.
- Prefer explicit `ValueError` messages for unsupported user combinations.
- Do not add CLI entry points unless project direction changes.

## Parser behavior

### GTF

- Keep the output shape expected by `visualizer.py`.
- Preserve lookup by gene ID, gene name, and transcript ID.
- Preserve transcript-coordinate projection semantics.
- Maintain Python fallback behavior when Cython helpers are unavailable.

### BED

- Preserve BED3 support and default optional fields.
- Use half-open interval overlap semantics for region filtering.
- Keep Python and Cython parsing behavior aligned.
- Keep transcript-coordinate projection dependent on a loaded `GTFParser`.

### BAM

- Use int64 accumulators for coverage sums.
- Keep `sum` and `mean` aggregation semantics explicit.
- Keep multi-BAM worker counts capped.
- Preserve serial fallback when parallel aggregation fails.
- Treat `pysam` as required only for BAM support.

## Rendering behavior

- `visualize_gene_transcripts(...)` should return a matplotlib figure.
- Preserve track order from user registration.
- Preserve numeric y-axis grouping rules.
- Keep split transcript labels aligned with prepared track order.
- Do not mix plotting side effects into parser classes.

## Concurrency

- Prepare only process-safe genomic BED tracks in process pools.
- Keep BAM and transcript-coordinate tracks out of process preparation unless their state transfer is made explicit and tested.
- Cap worker counts when adding parallel work.
- Preserve deterministic output order even when preparation or aggregation is concurrent.

## Error handling

- Validate external/user inputs at boundaries.
- Avoid broad fallbacks that hide unsupported states.
- Wrap file read failures with useful context.
- Let unexpected internal errors fail loudly during tests.

## Documentation

- Public docs belong in `docs/` or `README.md` inside this repository.
- Local-only analysis, benchmark, and validation notes belong in the workspace `workspace/` directory.
- Keep examples runnable with the public `from drvizer import DrViz` import.
- Do not document unsupported CLI workflows as current behavior.

## Testing standards

- Add a failing test before production changes for bugs or behavior changes.
- Keep tests small and behavior-focused.
- Prefer public API tests for user-visible behavior.
- Add lower-level parser tests when parser semantics change.
- Run the full pytest suite before reporting completion.

## Repository operations

- Software changes belong in the inner source repository.
- Local workflow material belongs outside the source repository.
- Use the DRS Python environment for validation.
- Before committing, inspect changed scope and run repository-specific change detection where required by maintainer guidance.
