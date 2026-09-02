# Roadmap

Living record of drVizer development phases. Each phase corresponds to a
discrete unit of work (audit remediation, infrastructure, packaging,
release) and links to the relevant CHANGELOG entries and release tags.

## Phase status

| Phase | Theme | Status |
|-------|-------|--------|
| 2.2   | LRU cache + lazy ProcessPool lifecycle | shipped |
| 3     | BAM kernel + chunk degradation + CRLF tolerance | shipped |
| 7     | Imap-unordered contract + first public release metadata | shipped |
| 8     | Tuple-return attribution + editable CI migration prep | shipped |
| 9     | Post-release improvements (CI wheel install, TestPyPI smoke test, version guard, dev docs) | in progress |

## Phase 9 detail

Concrete deliverables for the current phase:

1. **CI uses wheel install instead of editable install.** `.github/workflows/ci.yml`
   now builds a real wheel via `python -m build --wheel` and installs
   `dist/*.whl[test]`. This catches packaging regressions that the editable
   path silently tolerates.
2. **TestPyPI smoke test step.** `.github/workflows/release.yml` adds a
   post-upload step that pulls the just-published artifact back via pip and
   asserts the public API imports cleanly. Placed *after* the `twine upload`
   step because the smoke test depends on the upload having completed.
3. **Version consistency guard.** `tests/test_version.py` re-asserts the
   `__version__` / `importlib.metadata.version` / `pyproject.toml` triple
   without baking any literal version string into the test body.
4. **Development docs.** `docs/development/environment-setup.md` collects the
   editable-install gotcha and CI-vs-local workflow notes that previously
   lived inline in `CLAUDE.md`.
5. **Benchmark data download helper.** `benchmarks/download_benchmark_data.sh`
   fetches the (small) inputs the local benchmarks expect, with a graceful
   skip path when `samtools` is not on `PATH`.

## Cross-phase invariants

These have held across every shipped phase and are expected to keep
holding:

- `drvizer.__version__` and `importlib.metadata.version("drvizer")` agree.
- `pyproject.toml [project] version` matches the runtime attribute.
- The tuple-return contract from Phase 8 (`(bam_path, coverage_array)`)
  keeps aggregation order-independent -- see the `Note` docstring at the
  top of `src/drvizer/_parallel.py`.
