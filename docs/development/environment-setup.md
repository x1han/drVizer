# Environment setup

Day-to-day environment notes for drVizer contributors. The maintainer-oriented
summary lives in `repo/CLAUDE.md`; this document expands the gotchas and
recovery recipes that are too verbose to keep inline.

## Editable install gotcha

Switching between git branches in this repository -- especially after manual
ref manipulation such as `git update-ref refs/heads/main <sha>` or
`git symbolic-ref HEAD refs/heads/main` following a worktree cleanup -- can
leave the editable install pointing at the old `src/` tree. Symptoms range
from a stale `import` (the editable `.pth` file still references a worktree
that no longer exists) to `ImportError` for modules that have been moved,
renamed, or removed on the current branch.

### Diagnose

Check which `src/` the editable link resolves to:

```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/python -c "import drvizer; print(drvizer.__file__)"
```

If the printed path does not match the current branch's `src/drvizer/__init__.py`,
refresh the editable install:

```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/pip install --force-reinstall --no-deps -e .
```

Run this command after **every** merge or branch switch that involved manual
ref manipulation. The build/bytecode caches survive the branch switch but
the `.pth` entry that pip writes for editable installs points at a specific
filesystem location, so it must be regenerated when that location changes.

If files appear to be missing from the working tree after a manual ref update
(rather than just stale imports), the diagnosis is different --
`git reset --hard HEAD` restores the worktree contents from HEAD without
touching the editable install.

## CI vs. local install

The CI workflow (`ci.yml`) builds and installs a real wheel:

```bash
pip install build
python -m build --wheel
pip install dist/*.whl[test]
```

Local development typically uses the editable install instead, because the
edit-refresh cycle is faster than a wheel rebuild. Both paths exercise the
same `pyproject.toml` build backend, so behavior should agree modulo the
`.pth` resolution described above.

## Optional dependencies

`pysam` (BAM support) is not listed in `requirements.txt` -- install it
explicitly when working on BAM-related code:

```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/pip install pysam
```

Or via the optional extra defined in `pyproject.toml`:

```bash
/datf/hanxi/software/miniconda3/envs/DRS/bin/pip install -e .[bam]
```
