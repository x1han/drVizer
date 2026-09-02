"""Phase 9 -- version consistency guard (no hardcoded version string).

This module complements ``tests/test_phase_7_version_and_metadata.py``
by re-asserting the ``drvizer.__version__`` / metadata / pyproject
version triple without baking any literal version string into the test
file. The Phase 7 test guards the specific 0.1.1 release; this one
keeps working as the project bumps to 0.1.2, 0.2.0, etc., catching the
common drift where ``__version__`` is bumped but ``pyproject.toml`` is
forgotten (or vice versa).

Tests:

* ``test_runtime_metadata_pyproject_versions_agree``
    All three introspection paths return the same string. Uses
    ``importlib.metadata.version`` (installed wheel metadata),
    ``drvizer.__version__`` (runtime attribute), and ``pyproject.toml``
    ``[project] version`` (source-of-truth for the next release).
    No literal version string appears in the test body.

* ``test_pyproject_version_is_valid_pep440``
    The pyproject version parses via ``packaging.version.Version`` so
    a typo like ``"0.1..1"`` or ``"v0.1.1"`` fails this test before it
    reaches PyPI.
"""

from __future__ import annotations

import importlib.metadata
from pathlib import Path

try:
    import tomllib
except ImportError:  # Python < 3.11
    import tomli as tomllib

from packaging.version import Version, InvalidVersion


def _read_pyproject_version() -> str:
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    with pyproject_path.open("rb") as f:
        data = tomllib.load(f)
    return data["project"]["version"]


def test_runtime_metadata_pyproject_versions_agree():
    """``drvizer.__version__``, ``importlib.metadata.version('drvizer')``,
    and ``pyproject.toml`` must all return the same string. No literal
    version is hardcoded so the test stays useful across releases."""
    import drvizer

    runtime_version = drvizer.__version__
    metadata_version = importlib.metadata.version("drvizer")
    pyproject_version = _read_pyproject_version()

    assert runtime_version == metadata_version, (
        f"drvizer.__version__ ({runtime_version!r}) disagrees with "
        f"importlib.metadata.version ({metadata_version!r})"
    )
    assert runtime_version == pyproject_version, (
        f"drvizer.__version__ ({runtime_version!r}) disagrees with "
        f"pyproject.toml [project] version ({pyproject_version!r})"
    )


def test_pyproject_version_is_valid_pep440():
    """The pyproject version must parse as a PEP 440 Version.

    Catches typos (``0.1..1``) and stray prefixes (``v0.1.1``) that
    would otherwise sneak into a release. Uses ``packaging.version``
    rather than a hand-rolled regex so the rule set stays in sync
    with PyPI's actual parser."""
    pyproject_version = _read_pyproject_version()
    try:
        Version(pyproject_version)
    except InvalidVersion as exc:
        raise AssertionError(
            f"pyproject.toml [project] version {pyproject_version!r} "
            f"is not a valid PEP 440 version: {exc}"
        ) from exc
