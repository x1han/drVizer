"""Phase 7 -- version string and pyproject metadata.

The first public release is v0.1.0. Two assertions:

* ``drvizer.__version__`` equals the string ``"0.1.2"`` AND
  ``importlib.metadata.version("drvizer")`` returns the same string.
* The ``pyproject.toml [project]`` classifiers include Python 3.8
  through 3.12 and ``Development Status :: 4 - Beta``.

Tests:

* ``test_version_string_matches_0_1_2``
    Both introspection paths return ``"0.1.2"``.

* ``test_pyproject_classifiers_include_python_38_through_312``
    The PEP 621 classifiers list contains every Python 3.8-3.12
    classifier and the Development Status :: 4 - Beta classifier.

* ``test_drviz_reexport_round_trip``
    ``from drvizer import DrViz`` works (so a broken ``__init__.py``
    re-export is caught immediately).
"""

import importlib.metadata
try:
    import tomllib
except ImportError:  # Python < 3.11
    import tomli as tomllib
from pathlib import Path

import pytest


def test_version_string_matches_0_1_2():
    """Both ``drvizer.__version__`` and ``importlib.metadata.version``
    must return the v0.1.2 string."""
    import drvizer

    assert drvizer.__version__ == "0.1.2", (
        f"drvizer.__version__ must equal '0.1.2'; got {drvizer.__version__!r}"
    )
    metadata_version = importlib.metadata.version("drvizer")
    assert metadata_version == "0.1.2", (
        f"importlib.metadata.version('drvizer') must equal '0.1.2'; got {metadata_version!r}"
    )
    assert drvizer.__version__ == metadata_version, (
        f"drvizer.__version__ and importlib.metadata.version must agree; "
        f"got {drvizer.__version__!r} vs {metadata_version!r}"
    )


def test_pyproject_classifiers_include_python_38_through_312():
    """The PEP 621 classifiers list must include Python 3.8 through 3.12
    and the Development Status :: 4 - Beta classifier."""
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    with pyproject_path.open("rb") as f:
        data = tomllib.load(f)

    classifiers = data["project"]["classifiers"]
    assert isinstance(classifiers, list)

    expected_python_classifiers = [
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ]
    for expected in expected_python_classifiers:
        assert expected in classifiers, (
            f"pyproject classifiers must include {expected!r}; got {classifiers}"
        )

    assert "Development Status :: 4 - Beta" in classifiers, (
        f"pyproject classifiers must include 'Development Status :: 4 - Beta'; got {classifiers}"
    )


def test_drviz_reexport_round_trip():
    """``from drvizer import DrViz`` must succeed so a broken
    ``__init__.py`` re-export is caught immediately after the
    version bump."""
    import drvizer

    DrViz = getattr(drvizer, "DrViz", None)
    assert DrViz is not None, "drvizer must re-export DrViz"
    ParallelCoverageError = getattr(drvizer, "ParallelCoverageError", None)
    assert ParallelCoverageError is not None, "drvizer must re-export ParallelCoverageError"
    assert issubclass(ParallelCoverageError, RuntimeError)