import warnings
from pathlib import Path

from setuptools import Extension, find_packages, setup

ROOT = Path(__file__).parent

# P1-16 (Phase 1) + Phase 7: setup.py is a thin C-extension compile-only
# wrapper. PEP 621 metadata lives in pyproject.toml [project]; setup.py only
# carries what the build backend needs to compile the Cython extensions and
# locate the source layout (``src/``). All version/classifier/dependency
# truth is in pyproject.toml.
extensions = [
    Extension(
        "drvizer._cython_projection",
        ["src/drvizer/_cython_projection.pyx"],
        extra_compile_args=["-std=c99"],
    ),
    Extension(
        "drvizer._cython_gtf",
        ["src/drvizer/_cython_gtf.pyx"],
        extra_compile_args=["-std=c99"],
    ),
    Extension(
        "drvizer._cython_bed",
        ["src/drvizer/_cython_bed.pyx"],
        extra_compile_args=["-std=c99"],
    ),
]

try:
    from Cython.Build import cythonize
except ImportError:
    warnings.warn(
        "Cython is not installed; building drvizer without compiled "
        "extensions. Parsers will use Python fallbacks only.",
        ImportWarning,
        stacklevel=2,
    )
    ext_modules = []
else:
    ext_modules = cythonize(extensions, language_level="3")

setup(
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    ext_modules=ext_modules,
)