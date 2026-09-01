import warnings
from pathlib import Path

from setuptools import Extension, find_packages, setup

ROOT = Path(__file__).parent
long_description = (ROOT / "README.md").read_text(encoding="utf-8")
requirements = [
    line.strip()
    for line in (ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines()
    if line.strip() and not line.startswith("#")
]

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
    name="drvizer",
    version="1.0.0",
    author="x1han",
    author_email="han_xi@gzlab.ac.cn",
    description="A Python tool for parsing GTF/BED files and visualizing gene transcript structures",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/x1han/drVizer",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    ext_modules=ext_modules,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    python_requires=">=3.9",
    install_requires=requirements,
    extras_require={
        "bam": ["pysam"],
    },
    keywords="bioinformatics, genomics, transcriptomics, gtf, bed, visualization",
    project_urls={
        "Bug Reports": "https://github.com/x1han/drVizer/issues",
        "Source": "https://github.com/x1han/drVizer",
        "Documentation": "https://github.com/x1han/drVizer#readme",
    },
)
