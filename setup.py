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
]

try:
    from Cython.Build import cythonize
except ImportError:
    ext_modules = []
else:
    ext_modules = cythonize(extensions, language_level="3")

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
