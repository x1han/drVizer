from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="drvizer",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A Python tool for parsing GTF/BED files and visualizing gene transcript structures",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/drvizer",
    packages=["drvizer"],
    package_dir={"drvizer": "src"},
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