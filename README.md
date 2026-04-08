# drVizer 🧬

**drVizer** is a visualization tool for **direct RNA-Seq** related analyses. It helps turn transcript-centric genomic information into compact, publication-friendly figures, including:

- **isoform structures**
- **read-count / coverage-style tracks**
- **transposable element (TE) annotations**
- **RNA modification related annotations**

The name **drVizer** comes from **direct RNA visualization**.

## ✨ What drVizer is for

drVizer focuses on transcript structure visualization from **GTF annotations**, with optional overlay tracks from **BED** and **BAM**-derived data.

Typical use cases include:

- visualizing **isoform architecture** for one gene
- comparing transcript structures across genes or annotations
- adding **TE / repeat** annotations to transcript plots
- showing **coverage / read-count style** signal from BAM files
- displaying extra interval-based annotations such as **RNA modification** sites

## 📦 Installation

### Clone from GitHub and install locally

```bash
git clone https://github.com/x1han/drVizer.git
cd drVizer
pip install -e .
```

If you want a regular local install instead of editable mode:

```bash
pip install .
```

### Optional dependency for BAM tracks

Core dependencies are listed in `requirements.txt`.

If you want to use **BAM coverage tracks**, you may also need:

```bash
pip install pysam
```

## 🚀 Two ways to use drVizer

`drVizer` supports both:

1. **CLI usage** for straightforward command-line plotting
2. **Python API usage** for notebooks and scripted workflows

---

## 1. CLI usage 🖥️

After installation, the CLI entry point is:

```bash
drvizer
```

The supported CLI workflow is gene-centric: provide a GTF file, a target gene, and optionally one or more BED tracks.

### Visualize one gene from a GTF file

```bash
drvizer --gtf genes.gtf --gene ENSG00000136997 --output gene_plot.png
```

### Combine GTF and BED data in one figure

```bash
drvizer --gtf genes.gtf --bed repeats.bed --gene TP53 --output merged_plot.png
```

### Run the repository-local wrapper

If you are working directly inside a cloned repository checkout, you can also run:

```bash
python drvizer_cli.py --gtf genes.gtf --gene TP53 --output tp53_plot.png
```

---

## 2. Python API usage 📓

The high-level workflow is:

- `load_gtf(...)`
- `add_bed_track(...)`
- `add_bam_track(...)`
- `build()`
- `plot(...)`

### Basic workflow

```python
from drvizer import DrViz

parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed", label="TE")
    .build()
)

fig = parser.plot("TP53")
```

### Add multiple BED tracks

```python
from drvizer import DrViz

parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("te.bed", label="TE", color="tomato")
    .add_bed_track("rna_mod_sites.bed", label="m6A", color="royalblue")
    .build()
)

fig = parser.plot("TP53")
```

### Add BAM coverage

```python
from drvizer import DrViz

parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bam_track("reads.bam", label="Coverage", color="steelblue")
    .build()
)

fig = parser.plot("TP53")
```

### One-shot plotting

```python
from drvizer import DrViz

fig = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed", label="TE")
    .plot("TP53")
)
```

### Jupyter / notebook usage

```python
%matplotlib inline

from drvizer import DrViz

parser = (
    DrViz()
    .load_gtf("genes.gtf")
    .add_bed_track("repeats.bed", label="TE")
    .build()
)

fig = parser.plot("ENSG00000136997")
```

## 🧩 Main capabilities

- Parse one or multiple **GTF** files
- Access genes by **gene ID**, **gene name**, or **transcript ID**
- Add **BED** annotation tracks with `add_bed_track(...)`
- Add **BAM coverage** tracks with `add_bam_track(...)`
- Export figures as **PNG / PDF / SVG**
- Reuse parsed data efficiently across many genes with `build()`
- Use a single **high-level chainable API** centered on `DrViz`

## 📚 Main public API

```python
from drvizer import DrViz
```

## 🗂️ Project structure

- `src/drvizer/api.py` — high-level chainable API
- `src/drvizer/gtf_parser.py` — GTF parsing and transcript extraction
- `src/drvizer/bed_parser.py` — BED parsing and annotation grouping
- `src/drvizer/bam_parser.py` — BAM-based coverage extraction
- `src/drvizer/visualizer.py` — figure generation and track layout
- `src/drvizer/cli.py` — installed CLI entry point
- `drvizer_cli.py` — repository-local CLI wrapper

## ⚠️ Notes

- Python package import name is **`drvizer`**
- Repository name is **`drVizer`**
- Python imports are case-sensitive, so use:

```python
from drvizer import DrViz
```

not:

```python
from drVizer import DrViz
```
