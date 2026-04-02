# drVizer 🧬

**drVizer** is a visualization tool for **direct RNA-Seq** related analyses. It is designed to help inspect and present transcript-centric genomic information in a compact figure, including:

- **isoform structures**
- **read-count / coverage-style tracks**
- **transposable element (TE) annotations**
- **RNA modification related annotations**

The name **drVizer** comes from **direct RNA visualization**: a lightweight tool for turning transcript models and associated genomic annotations into publication-friendly plots.

## ✨ What drVizer is for

drVizer focuses on transcript structure visualization from **GTF annotations**, with optional overlay tracks from **BED** and **BAM**-derived data.

Typical use cases include:

- visualizing **isoform architecture** for one gene
- comparing transcript structures across genes or annotations
- adding **TE / repeat** annotations to transcript plots
- showing **coverage / read-count style** signal from BAM files
- displaying extra interval-based annotations such as **RNA modification** sites

## 📦 Installation

### Option 1: Clone from GitHub and install locally

```bash
git clone https://github.com/x1han/drVizer.git
cd drVizer
pip install -e .
```

If you only want a regular local install instead of editable mode:

```bash
pip install .
```

### Dependencies

Core dependencies are listed in `requirements.txt` and include:

- `pandas`
- `matplotlib`
- `numpy`

If you want to use **BAM coverage tracks**, you may also need:

```bash
pip install pysam
```

## 🚀 Two ways to use drVizer

drVizer supports both:

1. **CLI usage** for straightforward command-line plotting
2. **Python API usage** for notebooks and scripted workflows

---

## 1. CLI usage 🖥️

After installation, the CLI entry point is:

```bash
drvizer
```

### Example: visualize one gene from a GTF file

```bash
drvizer --gtf genes.gtf --gene ENSG00000136997 --output gene_plot.png
```

### Example: visualize BED annotations in a genomic region

```bash
drvizer --bed repeats.bed --region chr1:1000000-2000000 --output te_plot.png
```

### Example: combine GTF and BED data in one figure

```bash
drvizer --gtf genes.gtf --bed repeats.bed --gene TP53 --output merged_plot.png
```

### Repository-local wrapper

If you are working directly inside a cloned repository checkout, you can also run:

```bash
python drvizer_cli.py --gtf genes.gtf --gene TP53 --output tp53_plot.png
```

---

## 2. Python API usage 📓

### High-level API with `DrViz`

This is the most convenient notebook-friendly interface.

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

### Simple one-shot plotting

```python
from drvizer import DrViz

fig = (
    DrViz()
    .load_gtf("genes.gtf")
    .plot("TP53")
)
```

### Low-level API

If you want more control, you can use the parser and visualization functions directly.

```python
from drvizer import GTFParser, visualize_gene_transcripts

parser = GTFParser("genes.gtf")
parser.parse_gtf()
transcript_data = parser.get_transcript_data("TP53")

fig = visualize_gene_transcripts(transcript_data)
```

### Add BED annotations

```python
from drvizer import GTFParser, BEDParser, merge_parsers, visualize_gene_transcripts

gtf_parser = GTFParser("genes.gtf")
gtf_parser.parse_gtf()

bed_parser = BEDParser("repeats.bed")
bed_parser.parse_bed()

merged = merge_parsers(gtf_parser, bed_parser)
data = merged.get_transcript_data("TP53")
fig = visualize_gene_transcripts(data)
```

### Jupyter / notebook usage

```python
%matplotlib inline

from drvizer import DrViz

parser = DrViz().load_gtf("genes.gtf").build()
fig = parser.plot("ENSG00000136997")
```

## 🧩 Main capabilities

- Parse one or multiple **GTF** files
- Access genes by **gene ID**, **gene name**, or **transcript ID**
- Add **BED** annotation tracks
- Add **BAM coverage** tracks
- Export figures as **PNG / PDF / SVG**
- Reuse parsed data efficiently across many genes
- Use either a **high-level chainable API** or **lower-level parser classes**

## 📚 Main public API

Commonly used imports include:

```python
from drvizer import (
    DrViz,
    GTFParser,
    BEDParser,
    visualize_gene_transcripts,
    save_visualization,
    merge_parsers,
)
```

## 🗂️ Project structure

- `src/drvizer/gtf_parser.py` — GTF parsing and transcript extraction
- `src/drvizer/bed_parser.py` — BED parsing and annotation grouping
- `src/drvizer/bam_parser.py` — BAM-based coverage extraction
- `src/drvizer/visualizer.py` — figure generation and track layout
- `src/drvizer/api.py` — high-level chainable API
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
