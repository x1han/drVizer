# drVizer - Gene Transcript Visualization Tool
==============================================

This Python module provides functionality for parsing GTF files and visualizing gene transcript structures. It is designed to work as a standalone tool but can also be integrated into other applications.

## Features
----------

- Parse GTF files to extract transcript information
- Parse BED files to extract transposable element (TE) information
- Visualize gene transcript structures with improved formatting
- Export data in various formats (JSON, CSV)
- Filter transcripts based on various criteria
- Compare transcript structures from different sources
- Efficiently handle multiple genes with cached parsing
- Unified method for accessing genes by ID or name
- Transcript sorting by exon order
- Enhanced visualization with strand direction indicators
- Support for parsing single or multiple GTF files
- Merge GTF and BED data for comprehensive visualization

## Installation
--------------

Install the package:

```bash
pip install .
```

For development:

```bash
pip install -e .
pip install -r requirements.txt
```

## Usage
-----

### Command Line Interface
---------------------

Installed CLI usage:

```bash
# Parse and visualize a specific gene from GTF file
drvizer --gtf genes.gtf --gene ENSG00000136997 --output gene_plot.png

# Parse BED file for a genomic region
drvizer --bed repeats.bed --region chr1:1000000-2000000 --output te_plot.png

# Merge GTF and BED data for comprehensive visualization
drvizer --gtf genes.gtf --bed repeats.bed --gene TP53 --output merged_plot.png
```

Repository checkout usage:

```bash
python drvizer_cli.py --gtf genes.gtf --gene ENSG00000136997 --output gene_plot.png
```

### Python API
----------

#### Basic Usage
```python
from drvizer import parse_gtf_for_gene
from drvizer import visualize_gene_transcripts, save_visualization

transcript_data = parse_gtf_for_gene("path/to/file.gtf", "ENSG00000136997")
fig = visualize_gene_transcripts(transcript_data)
save_visualization(fig, "output.png")
```

#### Parse Multiple GTF Files
```python
from drvizer import parse_gtf_for_gene, parse_gtf_all_genes

gtf_files = ["path/to/file1.gtf", "path/to/file2.gtf", "path/to/file3.gtf"]
transcript_data = parse_gtf_for_gene(gtf_files, "ENSG00000136997")
all_genes_data = parse_gtf_all_genes(gtf_files)
```

#### Efficient Multiple Genes Usage
```python
from drvizer import GTFParser, visualize_gene_transcripts, save_visualization

parser = GTFParser(["path/to/file1.gtf", "path/to/file2.gtf"])
parser.parse_gtf()

for gene_id in ["ENSG00000136997", "ENSG00000236998", "ENSG00000188015.10"]:
    transcript_data = parser.get_transcript_data(gene_id)
    fig = visualize_gene_transcripts(transcript_data)
    save_visualization(fig, f"{gene_id}_transcripts.png")
```

#### Access Genes by Name or ID
```python
from drvizer import GTFParser

parser = GTFParser("path/to/file.gtf")
parser.parse_gtf()

transcript_data_by_id = parser.get_transcript_data("ENSG00000136997")
transcript_data_by_name = parser.get_transcript_data("TP53")
```

#### Parse BED Files for Transposable Elements
```python
from drvizer import BEDParser

bed_parser = BEDParser("path/to/file.bed")
bed_data = bed_parser.parse_bed()
anno_in_region = bed_parser.get_anno_in_region("chr1", 1000000, 2000000)
```

#### Merge GTF and BED Data for Comprehensive Visualization
```python
from drvizer import GTFParser, BEDParser, merge_parsers, visualize_gene_transcripts

gtf_parser = GTFParser("path/to/file.gtf")
bed_parser = BEDParser("path/to/file.bed")

gtf_parser.parse_gtf()
bed_parser.parse_bed()

merged_parser = merge_parsers(gtf_parser, bed_parser)
gene_data = merged_parser.get_transcript_data("ENSG00000136997")
fig = visualize_gene_transcripts(gene_data)
```

#### High-level API
```python
from drvizer import DrViz

parser = (
    DrViz()
    .load_gtf("path/to/file.gtf")
    .add_bed_track("path/to/file.bed", label="TE")
    .build()
)

fig = parser.plot("TP53")
```

#### Advanced Usage
```python
from drvizer import GTFParser, filter_transcripts, get_transcript_stats

parser = GTFParser("path/to/file.gtf")
parser.parse_gtf()

transcript_data = parser.get_transcript_data("ENSG00000136997")
filtered_data = filter_transcripts(transcript_data, min_exons=3)
stats = get_transcript_stats(transcript_data)
print(stats)
```

## Module Structure
----------------

- `src/gtf_parser.py`: Main parser for GTF files with efficient caching for multiple genes and unified gene access
- `src/bed_parser.py`: Parser for BED files to extract transposable element information
- `src/visualizer.py`: Functions for visualizing transcript structures with improved formatting and sorting
- `src/api.py`: High-level chainable Python API for notebook and interactive use
- `drvizer_cli.py`: Repository-local CLI wrapper

## Requirements
------------

- Python 3.7+
- pandas
- matplotlib
- numpy

## Integration with Jupyter
------------------------

```python
%matplotlib inline

from drvizer import DrViz

parser = DrViz().load_gtf("path/to/file.gtf").build()
fig = parser.plot("ENSG00000136997")
```

## Improvements
------------

The drVizer tool includes several improvements:

1. **Unified Method**: The `get_transcript_data` method now works with both gene IDs and gene names, eliminating the need for separate methods.
2. **Improved Visualization**: Better spacing prevents y-stick stacking and makes the visualizations clearer.
3. **Transcript Sorting**: Transcripts are now sorted by exon order for better organization.
4. **Enhanced Titles**: Titles now use 5' → 3' and 3' ← 5' indicators to clearly show strand direction.
5. **Better Formatting**: Overall improved formatting and layout for clearer visualizations.
6. **Multiple GTF Support**: Support for parsing single or multiple GTF files with unified API.
7. **BED File Support**: Added support for parsing BED files containing transposable element data.
8. **Merged Visualization**: Ability to merge GTF and BED data for comprehensive visualization of genes and overlapping TEs.
