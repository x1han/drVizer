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

1. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage
-----

### Command Line Interface
---------------------

The main entry point is `drvizer.py` in the root directory:

```bash
# Parse and visualize a specific gene from GTF file
python drvizer.py --gtf genes.gtf --gene ENSG00000136997 --output gene_plot.png

# Parse BED file for a genomic region
python drvizer.py --bed repeats.bed --region chr1:1000000-2000000 --output te_plot.png

# Merge GTF and BED data for comprehensive visualization
python drvizer.py --gtf genes.gtf --bed repeats.bed --gene TP53 --output merged_plot.png
```

### Python API
----------

#### Basic Usage
```python
from src.gtf_parser import parse_gtf_for_gene
from src.visualizer import visualize_gene_transcripts, save_visualization

# Parse GTF file for a specific gene
transcript_data = parse_gtf_for_gene("path/to/file.gtf", "ENSG00000136997")

# Visualize transcripts with improved formatting
fig = visualize_gene_transcripts(transcript_data)

# Save visualization
save_visualization(fig, "output.png")
```

#### Parse Multiple GTF Files
```python
from src.gtf_parser import parse_gtf_for_gene, parse_gtf_all_genes

# Parse multiple GTF files for a specific gene
gtf_files = ["path/to/file1.gtf", "path/to/file2.gtf", "path/to/file3.gtf"]
transcript_data = parse_gtf_for_gene(gtf_files, "ENSG00000136997")

# Parse all genes from multiple GTF files
all_genes_data = parse_gtf_all_genes(gtf_files)
```

#### Efficient Multiple Genes Usage
```python
from src.gtf_parser import GTFParser
from src.visualizer import visualize_gene_transcripts, save_visualization

# Initialize parser with single or multiple GTF files
parser = GTFParser(["path/to/file1.gtf", "path/to/file2.gtf"])

# Parse all genes at once (efficient)
parser.parse_gtf()  # This will parse all genes and cache the data

# Now we can access any gene efficiently using the unified method
gene_ids = ["ENSG00000136997", "ENSG00000236998", "ENSG00000188015.10"]
for gene_id in gene_ids:
    # This is now fast because the data is cached
    transcript_data = parser.get_transcript_data(gene_id)
    
    # Visualize transcripts
    fig = visualize_gene_transcripts(transcript_data)
    
    # Save visualization
    save_visualization(fig, f"{gene_id}_transcripts.png")
```

#### Access Genes by Name or ID
```python
from src.gtf_parser import GTFParser

# Initialize parser
parser = GTFParser("path/to/file.gtf")

# Parse all genes to build the name mapping
parser.parse_gtf()

# Access gene by ID
transcript_data_by_id = parser.get_transcript_data("ENSG00000136997")

# Access gene by name
transcript_data_by_name = parser.get_transcript_data("TP53")
```

#### Parse BED Files for Transposable Elements
```python
from src.bed_parser import BEDParser

# Initialize BED parser
bed_parser = BEDParser("path/to/file.bed")

# Parse all TEs
bed_data = bed_parser.parse_bed()

# Get TEs in a specific genomic region
te_in_region = bed_parser.get_te_in_region("chr1", 1000000, 2000000)
```

#### Merge GTF and BED Data for Comprehensive Visualization
```python
from src.gtf_parser import GTFParser
from src.visualizer import merge_parsers
from src.bed_parser import BEDParser
from src.visualizer import visualize_gene_transcripts

# Initialize parsers
gtf_parser = GTFParser("path/to/file.gtf")
bed_parser = BEDParser("path/to/file.bed")

# Parse the files
gtf_parser.parse_gtf()
bed_parser.parse_bed()

# Merge parsers
merged_parser = merge_parsers(gtf_parser, bed_parser)

# Visualize gene transcripts with overlapping TEs
gene_data = merged_parser.get_transcript_data("ENSG00000136997")
fig = visualize_gene_transcripts(gene_data)
```

#### Advanced Usage
```python
from src.gtf_parser import GTFParser
from src.utils import filter_transcripts, get_transcript_stats

# Initialize parser
parser = GTFParser("path/to/file.gtf")
parser.parse_gtf("ENSG00000136997")

# Get transcript data
transcript_data = parser.get_transcript_data("ENSG00000136997")

# Filter transcripts
filtered_data = filter_transcripts(transcript_data, min_exons=3)

# Get statistics
stats = get_transcript_stats(transcript_data)
print(stats)
```

## Module Structure
----------------

- `src/gtf_parser.py`: Main parser for GTF files with efficient caching for multiple genes and unified gene access
- `src/bed_parser.py`: Parser for BED files to extract transposable element information
- `src/visualizer.py`: Functions for visualizing transcript structures with improved formatting and sorting
- `src/utils.py`: Utility functions for data processing
- `drvizer.py`: Main command-line interface script

## Requirements
------------

- Python 3.7+
- pandas
- matplotlib
- numpy

## Integration with Jupyter
------------------------

This module can be easily used in Jupyter notebooks for interactive analysis:

```python
# In a Jupyter notebook
%matplotlib inline

from src.gtf_parser import parse_gtf_for_gene
from src.visualizer import visualize_gene_transcripts

# Parse and visualize
transcript_data = parse_gtf_for_gene("path/to/file.gtf", "ENSG00000136997")
fig = visualize_gene_transcripts(transcript_data)
fig.show()
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