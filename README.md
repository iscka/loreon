# LOREON - Long-Read ONT Metagenomic Pipeline

A comprehensive bioinformatics pipeline for processing Oxford Nanopore Technologies (ONT) metagenomic sequencing data, specifically designed for fungal ITS barcode analysis.

## Features

- **Sequence Quality Filtering**: Length-based filtering with configurable min/max thresholds
- **High-Performance Mapping**: Utilizes minimap2 for fast and accurate sequence alignment
- **Parallel Processing**: Multi-threaded execution for filtering, mapping, and analysis
- **OTU Table Generation**: Automated OTU (Operational Taxonomic Unit) aggregation with DuckDB
- **Interactive Reports**: HTML reports with Plotly visualizations
- **Multiple Database Formats**: Support for UNITE, SILVA, CBS, Eukariome, and custom formats
- **GUI Interface**: User-friendly PyQt5 graphical interface
- **Taxonomy Mapping**: Intelligent caching system for taxonomic information

## System Requirements

### Software Dependencies

**External Tools** (must be installed separately):
- [minimap2](https://github.com/lh3/minimap2) >= 2.17
- [samtools](http://www.htslib.org/) >= 1.10
- awk (pre-installed on most Unix systems)

**Python** >= 3.8

### Python Dependencies

```bash
pip install -r requirements.txt
```

Key packages:
- biopython >= 1.79
- pandas >= 1.3.0
- duckdb >= 0.8.0
- PyQt5 >= 5.15.0 (for GUI)
- plotly >= 5.0.0 (for reports)
- jinja2 >= 3.0.0

## Installation

### Using Conda (Recommended)

```bash
# Create environment
conda create -n loreon python=3.9
conda activate loreon

# Install external dependencies
conda install -c bioconda minimap2 samtools

# Install Python dependencies
pip install -r requirements.txt
```

### Manual Installation

1. Install minimap2 and samtools following their respective documentation
2. Install Python dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### GUI Mode (Recommended)

```bash
python pipeline_gui.py
```

The GUI provides:
1. **Path Selection**: Input folder, reference database, output folder
2. **Parameters Configuration**:
   - Database format (UNITE, SILVA, CBS, etc.)
   - Length filter (min/max)
   - Threading options
   - minimap2 parameters (k-mer, window size)
3. **Real-time Logging**: Progress tracking with status bar
4. **Automatic Report**: Opens HTML report upon completion

### Command-Line Mode

```bash
python metaGenomics_new.py \
  -i /path/to/input_data \
  -o /path/to/output \
  -d /path/to/database.fasta \
  -f unite \
  -min 200 \
  -max 300 \
  -T 8 \
  -t 2 \
  -k 15 \
  -w 10
```

#### Required Arguments

- `-i, --input_dir`: Input directory containing barcodeXX folders and unclassified sequences
- `-o, --output_dir`: Output directory for results
- `-d, --db_path`: Path to reference database FASTA file

#### Optional Arguments

**Filter Parameters:**
- `-min, --min_len`: Minimum sequence length (default: 200)
- `-max, --max_len`: Maximum sequence length (default: 300)

**Performance:**
- `-T, --total_threads`: Total CPU threads (default: 8)
- `-t, --threads_per_job`: Threads per single mapping job (default: 2)

**minimap2 Parameters:**
- `-k, --kmer_size`: K-mer size (default: 15)
- `-w, --window_size`: Minimizer window size (default: 10)

**Database:**
- `-f, --format`: Taxonomy format - choices: unite, silva, eukariome, cbs, none (default: unite)
- `--force-tax-map`: Force recreation of taxonomy map (ignores cache)

**Other:**
- `--debug`: Enable verbose debug output

## Pipeline Workflow

### Step 0: Taxonomy Parsing
- Parses reference database headers
- Creates taxonomy map (TSV format)
- Caches map for future runs

### Step 1: Sequence Filtering
- Scans input directories (barcodeXX, unclassified)
- Filters sequences by length
- Merges sequences per barcode
- Generates filtering statistics report

### Step 1.5: Filtering Report
- Creates Excel report with:
  - Mean/median lengths per barcode
  - Good/bad sequence counts
  - Total sequence counts

### Step 2: Mapping & Tabeling
- Maps sequences to reference database using minimap2
- Filters alignments with samtools
- Generates count tables per sample

### Step 2.5: Index Reformatting
- Reformats OTU IDs based on database format
- Handles format-specific ID parsing

### Step 3: Results Aggregation
- Aggregates all samples into single OTU table
- Joins with taxonomy information
- Sorts barcodes naturally (barcode01, barcode02, ..., barcode10, unclassified)
- Saves final Excel file

### Step 4: HTML Report Generation
- Interactive Plotly visualizations
- Filter statistics (good vs. bad counts, length distributions)
- Taxonomic composition (Family, Genus, Species)
- Top 25 OTU heatmap
- Full OTU table

## Input Directory Structure

```
input_data/
├── barcode01/
│   ├── reads_001.fastq.gz
│   ├── reads_002.fastq.gz
│   └── ...
├── barcode02/
│   └── ...
├── unclassified/
│   └── ...
└── ...
```

## Output Directory Structure

```
output_folder/
├── mapping/               # SAM files from minimap2
├── results/               # Per-sample count tables (*.tmp.txt)
├── filtered/              # Good sequences (merged per barcode)
├── unfiltered/            # Bad sequences (filtered out)
├── report_filtering_200_300.xlsx  # Filter statistics
├── OTU_Table_[name]_[db].xlsx     # Final OTU table with taxonomy
├── Report_[name]_[db].html        # Interactive HTML report
└── [db]_taxonomy_map.tsv          # Taxonomy cache
```

## Supported Database Formats

### UNITE (Fungi)
Header format: `>accession|GB_ID|...taxonomy_string`
- Example: `>SH1234567|KX123456|...|k__Fungi;p__Basidiomycota;...`

### SILVA (Bacteria/Archaea)
Header format: `>accession taxonomy_string`
- Example: `>AB123456 Bacteria;Proteobacteria;...`

### CBS (Fungi)
Custom format for Westerdijk Institute database
- Extracts genus/species from structured headers

### Eukariome
Custom eukaryotic database format

### None
No taxonomy parsing - uses raw sequence IDs

## Performance Optimization

### Threading Strategy

The pipeline uses a two-level threading model:

1. **Total Threads (`-T`)**: Overall parallelization level
   - Filter step: Processes multiple barcodes simultaneously
   - Mapping step: Determines parallel mapping jobs

2. **Threads per Job (`-t`)**: minimap2 threading
   - Each mapping job uses this many threads
   - Calculate parallel jobs as: `total_threads / threads_per_job`

**Example configurations:**

- **8-core machine, few large barcodes:**
  ```bash
  -T 8 -t 4  # 2 parallel jobs, 4 threads each
  ```

- **16-core machine, many small barcodes:**
  ```bash
  -T 16 -t 2  # 8 parallel jobs, 2 threads each
  ```

### Memory Considerations

- **Filtering**: Low memory (streaming processing)
- **Mapping**: Depends on database size (~2-4GB for UNITE)
- **DuckDB aggregation**: Efficient in-memory processing
- **Large datasets**: Consider reducing parallel jobs

## Report Generator (Standalone)

Generate HTML reports from existing results:

```bash
python report_generator.py \
  -f report_filtering_200_300.xlsx \
  -otu OTU_Table_final.xlsx \
  -o MyReport.html \
  -pn "My Project Name"
```

Arguments:
- `-f`: Filter report Excel file
- `-otu`: OTU table Excel file
- `-o`: Output HTML path
- `-pn`: Project name for report title

## Troubleshooting

### Common Issues

**1. "minimap2 not found"**
```bash
# Check installation
which minimap2
# Install with conda
conda install -c bioconda minimap2
```

**2. "Cannot create taxonomy map"**
- Check database file format matches `-f` argument
- Ensure write permissions in database or output folder
- Use `--force-tax-map` to regenerate

**3. "No 'good' sequences after filtering"**
- Check min/max length settings match your data
- Use `--debug` to see detailed filtering output
- Verify input FASTQ files are not empty/corrupt

**4. "Missing Python module"**
```bash
pip install biopython pandas duckdb openpyxl
```

**5. GUI crashes on macOS**
- Install PyQt5: `pip install PyQt5`
- May need: `brew install python-tk`

### Debug Mode

Enable verbose logging:
```bash
python metaGenomics_new.py ... --debug
```

Shows:
- Exact commands executed
- File paths being processed
- Detailed error messages
- Progress for each step

## File Organizer Utility

Organize raw .gz files into barcoded folders:

```bash
python fileorganizer.py /path/to/fastq_files
```

Transforms:
```
folder/
├── sample_barcode01.fastq.gz
├── sample_barcode02.fastq.gz
```

Into:
```
folder/
├── sample_barcode01/
│   └── sample_barcode01.fastq.gz
├── sample_barcode02/
│   └── sample_barcode02.fastq.gz
```

## Citation

If you use this pipeline in your research, please cite:

```
Scarponi, R. (2025). LOREON: Long-Read ONT Metagenomic Pipeline (v2.6).
University of Perugia. https://github.com/iscka/loreon
```

A `CITATION.cff` file is included for automated citation tools (Zenodo, GitHub).

## License

[License to be determined - see LICENSE file]

## Authors

- Roberto Scarponi - University of Perugia
- Developed in collaboration with the Microbiology Laboratory
- Based on research by Angela Conti, Gianluigi Cardinali, and colleagues

## Contact

For bug reports and feature requests, please open an issue on GitHub.

## Acknowledgments

- minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100.
- UNITE database: Abarenkov et al. (2024). UNITE: a database providing web-based methods for the molecular identification of ectomycorrhizal fungi.
- DuckDB: High-performance analytical database system
- Plotly: Interactive visualization library

## Version History

### v2.6 (Current)
- Merge-by-directory filtering
- Natural barcode sorting
- Force taxonomy map recreation option
- Improved GUI with persistent worker thread
- Enhanced error handling

### v2.x
- DuckDB-based OTU aggregation
- Interactive HTML reports
- PyQt5 GUI interface
- Multi-format database support

### v1.x
- Initial release
- Basic filtering and mapping
