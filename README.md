# LOREON - Long-Read ONT Metagenomic Pipeline

A comprehensive bioinformatics pipeline for processing Oxford Nanopore Technologies (ONT) metagenomic sequencing data, specifically designed for fungal ITS barcode analysis.

## Features

- **Sequence Quality Filtering**: Length-based filtering with configurable min/max thresholds
- **Alignment Quality Filtering**: Post-mapping filters on percent identity and reference coverage applied directly via `samtools view -e`
- **High-Performance Mapping**: Utilizes minimap2 for fast and accurate sequence alignment
- **Parallel Processing**: Multi-threaded execution for filtering, mapping, and analysis
- **OTU Table Generation**: Automated OTU (Operational Taxonomic Unit) aggregation with DuckDB
- **Interactive Reports**: HTML reports with Plotly visualizations, including a dedicated **Mapping QC** section
- **Multiple Database Formats**: Support for UNITE, SILVA, CBS, Eukariome, and custom formats
- **GUI Interface**: User-friendly PyQt5 graphical interface with splash screen and system tray icon
- **Taxonomy Mapping**: Intelligent caching system for taxonomic information
- **Windows Support**: Native Docker-based execution on Windows via standalone executable
- **Performance Profiling**: Optional per-step metrics with JSON and Word report output

## System Requirements

### Software Dependencies 

** Linux/MacOS**
**External Tools** (must be installed separately):
- [minimap2](https://github.com/lh3/minimap2) >= 2.17
- [samtools](http://www.htslib.org/) >= 1.10
- awk (pre-installed on most Unix systems)

> **Windows users**: 
- [python](https://www.python.org/downloads/) >= 3.8
- [docker](https://www.docker.com/products/docker-desktop)
- the other tools run automatically inside the Docker container. See the [Windows section](#windows-docker-mode) below.

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
- python-docx >= 1.1.0 (for Word profiling reports)

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

---

## Windows — Docker Mode

On Windows, LOREON runs in mixed mode: minimap2, samtools and the pipeline execute inside a docker container, while the GUI runs natively on the host.

### Prerequisites

- [Docker Desktop for Windows](https://www.docker.com/products/docker-desktop/) (with WSL 2 backend)
- Python >= 3.8 on the host (for the GUI)

### Quick Start — Standalone Executable

The easiest way to run LOREON on Windows is to run the installer. The installer will create a shortcut to the pre-built executable:

1. Download or build `windows_setupv/LOREON_Setup_2.6.exe` 
2. The installer will check for required dependencies and install them if missing
3. Double-click `LOREON_Setup_2.6.exe`
4. Open the loreon GUI, check **Use Docker** and click **Build Image** (first run only)
5. Select input folder, database file and output folder
6. Click **START PIPELINE**


### Building the Executable

```bat
build_exe.bat
```

The script:
1. Installs PyInstaller and Pillow if missing
2. Converts `loreon_app_icon_2.png` to `.ico`
3. Builds `dist/LOREON/LOREON.exe` with the icon embedded
4. Copies all Docker context files into `dist/LOREON/`



## Usage

### GUI Mode (Recommended)

```bash
python pipeline_gui.py
```

The GUI provides:
1. **Path Selection**: Input folder, reference database, output folder
2. **Parameters Configuration**:
   - Database format (UNITE, SILVA, CBS, etc.)
   - **Read & Alignment Filter** group:
     - Length filter (min/max bp)
     - Min Identity (%) — minimum alignment percent identity (default: **95.0%**)
     - Min Ref Cov (%) — minimum reference coverage (default: **90.0%**)
   - Threading options
   - minimap2 parameters (k-mer, window size)
3. **Real-time Logging**: Progress tracking with status bar
4. **Automatic Report**: Opens HTML report upon completion (includes Mapping QC section)
5. **Performance Profiling**: To profile the pipeline performance
6. **Debug Mode**: To enable verbose debug output

### Command-Line Mode

```bash
python metaGenomics_new.py \
  -i /path/to/input_data \
  -o /path/to/output \
  -d /path/to/database.fasta \
  -f unite \
  -min 200 \
  -max 300 \
  --min-percent-identity 95 \
  --min-ref-coverage 90 \
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

**Length Filter Parameters:**
- `-min, --min_len`: Minimum sequence length (default: 200)
- `-max, --max_len`: Maximum sequence length (default: 300)

**Alignment Quality Filters:**
- `--min-percent-identity PCT`: Minimum alignment percent identity to keep a mapped read (default: **95.0**, set to 0 to disable).
  Formula applied via `samtools view -e`: `(alen - NM) / alen × 100 ≥ PCT`
- `--min-ref-coverage PCT`: Minimum reference sequence coverage to keep a mapped read (default: **90.0**, set to 0 to disable).
  Formula applied via `samtools view -e`: `alen / rlen × 100 ≥ PCT`

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
- `--profile`: Enable performance profiling (see [Performance Profiling](#performance-profiling))

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
- Maps sequences to reference database using minimap2 (`map-ont` preset)
- Converts SAM → sorted BAM (samtools)
- Collects pre-filter mapping statistics via `samtools flagstat`
- Applies alignment quality filters via `samtools view -e`:
  - **Percent Identity**: `(alen - NM) × 100 ≥ min_percent_identity × alen`
  - **Reference Coverage**: `alen × 100 ≥ min_ref_coverage × rlen`
  - Flag-based filter: `-F 0x904` (removes unmapped, secondary, supplementary)
- Collects post-filter mapping statistics
- Generates count tables per sample
- Saves per-sample mapping statistics to `mapping_stats.json`

### Step 2.5: Index Reformatting
- Reformats OTU IDs based on database format
- Handles format-specific ID parsing

### Step 3: Results Aggregation
- Aggregates all samples into single OTU table
- Joins with taxonomy information
- Sorts barcodes naturally (barcode01, barcode02, ..., barcode10, unclassified)
- Saves final Excel file

### Step 4: HTML Report Generation

The HTML report is organised in numbered sections:

**Section 1 — Mapping QC Report** *(generated when `mapping_stats.json` is available)*
- KPI cards: Total Input Reads, Mapped Reads, Mapping Rate, Reads After Quality Filter, Quality Retention Rate
- Info badge showing the active quality filter thresholds (identity, reference coverage)
- Stacked bar chart: **Mapped vs. Unmapped reads** per sample
- Stacked bar chart: **Quality filter retention** per sample *(shown only when at least one quality filter is active)*
- Per-sample mapping statistics table

**Section 2 — Filter QC Analysis**
- Filter statistics: good vs. bad counts, mean/median length per barcode

**Section 3 — Metagenomic Analysis**
- Taxonomic composition (Family, Genus, Species) — interactive tabs
- Top 25 OTU heatmap (log-normalised)
- Full OTU table (searchable)

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
├── tabeling/              # Sorted/filtered BAM files and intermediate results
│   └── results/           # Per-sample count tables (*.tmp.txt)
├── good_seq/              # Good sequences (merged per barcode, post length-filter)
├── bad_seq/               # Sequences discarded by length filter
├── report_filtering_200_300.xlsx  # Length filter statistics (Excel)
├── mapping_stats.json             # Per-sample mapping & alignment quality stats
├── OTU_Table_[name]_[db].xlsx     # Final OTU table with taxonomy
├── Report_[name]_[db].html        # Interactive HTML report (incl. Mapping QC)
├── [db]_taxonomy_map.tsv          # Taxonomy cache
├── profile_report.json            # Profiling data (if --profile)
└── profile_report.docx            # Profiling Word report (if --profile)
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

## Performance Profiling

Enable with `--profile` flag (CLI) or the **Enable Performance Profiling** checkbox in the GUI.

At the end of the run, two files are written to the output folder:

- `profile_report.json` — raw structured data
- `profile_report.docx` — formatted Word table (blue header, alternating rows, TOTAL footer)

### Metrics collected per step

| Metric | Unit | Description |
|--------|------|-------------|
| `wall_time_s` | seconds | Elapsed clock time for the step |
| `cpu_time_s` | seconds | CPU time consumed (user + system) |
| `mem_mb` | MB | Resident memory (RSS) at step end |
| `taxonomy_entries` | count | Entries parsed from the reference DB *(Step 0)* |
| `sequences_total` | count | Total sequences submitted to filter *(Step 1)* |
| `sequences_good` | count | Sequences passing the length filter *(Step 1)* |
| `sequences_bad` | count | Sequences discarded by the filter *(Step 1)* |
| `throughput_seq_s` | seq/s | Filter throughput *(Step 1)* |
| `files_to_map` | count | FASTQ files submitted to mapping *(Step 2)* |
| `success` | count | Files mapped successfully *(Step 2)* |
| `failed` | count | Files that failed mapping *(Step 2)* |

The summary is also printed in the execution log at the end of the run.

### CLI example

```bash
python metaGenomics_new.py \
  -i /data/input -o /data/output -d /data/db/unite.fasta \
  -f unite --profile
```

---

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
  -pn "My Project Name" \
  -ms mapping_stats.json
```

Arguments:
- `-f, --filter_report`: Length filter report Excel file (required)
- `-otu, --otu_table`: OTU table Excel file (required)
- `-o, --output_html`: Output HTML path (required)
- `-pn, --project_name`: Project name for report title (optional)
- `-ms, --mapping_stats`: Path to `mapping_stats.json` for the Mapping QC section (optional — section is omitted if not provided)

> The `mapping_stats.json` file is generated automatically by the pipeline in the output folder. Pass it to `report_generator.py` to include the Mapping QC section in reports regenerated manually.

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

**5. "Zero reads after quality filter"**
- The `--min-percent-identity` and `--min-ref-coverage` defaults are 95% and 90% respectively. These are strict thresholds; if your reads are very noisy, lower them (e.g. `--min-percent-identity 85 --min-ref-coverage 70`) or set to `0` to disable entirely.
- Use `--debug` to see the exact `samtools view -e` expression being applied.
- The Mapping QC section in the HTML report shows per-sample retention rates, which help diagnose overly aggressive filtering.

**6. "samtools view: expression filter not supported"**
- The `-e` expression filter requires **samtools ≥ 1.13**. Update samtools:
  ```bash
  conda install -c bioconda "samtools>=1.13"
  ```
- If you need to use an older samtools version, set both `--min-percent-identity 0` and `--min-ref-coverage 0` to disable the expression-based filter.

**7. GUI crashes on macOS**
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
- Developed in collaboration with the Microbiology Laboratory, department of Biotechnology, University of Perugia
- Based on research by Angela Conti, Gianluigi Cardinali, and colleagues

## Contact

For bug reports and feature requests, please open an issue on GitHub.

## Acknowledgments

- minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100.
- UNITE database: Abarenkov et al. (2024). UNITE: a database providing web-based methods for the molecular identification of ectomycorrhizal fungi.
- DuckDB: High-performance analytical database system
- Plotly: Interactive visualization library