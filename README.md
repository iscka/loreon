# LOREON — Long-Read ONT MetageNOmic Pipeline

LOREON is a cross-platform bioinformatics application for processing **Oxford Nanopore Technologies (ONT)** long-read amplicon metagenomic data. It was originally designed for **fungal ITS** barcoding and now also supports **bacterial 16S**, and eukaryotic **18S (SSU)** / **LSU** markers.

The application takes demultiplexed FASTQ reads (one sub-folder or one file per sample), optionally filters them by length, maps them against a reference database with **minimap2**, counts the primary alignments per reference (optionally filtered by identity and read coverage), aggregates the results into an **OTU table** with **DuckDB**, and produces a self-explanatory **interactive HTML report**. Everything is driven either from a **PyQt5 graphical interface** or from the **command line**, and can run natively on Linux/macOS or inside a **Docker container** on Windows.

> **Name.** LOREON is not related to LoReAn (Long-Read Annotation; Cook et al., 2019, *Plant Physiology* 179:38–54), a pipeline for the annotation of eukaryotic genomes.

---

## Table of Contents

- [Key Features](#key-features)
- [How It Works (Pipeline Overview)](#how-it-works-pipeline-overview)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Example Dataset](#example-dataset)
- [Graphical Interface (GUI)](#graphical-interface-gui)
- [Command-Line Reference](#command-line-reference)
- [Amplicon Presets](#amplicon-presets)
- [Supported Database Formats](#supported-database-formats)
- [Input & Output Layout](#input--output-layout)
- [Interactive HTML Report](#interactive-html-report)
- [Standalone Report Generator](#standalone-report-generator)
- [Performance & Threading](#performance--threading)
- [Performance Profiling](#performance-profiling)
- [Windows / Docker Mode](#windows--docker-mode)
- [Building Distributables](#building-distributables)
- [Utilities](#utilities)
- [Troubleshooting](#troubleshooting)
- [Project Structure](#project-structure)
- [Citation](#citation)
- [License](#license)
- [Authors & Acknowledgments](#authors--acknowledgments)

---

## Key Features

- **Length filtering** — configurable minimum/maximum read length, with ready-made presets for common amplicon markers (ITS1, ITS2, full ITS, ITS–LSU, rRNA operon, 16S V3–V4, full 16S, 18S V4/V9). No window is applied unless one is chosen; the run reports the observed read lengths and warns when a window discards most of the reads.
- **High-throughput mapping** — minimap2 with the `map-ont` preset, piped directly into `samtools sort` to avoid intermediate SAM files on disk.
- **Alignment-quality filtering** — optional post-mapping filters on **percent identity** and **read (query) coverage**, computed by LOREON from the CIGAR string and the `NM` tag of each alignment (independent of the samtools version), plus removal of unmapped/secondary/supplementary alignments. Both filters are disabled by default.
- **Parallel execution** — a two-level threading model parallelises both the filtering stage (per barcode) and the mapping stage (multiple concurrent minimap2 jobs).
- **OTU aggregation with DuckDB** — per-sample count tables are pivoted, joined with taxonomy and written to a single Excel/TSV OTU table using an out-of-core analytical engine.
- **Interactive HTML report** — a single self-describing HTML file with Plotly charts: mapping QC, length-filter QC, taxonomic composition (Family/Genus/Species), a top-OTU heatmap and a searchable OTU table.
- **Multiple database formats** — UNITE, SILVA, CBS (Westerdijk), Eukariome, and a raw "none" mode.
- **Taxonomy caching** — the parsed taxonomy map is written once per database and reused on subsequent runs.
- **Cross-platform GUI** — PyQt5 interface with splash screen, system-tray icon, live execution log, progress bar, hardware-aware thread hints and an embedded license/credits viewer.
- **Windows Docker mode** — the GUI runs natively on the host while minimap2, samtools and the pipeline execute inside a reproducible Ubuntu container with pinned tool versions (minimap2 2.28, samtools 1.21).
- **Explicit failures** — the exit status of every external tool is checked; a failing sample is reported with the last lines of the tool's error output, progress is printed periodically during mapping, and an optional time limit per sample stops a stuck job.
- **Performance profiling** — optional per-step metrics (wall time, CPU time, memory, throughput) exported to JSON and to a formatted Word document.
- **Disk-space controls** — optional deletion of intermediate BAM/SAM files and gzip compression of merged FASTQ intermediates.

---

## How It Works (Pipeline Overview)

| Step | Name | What happens |
|------|------|--------------|
| **0** | Taxonomy parsing | Reference FASTA headers are parsed into a taxonomy map (`<db>_taxonomy_map.tsv`). The map is cached and reused unless `--force-tax-map` is set. |
| **1** | Sequence filtering | Every sample of the input folder (a sub-folder such as `barcodeXX`, or a single FASTQ file) is scanned; reads outside the `[min_len, max_len]` window, when one is set, are discarded. Kept reads are merged per sample into a single FASTQ (optionally gzipped). |
| **1.5** | Filter report | Per-barcode statistics (good/bad counts, mean & median length) are written to an Excel report. |
| **2** | Mapping & tabulation | Each merged FASTQ is mapped with minimap2 (`map-ont`), piped into `samtools sort`. The sorted BAM is then read once: primary alignments are counted per reference and, if requested, filtered by identity and read coverage. Pre- and post-filter statistics and a per-sample count table are produced; the alignments that pass are written to a filtered BAM. All figures, with the tool versions, are stored in `mapping_stats.json`. |
| **3** | Results aggregation | DuckDB pivots and joins all per-sample count tables with taxonomy into one OTU table, with barcodes sorted naturally (`barcode01 … barcode10 … unclassified`). Reference-ID reformatting per database format happens in-query. |
| **4** | HTML report | An interactive Plotly HTML report is generated from the OTU table, the filter report and `mapping_stats.json`. |

---

## System Requirements

### External tools (Linux / macOS, native mode)

Must be installed and available on `PATH`:

- [minimap2](https://github.com/lh3/minimap2) **≥ 2.24**
- [samtools](http://www.htslib.org/) **≥ 1.13**

These are the oldest versions LOREON was tested with; the analyses of the article used minimap2 2.28 and samtools 1.21, which are also the versions built into the Docker image. The versions found on `PATH` are printed at start-up and recorded in `mapping_stats.json`; a run with an older version stops with an explicit message.
- `awk` / `gawk` (pre-installed on most Unix systems)

### Python

Python **≥ 3.9**.

### Windows

- [Python](https://www.python.org/downloads/) ≥ 3.9, only to run the GUI from the source code: the application installed with the Windows installer does not need it
- [Docker Desktop](https://www.docker.com/products/docker-desktop) with the WSL 2 backend — minimap2, samtools and the pipeline run inside the container, so they do **not** need to be installed on the host.

### Python dependencies

Installed by the conda environment file or via `pip install -r requirements.txt`:

| Package | Purpose |
|---------|---------|
| `biopython` ≥ 1.79 | FASTA/FASTQ parsing |
| `pandas` ≥ 1.3, `numpy` ≥ 1.21 | Data handling |
| `duckdb` ≥ 0.8 | OTU aggregation engine |
| `openpyxl` ≥ 3.0, `xlsxwriter` ≥ 3.0 | Excel output |
| `isal` ≥ 1.0 *(optional)* | Faster gzip decompression (falls back to stdlib) |
| `PyQt5` ≥ 5.15 | Graphical interface |
| `plotly` ≥ 5.0, `jinja2` ≥ 3.0 | Interactive HTML report |
| `python-docx` ≥ 1.1 | Word profiling report |

---

## Installation

Download or clone the repository, open a terminal in its folder and follow one of the two routes. Then check the result with the [example dataset](#example-dataset).

### With Conda (recommended for Linux/macOS)

```bash
conda env create -f environment.yml
conda activate loreon
python pipeline_gui.py
```

`environment.yml` takes everything from the conda-forge and bioconda channels and pins Python 3.11, PyQt 5.15, minimap2 2.28 and samtools 1.21. The environment contains the Qt libraries, so no system package has to be installed: on clean Ubuntu 22.04 and 24.04 the display plugin of Qt finds all its libraries, the whole test suite passes and the interface opens. The commands were tested with [Miniforge](https://github.com/conda-forge/miniforge).

### With pip / virtualenv

```bash
sudo apt install python3-venv      # Ubuntu and Debian: the venv module is a separate package
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python pipeline_gui.py
```

With this route:

- minimap2 and samtools must be installed separately, unless Docker mode is used;
- PyQt5 is published on PyPI for x86-64 Linux, macOS and Windows; on ARM Linux use the conda route;
- on Linux PyQt5 uses the display libraries of the system. A desktop installation already has most of them; `libxcb-xinerama0` is the one most often missing on Ubuntu. If some are missing, LOREON does not open the window and prints their names with the `apt` command that installs them. On clean Ubuntu 22.04 and 24.04 that command is enough to open the interface.

### Checks at start-up

`python pipeline_gui.py` explains what is wrong, instead of failing with a Qt error, when PyQt5 is not installed, when system libraries are missing, and when there is no graphical session (for example over `ssh` without display forwarding; use the [command line](#command-line-reference) there). When the window opens, the log lists the components that are missing for an analysis on this computer (minimap2, samtools, Python packages, or Docker in Docker mode) and how to install them; **START PIPELINE** shows the same list instead of starting a run that cannot succeed.

The bundled `venv_manager.py` can also create and validate a dedicated virtual environment automatically; the packaged application for Linux and macOS uses it to locate a suitable Python interpreter for the pipeline subprocesses. On Windows the pipeline runs in Docker and no virtual environment is created.

---

## Quick Start

**Check the installation with the example dataset:**

```bash
python metaGenomics_new.py --example
```

**GUI (recommended):**

```bash
python pipeline_gui.py
```

**Command line:**

```bash
python metaGenomics_new.py \
  -i /path/to/input_data \
  -o /path/to/output \
  -d /path/to/database.fasta \
  -f unite \
  --preset its2 \
  -T 8 -t 2
```

---

## Example Dataset

The folder [`example_data`](example_data) holds three samples of 500 reads from yeast mock communities, a subset of the UNITE database and the OTU table expected for them (2.3 MB in all; sources and licences are in its [README](example_data/README.md)).

- Command line: `python metaGenomics_new.py --example` analyses it with the default settings, writes the results to `./loreon_example_output` (or to the folder given with `-o`) and compares the read counts with the expected table. The run takes less than a minute and ends with `The installation works.` and exit status `0`; if the tables differ, the differing cells are listed and the exit status is `1`.
- Interface: **RUN EXAMPLE DATASET** fills the form, runs the analysis in `LOREON_example_output` in the home folder and shows the outcome of the comparison. It works in Docker mode as well.
- Docker image: `docker run --rm loreon:latest python3 /app/metaGenomics_new.py --example -o /tmp/example`.

Any analysis can be checked in the same way against a table written by an earlier run: `--expected-table <OTU_Table.tsv>`.

---

## Graphical Interface (GUI)

Launch with `python pipeline_gui.py`. The window is organised into numbered groups:

1. **Paths (mandatory)** — pick the input folder, the reference-database FASTA and the output folder.
2. **Pipeline parameters**
   - Database format (UNITE, SILVA, Eukariome, CBS, none).
   - **Debug Mode** — verbose logging.
   - **Open HTML report when done** — auto-opens the report on completion.
   - **Force Taxonomy Map Recreation** — ignores the cached taxonomy map.
   - **Enable Performance Profiling** — writes the profiling report.
   - **Delete BAM/SAM after analysis** — frees disk space.
   - **Retain unmapped reads** — saves the reads that are not counted, one FASTQ file per sample (see [Reads that are not counted](#reads-that-are-not-counted)).
   - **Flat SH (UNITE only)** — collapses reference accessions sharing the same UNITE Species Hypothesis.
   - **Read & Alignment Filter** group — the **Length filter** menu offers *No length filter* (default), *Custom range* and the [amplicon presets](#amplicon-presets); **Scan read lengths...** reads a sample of the input folder and shows the length distribution, the share of reads each preset would keep and a window that can be applied with one click. **Min Identity (%)** and **Min Read Cov (%)** are disabled (0) by default: read the caution in [Alignment-quality filters](#alignment-quality-filters) before setting them.
3. **Performance (threads)** — the group header shows the number of detected CPU cores; set total threads and threads-per-job, with a live hint for the resulting number of parallel jobs, and optionally a maximum time per sample.
4. **Docker (required on Windows)** — enable Docker mode, **Check Image** and **Build Image** (first run only).
5. **START PIPELINE / RUN EXAMPLE DATASET / STOP** — run the analysis of the form, run the [example dataset](#example-dataset), or terminate the pipeline. In Docker mode **STOP** also stops the container.
6. **Execution Log** — real-time streaming log with a progress bar.

Additional conveniences: a **splash screen** on launch, a **system-tray icon** with show/quit actions, and an embedded **License & Credits** viewer.

---

## Command-Line Reference

`python metaGenomics_new.py [options]`

### Required

Not needed with `--example`.

| Flag | Description |
|------|-------------|
| `-i, --input_dir` | Input directory: one sub-folder per sample, or one FASTQ file per sample (see [Input & Output Layout](#input--output-layout)). |
| `-o, --output_dir` | Output directory for all results. |
| `-d, --db_path` | Path to the reference database FASTA. |

### Database

| Flag | Description |
|------|-------------|
| `-f, --format` | Taxonomy format: `unite` (default), `silva`, `eukariome`, `cbs`, `none`. |
| `--force-tax-map` | Rebuild the taxonomy map, ignoring the cache. |
| `--flat-sh` | UNITE only: aggregate OTU rows by Species Hypothesis. |

### Length filter

| Flag | Default | Description |
|------|---------|-------------|
| `-min, --min_len` | no limit | Minimum read length (bp). |
| `-max, --max_len` | no limit | Maximum read length (bp); `0` means no upper limit. |
| `--preset` | — | Amplicon preset that sets the length window of a common marker (see [presets](#amplicon-presets)). A bound given with `-min` or `-max` replaces the corresponding bound of the preset. |

Without `-min`, `-max` and `--preset` no length filter is applied. After the filtering step the run prints the observed read lengths (median and the range holding 90 % of the reads). If a window discards more than half of the reads, a warning suggests a window and the presets that would keep them.

### Alignment-quality filters

| Flag | Default | Description |
|------|---------|-------------|
| `--min-percent-identity PCT` | 0 (disabled) | Minimum alignment percent identity, `(L − NM) / L × 100 ≥ PCT`. |
| `--min-ref-coverage PCT` | 0 (disabled) | Minimum **read (query) coverage**: percentage of the read that must align, `L / read length × 100 ≥ PCT`. *(The flag name is retained for backward compatibility; the value refers to the read, not to the reference.)* |

`L` is the read length without soft-clipped bases and `NM` is the edit distance reported by minimap2. Both values are computed by LOREON while it counts the alignments, so the result is the same with every supported samtools version. Unmapped, secondary and supplementary alignments are always left out.

> **Use these two filters with caution.** With Oxford Nanopore reads, a threshold greater than 0 can discard a large part of the reads of a sample, or all of them. Both filters are disabled by default and the published analyses were run without them. Set a threshold only for a known reason, and after measuring its effect on a control.

How many reads a threshold removes depends on the data as much as on the organisms:

- **Identity** follows the accuracy of the reads, hence the flow cell, the kit and the basecalling model, and the distance between the organism and its closest reference. The thresholds in use for short reads, such as 97 %, remove a large share of nanopore reads even with recent chemistries.
- **Read coverage** counts every base of the read that is not aligned. Adapters and barcodes left on the reads, about 100 bases at the start of each read in the two runs below, and references shorter than the amplicon, as for ITS–LSU amplicons mapped on an ITS-only database, keep it low whatever the accuracy of the reads.

Share of the mapped reads that each threshold removes, measured on two sets of yeast mock communities:

| Threshold | MinION run of 2021, ITS–LSU amplicons on UNITE | R10.4.1 run of 2025, high-accuracy basecalling, on Eukaryome |
|-----------|-----------------------------------------------|---------------------------------------------------------------|
| Identity ≥ 90 % | 49 % | 4 % |
| Identity ≥ 95 % | 77 % | 15 % |
| Identity ≥ 97 % | 93 % | 29 % |
| Read coverage ≥ 50 % | 71 % | 46 % |
| Read coverage ≥ 80 % | 92 % | 99 % |

A sample whose organisms are all in the database, such as a pure culture, is not exempt. In the same runs the reads assigned to a species of the mock lost 12–32 % (2021) and 1–4 % (2025) at 90 % identity, 48–68 % and 6–13 % at 95 %, and 96–100 % at 80 % read coverage.

A pure culture or a mock community is the right control. Analyse it with the candidate threshold, sequenced and basecalled like the samples, and read the *Quality Retention (%)* column of the report: it is the share of mapped reads that the thresholds keep. A read-coverage threshold also needs reads without adapters and barcodes and references that span the whole amplicon. With `--retain-unmapped` the reads removed by the thresholds are saved in `below_threshold_seq/` and can be examined.

### Performance

| Flag | Default | Description |
|------|---------|-------------|
| `-T, --total_threads` | 8 | Total CPU threads. |
| `-t, --threads_per_job` | 2 | Threads per single mapping job (parallel jobs ≈ `T / t`). |
| `--sample-timeout MIN` | 0 (no limit) | Maximum time, in minutes, for mapping and counting one sample. A sample that exceeds it is stopped and reported as failed. |
| `--progress-interval SEC` | 60 | A progress message is printed during mapping whenever this many seconds pass without a sample finishing (`0` = none). |

### minimap2

| Flag | Default | Description |
|------|---------|-------------|
| `-k, --kmer_size` | 15 | K-mer size. |
| `-w, --window_size` | 10 | Minimizer window size. |

### Other

| Flag | Description |
|------|-------------|
| `--debug` | Verbose output (shows exact commands and per-step detail). |
| `--profile` | Enable performance profiling (JSON + Word report). |
| `--delete-bam` | Delete BAM/SAM files after tabulation to save disk space. |
| `--retain-unmapped` | Save the reads that are not counted, one FASTQ file per sample: reads that do not align in `unmapped_seq/` and, when an alignment-quality filter is set, aligned reads below the thresholds in `below_threshold_seq/` (see [Reads that are not counted](#reads-that-are-not-counted)). |
| `--gzip-intermediate` | Write merged FASTQ intermediates as `.fastq.gz` (minimap2 reads gzip natively). |
| `--example` | Analyse the [example dataset](#example-dataset) with the default settings and compare the result with its expected table. |
| `--expected-table TSV` | After the run, compare the read counts of the OTU table with this table (a `.tsv` written by LOREON); exit status `1` if they differ. |

### Exit status

| Code | Meaning |
|------|---------|
| `0` | All samples were processed. |
| `1` | The run stopped: invalid input, missing or outdated tool, no reads left after the length filter, every sample failed, the OTU table could not be written, or it differs from the table given with `--expected-table` or `--example`. |
| `3` | The OTU table was written, but one or more samples failed (unreadable input, mapping or counting error, time limit) and are missing from it. The failed samples and the cause of each failure are listed in the log. |

---

## Amplicon Presets

`--preset` (or the **Length filter** menu of the GUI) sets the length window of a common marker. A bound given with `-min` or `-max` replaces the corresponding bound of the preset.

| Preset | min–max (bp) | Marker |
|--------|--------------|--------|
| `its1` | 150–350 | Fungal ITS1 |
| `its2` | 150–400 | Fungal ITS2 |
| `its-full` | 400–900 | Full-length ITS |
| `its-lsu` | 900–1800 | ITS plus the 5′ end of the LSU (e.g. ITS1–NL4 amplicons) |
| `rrna-operon` | 3000–6500 | rRNA operon, SSU to LSU (e.g. NS1–RCA95m amplicons) |
| `16s-v3v4` | 380–530 | Bacterial 16S V3–V4 |
| `16s-full` | 1200–1700 | Full-length 16S |
| `18s-v9` | 100–220 | Eukaryotic 18S V9 |
| `18s-v4` | 350–550 | Eukaryotic 18S V4 |

The windows are indicative: check them against the read lengths of your own run, which are printed after the filtering step and shown by **Scan read lengths...** in the GUI.

---

## Supported Database Formats

| Format | Target | Header handling |
|--------|--------|-----------------|
| **UNITE** | Fungi (ITS) | Parses `>accession|GB_ID|SH|...|k__Fungi;p__…` headers; optional SH-level aggregation via `--flat-sh`. |
| **SILVA** | Bacteria / Archaea (16S/18S) | Parses `>accession taxonomy;string`. |
| **CBS** | Fungi | Westerdijk (CBS-KNAW) structured headers; extracts genus/species. |
| **Eukariome** | Eukaryotes | Custom eukaryotic header format. |
| **none** | Any | No taxonomy parsing — raw sequence IDs are used as OTU identifiers. |

---

## Input & Output Layout

**Input.** Two layouts are accepted.

One sub-folder per sample, as written by the ONT basecallers; all FASTQ files below a sub-folder are merged and the folder name is the sample name:

```
input_data/
├── barcode01/
│   ├── reads_001.fastq.gz
│   └── reads_002.fastq.gz
├── barcode02/
│   └── ...
└── unclassified/
    └── ...
```

One FASTQ file per sample placed directly in the folder, as downloaded from the Sequence Read Archive; the file name without its extensions is the sample name:

```
input_data/
├── SRR20631366.fastq
└── SRR20631367.fastq.gz
```

Accepted extensions are `.fastq`, `.fq`, `.fastq.gz` and `.fq.gz`. When the folder holds both sub-folders with reads and loose FASTQ files, the sub-folders are used and the loose files are listed in the log as ignored.

**Output:**

```
output_folder/
├── mapping/                          # minimap2 / sorted BAM intermediates
├── tabeling/                         # tabulation: filtered BAM (alignments that were counted) and index, per sample
│   └── results/                      # per-sample count tables (*.tmp.txt)
├── good_seq/                         # reads kept by the length filter (merged per sample)
├── bad_seq/                          # reads discarded by the length filter
├── unmapped_seq/                     # reads that do not align (with --retain-unmapped)
├── below_threshold_seq/              # aligned reads below the quality thresholds (with --retain-unmapped)
├── report_filtering_<min>_<max>.xlsx # length-filter statistics
├── mapping_stats.json                # per-sample mapping & quality-filter stats
├── OTU_Table_<name>_<db>.xlsx        # final OTU table with taxonomy
├── Report_<name>_<db>.html           # interactive HTML report
├── <db>_taxonomy_map.tsv             # cached taxonomy map
├── profile_report.json               # profiling data (with --profile)
└── profile_report.docx               # profiling Word report (with --profile)
```

An output folder can be reused: the per-sample count tables left by an earlier run are removed before mapping, so the OTU table always describes the current run only.

### Reads that are not counted

By default the reads that do not enter the OTU table are only counted. With **Retain unmapped reads** in the interface, or `--retain-unmapped` on the command line, they are saved as FASTQ, one file per sample, for later analysis:

| Folder | Reads | Written when |
|--------|-------|--------------|
| `unmapped_seq/` | Reads for which minimap2 finds no alignment. | The option is on. |
| `below_threshold_seq/` | Aligned reads removed by `--min-percent-identity` or `--min-ref-coverage`. | The option is on and at least one of the two thresholds is greater than 0. |

For every sample, the reads counted in the OTU table plus the reads of the two files are the reads that passed the length filter. Each read is written with its name, sequence and qualities as in the input, whatever the strand of its alignment; the comment after the read name is not kept. A sample with no such reads gets no file. The files are `.fastq`, or `.fastq.gz` with `--gzip-intermediate`, and they are kept when `--delete-bam` removes the alignment files. The HTML report lists, for each sample, the path of the two files and the number of reads they hold; the same numbers are in `mapping_stats.json`.

With the default settings almost every read aligns, so `unmapped_seq/` is small: the reads of taxa absent from the database are attributed to the closest reference and are counted there. An identity threshold moves them to `below_threshold_seq/`, together with the less accurate reads of the taxa that are present.

---

## Interactive HTML Report

A single self-contained HTML file with the following sections:

1. **Mapping QC** *(when `mapping_stats.json` is available)* — KPI cards (total input reads, mapped reads, mapping rate, reads after quality filter, retention rate), a badge showing the active thresholds, a stacked "mapped vs. unmapped" bar chart, a "quality-filter retention" chart and a per-sample statistics table. When the reads that are not counted were saved (`--retain-unmapped`), a second table gives the path of each file and the number of reads it holds.
2. **Filter QC** — good vs. bad read counts and mean/median length per barcode.
3. **Metagenomic analysis** — taxonomic composition at Family / Genus / Species level (interactive tabs), a log-normalised heatmap of the top OTUs, and a full searchable OTU table.

The report closes with a note on the reference database: each read is assigned to the reference of its best alignment, with no fallback to a higher taxonomic rank, so the reads of taxa absent from the database are attributed to the closest reference available and a high mapping rate is not evidence that the database is complete.

---

## Standalone Report Generator

Regenerate a report from existing results without re-running the pipeline:

```bash
python report_generator.py \
  -f report_filtering_200_300.xlsx \
  -otu OTU_Table_final.xlsx \
  -o MyReport.html \
  -pn "My Project" \
  -ms mapping_stats.json
```

| Flag | Description |
|------|-------------|
| `-f, --filter_report` | Length-filter report (TSV or XLSX). **Required.** |
| `-otu, --otu_table` | OTU table Excel file. **Required.** |
| `-o, --output_html` | Output HTML path. **Required.** |
| `-pn, --project_name` | Project name shown on the report. |
| `-ms, --mapping_stats` | `mapping_stats.json` to include the Mapping QC section. |
| `--min-len`, `--max-len` | Length-filter values used (for display; `0` disables). |

---

## Performance & Threading

LOREON uses a two-level threading model:

- **Total threads (`-T`)** — overall parallelism. The filter step processes multiple barcodes at once; the mapping step launches multiple concurrent minimap2 jobs.
- **Threads per job (`-t`)** — threads assigned to each minimap2 job. The number of parallel mapping jobs is approximately `T / t`.

**Example configurations:**

```bash
-T 8  -t 4   # 8-core machine, few large barcodes  → 2 jobs × 4 threads
-T 16 -t 2   # 16-core machine, many small barcodes → 8 jobs × 2 threads
```

Memory notes: filtering is streaming (low memory); mapping memory scales with the database (~2–4 GB for UNITE); DuckDB aggregation is efficient and spills to disk when needed. Use `--delete-bam` and `--gzip-intermediate` to reduce disk footprint on large ONT runs.

---

## Performance Profiling

Enable with `--profile` (CLI) or the **Enable Performance Profiling** checkbox (GUI). Two files are written to the output folder:

- `profile_report.json` — raw structured metrics.
- `profile_report.docx` — formatted Word table.

Metrics collected per step include `wall_time_s`, `cpu_time_s`, `mem_mb`, plus step-specific counters such as `taxonomy_entries`, `sequences_total/good/bad`, `throughput_seq_s`, `files_to_map`, `success` and `failed`. A summary is also printed to the execution log.

---

## Windows / Docker Mode

On Windows, LOREON runs in mixed mode: the GUI is native, while minimap2, samtools and the pipeline run inside a Docker container built from the included `Dockerfile` (Ubuntu 22.04 with minimap2 2.28 and samtools 1.21 compiled from the official source releases, whose checksums are verified during the build, plus the Python pipeline). The first build downloads and compiles the two tools and takes a few minutes.

1. Install Docker Desktop (WSL 2 backend) and, to run the GUI from the source code, Python ≥ 3.9.
2. Launch the GUI (`python pipeline_gui.py` or the installed executable).
3. In group **3. Docker**, tick **Use Docker** and click **Build Image** (first run only) — this builds the `loreon:latest` image.
4. Select the input folder, database FASTA and output folder, then click **START PIPELINE**.

Under the hood the host folders are bind-mounted into the container as `/data/input` (read-only), `/data/output` and `/data/db`. Each run uses a container with its own name (`loreon-<id>`); **STOP** kills that container.

Mapping needs memory in proportion to the reference database. If a tool is reported as *killed*, raise the memory assigned to Docker (Docker Desktop → Settings → Resources) or lower the number of parallel jobs.

The helper script `run_loreon.sh` runs the containerised pipeline directly from a shell, creating `~/Documents/loreon/{data,results,db}` and mounting them automatically.

---

## Building Distributables

Platform build scripts are included:

| Script | Target |
|--------|--------|
| `build_exe.bat` / `build_installer.bat` / `installer.iss` | Windows executable and Inno Setup installer |
| `build_mac.sh` / `build_mac.command` / `build_dmg.sh` | macOS app bundle and DMG |
| `build_linux.sh` / `build_deb.sh` | Linux binary and `.deb` package |

The Windows build converts `loreon_app_icon_2.png` to an `.ico`, bundles the pipeline with PyInstaller and copies the Docker context files alongside the executable.

The repository holds no pre-built installers: build them from a tagged version with these scripts, or take them from the release page when a release provides them. All scripts and the application report the same version number, defined once as `LOREON_VERSION` in `OtuUtils.py`.

---

## Utilities

**File organizer** — turn a flat folder of barcoded `.gz` files into the per-barcode folder layout. LOREON also reads such a folder directly, one file per sample, so this step is optional:

```bash
python fileorganizer.py /path/to/fastq_files
```

```
sample_barcode01.fastq.gz   →   sample_barcode01/sample_barcode01.fastq.gz
sample_barcode02.fastq.gz   →   sample_barcode02/sample_barcode02.fastq.gz
```

---

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| **`minimap2` not found** | Install it (`conda install -c bioconda minimap2`) or use Docker mode. |
| **`minimap2`/`samtools` is too old** | LOREON needs minimap2 ≥ 2.24 and samtools ≥ 1.13: `conda install -c bioconda minimap2=2.28 samtools=1.21`, or use Docker mode. |
| **No FASTQ files found in the input folder** | Use one sub-folder per sample or one FASTQ file per sample; see [Input & Output Layout](#input--output-layout). |
| **Cannot create taxonomy map** | Check that `-f` matches the database format; ensure write permission; retry with `--force-tax-map`. |
| **No reads left after the length filter** | The log reports the observed read lengths and suggests a window and the matching presets; set `-min`/`-max` or `--preset` accordingly (GUI: **Scan read lengths...**). |
| **Zero reads after quality filter** | The identity and read-coverage filters are disabled by default; if you enabled them, lower the thresholds or set them to `0`. Read coverage is low by construction when the reads are longer than the reference sequences. The Mapping QC section shows per-sample retention. |
| **Step 4 seems stuck** | A progress line (`Step 4 progress: n/N samples completed`) is printed whenever 60 s pass without a sample finishing, so a silent log is not expected. Use `--sample-timeout` (GUI: *Max min/sample*) to stop a sample that exceeds a time limit. |
| **A tool "was terminated by signal 9"**, or **"a worker process ended unexpectedly"** | The computer or the Docker container ran out of memory: assign more memory to Docker or use fewer parallel jobs (`-T`/`-t`). |
| **"n of N samples failed"** | The other samples are in the OTU table (exit status `3`); the log shows the failing tool, its exit status and the last lines of its error output next to each failed sample. |
| **Missing Python module** | `pip install -r requirements.txt`, or use the conda environment. |
| **The window does not open on Linux** | The terminal shows the reason: PyQt5 not installed, system libraries missing (with the `apt` command that installs them) or no graphical session. See [Installation](#installation). |
| **GUI crashes on macOS** | `pip install PyQt5`; you may also need `brew install python-tk`. |

Enable `--debug` for exact commands, processed paths and detailed per-step output.

---

## Project Structure

| File | Role |
|------|------|
| `pipeline_gui.py` | PyQt5 graphical interface (entry point). |
| `pipeline_worker.py` | Background thread that runs the pipeline / report as subprocesses (native or Docker). |
| `metaGenomics_new.py` | Core pipeline orchestrator and CLI. |
| `OtuUtils.py` | Filtering, mapping and tabulation utilities. |
| `ResultsReader.py` | DuckDB-based OTU aggregation. |
| `report_generator.py` | Interactive HTML report generator (uses `template.html`). |
| `pipeline_profiler.py` | Per-step performance profiling. |
| `PipelineLogger.py` | Logging helpers. |
| `venv_manager.py` | Virtual-environment discovery/creation. |
| `fileorganizer.py` | Barcode-folder organizer utility. |
| `Dockerfile` | Container image for the containerised pipeline. |
| `template.html` | Jinja2 template for the HTML report. |
| `requirements.txt` | Python dependencies. |
| `environment.yml` | Conda environment with pinned Python, PyQt, minimap2 and samtools. |
| `example_data/` | Example dataset with its expected OTU table. |
| `CHANGELOG.md` | Changes between versions. |
| `tests/` | Automated tests (`pytest`). |

### Running the tests

```bash
pip install pytest
pytest tests
```

The tests build a small synthetic dataset on the fly. Those that need minimap2 and samtools are skipped when the tools are not installed; the others replace the tools with stand-ins to check how failures are reported.

---

## Citation

If you use LOREON in your research, please cite:

> Scarponi, R. (2026). *LOREON: Long-Read ONT Metagenomic Pipeline* (v1.1.1). Zenodo. https://doi.org/10.5281/zenodo.22754395

```bibtex
@software{scarponi_loreon_2026,
  author    = {Scarponi, Roberto},
  title     = {{LOREON}: Long-Read {ONT} Metagenomic Pipeline},
  version   = {1.1.1},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22754395},
  url       = {https://doi.org/10.5281/zenodo.22754395}
}
```

A `CITATION.cff` file is included for automated citation tools (Zenodo, GitHub).

---

## License

LOREON is released under the **GNU General Public License v3.0** (GPL-3.0-only), the licence of the PyQt5 toolkit used by its interface — see the [`LICENSE`](LICENSE) file.

---

## Authors & Acknowledgments

- **Roberto Scarponi** — University of Perugia
- Developed in collaboration with the **Microbiology Laboratory**, Department of Biotechnology, University of Perugia
- Based on research by **Angela Conti**, **Gianluigi Cardinali**, and colleagues

**Acknowledgments**

- **minimap2** — Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34:3094–3100.
- **samtools / HTSlib** — Danecek et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2).
- **UNITE database** — Abarenkov et al. (2024).
- **SILVA database** — Quast et al. (2013); Yilmaz et al. (2014).
- **CBS database** — Westerdijk Fungal Biodiversity Institute (CBS-KNAW), Utrecht, The Netherlands.
- **DuckDB** — high-performance analytical database system.
- **Plotly** — interactive visualization library.
- **PyQt5** — cross-platform GUI toolkit.
