# Changelog

## Unreleased

- Layout of the repository: the application, that is the Python scripts and the files they use at run time, is in `sources/`; the scripts that build the installers are in `scripts/`. The commands of the README are run from `sources/`. The installed application and the published installers do not change.

## 1.1.1 (2026-11-24)

- New option **Retain unmapped reads** in the interface, `--retain-unmapped` on the command line: the reads that are not counted are saved as FASTQ, one file per sample. Reads that do not align go to `unmapped_seq/`; when an identity or read-coverage threshold is set, the aligned reads below the thresholds go to `below_threshold_seq/`. In the interface the option is ticked by default; on the command line the reads are saved only when it is given. It does not change the OTU table.
- The HTML report and `mapping_stats.json` give the path of each of these files and the number of reads it holds.
- Interface: **Enable Performance Profiling** and **Delete BAM/SAM after analysis** are ticked by default. With profiling on, the performance report `performance_report.html` opens when the analysis ends.
- Windows: the installer no longer installs Python packages or asks for Python, and the installed application no longer creates a virtual environment at first launch. Both were unused, because the interface is self-contained and the pipeline runs in Docker.
- Windows: the installer includes the icon used by its shortcuts, and `build_installer.bat` no longer ends with an error after a successful compilation.
- macOS: the application bundle is signed again after `Info.plist` and the data files are added. Its signature was not valid, and macOS treated a downloaded copy as modified.
- macOS: the packaged application adds `/opt/homebrew/bin` and `/usr/local/bin` to `PATH`. Opened from the Finder, it did not find minimap2, samtools or Docker.
- Linux: the commands `loreon` and `loreon-cli` installed by the Debian package work. The first was not executable; the second ran the pipeline with the system Python instead of the environment created by the application.
- Pre-built installers are attached to the release; the page `installers/` lists them, with their checksums and the steps to install them.

## 1.1.0 (not released)

Changes made in response to the review of the manuscript.

### Running the software

- The default settings no longer discard the reads: no read-length filter and no identity or read-coverage filter are applied unless requested (the defaults were 200–300 bp, 95 % and 90 %).
- Read-length presets are shared by the command line and the interface, with `its-lsu` and `rrna-operon` added. The run reports the observed read lengths and warns when a window discards more than half of the reads; the interface shows the length distribution of the input folder.
- The input folder can hold one FASTQ file per sample, as well as one sub-folder per sample.
- An example dataset with its expected table is included (`--example`, **RUN EXAMPLE DATASET**); any run can be compared with an expected table (`--expected-table`).
- A conda environment file pins Python, PyQt, minimap2 and samtools. At start-up the interface explains what is missing: PyQt5, system libraries, a graphical session, minimap2, samtools or Python packages.

### Results

- Identity and read coverage are computed from the CIGAR string and the `NM` tag of each alignment instead of a samtools filter expression, which samtools 1.13 rejected. Counts are unchanged where the expression was accepted.
- Per-sample counts left in the output folder by an earlier run no longer enter the OTU table.
- With the default settings the OTU tables are identical to those of version 1.0.0 run without filters.
- The HTML report closes with a note on the reference database: reads of taxa absent from the database are attributed to the closest reference available, so a high mapping rate is not evidence that the database is complete.

### Failures and interruption

- The analysis no longer waits forever when an external tool writes many error lines. Every tool failure is reported with its exit status and its last error lines.
- Samples that fail are reported as failed. The exit status is `1` when no result is produced and `3` when some samples failed.
- Optional time limit per sample, periodic progress messages, detection of a worker process that ended unexpectedly.
- **STOP** also stops the Docker container.

### Packaging and metadata

- The Docker image builds minimap2 2.28 and samtools 1.21 from the official sources. The versions of the two tools are checked at start-up (minimum 2.24 and 1.13) and recorded with the results.
- Licence changed to GPL-3.0-only, the licence of PyQt5.
- One version number for application, installers and citation; messages say "Tabulation".
- Pre-built installers are no longer stored in the repository.
- Automated tests (`pytest tests`).

## 1.0.0 (2026-09-14)

First public release.
