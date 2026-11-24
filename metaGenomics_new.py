#!/usr/bin/env python3

import sys
import os
import argparse
import re
import shutil
import subprocess
import multiprocessing
import json
import time
from pathlib import Path
from functools import partial
import pandas as pd
try:
    from OtuUtils import (
        BELOW_THRESHOLD_DIR,
        EXAMPLE_DATABASE,
        EXAMPLE_DIR,
        EXAMPLE_EXPECTED_TABLE,
        EXAMPLE_READS,
        INPUT_LAYOUT_HELP,
        LENGTH_PRESETS,
        LOREON_VERSION,
        NO_UPPER_LIMIT,
        UNMAPPED_DIR,
        compare_otu_tables,
        create_taxonomy_map_from_fasta,
        filter_and_merge_directory,
        find_input_samples,
        length_filter_advice,
        length_filter_is_active,
        mapping_improved,
        sample_name_of,
        tabeling_improved,
        _parse_fasta_headers
    )
    from ResultsReader import makeOtu_duckdb
except ImportError as e:
    print(f"ERROR: Cannot import modules. Details: {e}")
    sys.exit(1)


MIN_TOOL_VERSIONS = {'minimap2': (2, 24), 'samtools': (1, 13)}
PUBLISHED_TOOL_VERSIONS = {'minimap2': '2.28', 'samtools': '1.21'}
EXIT_PARTIAL_FAILURE = 3


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=f"LOREON Pipeline v{LOREON_VERSION}")
    io_group = parser.add_argument_group('Input/Output')
    io_group.add_argument(
        "-i", "--input_dir",
        help="Input directory: one sub-folder per sample (e.g. barcodeXX), or "
             "one FASTQ file per sample placed directly in the folder.",
        type=str
    )
    io_group.add_argument(
        "-o", "--output_dir", help="Output directory.",
        type=str
    )
    io_group.add_argument(
        "-d", "--db_path", help="Path to the FASTA database.",
        type=str
    )
    io_group.add_argument(
        "--example",
        help="Analyse the example dataset shipped with LOREON, with default "
             "settings, and compare the result with its expected table. Only "
             "-o (default: ./loreon_example_output) and the thread options are "
             "taken into account. Use it to check an installation.",
        action="store_true"
    )
    io_group.add_argument(
        "--expected-table",
        help="After the run, compare the read counts of the OTU table with "
             "this table (a .tsv written by LOREON). The run ends with status 1 "
             "if they differ.",
        type=str, default=None, metavar="TSV"
    )
    io_group.add_argument(
        "-f", "--format",
        help="Taxonomy format (for correct ID parsing). (Default: unite)",
        type=str,
        choices=['unite', 'silva', 'eukariome', 'cbs', 'none'],
        default='unite'
    )
    filter_group = parser.add_argument_group('Sequence Filter Parameters')
    filter_group.add_argument(
        "-min", "--min_len",
        help="Minimum read length in bp. Default: no lower limit, or the "
             "lower bound of --preset.",
        type=int, default=None
    )
    filter_group.add_argument(
        "-max", "--max_len",
        help="Maximum read length in bp (0 = no upper limit). Default: no "
             "upper limit, or the upper bound of --preset.",
        type=int, default=None
    )
    filter_group.add_argument(
        "--preset",
        help="Amplicon preset: sets the read-length window of a common "
             "marker. A bound given with -min or -max replaces the "
             "corresponding bound of the preset. Available: "
             + ", ".join(f"{name} ({low}-{high})" for name, (low, high) in LENGTH_PRESETS.items())
             + ". Without -min, -max and --preset no length filter is applied.",
        type=str,
        choices=list(LENGTH_PRESETS),
        default=None
    )
    perf_group = parser.add_argument_group('Performance Parameters')
    perf_group.add_argument(
        "-T", "--total_threads", help="Total CPU threads.",
        type=int, default=8
    )
    perf_group.add_argument(
        "-t", "--threads_per_job", help="Threads per single job.",
        type=int, default=2
    )
    perf_group.add_argument(
        "--sample-timeout",
        help="Maximum time, in minutes, for mapping and counting one sample "
             "(0 = no limit, default). A sample that exceeds it is stopped and "
             "reported as failed.",
        type=float, default=0.0, metavar="MIN"
    )
    perf_group.add_argument(
        "--progress-interval",
        help="During mapping, print a progress message whenever this many "
             "seconds pass without a sample finishing (0 = none). (Default: 60)",
        type=int, default=60, metavar="SEC"
    )
    map_group = parser.add_argument_group('minimap2 Parameters')
    map_group.add_argument(
        "-k", "--kmer_size", help="K-mer size.",
        type=int, default=15
    )
    map_group.add_argument(
        "-w", "--window_size", help="Minimizer window size.",
        type=int, default=10
    )
    quality_group = parser.add_argument_group('Alignment Quality Filters')
    quality_group.add_argument(
        "--min-percent-identity",
        help="Minimum percent identity for an alignment to be kept (0 = disabled, "
             "default). Identity is (L - NM) / L, where L is the read length "
             "without soft-clipped bases. E.g. 95 means ≥95%% identity.",
        type=float, default=0.0, metavar="PCT"
    )
    quality_group.add_argument(
        "--min-ref-coverage",
        help="Minimum read coverage: percentage of the read that must be aligned "
             "(0 = disabled, default). It is L / read length, where L is the read "
             "length without soft-clipped bases. E.g. 90 means ≥90%% of the read "
             "aligned. The flag name is kept for backward compatibility: the "
             "value refers to the read, not to the reference.",
        type=float, default=0.0, metavar="PCT"
    )
    other_group = parser.add_argument_group('Other')
    other_group.add_argument(
        "--debug", help="Enable debug output.",
        action="store_true"
    )
    other_group.add_argument(
        "--force-tax-map",
        help="Force recreation of the taxonomy map in the output folder.",
        action="store_true"
    )
    other_group.add_argument(
        "--profile",
        help="Enable performance profiling. Saves profile_report.json to the output folder.",
        action="store_true"
    )
    other_group.add_argument(
        "--delete-bam",
        help="Delete BAM and SAM files after tabulation to free disk space.",
        action="store_true"
    )
    other_group.add_argument(
        "--retain-unmapped",
        help="Save the reads that are not counted, one FASTQ file per sample: the "
             f"reads that do not align in {UNMAPPED_DIR}/ and, when "
             "--min-percent-identity or --min-ref-coverage is set, the aligned reads "
             f"below the thresholds in {BELOW_THRESHOLD_DIR}/.",
        action="store_true"
    )
    other_group.add_argument(
        "--gzip-intermediate",
        help="Write merged FASTQ files as .fastq.gz (halves disk usage on "
             "large ONT runs; minimap2 reads gz natively).",
        action="store_true"
    )
    other_group.add_argument(
        "--flat-sh",
        help="UNITE only: aggregate OTU rows by SH species hypothesis "
             "(collapse multiple reference accessions sharing the same SH).",
        action="store_true"
    )
    options = parser.parse_args(argv)
    if not options.example:
        missing = [flag for flag, value in (('-i/--input_dir', options.input_dir),
                                            ('-o/--output_dir', options.output_dir),
                                            ('-d/--db_path', options.db_path)) if not value]
        if missing:
            parser.error(f"the following arguments are required: {', '.join(missing)} "
                         f"(or use --example)")
    return options


def use_example_dataset(options):
    """Point the run to the example dataset and reset the settings that change its result."""
    if not (EXAMPLE_READS.is_dir() and EXAMPLE_DATABASE.is_file() and EXAMPLE_EXPECTED_TABLE.is_file()):
        print(f"CRITICAL ERROR: the example dataset was not found in {EXAMPLE_DIR}.")
        sys.exit(1)
    ignored = [flag for flag, given in (
        ('-i', options.input_dir is not None),
        ('-d', options.db_path is not None),
        ('-f', options.format != 'unite'),
        ('--expected-table', options.expected_table is not None),
        ('-min', options.min_len is not None),
        ('-max', options.max_len is not None),
        ('--preset', options.preset is not None),
        ('--min-percent-identity', options.min_percent_identity != 0.0),
        ('--min-ref-coverage', options.min_ref_coverage != 0.0),
        ('-k', options.kmer_size != 15),
        ('-w', options.window_size != 10),
        ('--flat-sh', options.flat_sh),
    ) if given]
    options.input_dir = str(EXAMPLE_READS)
    options.db_path = str(EXAMPLE_DATABASE)
    options.expected_table = str(EXAMPLE_EXPECTED_TABLE)
    if not options.output_dir:
        options.output_dir = str(Path.cwd() / 'loreon_example_output')
    options.format = 'unite'
    options.min_len = options.max_len = options.preset = None
    options.min_percent_identity = options.min_ref_coverage = 0.0
    options.kmer_size, options.window_size = 15, 10
    options.flat_sh = False
    options.force_tax_map = True
    print("Example dataset: three mock-community samples of 500 reads each and a subset of the "
          "UNITE database, analysed with the default settings.")
    print(f"  Output folder: {options.output_dir}")
    if ignored:
        print(f"  Options ignored with --example: {', '.join(ignored)}.")


def check_against_expected_table(expected_table, observed_table):
    """Print how the OTU table compares with an expected one. Returns True when the counts agree."""
    if not observed_table.exists():
        print(f"\nResult check FAILED: no OTU table was written to {observed_table.name}.")
        return False
    try:
        differences = compare_otu_tables(expected_table, observed_table)
    except (OSError, ValueError, IndexError) as e:
        print(f"\nResult check FAILED: the tables could not be compared ({type(e).__name__}: {e}).")
        return False
    if differences:
        print(f"\nResult check FAILED: the OTU table differs from the expected table "
              f"({expected_table}):")
        for line in differences:
            print(f"  {line}")
        print("  Different versions of minimap2 can assign a few reads differently: the versions")
        print("  used in this run are listed at the top of the log.")
        return False
    print(f"\nResult check: the OTU table matches the expected table ({expected_table.name}).")
    return True


def resolve_length_window(min_len, max_len, preset):
    """Return the (min, max) read-length window to apply.

    A bound given explicitly wins over the preset; a bound that is neither
    given nor set by a preset is 0, which means no limit on that side.
    """
    preset_min, preset_max = LENGTH_PRESETS[preset] if preset else (0, 0)
    return (
        preset_min if min_len is None else max(0, min_len),
        preset_max if max_len is None else max(0, max_len),
    )


def parse_tool_version(text):
    """Extract (major, minor) from the output of 'minimap2 --version' or 'samtools --version'."""
    first_line = text.strip().splitlines()[0] if text and text.strip() else ''
    match = re.search(r'(\d+)\.(\d+)', first_line)
    if not match:
        return None
    return int(match.group(1)), int(match.group(2))


def _tool_version_string(tool):
    try:
        result = subprocess.run([tool, '--version'], capture_output=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return ''
    output = (result.stdout or result.stderr or b'').decode('utf-8', errors='replace')
    lines = output.strip().splitlines()
    version = lines[0].strip() if lines else ''
    if version.lower().startswith(tool + ' '):
        version = version[len(tool) + 1:].strip()
    return version


def check_dependencies():
    """Check that minimap2 and samtools are installed and recent enough.

    Returns the version string of each tool. Stops the run when a tool is
    missing or older than the oldest version LOREON was tested with.
    """
    print("Checking dependencies...")
    dependencies = ["minimap2", "samtools"]
    missing_prog = [dep for dep in dependencies if not shutil.which(dep)]
    if missing_prog:
        print(f"CRITICAL ERROR: Dependencies not found: {', '.join(missing_prog)}")
        sys.exit(1)
    try:
        import pandas
    except ImportError as e:
        print(f"CRITICAL ERROR: Missing Python module: {e.name}")
        print("conda install pandas")
        sys.exit(1)

    versions = {}
    outdated = []
    for tool in dependencies:
        version = _tool_version_string(tool)
        versions[tool] = version
        parsed = parse_tool_version(version)
        if parsed is None:
            print(f"WARNING: could not read the version of {tool}; "
                  f"version {'.'.join(map(str, MIN_TOOL_VERSIONS[tool]))} or newer is required.")
        elif parsed < MIN_TOOL_VERSIONS[tool]:
            outdated.append(tool)

    if outdated:
        for tool in outdated:
            required = '.'.join(map(str, MIN_TOOL_VERSIONS[tool]))
            print(f"CRITICAL ERROR: {tool} {versions[tool]} is too old: "
                  f"LOREON requires {tool} {required} or newer.")
        print("  The published analyses used minimap2 "
              f"{PUBLISHED_TOOL_VERSIONS['minimap2']} and samtools {PUBLISHED_TOOL_VERSIONS['samtools']}.")
        print("  Install them, for example with "
              f"'conda install -c bioconda minimap2={PUBLISHED_TOOL_VERSIONS['minimap2']} "
              f"samtools={PUBLISHED_TOOL_VERSIONS['samtools']}', or run LOREON in Docker mode.")
        sys.exit(1)

    print("Dependencies found: " + ", ".join(f"{tool} {versions[tool]}" for tool in dependencies) + ".")
    return versions


def run_mapping_for_one_sample(fastq_file_path: Path, **kwargs):
    sample_name = fastq_file_path.name
    debug = kwargs.get('debug', False)
    min_percent_identity = kwargs.pop('min_percent_identity', 0.0)
    min_ref_coverage = kwargs.pop('min_ref_coverage', 0.0)
    delete_bam = kwargs.pop('delete_bam', False)
    retain_unmapped = kwargs.pop('retain_unmapped', False)
    gzip_output = kwargs.pop('gzip_output', False)
    sample_timeout_s = kwargs.pop('sample_timeout_s', None)
    started = time.monotonic()
    try:
        if debug:
            print(f"[Worker Map {sample_name}] Starting...")
        sam_file = mapping_improved(
            fastq_file=str(fastq_file_path), timeout_s=sample_timeout_s, **kwargs
        )
        if not sam_file:
            print(f"[Worker Map {sample_name}] ERROR: Mapping failed.", flush=True)
            return (sample_name, "Failure (Mapping)", None)
        remaining_s = None
        if sample_timeout_s:
            remaining_s = max(1.0, sample_timeout_s - (time.monotonic() - started))
        mapping_stats = tabeling_improved(
            samfile=sam_file,
            threads=kwargs.get('threads', 2),
            debug=debug,
            min_percent_identity=min_percent_identity,
            min_ref_coverage=min_ref_coverage,
            delete_bam=delete_bam,
            timeout_s=remaining_s,
            retain_unmapped=retain_unmapped,
            gzip_output=gzip_output,
        )
        if mapping_stats is None:
            print(f"[Worker Map {sample_name}] ERROR: Counting the alignments failed.", flush=True)
            return (sample_name, "Failure (Tabulation)", None)
        if debug:
            print(f"[Worker Map {sample_name}] Completed.")
        return (sample_name, "Success", mapping_stats)
    except (Exception, SystemExit) as e:
        print(f"[Worker Map {sample_name}] CRITICAL ERROR: {type(e).__name__}: {e}", flush=True)
        return (sample_name, f"Failure ({type(e).__name__}: {e})", None)


def main():
    os.environ.setdefault('PYTHONUNBUFFERED', '1')
    sys.stdout.reconfigure(line_buffering=True)

    if sys.platform.startswith('linux'):
        try:
            multiprocessing.set_start_method('fork', force=True)
        except RuntimeError:
            pass

    options = parse_arguments()
    print(f"--- Metagenomic Pipeline Start (v{LOREON_VERSION}) ---")
    if options.example:
        use_example_dataset(options)
    if options.expected_table and not Path(options.expected_table).is_file():
        print(f"CRITICAL ERROR: the expected table was not found: {options.expected_table}")
        sys.exit(1)

    requested_min, requested_max = options.min_len, options.max_len
    options.min_len, options.max_len = resolve_length_window(
        requested_min, requested_max, options.preset
    )
    if options.preset:
        overridden = [name for name, value in (('-min', requested_min), ('-max', requested_max))
                      if value is not None]
        print(f"[PRESET] Applied '{options.preset}': min_len={options.min_len}, "
              f"max_len={options.max_len}"
              + (f" ({' and '.join(overridden)} given explicitly)" if overridden else ""))
    if 0 < options.max_len < options.min_len:
        print(f"CRITICAL ERROR: the maximum read length ({options.max_len}) is smaller than "
              f"the minimum ({options.min_len}).")
        sys.exit(1)

    tool_versions = check_dependencies()

    profiler = None
    if options.profile:
        try:
            from pipeline_profiler import PipelineProfiler
            profiler = PipelineProfiler()
            print("[PROFILER] Performance profiling enabled.")
        except ImportError:
            print("[PROFILER] WARNING: pipeline_profiler.py not found. Profiling disabled.")

    input_dir = Path(options.input_dir)
    output_dir = Path(options.output_dir)
    db_path = Path(options.db_path)

    mapping_dir = output_dir / 'mapping'
    results_dir = output_dir / 'tabeling' / 'results'

    try:
        mapping_dir.mkdir(parents=True, exist_ok=True)
        results_dir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        print(f"CRITICAL ERROR: the output folder cannot be created or written: {output_dir} ({e})")
        sys.exit(1)

    analysis_name = output_dir.name
    db_name = db_path.stem

    otu_filename = f"OTU_Table_{analysis_name}_{db_name}.xlsx"
    final_output_file = output_dir / otu_filename

    print(f"OTU Table output file: {otu_filename}")

    if profiler:
        profiler.start_step('Step 1: Taxonomy')
    print(f"\n--- STEP 1: Taxonomy Parsing on Database (Format: {options.format}) ---")

    generated_tax_map_path = None
    tax_entries = 0

    if options.force_tax_map:
        print("Forcing taxonomy map creation (user request)...")
        tax_map_filename = f"{db_path.stem}_taxonomy_map.tsv"
        output_map_path = output_dir / tax_map_filename
        try:
            generated_tax_map_path, tax_entries = (
                create_taxonomy_map_from_fasta(
                    fasta_file=str(db_path),
                    output_dir=output_dir,
                    format=options.format,
                    force_create=True,
                    debug=options.debug
                )
            )
        except Exception as e:
            print(f"Critical error: cannot force tax map creation: {e}")
            sys.exit(1)
    if not generated_tax_map_path:
        generated_tax_map_path, tax_entries = (
            create_taxonomy_map_from_fasta(
                fasta_file=str(db_path),
                output_dir=output_dir,
                format=options.format,
                force_create=False,
                debug=options.debug
            )
        )

    if not generated_tax_map_path:
        print("CRITICAL ERROR: Taxonomy map not found and not created.")
        sys.exit(1)

    if profiler:
        profiler.end_step('Step 1: Taxonomy', taxonomy_entries=tax_entries)

    if profiler:
        profiler.start_step('Step 2: Filter')
    max_len_text = options.max_len if 0 < options.max_len < NO_UPPER_LIMIT else "no limit"
    print(f"\n--- STEP 2: Sequence Filter (Min: {options.min_len}, Max: {max_len_text}) ---")
    if not length_filter_is_active(options.min_len, options.max_len):
        print("  No read-length filter: reads of any length are used.")
    if options.min_percent_identity > 0:
        print(f"  Alignment quality filter: Min identity = {options.min_percent_identity}%")
    if options.min_ref_coverage > 0:
        print(f"  Alignment quality filter: Min read coverage = {options.min_ref_coverage}%")

    try:
        found_input = find_input_samples(input_dir)
    except (FileNotFoundError, NotADirectoryError):
        print(f"ERROR: Input directory not found: {input_dir}")
        sys.exit(1)

    directories_to_process = found_input['sources']
    if not directories_to_process:
        print(f"ERROR: No FASTQ files found in {input_dir}.")
        for line in INPUT_LAYOUT_HELP:
            print(line)
        sys.exit(1)

    if found_input['layout'] == 'files':
        names = {}
        for source in directories_to_process:
            names.setdefault(sample_name_of(source), []).append(source.name)
        clashes = {name: files for name, files in names.items() if len(files) > 1}
        if clashes:
            print("ERROR: these files would get the same sample name:")
            for name, files in sorted(clashes.items()):
                print(f"    {name}: {', '.join(files)}")
            print("  Keep one file per sample, or put the files of each sample in a sub-folder.")
            sys.exit(1)
        print(f"Found {len(directories_to_process)} FASTQ files in the input folder: "
              f"each file is one sample.")
    else:
        print(f"Found {len(directories_to_process)} directories to process...")
        ignored = found_input['ignored']
        if ignored:
            shown = ', '.join(path.name for path in ignored[:10])
            more = f" and {len(ignored) - 10} more" if len(ignored) > 10 else ""
            print(f"WARNING: {len(ignored)} FASTQ files placed directly in the input folder are "
                  f"ignored, because the folder also contains sub-folders with reads: {shown}{more}.")
            print("  Move them into a sub-folder to include them as a sample.")

    total_threads = max(1, options.total_threads)
    parallel_jobs = min(
        total_threads,
        len(directories_to_process),
        multiprocessing.cpu_count()
    )
    print(f"Job pool of {parallel_jobs} workers for filtering...")

    filter_worker_task = partial(
        filter_and_merge_directory,
        min_len=options.min_len,
        max_len=options.max_len,
        out_dir=output_dir,
        debug=options.debug,
        gzip_intermediate=options.gzip_intermediate
    )

    filter_stats_list = []
    failed_samples = []
    try:
        with multiprocessing.Pool(processes=parallel_jobs) as pool:
            results = pool.map(filter_worker_task, directories_to_process)
            filter_stats_list = [r for r in results if r is not None and not r.get('failed')]
            failed_samples = [
                (r['barcode'], f"Failure (Reading the input: {r.get('error', 'unknown error')})")
                for r in results if r is not None and r.get('failed')
            ]
    except Exception as e:
        print(f"CRITICAL ERROR during filtering pool: {e}")
        sys.exit(1)

    for sample, status in failed_samples:
        print(f"ERROR: the reads of '{sample}' could not be read; the sample is left out. {status}")

    if profiler:
        total_good = sum(s.get('good_count', 0) for s in filter_stats_list)
        total_bad = sum(s.get('bad_count', 0) for s in filter_stats_list)
        total_seq = sum(s.get('total_count', 0) for s in filter_stats_list)
        wall_s = time.perf_counter() - profiler._start_wall
        throughput = round(total_seq / wall_s) if wall_s > 0 else 0
        profiler.end_step('Step 2: Filter',
                          directories=len(directories_to_process),
                          sequences_total=total_seq,
                          sequences_good=total_good,
                          sequences_bad=total_bad,
                          throughput_seq_s=throughput)

    if profiler:
        profiler.start_step('Step 3: Filter Report')
    print("\n--- STEP 3: Filtering Report ---")
    if not filter_stats_list:
        print("No stats available. Skipping report generation.")
    else:
        try:
            report_df = pd.DataFrame(filter_stats_list)
            report_df = report_df.rename(columns={
                'mean_len': 'mean length', 'median_len': 'median length',
                'bad_count': 'bad count', 'good_count': 'good count',
                'total_count': 'total count'
            })
            final_report_cols = [
                'barcode', 'mean length', 'median length',
                'bad count', 'good count', 'total count'
            ]
            report_df = report_df[final_report_cols]
            report_name = f"report_filtering_{options.min_len}_{options.max_len}.tsv"
            report_path = output_dir / report_name
            report_df.to_csv(report_path, index=False, sep='\t')
            print(f"Filtering report saved to: {report_path}")
        except Exception as e:
            print(f"ERROR during filtering report creation: {e}")

    total_reads_in = sum(s.get('total_count', 0) for s in filter_stats_list)
    total_reads_kept = sum(s.get('good_count', 0) for s in filter_stats_list)
    sampled_lengths = []
    sampled_weights = []
    for s in filter_stats_list:
        sample = s.get('length_sample') or []
        if sample:
            sampled_lengths.extend(sample)
            sampled_weights.extend([s.get('total_count', 0) / len(sample)] * len(sample))
    length_advice = length_filter_advice(
        sampled_lengths, sampled_weights, total_reads_in, total_reads_kept,
        options.min_len, options.max_len
    )
    if length_advice:
        print()
        for line in length_advice:
            print(line)

    if profiler:
        profiler.end_step('Step 3: Filter Report')

    if profiler:
        profiler.start_step('Step 4: Mapping')
    print("\n--- STEP 4: Mapping and Tabulation ---")

    fastq_files_to_process = [
        Path(stats['good_file_out']) for stats in filter_stats_list
        if stats and stats['good_file_out']
    ]

    stale_counts = [f for f in results_dir.glob('*.tmp.txt')]
    if stale_counts:
        removed = 0
        for stale in stale_counts:
            try:
                stale.unlink()
                removed += 1
            except OSError as e:
                print(f"CRITICAL ERROR: counts left by an earlier run cannot be removed "
                      f"({stale.name}: {e}). Use an empty output folder.")
                sys.exit(1)
        print(f"Removed {removed} per-sample count files left in the output folder by an earlier run.")

    if not fastq_files_to_process:
        if total_reads_in == 0:
            print(f"ERROR: No reads found in {input_dir}. Cannot map.")
            for line in INPUT_LAYOUT_HELP:
                print(line)
        else:
            print(f"ERROR: No reads left after the read-length filter "
                  f"(Min: {options.min_len}, Max: {max_len_text}). Cannot map.")
            print("  Change the length window as suggested above and run again.")
        sys.exit(1)

    print(f"{len(fastq_files_to_process)} 'good' files (merged) ready for mapping.")

    threads_per_job_map = max(1, options.threads_per_job)
    total_threads_map = max(threads_per_job_map, options.total_threads)
    parallel_jobs_map = min(
        total_threads_map // threads_per_job_map,
        len(fastq_files_to_process),
        multiprocessing.cpu_count()
    )

    print(f"Starting mapping with {parallel_jobs_map} parallel jobs (each with {threads_per_job_map} threads)...")

    map_worker_kwargs = {
        'db_path': str(db_path),
        'threads': threads_per_job_map,
        'output_dir': str(mapping_dir),
        'kmer_size': options.kmer_size,
        'window_size': options.window_size,
        'debug': options.debug,
        'min_percent_identity': options.min_percent_identity,
        'min_ref_coverage': options.min_ref_coverage,
        'delete_bam': options.delete_bam,
        'retain_unmapped': options.retain_unmapped,
        'gzip_output': options.gzip_intermediate,
        'sample_timeout_s': options.sample_timeout * 60 if options.sample_timeout > 0 else None,
    }

    map_worker_task = partial(
        run_mapping_for_one_sample,
        **map_worker_kwargs
    )

    if options.sample_timeout > 0:
        print(f"Time limit per sample: {options.sample_timeout:g} min.")

    samples_to_map = len(fastq_files_to_process)
    total_samples = samples_to_map + len(failed_samples)
    progress_interval = options.progress_interval if options.progress_interval > 0 else None
    step_started = time.monotonic()
    last_progress = step_started
    try:
        with multiprocessing.Pool(processes=parallel_jobs_map) as pool:
            results = []
            worker_pids = {process.pid for process in multiprocessing.active_children()}
            pending = pool.imap_unordered(map_worker_task, fastq_files_to_process)
            while len(results) < samples_to_map:
                try:
                    result_tuple = pending.next(timeout=1.0)
                except multiprocessing.TimeoutError:
                    alive_pids = {process.pid for process in multiprocessing.active_children()}
                    if not worker_pids <= alive_pids:
                        completed_names = {r[0] for r in results}
                        print("ERROR: a worker process ended unexpectedly; the computer or the "
                              "Docker container may have run out of memory.", flush=True)
                        for fastq in fastq_files_to_process:
                            if fastq.name not in completed_names:
                                results.append((fastq.name,
                                                "Failure (a worker process ended unexpectedly)", None))
                                print(f"  - {fastq.name}: not completed", flush=True)
                        break
                    now = time.monotonic()
                    if progress_interval and now - last_progress >= progress_interval:
                        last_progress = now
                        elapsed = int(now - step_started)
                        print(f"  Step 4 progress: {len(results)}/{samples_to_map} samples completed "
                              f"(elapsed {elapsed // 3600:02d}:{elapsed % 3600 // 60:02d}:{elapsed % 60:02d})",
                              flush=True)
                    continue
                results.append(result_tuple)
                last_progress = time.monotonic()
                sample, status = result_tuple[0], result_tuple[1]
                print(f"  - {sample}: {status}", flush=True)

        print("\nParallel work completed. Summary:")
        success_count = 0
        mapping_stats_list = []
        for result_tuple in results:
            sample, status = result_tuple[0], result_tuple[1]
            stats = result_tuple[2] if len(result_tuple) > 2 else None
            if status == "Success":
                success_count += 1
                if stats:
                    mapping_stats_list.append(stats)
            else:
                failed_samples.append((sample, status))
        print(f"Completed successfully: {success_count} / {len(results)}")
    except Exception as e:
        print(f"CRITICAL ERROR during mapping pool: {e}")
        sys.exit(1)

    if failed_samples:
        print("\n" + "=" * 70)
        print(f"  ERROR: {len(failed_samples)} of {total_samples} samples failed "
              f"and are missing from the results:")
        for sample, status in failed_samples:
            print(f"    - {sample}: {status}")
        print("  The cause of each failure is reported above, next to the sample name.")
        print("=" * 70 + "\n")
    if success_count == 0:
        print("ERROR: mapping failed for every sample: no result can be produced.")
        sys.exit(1)

    if mapping_stats_list:
        total_reads_all = sum(s.get('total_reads', 0) for s in mapping_stats_list)
        total_mapped_all = sum(s.get('mapped_unfiltered', 0) for s in mapping_stats_list)
        total_filtered_all = sum(s.get('mapped_filtered', 0) for s in mapping_stats_list)
        overall_mapping_rate = round(total_mapped_all / total_reads_all * 100, 2) if total_reads_all > 0 else 0.0
        overall_retention = round(total_filtered_all / total_mapped_all * 100, 2) if total_mapped_all > 0 else 100.0

        if total_filtered_all == 0 and total_mapped_all > 0:
            print("\n" + "=" * 70)
            print("  WARNING: The alignment quality filter removed ALL mapped reads!")
            print(f"  → {total_mapped_all:,} reads mapped, 0 survived the filter.")
            print(f"  → min_percent_identity = {options.min_percent_identity}%")
            print(f"  → min_ref_coverage     = {options.min_ref_coverage}%")
            print("  The OTU table will be empty. Consider lowering or")
            print("  disabling (set to 0) the quality filter thresholds.")
            print("=" * 70 + "\n")

        mapping_json = {
            'samples': mapping_stats_list,
            'parameters': {
                'min_percent_identity': options.min_percent_identity,
                'min_ref_coverage': options.min_ref_coverage,
                'retain_unmapped': options.retain_unmapped,
            },
            'tools': tool_versions,
            'summary': {
                'total_reads_all': total_reads_all,
                'total_mapped_all': total_mapped_all,
                'total_after_quality_filter': total_filtered_all,
                'overall_mapping_rate_pct': overall_mapping_rate,
                'overall_quality_retention_pct': overall_retention,
            }
        }
        if options.retain_unmapped:
            unmapped_saved_all = sum(s.get('unmapped_saved', 0) for s in mapping_stats_list)
            below_saved_all = sum(s.get('below_threshold_saved', 0) for s in mapping_stats_list)
            mapping_json['summary']['unmapped_saved_all'] = unmapped_saved_all
            mapping_json['summary']['below_threshold_saved_all'] = below_saved_all
            print(f"Reads that do not align, saved in {UNMAPPED_DIR}/: {unmapped_saved_all:,}")
            if options.min_percent_identity > 0 or options.min_ref_coverage > 0:
                print(f"Aligned reads below the thresholds, saved in {BELOW_THRESHOLD_DIR}/: "
                      f"{below_saved_all:,}")
        mapping_stats_path = output_dir / 'mapping_stats.json'
        try:
            with open(mapping_stats_path, 'w', encoding='utf-8') as f:
                json.dump(mapping_json, f, indent=2)
            print(f"Mapping stats saved to: {mapping_stats_path.name}")
        except Exception as e:
            print(f"WARNING: Could not save mapping_stats.json: {e}")

    if profiler:
        success_count_p = sum(1 for r in results if r[1] == "Success")
        profiler.end_step('Step 4: Mapping',
                          files_to_map=len(fastq_files_to_process),
                          success=success_count_p,
                          failed=len(fastq_files_to_process) - success_count_p)

    if profiler:
        profiler.start_step('Step 5: OTU Aggregation')
    print("\n--- STEP 5: Results Aggregation ---")

    for earlier_table in (final_output_file, final_output_file.with_suffix('.tsv')):
        try:
            earlier_table.unlink()
        except FileNotFoundError:
            pass
        except OSError as e:
            print(f"WARNING: the table of an earlier run could not be removed "
                  f"({earlier_table.name}: {e}).")

    table_written = False
    try:
        table_written = makeOtu_duckdb(
            results_dir=str(results_dir),
            output_file=str(final_output_file),
            taxonomy_file=str(generated_tax_map_path),
            db_format=options.format,
            flat_sh=options.flat_sh,
        )
    except Exception as e:
        print(f"ERROR during final aggregation: {e}")

    if profiler:
        profiler.end_step('Step 5: OTU Aggregation')
        profiler.save(output_dir, db_name=db_name,
                      total_threads=options.total_threads,
                      threads_per_job=options.threads_per_job)

    if not table_written:
        print("\nERROR: the OTU table could not be written. See the messages above.")
        sys.exit(1)

    if failed_samples:
        print(f"\n--- Pipeline Complete with errors: {len(failed_samples)} of "
              f"{total_samples} samples failed ---")
    else:
        print(f"\n--- Pipeline Complete ---")
    print(f"OTU table saved to: {final_output_file.name}")

    result_matches = True
    if options.expected_table:
        result_matches = check_against_expected_table(
            Path(options.expected_table), final_output_file.with_suffix('.tsv')
        )
        if result_matches and options.example and not failed_samples:
            print("The installation works.")

    if failed_samples:
        sys.exit(EXIT_PARTIAL_FAILURE)
    if not result_matches:
        sys.exit(1)


if __name__ == "__main__":
    main()
