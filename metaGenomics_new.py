#!/usr/bin/env python3

import sys
import argparse
import shutil
import multiprocessing
from pathlib import Path
from functools import partial
import pandas as pd
try:
    from OtuUtils import (
        create_taxonomy_map_from_fasta,
        filter_and_merge_directory,
        mapping_improved,
        tabeling_improved,
        reformat_tmp_indices,
        _create_tax_map
    )
    from ResultsReader import makeOtu_duckdb
except ImportError as e:
    print(f"ERROR: Cannot import modules. Details: {e}")
    sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Pipeline v2.6 (Merge by Directory)")
    io_group = parser.add_argument_group('Input/Output')
    io_group.add_argument(
        "-i", "--input_dir", help="Input directory (containing barcodeXX folders and unclassified).",
        required=True, type=str
    )
    io_group.add_argument(
        "-o", "--output_dir", help="Output directory.",
        required=True, type=str
    )
    io_group.add_argument(
        "-d", "--db_path", help="Path to the FASTA database.",
        required=True, type=str
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
        "-min", "--min_len", help="Minimum length.",
        type=int, default=200
    )
    filter_group.add_argument(
        "-max", "--max_len", help="Maximum length.",
        type=int, default=300
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
    map_group = parser.add_argument_group('minimap2 Parameters')
    map_group.add_argument(
        "-k", "--kmer_size", help="K-mer size.",
        type=int, default=15
    )
    map_group.add_argument(
        "-w", "--window_size", help="Minimizer window size.",
        type=int, default=10
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
    return parser.parse_args()


def check_dependencies():
    print("Checking dependencies...")
    dependencies = ["minimap2", "samtools", "awk"]
    missing_prog = [dep for dep in dependencies if not shutil.which(dep)]
    if missing_prog:
        print(f"CRITICAL ERROR: Dependencies not found: {', '.join(missing_prog)}")
        sys.exit(1)
    try:
        import Bio, pandas, openpyxl
    except ImportError as e:
        print(f"CRITICAL ERROR: Missing Python module: {e.name}")
        print("conda install biopython pandas openpyxl")
        sys.exit(1)
    print("Dependencies found.")


def run_mapping_for_one_sample(fastq_file_path: Path, **kwargs):
    sample_name = fastq_file_path.name
    debug = kwargs.get('debug', False)
    try:
        if debug:
            print(f"[Worker Map {sample_name}] Starting...")
        sam_file = mapping_improved(
            fastq_file=str(fastq_file_path), **kwargs
        )
        if sam_file:
            tabeling_improved(
                samfile=sam_file,
                threads=kwargs.get('threads', 2),
                debug=debug
            )
            if debug:
                print(f"[Worker Map {sample_name}] Completed.")
            return (sample_name, "Success")
        else:
            print(f"[Worker Map {sample_name}] ERROR: Mapping failed.")
            return (sample_name, "Failure (Mapping)")
    except Exception as e:
        print(f"[Worker Map {sample_name}] CRITICAL ERROR: {e}")
        return (sample_name, f"Failure ({e})")


def main():
    options = parse_arguments()
    print("--- Metagenomic Pipeline Start (v2.6) ---")
    check_dependencies()

    input_dir = Path(options.input_dir)
    output_dir = Path(options.output_dir)
    db_path = Path(options.db_path)

    mapping_dir = output_dir / 'mapping'
    results_dir = output_dir / 'tabeling' / 'results'

    mapping_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    analysis_name = output_dir.name
    db_name = db_path.stem

    otu_filename = f"OTU_Table_{analysis_name}_{db_name}.xlsx"
    final_output_file = output_dir / otu_filename

    print(f"OTU Table output file: {otu_filename}")

    print(f"\n--- STEP 0: Taxonomy Parsing on Database (Format: {options.format}) ---")
    tax_map_filename = f"{db_path.stem}_taxonomy_map.tsv"
    db_map_path = db_path.parent / tax_map_filename
    output_map_path = output_dir / tax_map_filename

    generated_tax_map_path = None

    if options.force_tax_map:
        print("Forcing taxonomy map creation (user request)...")
        try:
            generated_tax_map_path = _create_tax_map(fasta_path=db_path, tax_map_file=output_map_path,
                                                     format=options.format, debug=options.debug)
        except Exception as e:
            print(f"Critical error: cannot force tax map creation: {e}")
            sys.exit(1)
    if not generated_tax_map_path:
        if db_map_path.exists():
            print(f"Taxonomy map found (Cache DB): {db_map_path.name}")
            generated_tax_map_path = db_map_path
        elif output_map_path.exists():
            print(f"Taxonomy map found (Cache Output): {output_map_path.name}")
            generated_tax_map_path = output_map_path
        else:
            print("Taxonomy map not found. Creating...")
            try:
                print(f"Attempting to write to: {db_map_path}")
                generated_tax_map_path = _create_tax_map(fasta_path=db_path, tax_map_file=db_map_path,
                                                         format=options.format, debug=options.debug)
            except (OSError, PermissionError) as e:
                print(f"  [Info] DB folder not writable ({e}). Writing map to output folder.")
                try:
                    generated_tax_map_path = _create_tax_map(fasta_path=db_path, tax_map_file=output_map_path,
                                                             format=options.format, debug=options.debug)
                except Exception as e_out:
                    print(f"Critical error: cannot create tax map: {e_out}")
                    sys.exit(1)
            except Exception as e_parse:
                print(f"Critical error during parsing: {e_parse}")
                sys.exit(1)
    if not generated_tax_map_path:
        print("CRITICAL ERROR: Taxonomy map not found and not created.")
        sys.exit(1)

    print(f"\n--- STEP 1: Sequence Filter (Min: {options.min_len}, Max: {options.max_len}) ---")

    try:
        directories_to_process = [
            d for d in input_dir.iterdir()
            if d.is_dir() and not d.name.startswith('.')
        ]
    except FileNotFoundError:
        print(f"ERROR: Input directory not found: {input_dir}")
        sys.exit(1)

    if not directories_to_process:
        print(f"ERROR: No subdirectories (e.g. barcodeXX, unclassified) found in {input_dir}.")
        sys.exit(1)

    print(f"Found {len(directories_to_process)} directories to process...")

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
        debug=options.debug
    )

    filter_stats_list = []
    try:
        with multiprocessing.Pool(processes=parallel_jobs) as pool:
            results = pool.map(filter_worker_task, directories_to_process)
            filter_stats_list = [r for r in results if r is not None]
    except Exception as e:
        print(f"CRITICAL ERROR during filtering pool: {e}")
        sys.exit(1)

    print("\n--- STEP 1.5: Filtering Report ---")
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
            report_name = f"report_filtering_{options.min_len}_{options.max_len}.xlsx"
            report_path = output_dir / report_name
            report_df.to_excel(report_path, index=False, engine='openpyxl')
            print(f"Filtering report saved to: {report_path}")
        except Exception as e:
            print(f"ERROR during filtering report creation: {e}")

    print("\n--- Step 2: Mapping and Tabeling ---")

    fastq_files_to_process = [
        Path(stats['good_file_out']) for stats in filter_stats_list
        if stats and stats['good_file_out']
    ]

    if not fastq_files_to_process:
        print("ERROR: No 'good' sequences after filtering. Cannot map.")
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
        'debug': options.debug
    }

    map_worker_task = partial(
        run_mapping_for_one_sample,
        **map_worker_kwargs
    )

    try:
        with multiprocessing.Pool(processes=parallel_jobs_map) as pool:
            results = pool.map(map_worker_task, fastq_files_to_process)

        print("\nParallel work completed. Summary:")
        success_count = 0
        for sample, status in results:
            print(f"  - {sample}: {status}")
            if status == "Success":
                success_count += 1
        print(f"Completed successfully: {success_count} / {len(results)}")
    except Exception as e:
        print(f"CRITICAL ERROR during mapping pool: {e}")
        sys.exit(1)

    print("\n--- STEP 2.5: Index Reformatting ---")
    try:
        reformat_tmp_indices(
            results_dir=str(results_dir),
            format=options.format,
            debug=options.debug
        )
    except Exception as e:
        print(f"ERROR during index reformatting: {e}")

    print("\n--- Step 3: Results Aggregation ---")

    try:
        makeOtu_duckdb(
            results_dir=str(results_dir),
            output_file=str(final_output_file),
            taxonomy_file=str(generated_tax_map_path),
        )
    except Exception as e:
        print(f"ERROR during final aggregation: {e}")

    print(f"\n--- Pipeline Complete ---")
    print(f"OTU table saved to: {final_output_file.name}")


if __name__ == "__main__":
    main()
