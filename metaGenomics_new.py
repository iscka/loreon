#!/usr/bin/env python3

import sys
import os
import argparse
import shutil
import multiprocessing
import json
import time
from pathlib import Path
from functools import partial
import pandas as pd
try:
    from OtuUtils import (
        create_taxonomy_map_from_fasta,
        filter_and_merge_directory,
        mapping_improved,
        tabeling_improved,
        _parse_fasta_headers
    )
    from ResultsReader import makeOtu_duckdb
except ImportError as e:
    print(f"ERROR: Cannot import modules. Details: {e}")
    sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(description="LOREON Pipeline v1.0")
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
    filter_group.add_argument(
        "--preset",
        help="Amplicon preset: sets min/max length to reasonable defaults "
             "for common markers. Overrides -min/-max if they are at their "
             "default values. Available: its1, its2, its-full, 16s-v3v4, "
             "16s-full, 18s-v9, 18s-v4.",
        type=str,
        choices=['its1', 'its2', 'its-full', '16s-v3v4', '16s-full',
                 '18s-v9', '18s-v4'],
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
        help="Minimum percent identity for an alignment to be kept (0 = disabled). "
             "E.g. 95 means ≥95%% identity.",
        type=float, default=95.0, metavar="PCT"
    )
    quality_group.add_argument(
        "--min-ref-coverage",
        help="Minimum query coverage: percentage of the read that must be aligned "
             "(0 = disabled). E.g. 90 means ≥90%% of the read's length aligned "
             "(alen/qlen). Note: flag name kept for backward compatibility; the "
             "previous 'ref coverage' semantics were incorrect for ONT reads.",
        type=float, default=90.0, metavar="PCT"
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
        help="Delete BAM and SAM files after tabeling to free disk space.",
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
    return parser.parse_args()


def check_dependencies():
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
    print("Dependencies found.")


def run_mapping_for_one_sample(fastq_file_path: Path, **kwargs):
    sample_name = fastq_file_path.name
    debug = kwargs.get('debug', False)
    min_percent_identity = kwargs.pop('min_percent_identity', 0.0)
    min_ref_coverage = kwargs.pop('min_ref_coverage', 0.0)
    delete_bam = kwargs.pop('delete_bam', False)
    try:
        if debug:
            print(f"[Worker Map {sample_name}] Starting...")
        sam_file = mapping_improved(
            fastq_file=str(fastq_file_path), **kwargs
        )
        if sam_file:
            mapping_stats = tabeling_improved(
                samfile=sam_file,
                threads=kwargs.get('threads', 2),
                debug=debug,
                min_percent_identity=min_percent_identity,
                min_ref_coverage=min_ref_coverage,
                delete_bam=delete_bam,
            )
            if debug:
                print(f"[Worker Map {sample_name}] Completed.")
            return (sample_name, "Success", mapping_stats)
        else:
            print(f"[Worker Map {sample_name}] ERROR: Mapping failed.")
            return (sample_name, "Failure (Mapping)", None)
    except Exception as e:
        print(f"[Worker Map {sample_name}] CRITICAL ERROR: {e}")
        return (sample_name, f"Failure ({e})", None)


def main():
    os.environ.setdefault('PYTHONUNBUFFERED', '1')
    sys.stdout.reconfigure(line_buffering=True)

    if sys.platform.startswith('linux'):
        try:
            multiprocessing.set_start_method('fork', force=True)
        except RuntimeError:
            pass

    options = parse_arguments()
    print("--- Metagenomic Pipeline Start (v1.0) ---")

    _PRESETS = {
        'its1':     (150, 350),
        'its2':     (150, 400),
        'its-full': (400, 900),
        '16s-v3v4': (380, 530),
        '16s-full': (1200, 1700),
        '18s-v9':   (100, 220),
        '18s-v4':   (350, 550),
    }
    if options.preset:
        p_min, p_max = _PRESETS[options.preset]
        if options.min_len == 200 and options.max_len == 300:
            options.min_len, options.max_len = p_min, p_max
            print(f"[PRESET] Applied '{options.preset}': min_len={p_min}, max_len={p_max}")
        else:
            print(f"[PRESET] '{options.preset}' requested but min_len/max_len "
                  f"already customised — keeping user values.")

    check_dependencies()

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

    mapping_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

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
    print(f"\n--- STEP 2: Sequence Filter (Min: {options.min_len}, Max: {options.max_len}) ---")
    if options.min_percent_identity > 0:
        print(f"  Alignment quality filter: Min identity = {options.min_percent_identity}%")
    if options.min_ref_coverage > 0:
        print(f"  Alignment quality filter: Min ref coverage = {options.min_ref_coverage}%")

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
        debug=options.debug,
        gzip_intermediate=options.gzip_intermediate
    )

    filter_stats_list = []
    try:
        with multiprocessing.Pool(processes=parallel_jobs) as pool:
            results = pool.map(filter_worker_task, directories_to_process)
            filter_stats_list = [r for r in results if r is not None]
    except Exception as e:
        print(f"CRITICAL ERROR during filtering pool: {e}")
        sys.exit(1)

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

    if profiler:
        profiler.end_step('Step 3: Filter Report')

    if profiler:
        profiler.start_step('Step 4: Mapping')
    print("\n--- STEP 4: Mapping and Tabeling ---")

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
        'debug': options.debug,
        'min_percent_identity': options.min_percent_identity,
        'min_ref_coverage': options.min_ref_coverage,
        'delete_bam': options.delete_bam,
    }

    map_worker_task = partial(
        run_mapping_for_one_sample,
        **map_worker_kwargs
    )

    try:
        with multiprocessing.Pool(processes=parallel_jobs_map) as pool:
            results = []
            for result_tuple in pool.imap_unordered(map_worker_task, fastq_files_to_process):
                results.append(result_tuple)
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
        print(f"Completed successfully: {success_count} / {len(results)}")
    except Exception as e:
        print(f"CRITICAL ERROR during mapping pool: {e}")
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
            },
            'summary': {
                'total_reads_all': total_reads_all,
                'total_mapped_all': total_mapped_all,
                'total_after_quality_filter': total_filtered_all,
                'overall_mapping_rate_pct': overall_mapping_rate,
                'overall_quality_retention_pct': overall_retention,
            }
        }
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

    try:
        makeOtu_duckdb(
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

    print(f"\n--- Pipeline Complete ---")
    print(f"OTU table saved to: {final_output_file.name}")


if __name__ == "__main__":
    main()
