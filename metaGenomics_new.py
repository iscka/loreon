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
    parser = argparse.ArgumentParser(description="Pipeline v2.7 (Optimized)")
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
    taxq_group = parser.add_argument_group('Taxonomic Quality / EM')
    taxq_group.add_argument(
        "--tax-quality", dest="tax_quality", action="store_true",
        help="Enable taxonomic-quality metrics via a read-level BAM pass "
             "(re-enables minimap2 secondary alignments; adds per-OTU "
             "resolution profile + confusion network). Opt-in; inflates the "
             "BAM ~20-30%%."
    )
    taxq_group.add_argument(
        "--em", dest="em", action="store_true",
        help="Enable fractional (EM) OTU assignment. Writes a second "
             "OTU_Table_..._EM table alongside the primary-count table. "
             "Implies the read-level pass (and secondary alignments)."
    )
    taxq_group.add_argument(
        "--consensus-tau", dest="consensus_tau", type=float, default=0.9,
        help="Stringency of the weighted taxonomic consensus (default 0.9; "
             "1.0 = strict LCA)."
    )
    taxq_group.add_argument(
        "--em-rank", dest="em_rank", choices=['reference', 'species', 'genus'],
        default='reference',
        help="Abundance unit for the EM table. Only 'reference' is implemented "
             "(the EM's latent unit); the EM table carries full taxonomy, so "
             "genus/species abundances are obtained by aggregating it per rank. "
             "Passing another value is refused rather than silently ignored."
    )
    taxq_group.add_argument(
        "--sec-N", dest="sec_N", type=int, default=20,
        help="minimap2 -N (max secondary) in tax-quality mode (default 20)."
    )
    taxq_group.add_argument(
        "--sec-p", dest="sec_p", type=float, default=0.8,
        help="minimap2 -p (secondary-to-primary score ratio) in tax-quality "
             "mode (default 0.8 = minimap2's own default). Lowering it exposes "
             "a longer ambiguity tail but PERTURBS primary selection: measured "
             "on a mock ITS community vs --secondary=no, -p 0.8 keeps 94.2%% of "
             "primaries identical while -p 0.5 keeps only 86.2%%, at almost no "
             "cost in secondaries retained (10.0 vs 19.3 per read, both with "
             "~98%% of reads carrying at least one)."
    )
    taxq_group.add_argument(
        "--em-max-iter", dest="em_max_iter", type=int, default=200,
        help="EM maximum iterations (default 200)."
    )
    taxq_group.add_argument(
        "--em-tol", dest="em_tol", type=float, default=1e-6,
        help="EM convergence tolerance on log-likelihood delta (default 1e-6)."
    )
    taxq_group.add_argument(
        "--em-prior", dest="em_prior", type=float, default=0.0,
        help="Weak Dirichlet prior for the EM (default 0 = off; small values "
             "stabilise rare references)."
    )
    taxq_group.add_argument(
        "--score-temperature", dest="score_temperature", type=float, default=6.0,
        help="Softmax temperature over minimap2 AS, in raw score units "
             "(default 6.0 = one mismatch). Lower is more decisive; a "
             "likelihood-calibrated value is nearer 2. Shared by the consensus "
             "and the EM."
    )
    taxq_group.add_argument(
        "--em-init", dest="em_init", choices=['uniform', 'primary'],
        default='uniform',
        help="EM initialisation (default: uniform, as kallisto/RSEM/Bracken). "
             "'primary' warm-starts from primary counts but cannot ever assign "
             "mass to a reference that never wins a primary alignment."
    )
    taxq_group.add_argument(
        "--em-model", dest="em_model", choices=['heuristic', 'generative'],
        default='heuristic',
        help="Which model computes the EM/quality output. 'heuristic' (DEFAULT, "
             "unchanged behaviour) uses softmax weights over minimap2 AS. "
             "'generative' uses the unified probabilistic model of the formal "
             "specification: indel-aware error model estimated from the data, "
             "chimeras inside the likelihood, true posterior responsibilities. "
             "Entirely opt-in — nothing changes unless you ask for it."
    )
    taxq_group.add_argument(
        "--em-bootstrap", dest="em_bootstrap", type=int, default=0, metavar="B",
        help="Generative model only: bootstrap replicates for credible intervals "
             "on abundances and diversity (0 = off; 200-1000 recommended). "
             "Writes <sample>.ci.tsv."
    )
    taxq_group.add_argument(
        "--em-ci-alpha", dest="em_ci_alpha", type=float, default=0.05,
        help="Generative model only: 1-alpha credible level (default 0.05 = 95%%)."
    )
    taxq_group.add_argument(
        "--em-abundance", dest="em_abundance", choices=['all', 'clean'],
        default='all',
        help="Generative model only: 'clean' excludes chimeric contributions "
             "from the abundance estimate."
    )
    taxq_group.add_argument(
        "--em-no-estimate-error", dest="em_no_estimate_error", action="store_true",
        help="Generative model only: keep the error model fixed instead of "
             "estimating it from the data."
    )
    taxq_group.add_argument(
        "--chimera-refined", dest="chimera_refined", action="store_true",
        help="Enable the refined cross-taxon chimera tier (compares primary vs "
             "supplementary taxa). Tier-1 split-rate is always computed."
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
    # --- Parallel execution setup ---
    os.environ.setdefault('PYTHONUNBUFFERED', '1')
    sys.stdout.reconfigure(line_buffering=True)

    if sys.platform.startswith('linux'):
        try:
            multiprocessing.set_start_method('fork', force=True)
        except RuntimeError:
            pass

    options = parse_arguments()
    print("--- Metagenomic Pipeline Start (v2.7 Optimized) ---")

    # Apply amplicon preset — only when the user left min/max at defaults,
    # so that explicit -min/-max always win.
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

    if options.em and options.em_rank != 'reference':
        print(f"ERROR: --em-rank '{options.em_rank}' is not implemented. The EM "
              f"runs at the 'reference' latent unit; aggregate the EM OTU table "
              f"per rank instead (it carries the full taxonomy).")
        sys.exit(2)

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

    # =======================================================================
    # STEP 1: Taxonomy Parsing
    # =======================================================================
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

    # =======================================================================
    # STEP 2: Sequence Filter
    # =======================================================================
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

    # =======================================================================
    # STEP 3: Filter Report  (OPT-7: TSV instead of XLSX)
    # =======================================================================
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
            # OPT-7: TSV is ~100x faster than XLSX round-trip
            report_name = f"report_filtering_{options.min_len}_{options.max_len}.tsv"
            report_path = output_dir / report_name
            report_df.to_csv(report_path, index=False, sep='\t')
            print(f"Filtering report saved to: {report_path}")
        except Exception as e:
            print(f"ERROR during filtering report creation: {e}")

    if profiler:
        profiler.end_step('Step 3: Filter Report')

    # =======================================================================
    # STEP 4: Mapping and Tabeling  (OPT-6: imap_unordered for overlap)
    # =======================================================================
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

    # Tax-quality / EM mode re-enables secondary alignments in the full BAM and
    # requires the full BAM to survive tabeling for the read-level pass.
    # --chimera-refined also needs the read-level pass (plan §3.4): previously it
    # was accepted, echoed into mapping_stats.json as enabled, and silently did
    # nothing while the report still drew a "0.0% chimeras" card.
    tax_mode = options.tax_quality or options.em or options.chimera_refined
    worker_delete_bam = options.delete_bam and not tax_mode
    if tax_mode:
        print(f"[TaxQuality] Read-level pass enabled "
              f"(tax_quality={options.tax_quality}, em={options.em}, "
              f"secondary -N {options.sec_N} -p {options.sec_p}).")

    map_worker_kwargs = {
        'db_path': str(db_path),
        'threads': threads_per_job_map,
        'output_dir': str(mapping_dir),
        'kmer_size': options.kmer_size,
        'window_size': options.window_size,
        'debug': options.debug,
        'min_percent_identity': options.min_percent_identity,
        'min_ref_coverage': options.min_ref_coverage,
        'delete_bam': worker_delete_bam,
        'tax_quality_mode': tax_mode,
        'sec_N': options.sec_N,
        'sec_p': options.sec_p,
    }

    map_worker_task = partial(
        run_mapping_for_one_sample,
        **map_worker_kwargs
    )

    try:
        # OPT-6: imap_unordered allows results to stream as they finish
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

    # Close Step 4 here: Step 4b is a sibling, not a child.  Opening 4b inside
    # Step 4's window made the profiler count its time twice, so the per-step
    # column summed past the true elapsed time.
    if profiler:
        success_count_p = sum(1 for r in results if r[1] == "Success")
        profiler.end_step('Step 4: Mapping',
                          files_to_map=len(fastq_files_to_process),
                          success=success_count_p,
                          failed=len(fastq_files_to_process) - success_count_p)

    # =======================================================================
    # STEP 4b: Taxonomic Quality / Chimera / EM read-level pass (opt-in)
    # =======================================================================
    quality_files = []
    em_files_present = False
    if tax_mode and mapping_stats_list:
        if profiler:
            profiler.start_step('Step 4b: TaxQuality/EM')
        model_name = ('generative (unified probabilistic model)'
                      if options.em_model == 'generative'
                      else 'heuristic (softmax over AS)')
        print(f"\n--- STEP 4b: Taxonomic Quality / EM (read-level pass) ---")
        print(f"    model: {model_name}")
        if options.em_model == 'generative' and options.em_bootstrap:
            print(f"    bootstrap: B={options.em_bootstrap}, "
                  f"credible level {(1-options.em_ci_alpha)*100:.0f}%")
        try:
            import TaxQuality
        except ImportError as e:
            print(f"WARNING: TaxQuality module unavailable ({e}); skipping.")
            TaxQuality = None
        TaxQualityModel = None
        if options.em_model == 'generative':
            try:
                import TaxQualityModel
            except ImportError as e:
                print(f"WARNING: TaxQualityModel unavailable ({e}); "
                      f"falling back to the heuristic model.")
                options.em_model = 'heuristic'

        if TaxQuality is not None:
            for stats in mapping_stats_list:
                sample = stats.get('sample')
                if not sample:
                    continue
                bam_path = mapping_dir / f"{sample}_sorted_full.bam"
                if not bam_path.exists():
                    print(f"[TaxQuality] {sample}: full BAM missing ({bam_path.name}); skipping.")
                    continue
                q_out = results_dir / f"{sample}.qual.tsv" if options.tax_quality else None
                em_out = results_dir / f"{sample}.em.tsv" if options.em else None
                ci_out = (results_dir / f"{sample}.ci.tsv"
                          if (options.em_model == 'generative' and options.em_bootstrap)
                          else None)
                try:
                    if options.em_model == 'generative' and TaxQualityModel is not None:
                        res = TaxQualityModel.run_sample(
                            bam_path=bam_path,
                            tax_map_path=str(generated_tax_map_path),
                            db_format=options.format,
                            flat_sh=(options.flat_sh and options.format == 'unite'),
                            threads=threads_per_job_map,
                            min_percent_identity=options.min_percent_identity,
                            min_ref_coverage=options.min_ref_coverage,
                            tau=options.consensus_tau,
                            max_iter=options.em_max_iter,
                            tol=options.em_tol,
                            prior_alpha=options.em_prior,
                            init=options.em_init,
                            estimate_error=not options.em_no_estimate_error,
                            bootstrap_B=options.em_bootstrap,
                            alpha=options.em_ci_alpha,
                            abundance=options.em_abundance,
                            quality_out=q_out,
                            em_out=em_out,
                            ci_out=ci_out,
                            debug=options.debug,
                        )
                    else:
                        res = TaxQuality.run_sample(
                            bam_path=bam_path,
                            tax_map_path=str(generated_tax_map_path),
                            db_format=options.format,
                            # Flat-SH must reach the read-level pass, otherwise the
                            # quality keys ("acc|SH") never match the OTU table
                            # index ("SH") and every quality column joins to NaN.
                            flat_sh=(options.flat_sh and options.format == 'unite'),
                            threads=threads_per_job_map,
                            # Both gates are threaded through so the EM table and
                            # the quality columns describe the SAME read set as
                            # the primary OTU table.
                            min_percent_identity=options.min_percent_identity,
                            min_ref_coverage=options.min_ref_coverage,
                            tau=options.consensus_tau,
                            temperature=options.score_temperature,
                            do_quality=options.tax_quality,
                            do_em=options.em,
                            chimera_refined=options.chimera_refined,
                            em_max_iter=options.em_max_iter,
                            em_tol=options.em_tol,
                            em_prior=options.em_prior,
                            em_init=options.em_init,
                            quality_out=q_out,
                            em_out=em_out,
                            debug=options.debug,
                        )
                except Exception as e:
                    print(f"[TaxQuality] {sample}: ERROR during read-level pass: {e}")
                    continue

                # Merge the new metric blocks into this sample's stats entry
                for key in ('taxonomic_quality', 'noise', 'em', 'model',
                            'uncertainty', 'reads_excluded_by_quality_filter'):
                    if key in res:
                        stats[key] = res[key]
                if q_out is not None and q_out.exists():
                    quality_files.append(str(q_out))
                if em_out is not None and em_out.exists():
                    em_files_present = True
                info = res.get('model') or res.get('em') or {}
                extra = ""
                if info:
                    extra = f" ({info.get('iterations','?')} iters"
                    if 'pi' in info:
                        extra += f", pi={info['pi']:.4f}"
                    if not info.get('converged', True):
                        extra += ", NOT CONVERGED"
                    extra += ")"
                print(f"[TaxQuality] {sample}: done{extra}")

        # Honour --delete-bam now that the read-level pass is finished.
        # Deliberately OUTSIDE `if TaxQuality is not None`: when the import
        # fails the pass is skipped, but the user still asked for the BAMs to be
        # deleted — previously that combination deleted nothing at all.  Both
        # the full BAM (mapping/) and the filtered BAM (tabeling/) are removed;
        # tabeling_improved's own cleanup was disabled for this run, so leaving
        # the filtered BAM behind leaked the larger of the two files.
        if options.delete_bam:
            tab_dir = output_dir / 'tabeling'
            removed = 0
            for stats in mapping_stats_list:
                sample = stats.get('sample')
                if not sample:
                    continue
                for p in (mapping_dir / f"{sample}_sorted_full.bam",
                          mapping_dir / f"{sample}_sorted_full.bam.bai",
                          tab_dir / f"{sample}_sorted_filtered.bam",
                          tab_dir / f"{sample}_sorted_filtered.bam.bai"):
                    if p.exists():
                        try:
                            p.unlink()
                            removed += 1
                        except Exception:
                            pass
            print(f"[--delete-bam] Removed {removed} intermediate BAM/index files.")

        if profiler:
            profiler.end_step('Step 4b: TaxQuality/EM')

    # Save mapping stats to JSON for the HTML report
    if mapping_stats_list:
        total_reads_all = sum(s.get('total_reads', 0) for s in mapping_stats_list)
        total_mapped_all = sum(s.get('mapped_unfiltered', 0) for s in mapping_stats_list)
        total_filtered_all = sum(s.get('mapped_filtered', 0) for s in mapping_stats_list)
        overall_mapping_rate = round(total_mapped_all / total_reads_all * 100, 2) if total_reads_all > 0 else 0.0
        overall_retention = round(total_filtered_all / total_mapped_all * 100, 2) if total_mapped_all > 0 else 100.0

        # Warn the user when the quality filter removes everything
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
                'tax_quality': options.tax_quality,
                'em': options.em,
                'consensus_tau': options.consensus_tau,
                'chimera_refined': options.chimera_refined,
                'score_temperature': options.score_temperature,
                'em_init': options.em_init,
                'em_model': options.em_model,
                'em_bootstrap': options.em_bootstrap,
                'em_prior': options.em_prior,
                'flat_sh': options.flat_sh,
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



    # =======================================================================
    # STEP 5: OTU Aggregation
    # (OPT-A1: Step 5 "Reformat" eliminated — DuckDB handles OTU ID
    #  reformatting via _build_reformat_sql() inside the PIVOT query)
    # =======================================================================
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
            quality_files=quality_files if quality_files else None,
        )
    except Exception as e:
        print(f"ERROR during final aggregation: {e}")

    # Feature B: second, non-destructive OTU table with fractional EM counts.
    if options.em and em_files_present:
        em_output_file = output_dir / f"OTU_Table_{analysis_name}_{db_name}_EM.xlsx"
        print(f"\n--- STEP 5b: EM OTU Table ({em_output_file.name}) ---")
        try:
            makeOtu_duckdb(
                results_dir=str(results_dir),
                output_file=str(em_output_file),
                taxonomy_file=str(generated_tax_map_path),
                db_format=options.format,
                flat_sh=options.flat_sh,
                file_suffix='.em.tsv',
                # EM counts are genuinely fractional; rounding them to int
                # silently destroyed read mass across a long tail of small
                # assignments.
                fractional_counts=True,
            )
            print(f"EM OTU table saved to: {em_output_file.name}")
        except Exception as e:
            print(f"ERROR during EM aggregation: {e}")

    if profiler:
        profiler.end_step('Step 5: OTU Aggregation')
        profiler.save(output_dir, db_name=db_name,
                      total_threads=options.total_threads,
                      threads_per_job=options.threads_per_job)

    print(f"\n--- Pipeline Complete ---")
    print(f"OTU table saved to: {final_output_file.name}")


if __name__ == "__main__":
    main()
