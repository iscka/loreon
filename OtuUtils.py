import statistics
import random
import re
import subprocess
from pathlib import Path
import sys
import multiprocessing
from functools import partial

# OPT-A3: isal (Intel ISA-L) is 4-5x faster than stdlib gzip for decompression
try:
    from isal import igzip as gzip
except ImportError:
    import gzip


# ---------------------------------------------------------------------------
# OPT-1: Raw FASTQ/FASTA iterators — replace BioPython SeqIO for speed
# ---------------------------------------------------------------------------

def _iter_fastq_raw(handle):
    """Yield (header, seq, plus, qual, seq_len) tuples from a FASTQ handle.
    ~5-10x faster than Bio.SeqIO.parse for simple length filtering."""
    while True:
        header = handle.readline()
        if not header:
            return
        seq = handle.readline().rstrip('\n\r')
        plus = handle.readline()
        qual = handle.readline().rstrip('\n\r')
        yield header, seq, plus, qual, len(seq)


def _iter_fasta_raw(handle):
    """Yield (header, seq, seq_len) tuples from a FASTA handle.
    Handles multi-line sequences."""
    header = None
    seq_parts = []
    for line in handle:
        line = line.rstrip('\n\r')
        if line.startswith('>'):
            if header is not None:
                seq = ''.join(seq_parts)
                yield header, seq, len(seq)
            header = line
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        seq = ''.join(seq_parts)
        yield header, seq, len(seq)


# ---------------------------------------------------------------------------
# OPT-8: Reservoir sampler for median estimation
# ---------------------------------------------------------------------------

_RESERVOIR_SIZE = 10_000


def _reservoir_add(reservoir, n, value):
    """Vitter's Algorithm R — maintains a random sample of fixed size."""
    n += 1
    if n <= _RESERVOIR_SIZE:
        reservoir.append(value)
    else:
        j = random.randint(0, n - 1)
        if j < _RESERVOIR_SIZE:
            reservoir[j] = value
    return n


# ---------------------------------------------------------------------------
# Filter & merge
# ---------------------------------------------------------------------------

def filter_and_merge_directory(
    input_dir: Path,
    min_len: int,
    max_len: int,
    out_dir: Path,
    debug: bool = False,
    gzip_intermediate: bool = False
):
    dir_name = input_dir.name

    # Treat 0 as "no limit" (consistent with quality filter logic)
    if max_len <= 0:
        max_len = 999_999_999

    file_list = list(input_dir.glob('**/*.f*q.gz')) + \
                list(input_dir.glob('**/*.f*q'))

    file_list = [f for f in file_list if not f.name.startswith('._')]

    if not file_list:
        if debug:
            print(f"  [Filter {dir_name}] No source files found. Skipping.")
        return None

    # OPT-B5: optional on-the-fly gzip of merged intermediates — halves disk
    # usage for large ONT runs. minimap2 reads .gz natively, so downstream
    # mapping is unaffected.
    out_suffix = '.fastq.gz' if gzip_intermediate else '.fastq'

    print(f"  [Filter {dir_name}] Merging {len(file_list)} files into '{dir_name}{out_suffix}'...", flush=True)

    try:
        good_dir = out_dir / 'good_seq'
        bad_dir = out_dir / 'bad_seq'
        good_dir.mkdir(parents=True, exist_ok=True)
        bad_dir.mkdir(parents=True, exist_ok=True)

        base_name = dir_name + out_suffix
        good_file_path = good_dir / base_name
        bad_file_path = bad_dir / base_name

    except Exception as e:
        print(f"Error setting up paths for '{dir_name}': {e}")
        return None

    total_count, good_count, bad_count, sum_of_lengths = 0, 0, 0, 0
    # OPT-8: reservoir sampler instead of full list
    reservoir = []
    reservoir_n = 0

    def _open_out(path):
        # isal.igzip.open is used when available (imported as `gzip`)
        if gzip_intermediate:
            return gzip.open(path, 'wt', encoding='utf-8', compresslevel=3)
        return open(path, 'wt', encoding='utf-8', buffering=1_048_576)

    try:
        # OPT-A7: 1 MB write buffers for sequential I/O (uncompressed only)
        with _open_out(good_file_path) as f_good, \
             _open_out(bad_file_path) as f_bad:

            for file_path in file_list:
                if debug:
                    print(f"    -> Reading {file_path.name}...")

                in_is_gzipped = file_path.name.endswith('.gz')
                stem_name = file_path.name.removesuffix('.gz')
                is_fastq = stem_name.endswith(('.fastq', '.fq'))

                read_open_func = gzip.open if in_is_gzipped else open

                # OPT-A7: Larger read buffer (1 MB) for sequential I/O
                with read_open_func(file_path, 'rt', encoding='utf-8',
                                    **({'buffering': 1_048_576} if not in_is_gzipped else {})) as f_in:
                    if is_fastq:
                        for header, seq, plus, qual, seq_len in _iter_fastq_raw(f_in):
                            total_count += 1
                            sum_of_lengths += seq_len
                            reservoir_n = _reservoir_add(reservoir, reservoir_n, seq_len)

                            if min_len <= seq_len <= max_len:
                                f_good.write(header)
                                f_good.write(seq + '\n')
                                f_good.write(plus)
                                f_good.write(qual + '\n')
                                good_count += 1
                            else:
                                f_bad.write(header)
                                f_bad.write(seq + '\n')
                                f_bad.write(plus)
                                f_bad.write(qual + '\n')
                                bad_count += 1
                    else:
                        # FASTA input → write as FASTQ with dummy quality
                        for header_line, seq, seq_len in _iter_fasta_raw(f_in):
                            total_count += 1
                            sum_of_lengths += seq_len
                            reservoir_n = _reservoir_add(reservoir, reservoir_n, seq_len)

                            rec_id = header_line[1:].split()[0] if header_line.startswith('>') else 'unknown'
                            fq_header = f"@{rec_id}\n"
                            fq_qual = 'I' * seq_len + '\n'

                            if min_len <= seq_len <= max_len:
                                f_good.write(fq_header)
                                f_good.write(seq + '\n')
                                f_good.write('+\n')
                                f_good.write(fq_qual)
                                good_count += 1
                            else:
                                f_bad.write(fq_header)
                                f_bad.write(seq + '\n')
                                f_bad.write('+\n')
                                f_bad.write(fq_qual)
                                bad_count += 1

    except Exception as e:
        print(f"ERROR during processing of {file_path.name}: {e}")
        return None

    mean_len = (sum_of_lengths / total_count) if total_count > 0 else 0
    median_len = statistics.median(reservoir) if reservoir else 0

    print(f"  [Filter {dir_name}] Merge complete.", flush=True)
    print(f"     Total: {total_count}, Good: {good_count}, Discarded: {bad_count}", flush=True)

    if good_count == 0 and good_file_path.exists():
        good_file_path.unlink()
    if bad_count == 0 and bad_file_path.exists():
        bad_file_path.unlink()

    return {
        'barcode': dir_name,
        'mean_len': mean_len,
        'median_len': median_len,
        'bad_count': bad_count,
        'good_count': good_count,
        'total_count': total_count,
        'good_file_out': str(good_file_path) if good_count > 0 else None,
    }


# ---------------------------------------------------------------------------
# Mapping (minimap2)
# ---------------------------------------------------------------------------

def mapping_improved(fastq_file: str, db_path: str, threads: int, output_dir: str,
                     kmer_size: int = 15, window_size: int = 10, debug: bool = False):
    """OPT-A2: Pipe minimap2 → samtools sort → BAM (eliminates SAM file on disk).

    Returns the path to the sorted BAM file, or None on failure.
    """
    try:
        fastq_path = Path(fastq_file)
        db_path_obj = Path(db_path)
        output_dir_path = Path(output_dir)

        # OPT-B5: strip both `.gz` and `.fastq`/`.fq` so that BAM names
        # (and downstream sample names) remain `barcodeXX_sorted_full.bam`
        # regardless of whether the input is gzipped.
        fastq_stem = fastq_path.name
        for suffix in ('.gz', '.fastq', '.fq'):
            if fastq_stem.lower().endswith(suffix):
                fastq_stem = fastq_stem[: -len(suffix)]
        bam_filename = fastq_stem + '_sorted_full.bam'
        bam_output_path = output_dir_path / bam_filename
        thread_str = str(threads)

        if debug:
            print(f"[{fastq_path.name}] Starting mapping (minimap2 → BAM pipe)...", flush=True)
            print(f"  DB: {db_path_obj.name}", flush=True)
            print(f"  Output: {bam_output_path}", flush=True)

    except Exception as e:
        print(f"ERROR in path setup for {fastq_file}: {e}")
        return None

    # OPT-B2: `--secondary=no` suppresses secondary alignments upstream.
    # Previously minimap2 emitted up to 5 secondary per read and they were
    # all discarded later by `samtools view -F 0x904`, wasting I/O and
    # BAM size (~20-30%). Disable them at the source.
    cmd_mm2 = [
        'minimap2',
        '-ax', 'map-ont',
        '--secondary=no',
        '-t', thread_str,
        '-k', str(kmer_size),
        '-w', str(window_size),
        str(db_path_obj),
        str(fastq_path)
    ]

    cmd_sort = [
        'samtools', 'sort',
        f'-@{thread_str}',
        '-o', str(bam_output_path),
        '-'
    ]

    try:
        # OPT-A2: Pipe minimap2 stdout directly to samtools sort
        # No intermediate SAM file — saves GBs of disk I/O
        proc_mm2 = subprocess.Popen(
            cmd_mm2, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        proc_sort = subprocess.Popen(
            cmd_sort, stdin=proc_mm2.stdout,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        proc_mm2.stdout.close()  # allow minimap2 to receive SIGPIPE if sort exits

        _, sort_stderr = proc_sort.communicate()
        mm2_stderr = proc_mm2.stderr.read().decode('utf-8', errors='replace')
        mm2_rc = proc_mm2.wait()

        if mm2_rc != 0:
            print(f"ERROR running minimap2 for {fastq_path.name} (rc={mm2_rc}).", flush=True)
            print(f"Command: {' '.join(cmd_mm2)}", flush=True)
            if mm2_stderr:
                print("Stderr:", mm2_stderr.strip(), flush=True)
            return None

        if proc_sort.returncode != 0:
            raise subprocess.CalledProcessError(
                proc_sort.returncode, cmd_sort,
                stderr=sort_stderr.decode('utf-8', errors='replace')
            )

        if debug:
            print(f"[{fastq_path.name}] Mapping + sort complete (BAM pipe).", flush=True)
            if mm2_stderr:
                print("  Info (from minimap2 stderr):", flush=True)
                print(mm2_stderr.strip(), flush=True)

        return str(bam_output_path)

    except subprocess.CalledProcessError as e:
        print(f"ERROR running pipeline for {fastq_path.name}.", flush=True)
        print(f"Command: {e.cmd}", flush=True)
        if e.stderr:
            print("Stderr:", e.stderr, flush=True)
        return None

    except FileNotFoundError:
        print("CRITICAL ERROR: 'minimap2' or 'samtools' command not found.", flush=True)
        print("Ensure both are installed and on your PATH.", flush=True)
        sys.exit(1)

    except Exception as e:
        print(f"Unexpected error during mapping of {fastq_path.name}: {e}", flush=True)
        return None


# ---------------------------------------------------------------------------
# Tabeling — OPT-2: reduced from 9 to 4 samtools subprocess calls
# ---------------------------------------------------------------------------

def _parse_flagstat(flagstat_output: str) -> dict:
    """Parse samtools flagstat output and return primary read counts."""
    total = 0
    mapped = 0
    for line in flagstat_output.splitlines():
        parts = line.split()
        if not parts:
            continue
        try:
            count = int(parts[0])
        except ValueError:
            continue
        if 'primary' in line and 'mapped' in line and 'secondary' not in line and 'supplementary' not in line:
            mapped = count
        elif 'primary' in line and 'mapped' not in line and 'duplicate' not in line:
            total = count
    return {'total': total, 'mapped': mapped}


def tabeling_improved(samfile: str, threads: int = 4, debug: bool = False,
                      min_percent_identity: float = 0.0, min_ref_coverage: float = 0.0,
                      delete_bam: bool = False):
    """Process alignment file (BAM or SAM) into OTU count results.

    OPT-A2: Now accepts sorted BAM from piped mapping (no intermediate SAM).
    OPT-A4: Post-filter stats derived from idxstats (eliminates second flagstat).
    """
    fname = "unknown"
    filtered_sorted_bam = None
    input_bam = samfile  # may be BAM from piped mapping or legacy SAM

    try:
        input_path = Path(input_bam)
        # Strip _sorted_full suffix if present (from piped mapping)
        fname = input_path.stem.replace('_sorted_full', '')

        base_dir = input_path.parent.parent

        results_dir = base_dir / 'tabeling' / 'results'
        tab_dir = base_dir / 'tabeling'

        results_dir.mkdir(parents=True, exist_ok=True)
        tab_dir.mkdir(parents=True, exist_ok=True)

        filtered_sorted_bam = tab_dir / f"{fname}_sorted_filtered.bam"
        results_file = results_dir / f"{fname}.tmp.txt"

        thread_str = str(threads)

        if debug:
            print(f"--- [Tabeling: {fname}] with {thread_str} threads ---", flush=True)

    except Exception as e:
        print(f"[{fname}] Error in path setup: {e}")
        return None

    # Build optional samtools -e expression for quality filters
    # NOTE (OPT-B1): coverage is now computed as QUERY coverage (alen/qlen).
    # The previous formula `alen*100>=X*rlen` was semantically wrong: in
    # htslib filter expressions `rlen` is the reference consumption of THIS
    # alignment (M+D from CIGAR), not the full reference length — so the
    # ratio measured insertion/deletion balance rather than "% of reference
    # covered". Query coverage (fraction of the read that actually aligned,
    # excluding soft clips) is the meaningful quantity for metabarcoding
    # amplicons: reads should align ≥90% of their own length.
    quality_conditions = []
    if min_percent_identity > 0.0:
        quality_conditions.append(
            f"(alen-[NM])*100>={min_percent_identity}*alen"
        )
    if min_ref_coverage > 0.0:
        quality_conditions.append(
            f"alen*100>={min_ref_coverage}*qlen"
        )
    e_filter_args = []          # list form for Popen (cross-platform safe)
    if quality_conditions:
        combined_expr = " && ".join(quality_conditions)
        e_filter_args = ['-e', combined_expr]
        print(f"[{fname}] Quality filter expression: {combined_expr}", flush=True)

    mapping_stats = {
        'sample': fname,
        'total_reads': 0,
        'mapped_unfiltered': 0,
        'mapped_filtered': 0,
        'mapping_rate_pct': 0.0,
        'quality_retention_pct': 100.0,
        'min_percent_identity': min_percent_identity,
        'min_ref_coverage': min_ref_coverage,
    }

    try:
        # --- Step 1: Pre-filter stats via flagstat on input BAM ---
        # (Fast on BAM; much faster than old approach on SAM text)
        result_flagstat_raw = subprocess.run(
            ['samtools', 'flagstat', f'-@{thread_str}', input_bam],
            check=True, capture_output=True, text=True, encoding='utf-8'
        )
        stats_raw = _parse_flagstat(result_flagstat_raw.stdout)
        mapping_stats['total_reads'] = stats_raw['total']
        mapping_stats['mapped_unfiltered'] = stats_raw['mapped']
        if stats_raw['total'] > 0:
            mapping_stats['mapping_rate_pct'] = round(stats_raw['mapped'] / stats_raw['total'] * 100, 2)

        # Compute ll1 (unmapped count before filter) from flagstat
        ll1_unmapped = stats_raw['total'] - stats_raw['mapped']
        ll1 = f"*\t0\t0\t{ll1_unmapped}"

        if debug:
            print(f"[{fname}] Pre-filter: total={stats_raw['total']}, mapped={stats_raw['mapped']}", flush=True)

        # --- Step 2: Filter -F 0x904 + quality → sort → filtered BAM ---
        cmd_view = [
            'samtools', 'view', f'-@{thread_str}', '-b', '-F', '0x904',
        ] + e_filter_args + [input_bam]
        cmd_sort = [
            'samtools', 'sort', f'-@{thread_str}', '-o', str(filtered_sorted_bam), '-',
        ]
        if debug:
            print(f"[{fname}] CMD (view): {' '.join(cmd_view)}", flush=True)
            print(f"[{fname}] CMD (sort): {' '.join(cmd_sort)}", flush=True)
        proc_view = subprocess.Popen(cmd_view, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc_sort = subprocess.Popen(cmd_sort, stdin=proc_view.stdout,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc_view.stdout.close()  # allow view to receive SIGPIPE if sort exits
        _, sort_stderr = proc_sort.communicate()
        view_rc = proc_view.wait()
        if view_rc != 0:
            raise subprocess.CalledProcessError(view_rc, cmd_view)
        if proc_sort.returncode != 0:
            raise subprocess.CalledProcessError(proc_sort.returncode, cmd_sort, stderr=sort_stderr)

        # --- Step 3: Index the filtered BAM ---
        subprocess.run(['samtools', 'index', f'-@{thread_str}', str(filtered_sorted_bam)],
                       check=True, capture_output=True)

        # --- Step 4: idxstats → results file + derive post-filter stats ---
        # OPT-A4: Replaces second flagstat call — derive mapped_filtered from
        # the sum of idxstats column 3 (mapped reads per reference)
        result_idxstats = subprocess.run(
            ['samtools', 'idxstats', f'-@{thread_str}', str(filtered_sorted_bam)],
            check=True, capture_output=True, text=True, encoding='utf-8'
        )

        total_mapped_filtered = 0
        unmapped_in_filtered = 0
        # OPT-B6: collect per-reference counts for QC stats in mapping_stats.json
        per_ref_counts = []
        with open(results_file, 'w', encoding='utf-8') as f:
            for line in result_idxstats.stdout.splitlines():
                parts = line.strip().split('\t')
                if not parts or not parts[0]:
                    continue
                if parts[0] == '*':
                    # Capture unmapped count from idxstats * line
                    unmapped_in_filtered = int(parts[3]) if len(parts) >= 4 else 0
                    continue
                if len(parts) >= 3 and parts[2] != '0':
                    ref_count = int(parts[2])
                    total_mapped_filtered += ref_count
                    per_ref_counts.append((parts[0], ref_count))
                    f.write(line.strip() + '\n')
            # Append unmapped counts for backward compatibility
            f.write(ll1 + '\n')
            f.write(f"*\t0\t0\t{unmapped_in_filtered}\n")

        mapping_stats['mapped_filtered'] = total_mapped_filtered
        if stats_raw['mapped'] > 0:
            mapping_stats['quality_retention_pct'] = round(
                total_mapped_filtered / stats_raw['mapped'] * 100, 2
            )

        # OPT-B6: per-reference QC summary — helps diagnose over/under-represented
        # references in metabarcoding runs (e.g. a single SH capturing 80% of reads)
        mapping_stats['references_detected'] = len(per_ref_counts)
        if per_ref_counts:
            counts_only = sorted((c for _, c in per_ref_counts), reverse=True)
            top = sorted(per_ref_counts, key=lambda x: -x[1])[:10]
            # Percentiles computed from sorted-descending list
            n = len(counts_only)

            def _pct(p):
                # p in [0,100]; 50 = median, 95 = top-5% threshold
                idx = min(n - 1, max(0, int(round((1 - p / 100.0) * (n - 1)))))
                return counts_only[idx]

            top_count = counts_only[0]
            top_share = round(top_count / total_mapped_filtered * 100, 2) if total_mapped_filtered else 0.0
            mapping_stats['per_ref'] = {
                'top10': [{'ref': r, 'count': c} for r, c in top],
                'top_ref_share_pct': top_share,
                'median_count_per_ref': _pct(50),
                'p95_count_per_ref': _pct(95),
                'singleton_refs': sum(1 for c in counts_only if c == 1),
            }
        else:
            mapping_stats['per_ref'] = {
                'top10': [], 'top_ref_share_pct': 0.0,
                'median_count_per_ref': 0, 'p95_count_per_ref': 0,
                'singleton_refs': 0,
            }

        if debug:
            print(f"[{fname}] Post-filter: mapped={total_mapped_filtered} (from idxstats)", flush=True)
            print(f"--- [Tabeling: {fname}] Complete ---", flush=True)

        return mapping_stats

    except subprocess.CalledProcessError as e:
        print(f"ERROR [Tabeling: {fname}].", flush=True)
        print(f"Failed command: {e.cmd}", flush=True)
        print("Stderr:", flush=True)
        print(e.stderr, flush=True)
        return None

    except FileNotFoundError as e:
        print(f"CRITICAL ERROR: Command not found. {e}", flush=True)
        print("Ensure 'samtools' is installed and on your PATH.", flush=True)
        sys.exit(1)

    finally:
        # OPT-14 / Delete BAM: clean up intermediate BAM files
        if delete_bam:
            for p in [filtered_sorted_bam,
                      Path(f"{filtered_sorted_bam}.bai") if filtered_sorted_bam else None,
                      Path(input_bam),
                      Path(f"{input_bam}.bai") if Path(f"{input_bam}.bai").exists() else None]:
                if p and Path(p).exists():
                    try:
                        Path(p).unlink()
                    except Exception:
                        pass


# ---------------------------------------------------------------------------
# OPT-11: Parallel reformat_tmp_indices
# ---------------------------------------------------------------------------

def _reformat_single_file(file_path: Path, format: str, debug: bool):
    """Reformat a single .tmp.txt file (worker for parallel pool)."""
    if debug:
        print(f"  Reformatting: {file_path.name}", flush=True)

    processed_lines = []
    try:
        with open(file_path, 'r', encoding='utf-8') as f_in:
            for line in f_in:
                line_stripped = line.strip()
                if not line_stripped:
                    continue
                parts = line_stripped.split('\t')
                if not parts:
                    continue

                index_raw = parts[0]

                if index_raw == '*':
                    processed_lines.append(line_stripped)
                    continue

                if format == 'unite':
                    index_parts = index_raw.split('|')
                    if len(index_parts) >= 5:
                        parts[0] = f"{index_parts[1]}|{index_parts[2]}"
                        processed_lines.append('\t'.join(parts))
                    else:
                        if debug:
                            print(f"  [reformat Warning] Skipping invalid UNITE line: {index_raw[:50]}...")
                        processed_lines.append(line_stripped)

                elif format == 'eukariome':
                    index_parts = index_raw.split(';', 1)
                    if len(index_parts) >= 1:
                        parts[0] = index_parts[0]
                        processed_lines.append('\t'.join(parts))
                    else:
                        processed_lines.append(line_stripped)

                elif format == 'cbs' or format == 'none':
                    index_parts = index_raw.split('|', 1)
                    if len(index_parts) >= 1:
                        parts[0] = index_parts[0]
                        processed_lines.append('\t'.join(parts))
                    else:
                        processed_lines.append(line_stripped)

                elif format == 'silva':
                    processed_lines.append(line_stripped)

        with open(file_path, 'w', encoding='utf-8') as f_out:
            for pl in processed_lines:
                f_out.write(pl + '\n')

        return file_path.name

    except Exception as e:
        print(f"ERROR [reformat_tmp] processing {file_path.name}: {e}")
        return None


def reformat_tmp_indices(
        results_dir: str,
        format: str = 'unite',
        debug: bool = False
):
    results_path = Path(results_dir)
    if not results_path.exists():
        print(f"ERROR [reformat_tmp]: Directory not found: {results_path}")
        return

    print(f"\n--- Step 5: Reformatting Taxonomic Indices (Format: {format}) ---")

    try:
        file_list = [f for f in results_path.glob('*.tmp.txt') if not f.name.startswith('._')]
        if not file_list:
            print("WARNING [reformat_tmp]: No valid .tmp.txt files found.")
            return
    except Exception as e:
        print(f"ERROR [reformat_tmp]: Cannot search files in {results_path}: {e}")
        return

    # Parallel processing for multiple files, sequential for single file
    if len(file_list) <= 1:
        for fp in file_list:
            _reformat_single_file(fp, format, debug)
    else:
        pool_size = min(len(file_list), multiprocessing.cpu_count(), 8)
        worker = partial(_reformat_single_file, format=format, debug=debug)
        with multiprocessing.Pool(processes=pool_size) as pool:
            results = pool.map(worker, file_list)
        processed_count = sum(1 for r in results if r is not None)
        print(f"Reformatting complete. {processed_count} files updated.")
        return

    print(f"Reformatting complete. {len(file_list)} files updated.")


# ---------------------------------------------------------------------------
# OPT-4: Unified taxonomy map creation — raw header parser, no BioPython
# ---------------------------------------------------------------------------

def _parse_fasta_headers(fasta_path: Path, tax_map_file: Path, format: str, debug: bool):
    """Parse FASTA headers (without reading sequences) and write taxonomy map.
    Returns the number of entries written. Replaces both _create_tax_map and
    _run_fasta_parsing with a single, faster implementation."""
    # Pre-compile the regex for formats that need it
    rank_re = re.compile(r'[kpcofgs]__') if format in ('unite', 'eukariome', 'cbs') else None

    try:
        count = 0
        with open(fasta_path, 'r', encoding='utf-8') as f_in, \
             open(tax_map_file, 'w', encoding='utf-8') as f_out:
            f_out.write("OTU_ID\tTaxonomy\n")

            for line in f_in:
                if not line.startswith('>'):
                    continue  # skip sequence lines — we only need headers

                header = line[1:].rstrip('\n\r')
                final_id = ""
                taxonomy = ""

                if format == 'unite':
                    parts = header.split('|')
                    if len(parts) >= 5:
                        final_id = f"{parts[1]}|{parts[2]}"
                        taxonomy = rank_re.sub('', parts[-1])
                    else:
                        continue

                elif format == 'eukariome':
                    parts = header.split(';', 1)
                    if len(parts) == 2:
                        final_id = parts[0]
                        taxonomy = rank_re.sub('', parts[1])
                    else:
                        continue

                elif format == 'cbs' or format == 'none':
                    parts = header.split('|', 1)
                    if len(parts) == 2:
                        final_id = parts[0]
                        taxonomy = rank_re.sub('', parts[1]) if format == 'cbs' else parts[1]
                    else:
                        continue

                elif format == 'silva':
                    parts = header.split(' ', 1)
                    if len(parts) == 2:
                        final_id = parts[0]
                        taxonomy = parts[1]
                    else:
                        continue

                if final_id:
                    f_out.write(f"{final_id}\t{taxonomy}\n")
                    count += 1

        return count

    except Exception as e:
        print(f"CRITICAL ERROR during parsing of {fasta_path.name}: {e}")
        if tax_map_file.exists():
            tax_map_file.unlink()
        raise


def create_taxonomy_map_from_fasta(
        fasta_file: str,
        output_dir: Path,
        format: str = 'unite',
        force_create: bool = False,
        debug: bool = False
):
    try:
        fasta_path = Path(fasta_file)
        if not fasta_path.exists():
            print(f"CRITICAL ERROR: FASTA database file does not exist: {fasta_path}")
            sys.exit(1)

        db_folder = fasta_path.parent
        fasta_stem = fasta_path.stem
        tax_map_name = f"{fasta_stem}_taxonomy_map.tsv"

        db_map_path = db_folder / tax_map_name
        output_map_path = output_dir / tax_map_name

    except Exception as e:
        print(f"CRITICAL ERROR in taxonomy map path setup: {e}")
        sys.exit(1)

    if force_create:
        print("Forcing taxonomy map creation (user request)...")
        try:
            count = _parse_fasta_headers(fasta_path, output_map_path, format, debug)
            print(f"Forced map created at: {output_map_path} ({count} entries)")
            return output_map_path, count
        except Exception as e:
            print(f"CRITICAL ERROR: Cannot create forced map: {e}")
            sys.exit(1)

    # OPT-B3: invalidate cache if the source FASTA is newer than the cached
    # taxonomy map — previously, updating UNITE/SILVA while keeping the same
    # filename silently reused the stale map, producing wrong taxonomy.
    def _cache_is_fresh(cache_path: Path) -> bool:
        try:
            return cache_path.stat().st_mtime >= fasta_path.stat().st_mtime
        except OSError:
            return False

    if db_map_path.exists():
        if _cache_is_fresh(db_map_path):
            print(f"Taxonomy map found (Cache DB): {db_map_path}")
            count = sum(1 for _ in open(db_map_path, encoding='utf-8')) - 1
            return db_map_path, count
        else:
            print(f"Taxonomy map stale (FASTA newer than cache): {db_map_path} — regenerating")

    if output_map_path.exists():
        if _cache_is_fresh(output_map_path):
            print(f"Taxonomy map found (Cache Output): {output_map_path}")
            count = sum(1 for _ in open(output_map_path, encoding='utf-8')) - 1
            return output_map_path, count
        else:
            print(f"Taxonomy map stale (FASTA newer than cache): {output_map_path} — regenerating")

    print("Taxonomy map not found in cache. Creating (may take time)...")

    try:
        print(f"Attempting to write to: {db_map_path}")
        count = _parse_fasta_headers(fasta_path, db_map_path, format, debug)
        print(f"Taxonomy map created: {db_map_path.name} ({count} entries).")
        return db_map_path, count
    except (OSError, PermissionError) as e:
        print(f"  [Info] DB folder not writable ({e}). Writing to output folder.")
    except Exception as e:
        print(f"CRITICAL ERROR during parsing: {e}")
        sys.exit(1)

    try:
        print(f"Attempting to write to: {output_map_path}")
        count = _parse_fasta_headers(fasta_path, output_map_path, format, debug)
        print(f"Taxonomy map created: {output_map_path.name} ({count} entries).")
        return output_map_path, count
    except Exception as e:
        print(f"CRITICAL ERROR: Cannot create taxonomy map in any folder: {e}")
        sys.exit(1)
