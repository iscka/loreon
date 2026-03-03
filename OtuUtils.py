import statistics
import re
import subprocess
from pathlib import Path
import sys
import gzip

from Bio import SeqIO


def filter_and_merge_directory(
    input_dir: Path,
    min_len: int,
    max_len: int,
    out_dir: Path,
    debug: bool = False
):
    dir_name = input_dir.name

    file_list = list(input_dir.glob('**/*.f*q.gz')) + \
                list(input_dir.glob('**/*.f*q'))

    file_list = [f for f in file_list if not f.name.startswith('._')]

    if not file_list:
        if debug:
            print(f"  [Filter {dir_name}] No source files found. Skipping.")
        return None

    print(f"  [Filter {dir_name}] Merging {len(file_list)} files into '{dir_name}.fastq'...")

    try:
        good_dir = out_dir / 'good_seq'
        bad_dir = out_dir / 'bad_seq'
        good_dir.mkdir(parents=True, exist_ok=True)
        bad_dir.mkdir(parents=True, exist_ok=True)

        base_name = dir_name + '.fastq'
        good_file_path = good_dir / base_name
        bad_file_path = bad_dir / base_name

        out_file_format = 'fastq'
        write_mode = 'wt'

    except Exception as e:
        print(f"Error setting up paths for '{dir_name}': {e}")
        return None

    total_count, good_count, bad_count, sum_of_lengths = 0, 0, 0, 0
    all_lengths = []

    try:
        with open(good_file_path, write_mode, encoding='utf-8') as f_good, \
             open(bad_file_path, write_mode, encoding='utf-8') as f_bad:

            for file_path in file_list:
                if debug:
                    print(f"    -> Reading {file_path.name}...")

                in_is_gzipped = file_path.name.endswith('.gz')
                stem_name = file_path.name.removesuffix('.gz')
                in_file_format = 'fastq' if stem_name.endswith(('.fastq', '.fq')) else 'fasta'

                read_open_func = gzip.open if in_is_gzipped else open
                read_mode = 'rt'

                with read_open_func(file_path, read_mode, encoding='utf-8') as f_in:
                    for record in SeqIO.parse(f_in, in_file_format):
                        total_count += 1
                        seq_len = len(record.seq)
                        all_lengths.append(seq_len)
                        sum_of_lengths += seq_len

                        if min_len <= seq_len <= max_len:
                            SeqIO.write(record, f_good, out_file_format)
                            good_count += 1
                        else:
                            SeqIO.write(record, f_bad, out_file_format)
                            bad_count += 1

    except Exception as e:
        print(f"ERROR during processing of {file_path.name}: {e}")
        return None

    mean_len = (sum_of_lengths / total_count) if total_count > 0 else 0
    median_len = statistics.median(all_lengths) if total_count > 0 else 0

    print(f"  [Filter {dir_name}] Merge complete.")
    print(f"     Total: {total_count}, Good: {good_count}, Discarded: {bad_count}")

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


def mapping_improved(fastq_file: str, db_path: str, threads: int, output_dir: str,
                     kmer_size: int = 15, window_size: int = 10, debug: bool = False):
    try:
        fastq_path = Path(fastq_file)
        db_path_obj = Path(db_path)
        output_dir_path = Path(output_dir)

        sam_filename = fastq_path.stem + '.sam'
        sam_output_path = output_dir_path / sam_filename
        thread_str = str(threads)

        if debug:
            print(f"[{fastq_path.name}] Starting mapping (minimap2)...")
            print(f"  DB: {db_path_obj.name}")
            print(f"  Output: {sam_output_path}")

    except Exception as e:
        print(f"ERROR in path setup for {fastq_file}: {e}")
        return None

    cmd_list = [
        'minimap2',
        '-ax', 'map-ont',
        '-t', thread_str,
        '-k', str(kmer_size),
        '-w', str(window_size),
        str(db_path_obj),
        str(fastq_path)
    ]

    try:
        with open(sam_output_path, 'w') as f_sam_out:
            result = subprocess.run(
                cmd_list,
                stdout=f_sam_out,
                stderr=subprocess.PIPE,
                check=True,
                text=True,
                encoding='utf-8'
            )

        if debug:
            print(f"[{fastq_path.name}] Mapping complete.")
            if result.stderr:
                print("  Info (from minimap2 stderr):")
                print(result.stderr.strip())

        return str(sam_output_path)

    except subprocess.CalledProcessError as e:
        print(f"ERROR running minimap2 for {fastq_path.name}.")
        print(f"Command: {' '.join(cmd_list)}")
        print("Stderr:")
        print(e.stderr)
        return None

    except FileNotFoundError:
        print("CRITICAL ERROR: 'minimap2' command not found.")
        print("Ensure minimap2 is installed and on your PATH.")
        sys.exit(1)

    except Exception as e:
        print(f"Unexpected error during mapping of {fastq_path.name}: {e}")
        return None


def tabeling_improved(samfile: str, threads: int = 4, debug: bool = False):
    fname = "unknown"
    unfiltered_sorted_bam = None

    try:
        sam_path = Path(samfile)
        fname = sam_path.stem

        base_dir = sam_path.parent.parent

        results_dir = base_dir / 'tabeling' / 'results'
        tab_dir = base_dir / 'tabeling'

        results_dir.mkdir(parents=True, exist_ok=True)
        tab_dir.mkdir(parents=True, exist_ok=True)

        unfiltered_sorted_bam = tab_dir / f"{fname}.sorted.bam"
        filtered_sorted_bam = tab_dir / f"{fname}_sorted_filtered.bam"
        results_file = results_dir / f"{fname}.tmp.txt"

        thread_str = str(threads)

        if debug:
            print(f"--- [Tabeling: {fname}] with {thread_str} threads ---")

    except Exception as e:
        print(f"[{fname}] Error in path setup: {e}")
        return

    try:
        cmd_unfiltered = (
            f"samtools view -@{thread_str} -bS '{samfile}' | "
            f"samtools sort -@{thread_str} -o '{unfiltered_sorted_bam}' -"
        )
        if debug:
            print(f"[{fname}] CMD (Unfiltered): {cmd_unfiltered}")
        subprocess.run(cmd_unfiltered, shell=True, check=True, capture_output=True, text=True, encoding='utf-8')

        subprocess.run(['samtools', 'index', f'-@{thread_str}', str(unfiltered_sorted_bam)],
                       check=True, capture_output=True)

        result_stats1 = subprocess.run(
            ['samtools', 'idxstats', f'-@{thread_str}', str(unfiltered_sorted_bam)],
            check=True, capture_output=True, text=True, encoding='utf-8'
        )
        ll1 = ""
        for line in result_stats1.stdout.splitlines():
            if line.startswith('*'):
                ll1 = line.strip()
                break
        if debug:
            print(f"[{fname}] ll1 (unmapped before filter): {ll1}")

        cmd_filtered_pipe = (
            f"samtools view -@{thread_str} -b '{unfiltered_sorted_bam}' | "
            f"samtools view -@{thread_str} -b -F 0x904 - | "
            f"samtools sort -@{thread_str} -o '{filtered_sorted_bam}' -"
        )
        if debug:
            print(f"[{fname}] CMD (Filtered Pipe): {cmd_filtered_pipe}")
        subprocess.run(cmd_filtered_pipe, shell=True, check=True, capture_output=True, text=True, encoding='utf-8')

        cmd_index = ['samtools', 'index', f'-@{thread_str}', str(filtered_sorted_bam)]
        if debug:
            print(f"[{fname}] CMD (Index Filtered): {' '.join(cmd_index)}")
        subprocess.run(cmd_index, check=True, capture_output=True, text=True, encoding='utf-8')

        result_stats2 = subprocess.run(
            ['samtools', 'idxstats', f'-@{thread_str}', str(filtered_sorted_bam)],
            check=True, capture_output=True, text=True, encoding='utf-8'
        )
        ll2 = ""
        for line in result_stats2.stdout.splitlines():
            if line.startswith('*'):
                ll2 = line.strip()
                break
        if debug:
            print(f"[{fname}] ll2 (unmapped after filter): {ll2}")

        cmd_results = (
            f"samtools idxstats -@{thread_str} '{filtered_sorted_bam}' | "
            f"awk '($3 != \"0\")' > '{results_file}'"
        )
        if debug:
            print(f"[{fname}] CMD (Results File): {cmd_results}")
        subprocess.run(cmd_results, shell=True, check=True, capture_output=True, text=True, encoding='utf-8')

        if debug:
            print(f"[{fname}] Appending stats to {results_file}")
        with open(results_file, 'a') as f:
            if ll1:
                f.write(f"{ll1}\n")
            if ll2:
                f.write(f"{ll2}\n")

        if debug:
            print(f"--- [Tabeling: {fname}] Complete ---")

    except subprocess.CalledProcessError as e:
        print(f"ERROR [Tabeling: {fname}].")
        print(f"Failed command: {e.cmd}")
        print("Stderr:")
        print(e.stderr)
        return

    except FileNotFoundError as e:
        print(f"CRITICAL ERROR: Command not found. {e}")
        print("Ensure 'samtools' and 'awk' are installed and on your PATH.")
        sys.exit(1)

    finally:
        if unfiltered_sorted_bam and unfiltered_sorted_bam.exists():
            unfiltered_sorted_bam.unlink()
        bai_path = Path(f"{unfiltered_sorted_bam}.bai") if unfiltered_sorted_bam else None
        if bai_path and bai_path.exists():
            bai_path.unlink()


def reformat_tmp_indices(
        results_dir: str,
        format: str = 'unite',
        debug: bool = False
):
    results_path = Path(results_dir)
    if not results_path.exists():
        print(f"ERROR [reformat_tmp]: Directory not found: {results_path}")
        return

    print(f"\n--- Phase 2.5: Reformatting Taxonomic Indices (Format: {format}) ---")

    try:
        file_list = [f for f in results_path.glob('*.tmp.txt') if not f.name.startswith('._')]
        if not file_list:
            print("WARNING [reformat_tmp]: No valid .tmp.txt files found.")
            return
    except Exception as e:
        print(f"ERROR [reformat_tmp]: Cannot search files in {results_path}: {e}")
        return

    processed_count = 0
    for file_path in file_list:
        if debug:
            print(f"  Reformatting: {file_path.name}")

        processed_lines = []
        try:
            with open(file_path, 'r', encoding='utf-8') as f_in:
                lines = f_in.readlines()

            for line in lines:
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
                        new_index = f"{index_parts[1]}|{index_parts[2]}"
                        parts[0] = new_index
                        processed_lines.append('\t'.join(parts))
                    else:
                        if debug:
                            print(f"  [reformat Warning] Skipping invalid UNITE line: {index_raw[:50]}...")
                        processed_lines.append(line_stripped)

                elif format == 'eukariome':
                    index_parts = index_raw.split(';', 1)
                    if len(index_parts) >= 1:
                        new_index = index_parts[0]
                        parts[0] = new_index
                        processed_lines.append('\t'.join(parts))
                    else:
                        processed_lines.append(line_stripped)

                elif format == 'cbs' or format == 'none':
                    index_parts = index_raw.split('|', 1)
                    if len(index_parts) >= 1:
                        new_index = index_parts[0]
                        parts[0] = new_index
                        processed_lines.append('\t'.join(parts))
                    else:
                        processed_lines.append(line_stripped)

                elif format == 'silva':
                    processed_lines.append(line_stripped)

            with open(file_path, 'w', encoding='utf-8') as f_out:
                f_out.write('\n'.join(processed_lines) + '\n')

            processed_count += 1

        except Exception as e:
            print(f"ERROR [reformat_tmp] processing {file_path.name}: {e}")

    print(f"Reformatting complete. {processed_count} files updated.")


def _run_fasta_parsing(fasta_path: Path, tax_map_file: Path, format: str, debug: bool):
    try:
        with open(fasta_path, "r") as f_in, open(tax_map_file, "w") as f_out:
            f_out.write("OTU_ID\tTaxonomy\n")
            count = 0

            for record in SeqIO.parse(f_in, "fasta"):
                header = record.description
                final_id = ""
                taxonomy = ""

                if format == 'unite':
                    parts = header.split('|')
                    if len(parts) >= 5:
                        final_id = f"{parts[1]}|{parts[2]}"
                        tax_raw = parts[-1]
                        taxonomy = re.sub(r'[kpcofgs]__', '', tax_raw)
                    else:
                        continue

                elif format == 'eukariome':
                    parts = header.split(';', 1)
                    if len(parts) == 2:
                        final_id = parts[0]
                        tax_raw = parts[1]
                        taxonomy = re.sub(r'[kpcofgs]__', '', tax_raw)
                    else:
                        continue

                elif format == 'cbs' or format == 'none':
                    parts = header.split('|', 1)
                    if len(parts) == 2:
                        final_id = parts[0]
                        taxonomy = parts[1]
                        if format == 'cbs':
                            taxonomy = re.sub(r'[kpcofgs]__', '', taxonomy)
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


def _create_tax_map(fasta_path: Path, tax_map_file: Path, format: str, debug: bool):
    try:
        with open(fasta_path, "r") as f_in, open(tax_map_file, "w") as f_out:
            f_out.write("OTU_ID\tTaxonomy\n")
            count = 0

            for record in SeqIO.parse(f_in, "fasta"):
                header = record.description
                final_id, taxonomy = "", ""

                if format == 'unite':
                    parts = header.split('|')
                    if len(parts) >= 5:
                        final_id = f"{parts[1]}|{parts[2]}"
                        tax_raw = parts[-1]
                        taxonomy = re.sub(r'[kpcofgs]__', '', tax_raw)
                    else:
                        continue

                elif format == 'eukariome':
                    parts = header.split(';', 1)
                    if len(parts) == 2:
                        final_id, tax_raw = parts
                        taxonomy = re.sub(r'[kpcofgs]__', '', tax_raw)
                    else:
                        continue

                elif format == 'cbs' or format == 'none':
                    parts = header.split('|', 1)
                    if len(parts) == 2:
                        final_id, taxonomy = parts
                        if format == 'cbs':
                            taxonomy = re.sub(r'[kpcofgs]__', '', taxonomy)
                    else:
                        continue

                elif format == 'silva':
                    parts = header.split(' ', 1)
                    if len(parts) == 2:
                        final_id, taxonomy = parts
                    else:
                        continue

                if final_id:
                    f_out.write(f"{final_id}\t{taxonomy}\n")
                    count += 1

        print(f"Taxonomy map created: {tax_map_file.name} ({count} entries).")
        return tax_map_file

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
            count = _run_fasta_parsing(fasta_path, output_map_path, format, debug)
            print(f"Forced map created at: {output_map_path} ({count} entries)")
            return output_map_path
        except Exception as e:
            print(f"CRITICAL ERROR: Cannot create forced map: {e}")
            sys.exit(1)

    if db_map_path.exists():
        print(f"Taxonomy map found (Cache DB): {db_map_path}")
        return db_map_path

    if output_map_path.exists():
        print(f"Taxonomy map found (Cache Output): {output_map_path}")
        return output_map_path

    print("Taxonomy map not found in cache. Creating (may take time)...")

    try:
        print(f"Attempting to write to: {db_map_path}")
        count = _run_fasta_parsing(fasta_path, db_map_path, format, debug)
        print(f"Map created in DB directory ({count} entries).")
        return db_map_path
    except (OSError, PermissionError) as e:
        print(f"  [Info] DB folder not writable ({e}). Writing to output folder.")
    except Exception as e:
        print(f"CRITICAL ERROR during parsing: {e}")
        sys.exit(1)

    try:
        print(f"Attempting to write to: {output_map_path}")
        count = _run_fasta_parsing(fasta_path, output_map_path, format, debug)
        print(f"Map created in output directory ({count} entries).")
        return output_map_path
    except Exception as e:
        print(f"CRITICAL ERROR: Cannot create taxonomy map in any folder: {e}")
        sys.exit(1)
