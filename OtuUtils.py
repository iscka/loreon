import statistics
import random
import re
import signal
import subprocess
import tempfile
import threading
import time
from pathlib import Path
import sys
import multiprocessing
from functools import partial

try:
    from isal import igzip as gzip
except ImportError:
    import gzip


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


LENGTH_PRESETS = {
    'its1': (150, 350),
    'its2': (150, 400),
    'its-full': (400, 900),
    'its-lsu': (900, 1800),
    'rrna-operon': (3000, 6500),
    '16s-v3v4': (380, 530),
    '16s-full': (1200, 1700),
    '18s-v9': (100, 220),
    '18s-v4': (350, 550),
}

LENGTH_PRESET_LABELS = {
    'its1': 'Fungal ITS1',
    'its2': 'Fungal ITS2',
    'its-full': 'Full-length ITS',
    'its-lsu': 'ITS-LSU',
    'rrna-operon': 'rRNA operon',
    '16s-v3v4': 'Bacterial 16S V3-V4',
    '16s-full': 'Full-length 16S',
    '18s-v9': 'Eukaryotic 18S V9',
    '18s-v4': 'Eukaryotic 18S V4',
}

NO_UPPER_LIMIT = 999_999
_LENGTH_SAMPLE_SIZE = 500
_HISTOGRAM_BIN_WIDTHS = (10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000)


def length_filter_is_active(min_len, max_len):
    """True when the length window excludes reads of some length."""
    return min_len > 0 or 0 < max_len < NO_UPPER_LIMIT


def _weighted_lengths(lengths, weights=None):
    if weights is None:
        pairs = sorted((int(length), 1.0) for length in lengths)
    else:
        pairs = sorted((int(length), float(weight)) for length, weight in zip(lengths, weights))
    return pairs, sum(weight for _, weight in pairs)


def _quantile(pairs, total, fraction):
    target = fraction * total
    cumulative = 0.0
    for length, weight in pairs:
        cumulative += weight
        if cumulative >= target:
            return length
    return pairs[-1][0]


def _share_within(pairs, total, low, high):
    inside = sum(weight for length, weight in pairs if low <= length <= high)
    return inside / total * 100 if total > 0 else 0.0


def summarize_lengths(lengths, weights=None):
    """Summarise read lengths: count, minimum, 5th percentile, median, 95th percentile, maximum.

    weights gives the number of reads each sampled length stands for.
    Returns None when no length is given.
    """
    pairs, total = _weighted_lengths(lengths, weights)
    if not pairs or total <= 0:
        return None
    return {
        'count': int(round(total)),
        'min': pairs[0][0],
        'p05': _quantile(pairs, total, 0.05),
        'median': _quantile(pairs, total, 0.5),
        'p95': _quantile(pairs, total, 0.95),
        'max': pairs[-1][0],
    }


def suggest_length_window(summary):
    """Window, rounded outwards to 50 bp, that holds the central 90 % of the reads."""
    low = (summary['p05'] // 50) * 50
    high = -(-summary['p95'] // 50) * 50
    return low, high


def preset_shares(lengths, weights=None):
    """Percentage of reads whose length falls inside each preset window."""
    pairs, total = _weighted_lengths(lengths, weights)
    return {name: _share_within(pairs, total, low, high)
            for name, (low, high) in LENGTH_PRESETS.items()}


def _describe_window(min_len, max_len):
    if 0 < max_len < NO_UPPER_LIMIT:
        return f"{min_len}-{max_len} bp"
    return f"{min_len} bp or longer"


def length_filter_advice(lengths, weights, total, kept, min_len, max_len):
    """Lines to print about the outcome of the read-length filter.

    lengths (with optional weights) is a sample of the read lengths before
    filtering; total and kept are the read counts before and after it.
    A warning is returned when the filter discards more than half of the
    reads, a note when no filter is active, and nothing otherwise.
    """
    summary = summarize_lengths(lengths, weights) if lengths else None
    if total <= 0 or summary is None:
        return []
    observed = (f"  Observed read lengths: median {summary['median']} bp; "
                f"90 % of the reads between {summary['p05']} and {summary['p95']} bp.")
    if not length_filter_is_active(min_len, max_len):
        return [
            f"NOTE: no read-length filter is active: all {total} reads are used.",
            observed,
            "  If off-target products of a different length are present, restrict the analysis",
            "  to the expected amplicon length (-min/-max, --preset, or the length menu of the GUI).",
        ]
    discarded = total - kept
    if discarded * 2 <= total:
        return []
    low, high = suggest_length_window(summary)
    lines = [
        f"WARNING: the read-length filter ({_describe_window(min_len, max_len)}) discarded "
        f"{discarded} of {total} reads ({discarded / total * 100:.1f} %).",
        observed,
        "  If the window does not match your amplicon, change it:",
        f"    a window of {low}-{high} bp would keep the central 90 % of the reads;",
    ]
    shares = preset_shares(lengths, weights)
    matching = [name for name in LENGTH_PRESETS if shares[name] >= 50.0]
    matching.sort(key=lambda name: -shares[name])
    for name in matching[:3]:
        preset_low, preset_high = LENGTH_PRESETS[name]
        lines.append(f"    preset '{name}' ({preset_low}-{preset_high} bp) would keep "
                     f"{shares[name]:.1f} % of the reads.")
    lines.append(f"  If the window is right, most reads of this run are products of a different "
                 f"length and the analysis uses the remaining {kept} reads.")
    return lines


def format_length_report(lengths, weights=None):
    """Text report of a read-length distribution: summary, histogram and share per preset."""
    summary = summarize_lengths(lengths, weights)
    if summary is None:
        return "No reads found."
    pairs, total = _weighted_lengths(lengths, weights)
    low, high = suggest_length_window(summary)
    lines = [
        f"Reads examined: {summary['count']}",
        "",
        "Read length (bp):",
        f"  minimum          {summary['min']:>7}",
        f"  5th percentile   {summary['p05']:>7}",
        f"  median           {summary['median']:>7}",
        f"  95th percentile  {summary['p95']:>7}",
        f"  maximum          {summary['max']:>7}",
        f"Window holding the central 90 % of the reads: {low}-{high} bp",
        "",
        "Length distribution:",
    ]
    first = _quantile(pairs, total, 0.005)
    last = _quantile(pairs, total, 0.995)
    width = next((w for w in _HISTOGRAM_BIN_WIDTHS if w * 24 >= last - first + 1), _HISTOGRAM_BIN_WIDTHS[-1])
    start = (first // width) * width
    stop = (last // width + 1) * width
    bins = {}
    below = above = 0.0
    for length, weight in pairs:
        if length < start:
            below += weight
        elif length >= stop:
            above += weight
        else:
            index = (length - start) // width
            bins[index] = bins.get(index, 0.0) + weight
    peak = max(list(bins.values()) + [below, above, 1e-9])

    def _row(label, weight):
        bar = '#' * max(1, int(round(weight / peak * 30))) if weight > 0 else ''
        return f"  {label:>13} | {bar} {weight / total * 100:.1f} %"

    if below > 0:
        lines.append(_row(f"< {start}", below))
    gap = False
    for index in range((stop - start) // width):
        weight = bins.get(index, 0.0)
        if weight <= 0:
            if not gap:
                lines.append(f"  {'...':>13} |")
                gap = True
            continue
        gap = False
        lines.append(_row(f"{start + index * width}-{start + (index + 1) * width - 1}", weight))
    if above > 0:
        lines.append(_row(f">= {stop}", above))
    lines.append("")
    lines.append("Share of the reads inside each preset window:")
    for name, (preset_low, preset_high) in LENGTH_PRESETS.items():
        share = _share_within(pairs, total, preset_low, preset_high)
        lines.append(f"  {name:<12} {preset_low:>5}-{preset_high:<5} bp  {share:5.1f} %")
    return '\n'.join(lines)


def sample_read_lengths(input_dir, reads_per_file=1000, max_files=200):
    """Lengths of the first reads of the FASTQ files below a folder.

    At most max_files files, evenly spread over the sorted file list, are
    opened and at most reads_per_file reads are taken from each.
    """
    input_dir = Path(input_dir)
    files = sorted(
        f for pattern in ('**/*.f*q.gz', '**/*.f*q') for f in input_dir.glob(pattern)
        if not f.name.startswith('._')
    )
    if len(files) > max_files:
        step = len(files) / max_files
        files = [files[int(i * step)] for i in range(max_files)]
    lengths = []
    for file_path in files:
        opener = gzip.open if file_path.name.endswith('.gz') else open
        try:
            with opener(file_path, 'rt', encoding='utf-8', errors='replace') as handle:
                for count, record in enumerate(_iter_fastq_raw(handle)):
                    if count >= reads_per_file:
                        break
                    lengths.append(record[4])
        except (OSError, EOFError):
            continue
    return lengths


def _thin_length_sample(reservoir):
    """Evenly spaced order statistics of a length sample, at most _LENGTH_SAMPLE_SIZE values."""
    ordered = sorted(reservoir)
    if len(ordered) <= _LENGTH_SAMPLE_SIZE:
        return ordered
    step = (len(ordered) - 1) / (_LENGTH_SAMPLE_SIZE - 1)
    return [ordered[int(round(i * step))] for i in range(_LENGTH_SAMPLE_SIZE)]


def filter_and_merge_directory(
    input_dir: Path,
    min_len: int,
    max_len: int,
    out_dir: Path,
    debug: bool = False,
    gzip_intermediate: bool = False
):
    dir_name = input_dir.name

    if max_len <= 0:
        max_len = 999_999_999

    file_list = list(input_dir.glob('**/*.f*q.gz')) + \
                list(input_dir.glob('**/*.f*q'))

    file_list = [f for f in file_list if not f.name.startswith('._')]

    if not file_list:
        if debug:
            print(f"  [Filter {dir_name}] No source files found. Skipping.")
        return None

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
        print(f"Error setting up paths for '{dir_name}': {e}", flush=True)
        return {'barcode': dir_name, 'failed': True, 'error': f"{type(e).__name__}: {e}"}

    total_count, good_count, bad_count, sum_of_lengths = 0, 0, 0, 0
    reservoir = []
    reservoir_n = 0
    file_path = None

    def _open_out(path):
        if gzip_intermediate:
            return gzip.open(path, 'wt', encoding='utf-8', compresslevel=3)
        return open(path, 'wt', encoding='utf-8', buffering=1_048_576)

    try:
        with _open_out(good_file_path) as f_good, \
             _open_out(bad_file_path) as f_bad:

            for file_path in file_list:
                if debug:
                    print(f"    -> Reading {file_path.name}...")

                in_is_gzipped = file_path.name.endswith('.gz')
                stem_name = file_path.name.removesuffix('.gz')
                is_fastq = stem_name.endswith(('.fastq', '.fq'))

                read_open_func = gzip.open if in_is_gzipped else open

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
        source = file_path.name if file_path is not None else dir_name
        print(f"ERROR during processing of {source}: {e}", flush=True)
        _remove_quietly(good_file_path)
        _remove_quietly(bad_file_path)
        return {'barcode': dir_name, 'failed': True, 'error': f"{type(e).__name__}: {e}"}

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
        'length_sample': _thin_length_sample(reservoir),
    }


_STDERR_TAIL_LINES = 15
_BROKEN_PIPE_EXIT = -getattr(signal, 'SIGPIPE', 13)
_CIGAR_OPERATION = re.compile(rb'(\d+)([MIDNSHP=X])')
_SOFT_CLIP_HEAD = re.compile(rb'(?:\d+H)?(\d+)S')
_SOFT_CLIP_TAIL = re.compile(rb'(\d+)S(?:\d+H)?$')


class _TimeLimit:
    """Kill a group of processes when a time limit expires."""

    def __init__(self, seconds, processes):
        self.expired = False
        self.seconds = seconds
        self._processes = processes
        self._timer = None
        if seconds and seconds > 0:
            self._timer = threading.Timer(seconds, self._expire)
            self._timer.daemon = True
            self._timer.start()

    def _expire(self):
        if all(process.poll() is not None for process in self._processes):
            return
        self.expired = True
        _kill_all(self._processes)

    def cancel(self):
        if self._timer is not None:
            self._timer.cancel()


def _kill_all(processes):
    for process in list(processes):
        if process.poll() is None:
            try:
                process.kill()
            except OSError:
                pass


def _close_quietly(stream):
    if stream is None:
        return
    try:
        stream.close()
    except OSError:
        pass


def _read_tail(handle, max_lines=_STDERR_TAIL_LINES):
    """Return the last non-empty lines written to a temporary file."""
    try:
        handle.seek(0, 2)
        size = handle.tell()
        handle.seek(max(0, size - 65536))
        text = handle.read().decode('utf-8', errors='replace')
    except (OSError, ValueError):
        return []
    lines = [line for line in text.splitlines() if line.strip()]
    return lines[-max_lines:]


def _describe_exit(tool, returncode):
    killed_hint = " (killed: the computer or the Docker container may have run out of memory)"
    if returncode < 0:
        text = f"{tool} was terminated by signal {-returncode}"
        return text + killed_hint if returncode == -9 else text
    text = f"{tool} failed with exit status {returncode}"
    return text + killed_hint if returncode == 137 else text


def _report_failures(label, failures):
    """Print which tools of a pipeline failed, with the end of their error output."""
    for tool, returncode, command, error_lines in failures:
        print(f"ERROR [{label}]: {_describe_exit(tool, returncode)}.", flush=True)
        print(f"  Command: {' '.join(command)}", flush=True)
        if error_lines:
            print("  Last lines of its error output:", flush=True)
            for line in error_lines:
                print(f"    {line}", flush=True)


def _report_time_limit(label, step, limit):
    print(f"ERROR [{label}]: {step} exceeded the time limit of {limit.seconds:g} s per sample "
          f"and was stopped.", flush=True)
    print("  Raise or disable the limit (--sample-timeout) if the sample is expected to take longer.",
          flush=True)


def _remove_quietly(path):
    try:
        Path(path).unlink()
    except OSError:
        pass


def _remove_partial_bam(bam_path):
    """Delete an incomplete BAM file and the temporary files samtools sort left next to it."""
    bam_path = Path(bam_path)
    _remove_quietly(bam_path)
    try:
        leftovers = list(bam_path.parent.glob(bam_path.name + '.tmp.*'))
    except OSError:
        leftovers = []
    for leftover in leftovers:
        _remove_quietly(leftover)


def mapping_improved(fastq_file: str, db_path: str, threads: int, output_dir: str,
                     kmer_size: int = 15, window_size: int = 10, debug: bool = False,
                     timeout_s: float = None):
    """Map one FASTQ file with minimap2 and write a coordinate-sorted BAM.

    minimap2 is piped into samtools sort, so no SAM file is written to disk.
    The error output of both tools goes to temporary files and the exit
    status of both is checked. Returns the path of the BAM file, or None
    when a tool fails, is missing, or exceeds timeout_s seconds.
    """
    try:
        fastq_path = Path(fastq_file)
        db_path_obj = Path(db_path)
        output_dir_path = Path(output_dir)

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

    label = fastq_path.name

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

    processes = []
    limit = None
    try:
        with tempfile.TemporaryFile() as mm2_err, tempfile.TemporaryFile() as sort_err:
            try:
                proc_mm2 = subprocess.Popen(cmd_mm2, stdout=subprocess.PIPE, stderr=mm2_err)
                processes.append(proc_mm2)
                proc_sort = subprocess.Popen(cmd_sort, stdin=proc_mm2.stdout,
                                             stdout=subprocess.DEVNULL, stderr=sort_err)
                processes.append(proc_sort)
                proc_mm2.stdout.close()
                limit = _TimeLimit(timeout_s, processes)
                sort_rc = proc_sort.wait()
                mm2_rc = proc_mm2.wait()
            finally:
                if limit is not None:
                    limit.cancel()
                _kill_all(processes)
                for process in processes:
                    process.wait()
                    _close_quietly(process.stdout)

            if limit.expired:
                _report_time_limit(label, 'mapping', limit)
                _remove_partial_bam(bam_output_path)
                return None

            failures = []
            if mm2_rc != 0 and not (mm2_rc == _BROKEN_PIPE_EXIT and sort_rc != 0):
                failures.append(('minimap2', mm2_rc, cmd_mm2, _read_tail(mm2_err)))
            if sort_rc != 0:
                failures.append(('samtools sort', sort_rc, cmd_sort, _read_tail(sort_err)))
            if failures:
                _report_failures(label, failures)
                _remove_partial_bam(bam_output_path)
                return None

            if not bam_output_path.exists():
                print(f"ERROR [{label}]: mapping ended without errors but {bam_output_path.name} "
                      f"was not written.", flush=True)
                return None

            if debug:
                print(f"[{label}] Mapping + sort complete (BAM pipe).", flush=True)
                mm2_lines = _read_tail(mm2_err, max_lines=40)
                if mm2_lines:
                    print("  Info (from minimap2 stderr):", flush=True)
                    for line in mm2_lines:
                        print(f"  {line}", flush=True)

        return str(bam_output_path)

    except FileNotFoundError as e:
        print(f"ERROR [{label}]: command not found: {e.filename}.", flush=True)
        print("  Ensure 'minimap2' and 'samtools' are installed and on your PATH.", flush=True)
        return None

    except Exception as e:
        print(f"ERROR [{label}]: unexpected error during mapping: {e}", flush=True)
        _remove_partial_bam(bam_output_path)
        return None


def _query_length_from_cigar(cigar):
    return sum(int(count) for count, operation in _CIGAR_OPERATION.findall(cigar)
               if operation in b'MIS=X')


def _edit_distance(tags):
    """Return the value of the NM tag from the optional fields of a SAM record, or None."""
    if tags.startswith(b'NM:i:'):
        start = 5
    else:
        start = tags.find(b'\tNM:i:')
        if start < 0:
            return None
        start += 6
    end = tags.find(b'\t', start)
    try:
        return int(tags[start:end] if end >= 0 else tags[start:])
    except ValueError:
        return None


def aligned_query_length(cigar, query_length):
    """Number of read bases inside the alignment: read length minus soft-clipped bases."""
    clipped = 0
    head = _SOFT_CLIP_HEAD.match(cigar)
    if head:
        clipped += int(head.group(1))
    if not (head and head.end() == len(cigar)):
        tail = _SOFT_CLIP_TAIL.search(cigar[-24:])
        if tail:
            clipped += int(tail.group(1))
    return query_length - clipped


def passes_quality_filter(aligned_length, query_length, edit_distance,
                          min_percent_identity, min_ref_coverage):
    """Apply the identity and read-coverage thresholds to one alignment.

    aligned_length is the read length without soft-clipped bases.
    Identity is (aligned_length - NM) / aligned_length; read coverage is
    aligned_length / read length. Both thresholds are percentages, the
    comparison is inclusive, and a threshold of 0 is not applied.
    """
    if aligned_length <= 0:
        return False
    if min_percent_identity > 0.0:
        if edit_distance is None:
            return False
        if (aligned_length - edit_distance) * 100 < min_percent_identity * aligned_length:
            return False
    if min_ref_coverage > 0.0 and aligned_length * 100 < min_ref_coverage * query_length:
        return False
    return True


def tabulate_sam_stream(lines, min_percent_identity=0.0, min_ref_coverage=0.0, sink=None):
    """Count the primary alignments of a SAM stream per reference.

    lines yields SAM lines as bytes, header first. Secondary and
    supplementary records are skipped. Mapped reads that pass the identity
    and read-coverage thresholds are counted; when sink is given, the header
    and the records that pass are written to it unchanged.

    Returns a dict with the references in header order, the counts per
    reference, and the numbers of primary, mapped and passing reads.
    """
    references = []
    counts = {}
    total = mapped = passed = 0
    filtering = min_percent_identity > 0.0 or min_ref_coverage > 0.0
    write = sink.write if sink is not None else None
    in_header = True

    for line in lines:
        if in_header:
            if line.startswith(b'@'):
                if line.startswith(b'@SQ\t'):
                    name = None
                    length = 0
                    for field in line.rstrip(b'\r\n').split(b'\t')[1:]:
                        if field.startswith(b'SN:'):
                            name = field[3:]
                        elif field.startswith(b'LN:'):
                            length = int(field[3:])
                    if name is not None:
                        references.append((name, length))
                if write is not None:
                    write(line)
                continue
            in_header = False

        fields = line.split(b'\t', 11)
        if len(fields) < 11:
            continue
        flag = int(fields[1])
        if flag & 0x900:
            continue
        total += 1
        if flag & 0x4:
            continue
        mapped += 1

        if filtering:
            cigar = fields[5]
            sequence = fields[9]
            query_length = len(sequence) if sequence != b'*' else _query_length_from_cigar(cigar)
            edit_distance = _edit_distance(fields[11]) if len(fields) > 11 else None
            if not passes_quality_filter(aligned_query_length(cigar, query_length), query_length,
                                         edit_distance, min_percent_identity, min_ref_coverage):
                continue

        passed += 1
        reference = fields[2]
        counts[reference] = counts.get(reference, 0) + 1
        if write is not None:
            write(line)

    def _text(value):
        return value.decode('utf-8', errors='surrogateescape')

    return {
        'references': [(_text(name), length) for name, length in references],
        'counts': {_text(name): count for name, count in counts.items()},
        'total_reads': total,
        'mapped_reads': mapped,
        'passed_reads': passed,
    }


def tabeling_improved(samfile: str, threads: int = 4, debug: bool = False,
                      min_percent_identity: float = 0.0, min_ref_coverage: float = 0.0,
                      delete_bam: bool = False, timeout_s: float = None):
    """Count the primary alignments of one sample per reference.

    The sorted BAM written by mapping_improved is read once. Unmapped,
    secondary and supplementary records are left out; the identity and
    read-coverage thresholds are evaluated on the CIGAR string and the NM
    tag of each alignment, so the result does not depend on the
    filter-expression support of the installed samtools. The alignments
    that pass are also written to <sample>_sorted_filtered.bam.

    Returns the mapping statistics of the sample, or None when a tool
    fails, is missing, or exceeds timeout_s seconds.
    """
    fname = "unknown"
    filtered_sorted_bam = None
    input_bam = samfile

    try:
        input_path = Path(input_bam)
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

    if min_percent_identity > 0.0 or min_ref_coverage > 0.0:
        print(f"[{fname}] Quality filter: identity >= {min_percent_identity}%, "
              f"read coverage >= {min_ref_coverage}%", flush=True)

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

    cmd_read = ['samtools', 'view', f'-@{thread_str}', '-h', '-F', '0x900', str(input_bam)]
    cmd_write = ['samtools', 'view', f'-@{thread_str}', '-b', '-o', str(filtered_sorted_bam), '-']
    cmd_index = ['samtools', 'index', f'-@{thread_str}', str(filtered_sorted_bam)]

    if debug:
        print(f"[{fname}] CMD (read): {' '.join(cmd_read)}", flush=True)
        print(f"[{fname}] CMD (write): {' '.join(cmd_write)}", flush=True)

    processes = []
    limit = None
    completed = False
    started = time.monotonic()
    try:
        with tempfile.TemporaryFile() as read_err, tempfile.TemporaryFile() as write_err:
            counted = None
            stream_error = None
            try:
                proc_read = subprocess.Popen(cmd_read, stdout=subprocess.PIPE, stderr=read_err,
                                             bufsize=1_048_576)
                processes.append(proc_read)
                proc_write = subprocess.Popen(cmd_write, stdin=subprocess.PIPE,
                                              stdout=subprocess.DEVNULL, stderr=write_err,
                                              bufsize=1_048_576)
                processes.append(proc_write)
                limit = _TimeLimit(timeout_s, processes)
                try:
                    counted = tabulate_sam_stream(proc_read.stdout, min_percent_identity,
                                                  min_ref_coverage, sink=proc_write.stdin)
                    proc_write.stdin.close()
                except BrokenPipeError:
                    counted = None
                except ValueError as e:
                    counted = None
                    stream_error = e
                    _kill_all(processes)
                _close_quietly(proc_write.stdin)
                _close_quietly(proc_read.stdout)
                read_rc = proc_read.wait()
                write_rc = proc_write.wait()
            finally:
                if limit is not None:
                    limit.cancel()
                _kill_all(processes)
                for process in processes:
                    process.wait()
                    _close_quietly(process.stdin)
                    _close_quietly(process.stdout)

            if limit.expired:
                _report_time_limit(fname, 'counting the alignments', limit)
                return None

            if stream_error is not None:
                print(f"ERROR [{fname}]: the alignment records could not be read: {stream_error}",
                      flush=True)
                return None

            failures = []
            if read_rc != 0 and not (read_rc == _BROKEN_PIPE_EXIT and write_rc != 0):
                failures.append(('samtools view', read_rc, cmd_read, _read_tail(read_err)))
            if write_rc != 0:
                failures.append(('samtools view', write_rc, cmd_write, _read_tail(write_err)))
            if failures or counted is None:
                if failures:
                    _report_failures(fname, failures)
                else:
                    print(f"ERROR [{fname}]: writing the filtered alignments was interrupted.",
                          flush=True)
                return None

        index_timeout = None
        if timeout_s and timeout_s > 0:
            index_timeout = max(1.0, timeout_s - (time.monotonic() - started))
        try:
            index_run = subprocess.run(cmd_index, capture_output=True, timeout=index_timeout)
        except subprocess.TimeoutExpired:
            _report_time_limit(fname, 'indexing the alignments', limit)
            return None
        if index_run.returncode != 0:
            error_lines = index_run.stderr.decode('utf-8', errors='replace').splitlines()
            _report_failures(fname, [('samtools index', index_run.returncode, cmd_index,
                                      [line for line in error_lines if line.strip()][-_STDERR_TAIL_LINES:])])
            return None

        mapping_stats['total_reads'] = counted['total_reads']
        mapping_stats['mapped_unfiltered'] = counted['mapped_reads']
        if counted['total_reads'] > 0:
            mapping_stats['mapping_rate_pct'] = round(
                counted['mapped_reads'] / counted['total_reads'] * 100, 2
            )
        unmapped = counted['total_reads'] - counted['mapped_reads']

        if debug:
            print(f"[{fname}] Pre-filter: total={counted['total_reads']}, "
                  f"mapped={counted['mapped_reads']}", flush=True)

        total_mapped_filtered = 0
        per_ref_counts = []
        with open(results_file, 'w', encoding='utf-8', errors='surrogateescape') as f:
            for reference, length in counted['references']:
                ref_count = counted['counts'].get(reference, 0)
                if ref_count:
                    total_mapped_filtered += ref_count
                    per_ref_counts.append((reference, ref_count))
                    f.write(f"{reference}\t{length}\t{ref_count}\t0\n")
            f.write(f"*\t0\t0\t{unmapped}\n")
            f.write("*\t0\t0\t0\n")

        mapping_stats['mapped_filtered'] = total_mapped_filtered
        if counted['mapped_reads'] > 0:
            mapping_stats['quality_retention_pct'] = round(
                total_mapped_filtered / counted['mapped_reads'] * 100, 2
            )

        mapping_stats['references_detected'] = len(per_ref_counts)
        if per_ref_counts:
            counts_only = sorted((c for _, c in per_ref_counts), reverse=True)
            top = sorted(per_ref_counts, key=lambda x: -x[1])[:10]
            n = len(counts_only)

            def _pct(p):
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
            print(f"[{fname}] Post-filter: mapped={total_mapped_filtered}", flush=True)
            print(f"--- [Tabeling: {fname}] Complete ---", flush=True)

        completed = True
        return mapping_stats

    except FileNotFoundError as e:
        print(f"ERROR [{fname}]: command not found: {e.filename}.", flush=True)
        print("  Ensure 'samtools' is installed and on your PATH.", flush=True)
        return None

    except Exception as e:
        print(f"ERROR [{fname}]: unexpected error while counting the alignments: {e}", flush=True)
        return None

    finally:
        if not completed:
            _remove_quietly(results_file)
        if not completed and filtered_sorted_bam is not None and not delete_bam:
            _remove_quietly(filtered_sorted_bam)
            _remove_quietly(f"{filtered_sorted_bam}.bai")
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


def _parse_fasta_headers(fasta_path: Path, tax_map_file: Path, format: str, debug: bool):
    """Parse FASTA headers (without reading sequences) and write taxonomy map.
    Returns the number of entries written. Replaces both _create_tax_map and
    _run_fasta_parsing with a single, faster implementation."""
    rank_re = re.compile(r'[kpcofgs]__') if format in ('unite', 'eukariome', 'cbs') else None

    try:
        count = 0
        with open(fasta_path, 'r', encoding='utf-8') as f_in, \
             open(tax_map_file, 'w', encoding='utf-8') as f_out:
            f_out.write("OTU_ID\tTaxonomy\n")

            for line in f_in:
                if not line.startswith('>'):
                    continue

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
