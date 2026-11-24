"""Shared helpers for the LOREON test suite: synthetic data, stand-in tools, guarded execution."""
import os
import random
import re
import shutil
import signal
import subprocess
import sys
import textwrap
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SOURCES_DIR = REPO_ROOT / 'sources'
SCRIPTS_DIR = REPO_ROOT / 'scripts'
HAVE_TOOLS = shutil.which('minimap2') is not None and shutil.which('samtools') is not None
SYSTEM_PATH = '/usr/bin:/bin:/usr/sbin:/sbin'

TAXA = [
    ('Alpha_one', 'AC000001', 'SH0000001.10FU',
     'k__Fungi;p__Ascomycota;c__Classa;o__Ordera;f__Familya;g__Alpha;s__Alpha_one'),
    ('Beta_two', 'AC000002', 'SH0000002.10FU',
     'k__Fungi;p__Ascomycota;c__Classb;o__Orderb;f__Familyb;g__Beta;s__Beta_two'),
    ('Gamma_three', 'AC000003', 'SH0000003.10FU',
     'k__Fungi;p__Basidiomycota;c__Classc;o__Orderc;f__Familyc;g__Gamma;s__Gamma_three'),
]


def random_dna(rng, length):
    return ''.join(rng.choice('ACGT') for _ in range(length))


def mutate(rng, sequence, rate):
    """Return a copy of the sequence with substitutions, insertions and deletions at the given rate."""
    out = []
    for base in sequence:
        if rng.random() >= rate:
            out.append(base)
            continue
        kind = rng.random()
        if kind < 0.6:
            out.append(rng.choice([b for b in 'ACGT' if b != base]))
        elif kind < 0.8:
            out.append(base)
            out.append(rng.choice('ACGT'))
    return ''.join(out)


def write_fastq(path, reads):
    with open(path, 'w', encoding='utf-8') as handle:
        for name, sequence in reads:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def build_dataset(root, seed=7, samples=('barcode01', 'barcode02')):
    """Create a reference FASTA and one FASTQ folder per sample.

    Each sample holds reads that are close to a reference, reads that are
    divergent, reads that align over half of their length only, and reads
    of random sequence that do not map.
    """
    rng = random.Random(seed)
    root = Path(root)
    input_dir = root / 'input'
    db_dir = root / 'db'
    db_dir.mkdir(parents=True, exist_ok=True)
    references = [random_dna(rng, 1500) for _ in TAXA]
    db_path = db_dir / 'reference.fasta'
    with open(db_path, 'w', encoding='utf-8') as handle:
        for (name, accession, sh, taxonomy), sequence in zip(TAXA, references):
            handle.write(f">{name}|{accession}|{sh}|refs|{taxonomy}\n{sequence}\n")
    for sample in samples:
        reads = []
        for i in range(24):
            start = rng.randrange(0, 300)
            reads.append((f"{sample}_close_{i}", mutate(rng, references[i % 3][start:start + 1100], 0.02)))
        for i in range(12):
            start = rng.randrange(0, 300)
            reads.append((f"{sample}_divergent_{i}", mutate(rng, references[i % 3][start:start + 1100], 0.10)))
        for i in range(8):
            start = rng.randrange(0, 900)
            fragment = mutate(rng, references[i % 3][start:start + 500], 0.02)
            reads.append((f"{sample}_partial_{i}", fragment + random_dna(rng, 600)))
        for i in range(4):
            reads.append((f"{sample}_random_{i}", random_dna(rng, 1000)))
        sample_dir = input_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        write_fastq(sample_dir / 'reads.fastq', reads)
    return {'input_dir': input_dir, 'db_path': db_path, 'samples': list(samples)}


def make_tool(bin_dir, name, script):
    """Write an executable stand-in for an external tool."""
    bin_dir = Path(bin_dir)
    bin_dir.mkdir(parents=True, exist_ok=True)
    path = bin_dir / name
    path.write_text('#!/bin/bash\n' + textwrap.dedent(script), encoding='utf-8')
    path.chmod(0o755)
    return path


def run_guarded(args, timeout, env=None, cwd=None):
    """Run a command and return (returncode, output, hung).

    When the time limit is exceeded the whole process group is killed and
    hung is True, so a blocked pipeline cannot outlive the test.
    """
    process = subprocess.Popen(
        args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        env=env, cwd=cwd, start_new_session=True
    )
    try:
        output, _ = process.communicate(timeout=timeout)
        return process.returncode, output, False
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        output, _ = process.communicate()
        return None, output, True


def run_python(code, timeout, env=None):
    """Run a Python snippet with the application folder on the import path."""
    full_env = dict(os.environ if env is None else env)
    full_env['PYTHONPATH'] = str(SOURCES_DIR)
    return run_guarded([sys.executable, '-c', textwrap.dedent(code)], timeout=timeout, env=full_env)


def run_cli(arguments, timeout, env=None):
    """Run the command-line pipeline and return (returncode, output, hung)."""
    command = [sys.executable, str(SOURCES_DIR / 'metaGenomics_new.py')] + [str(a) for a in arguments]
    return run_guarded(command, timeout=timeout, env=env)


def read_counts(tmp_file):
    """Return {reference: count} from a per-sample count table, skipping the unmapped rows."""
    counts = {}
    with open(tmp_file, encoding='utf-8') as handle:
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if parts[0] == '*':
                continue
            counts[parts[0]] = int(parts[2])
    return counts


def expected_counts(bam_path, min_percent_identity, min_ref_coverage):
    """Count primary alignments per reference with an independent implementation of the filter."""
    output = subprocess.run(
        ['samtools', 'view', '-F', '0x904', str(bam_path)],
        check=True, capture_output=True, text=True
    ).stdout
    counts = {}
    for line in output.splitlines():
        fields = line.split('\t')
        operations = re.findall(r'(\d+)([MIDNSHP=X])', fields[5])
        query_length = sum(int(n) for n, op in operations if op in 'MIS=X')
        soft_clipped = sum(int(n) for n, op in operations if op == 'S')
        edit_distance = next(int(tag[5:]) for tag in fields[11:] if tag.startswith('NM:i:'))
        aligned = query_length - soft_clipped
        if (aligned - edit_distance) * 100 < min_percent_identity * aligned:
            continue
        if aligned * 100 < min_ref_coverage * query_length:
            continue
        counts[fields[2]] = counts.get(fields[2], 0) + 1
    return counts
