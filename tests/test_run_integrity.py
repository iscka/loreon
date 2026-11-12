"""A run must report every sample it could not process and must not reuse results of an earlier run."""
import gzip
import os
import shutil
import subprocess
import time

import pytest

import OtuUtils
from helpers import SYSTEM_PATH, build_dataset, make_tool, run_cli


def table_header(out_dir):
    return (out_dir / 'OTU_Table_out_reference.tsv').read_text().splitlines()[0].split('\t')


@pytest.mark.real_tools
def test_counts_of_an_earlier_run_do_not_enter_the_table(tmp_path, dataset):
    out_dir = tmp_path / 'out'
    arguments = ['-i', dataset['input_dir'], '-o', out_dir, '-d', dataset['db_path'], '-T', '2', '-t', '1']
    returncode, output, hung = run_cli(arguments, timeout=120)
    assert returncode == 0, output
    assert 'barcode02' in table_header(out_dir)
    shutil.rmtree(dataset['input_dir'] / 'barcode02')
    returncode, output, hung = run_cli(arguments, timeout=120)
    assert returncode == 0, output
    header = table_header(out_dir)
    assert 'barcode01' in header
    assert 'barcode02' not in header


@pytest.mark.real_tools
def test_run_fails_when_the_table_cannot_be_written(tmp_path, dataset):
    out_dir = tmp_path / 'out'
    (out_dir / 'OTU_Table_out_reference.xlsx').mkdir(parents=True)
    (out_dir / 'OTU_Table_out_reference.tsv').mkdir()
    returncode, output, hung = run_cli(
        ['-i', dataset['input_dir'], '-o', out_dir, '-d', dataset['db_path'], '-T', '2', '-t', '1'],
        timeout=120)
    assert not hung
    assert returncode == 1, output
    assert 'Pipeline Complete' not in output


@pytest.mark.real_tools
def test_sample_that_cannot_be_read_is_reported_as_failed(tmp_path, dataset):
    broken = dataset['input_dir'] / 'barcode03'
    broken.mkdir()
    records = ''.join(f"@r{i}\n{'ACGT' * 300}\n+\n{'I' * 1200}\n" for i in range(200)).encode()
    compressed = gzip.compress(records)
    (broken / 'reads.fastq.gz').write_bytes(compressed[: len(compressed) // 2])
    out_dir = tmp_path / 'out'
    returncode, output, hung = run_cli(
        ['-i', dataset['input_dir'], '-o', out_dir, '-d', dataset['db_path'], '-T', '2', '-t', '1'],
        timeout=120)
    assert not hung
    assert returncode == 3, output
    assert '1 of 3 samples failed' in output
    header = table_header(out_dir)
    assert 'barcode01' in header and 'barcode02' in header
    assert 'barcode03' not in header
    assert not (out_dir / 'good_seq' / 'barcode03.fastq').exists()


@pytest.mark.fake_tools
def test_run_ends_when_a_worker_process_dies(tmp_path):
    bin_dir = tmp_path / 'bin'
    make_tool(bin_dir, 'minimap2', """
        if [ "$1" = "--version" ]; then echo "2.28-r1209"; exit 0; fi
        kill -9 $PPID
        sleep 1
        exit 1
    """)
    make_tool(bin_dir, 'samtools', """
        if [ "$1" = "--version" ]; then echo "samtools 1.21"; exit 0; fi
        cat > /dev/null
        exit 0
    """)
    data = build_dataset(tmp_path / 'data')
    env = dict(os.environ, PATH=f"{bin_dir}:{SYSTEM_PATH}")
    returncode, output, hung = run_cli(
        ['-i', data['input_dir'], '-o', tmp_path / 'out', '-d', data['db_path'],
         '-T', '2', '-t', '1', '--progress-interval', '1'],
        timeout=40, env=env)
    assert not hung, 'the run never ended'
    assert returncode == 1
    assert 'ended unexpectedly' in output


def test_time_limit_does_not_flag_processes_that_already_ended():
    process = subprocess.Popen(['true'])
    process.wait()
    limit = OtuUtils._TimeLimit(0.2, [process])
    time.sleep(0.6)
    assert limit.expired is False
