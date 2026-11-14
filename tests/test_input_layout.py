"""The input folder may hold one sub-folder per sample or one FASTQ file per sample."""
import gzip
import shutil

import pytest

import OtuUtils
from helpers import run_cli, write_fastq


def flatten(dataset, names=('SRR20631366.fastq', 'SRR20631367.fastq.gz')):
    """Turn the per-sample folders of the dataset into FASTQ files placed directly in the folder."""
    input_dir = dataset['input_dir']
    for sample, name in zip(dataset['samples'], names):
        source = input_dir / sample / 'reads.fastq'
        if name.endswith('.gz'):
            (input_dir / name).write_bytes(gzip.compress(source.read_bytes()))
        else:
            shutil.copy(source, input_dir / name)
        shutil.rmtree(input_dir / sample)
    return input_dir


def table_header(out_dir):
    return (out_dir / 'OTU_Table_out_reference.tsv').read_text().splitlines()[0].split('\t')


def test_sub_folders_are_found_as_samples(dataset):
    found = OtuUtils.find_input_samples(dataset['input_dir'])
    assert found['layout'] == 'folders'
    assert [source.name for source in found['sources']] == ['barcode01', 'barcode02']
    assert found['ignored'] == []


def test_fastq_files_placed_in_the_folder_are_found_as_samples(dataset):
    input_dir = flatten(dataset)
    found = OtuUtils.find_input_samples(input_dir)
    assert found['layout'] == 'files'
    assert [source.name for source in found['sources']] == ['SRR20631366.fastq', 'SRR20631367.fastq.gz']


def test_loose_files_next_to_sub_folders_are_listed_as_ignored(dataset):
    write_fastq(dataset['input_dir'] / 'extra.fastq', [('r1', 'ACGT' * 300)])
    found = OtuUtils.find_input_samples(dataset['input_dir'])
    assert found['layout'] == 'folders'
    assert [path.name for path in found['ignored']] == ['extra.fastq']


def test_hidden_and_resource_files_are_not_samples(dataset):
    input_dir = flatten(dataset)
    write_fastq(input_dir / '._SRR20631366.fastq', [('r1', 'ACGT' * 300)])
    (input_dir / 'notes.txt').write_text('not a read file')
    found = OtuUtils.find_input_samples(input_dir)
    assert [source.name for source in found['sources']] == ['SRR20631366.fastq', 'SRR20631367.fastq.gz']


def test_sample_name_of_a_file_drops_the_read_file_extensions():
    assert OtuUtils.sample_name_of(OtuUtils.Path('SRR20631367.fastq.gz')) == 'SRR20631367'
    assert OtuUtils.sample_name_of(OtuUtils.Path('run.1.fq')) == 'run.1'


@pytest.mark.real_tools
def test_flat_folder_gives_one_column_per_file(tmp_path, dataset):
    input_dir = flatten(dataset)
    out_dir = tmp_path / 'out'
    returncode, output, hung = run_cli(
        ['-i', input_dir, '-o', out_dir, '-d', dataset['db_path'], '-T', '2', '-t', '1'], timeout=120)
    assert not hung
    assert returncode == 0, output
    header = table_header(out_dir)
    assert 'SRR20631366' in header and 'SRR20631367' in header
    rows = (out_dir / 'OTU_Table_out_reference.tsv').read_text().splitlines()
    assert len(rows) == 4


@pytest.mark.real_tools
def test_flat_folder_and_sub_folders_give_the_same_counts(tmp_path, dataset):
    import pandas as pd
    arguments = ['-d', dataset['db_path'], '-T', '2', '-t', '1']
    returncode, output, hung = run_cli(['-i', dataset['input_dir'], '-o', tmp_path / 'folders'] + arguments,
                                       timeout=120)
    assert returncode == 0, output
    input_dir = flatten(dataset, names=('barcode01.fastq', 'barcode02.fastq'))
    returncode, output, hung = run_cli(['-i', input_dir, '-o', tmp_path / 'files'] + arguments, timeout=120)
    assert returncode == 0, output
    by_folder = pd.read_csv(tmp_path / 'folders' / 'OTU_Table_folders_reference.tsv', sep='\t', index_col=0)
    by_file = pd.read_csv(tmp_path / 'files' / 'OTU_Table_files_reference.tsv', sep='\t', index_col=0)
    assert by_folder[['barcode01', 'barcode02']].equals(by_file[['barcode01', 'barcode02']])


@pytest.mark.real_tools
def test_loose_files_next_to_sub_folders_are_reported(tmp_path, dataset):
    write_fastq(dataset['input_dir'] / 'extra.fastq', [('r1', 'ACGT' * 300)])
    returncode, output, hung = run_cli(
        ['-i', dataset['input_dir'], '-o', tmp_path / 'out', '-d', dataset['db_path'], '-T', '2', '-t', '1'],
        timeout=120)
    assert returncode == 0, output
    assert 'extra.fastq' in output
    assert 'ignored' in output.lower()


def test_folder_without_reads_explains_the_two_layouts(tmp_path, dataset):
    empty = tmp_path / 'empty'
    empty.mkdir()
    returncode, output, hung = run_cli(
        ['-i', empty, '-o', tmp_path / 'out', '-d', dataset['db_path']], timeout=120)
    assert returncode == 1
    assert 'one sub-folder per sample' in output
    assert 'one FASTQ file per sample' in output


def test_files_that_would_share_a_sample_name_are_rejected(tmp_path, dataset):
    input_dir = flatten(dataset, names=('same.fastq', 'same.fastq.gz'))
    returncode, output, hung = run_cli(
        ['-i', input_dir, '-o', tmp_path / 'out', '-d', dataset['db_path']], timeout=120)
    assert returncode == 1
    assert 'same sample name' in output
    assert 'same.fastq' in output
