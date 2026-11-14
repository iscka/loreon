"""Example dataset shipped with LOREON and the comparison with an expected table."""
import os
import shutil
import sys

import pytest

import OtuUtils
import metaGenomics_new
from helpers import REPO_ROOT, run_cli, run_guarded

HEADER = "OTU_ID\tTaxonomy\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tDetails\tS1\tS2\n"


def write_table(path, rows, header=HEADER):
    path.write_text(header + ''.join(rows))
    return path


def row(otu, first, second):
    return f"{otu}\tFungi;A\tFungi\tA\t\t\t\t\t\t\t{first}\t{second}\n"


def test_example_folder_is_complete_and_small():
    assert sorted(p.name for p in OtuUtils.EXAMPLE_READS.iterdir()) == \
        ['MockA.fastq.gz', 'MockB.fastq.gz', 'MockC.fastq.gz']
    assert OtuUtils.EXAMPLE_DATABASE.is_file()
    assert OtuUtils.EXAMPLE_EXPECTED_TABLE.is_file()
    assert (OtuUtils.EXAMPLE_DIR / 'README.md').is_file()
    size = sum(p.stat().st_size for p in OtuUtils.EXAMPLE_DIR.rglob('*') if p.is_file())
    assert size < 3_000_000


def test_identical_tables_have_no_differences(tmp_path):
    rows = [row('A|1', 5, 0), row('B|2', 1, 7)]
    first = write_table(tmp_path / 'a.tsv', rows)
    second = write_table(tmp_path / 'b.tsv', rows)
    assert OtuUtils.compare_otu_tables(first, second) == []


def test_row_and_column_order_do_not_matter(tmp_path):
    first = write_table(tmp_path / 'a.tsv', [row('A|1', 5, 0), row('B|2', 1, 7)])
    swapped_header = HEADER.replace('S1\tS2', 'S2\tS1')
    second = write_table(tmp_path / 'b.tsv', [row('B|2', 7, 1), row('A|1', 0, 5)], header=swapped_header)
    assert OtuUtils.compare_otu_tables(first, second) == []


def test_changed_count_is_reported(tmp_path):
    first = write_table(tmp_path / 'a.tsv', [row('A|1', 5, 0), row('B|2', 1, 7)])
    second = write_table(tmp_path / 'b.tsv', [row('A|1', 5, 0), row('B|2', 1, 9)])
    differences = OtuUtils.compare_otu_tables(first, second)
    assert len(differences) == 1
    assert 'B|2' in differences[0] and 'S2' in differences[0]
    assert '7' in differences[0] and '9' in differences[0]


def test_missing_and_extra_rows_are_reported(tmp_path):
    first = write_table(tmp_path / 'a.tsv', [row('A|1', 5, 0), row('B|2', 1, 7)])
    second = write_table(tmp_path / 'b.tsv', [row('A|1', 5, 0), row('C|3', 2, 2)])
    text = '\n'.join(OtuUtils.compare_otu_tables(first, second))
    assert 'B|2' in text and 'missing' in text
    assert 'C|3' in text and 'not expected' in text


def test_missing_sample_is_reported(tmp_path):
    first = write_table(tmp_path / 'a.tsv', [row('A|1', 5, 0)])
    second = write_table(tmp_path / 'b.tsv', ["A|1\tFungi;A\tFungi\tA\t\t\t\t\t\t\t5\n"],
                         header=HEADER.replace('\tS2', ''))
    text = '\n'.join(OtuUtils.compare_otu_tables(first, second))
    assert 'S2' in text and 'missing' in text


def test_long_lists_of_differences_are_shortened(tmp_path):
    first = write_table(tmp_path / 'a.tsv', [row(f'X|{i}', i, 0) for i in range(40)])
    second = write_table(tmp_path / 'b.tsv', [row(f'X|{i}', i + 1, 0) for i in range(40)])
    differences = OtuUtils.compare_otu_tables(first, second, max_lines=10)
    assert len(differences) == 11
    assert '30 more' in differences[-1]


def test_example_needs_no_paths():
    options = metaGenomics_new.parse_arguments(['--example'])
    assert options.example is True


def test_paths_are_required_without_the_example(capsys):
    with pytest.raises(SystemExit) as stop:
        metaGenomics_new.parse_arguments(['-o', 'out'])
    assert stop.value.code == 2
    assert '-i' in capsys.readouterr().err


@pytest.mark.real_tools
def test_example_run_matches_the_expected_table(tmp_path):
    returncode, output, hung = run_cli(['--example', '-o', tmp_path / 'out', '-T', '2', '-t', '1'],
                                       timeout=180)
    assert not hung
    assert returncode == 0, output
    assert 'matches the expected table' in output
    assert 'The installation works' in output


@pytest.mark.real_tools
def test_example_writes_its_taxonomy_map_to_the_output_folder(tmp_path):
    database_folder = OtuUtils.EXAMPLE_DATABASE.parent
    before = sorted(entry.name for entry in database_folder.iterdir())
    out = tmp_path / 'out'
    returncode, output, hung = run_cli(['--example', '-o', out, '-T', '2', '-t', '1'], timeout=180)
    assert not hung
    assert returncode == 0, output
    assert (out / f'{OtuUtils.EXAMPLE_DATABASE.stem}_taxonomy_map.tsv').is_file()
    assert sorted(entry.name for entry in database_folder.iterdir()) == before


@pytest.mark.real_tools
@pytest.mark.skipif(os.name == 'nt' or (hasattr(os, 'geteuid') and os.geteuid() == 0),
                    reason='needs a folder that the current user cannot write to')
def test_example_runs_cleanly_from_a_read_only_application_folder(tmp_path):
    app = tmp_path / 'app'
    app.mkdir()
    for source in REPO_ROOT.glob('*.py'):
        shutil.copy2(source, app / source.name)
    shutil.copy2(REPO_ROOT / 'template.html', app / 'template.html')
    shutil.copytree(OtuUtils.EXAMPLE_DIR, app / 'example_data',
                    ignore=shutil.ignore_patterns('*_taxonomy_map.tsv'))
    locked = [app / 'example_data'] + sorted((app / 'example_data').rglob('*'))
    try:
        for entry in locked:
            entry.chmod(0o555 if entry.is_dir() else 0o444)
        returncode, output, hung = run_guarded(
            [sys.executable, str(app / 'metaGenomics_new.py'), '--example', '-o', str(tmp_path / 'out'),
             '-T', '2', '-t', '1'], timeout=180)
    finally:
        for entry in locked:
            entry.chmod(0o755 if entry.is_dir() else 0o644)
    assert not hung
    assert returncode == 0, output
    assert 'CRITICAL ERROR' not in output
    assert 'Permission denied' not in output
    assert 'The installation works' in output


@pytest.mark.real_tools
def test_run_that_differs_from_the_expected_table_fails(tmp_path):
    expected = tmp_path / 'expected.tsv'
    lines = OtuUtils.EXAMPLE_EXPECTED_TABLE.read_text().splitlines()
    cells = lines[1].split('\t')
    cells[-1] = str(int(cells[-1]) + 3)
    lines[1] = '\t'.join(cells)
    expected.write_text('\n'.join(lines) + '\n')
    returncode, output, hung = run_cli(
        ['-i', OtuUtils.EXAMPLE_READS, '-d', OtuUtils.EXAMPLE_DATABASE, '-o', tmp_path / 'out',
         '--expected-table', expected, '-T', '2', '-t', '1'], timeout=180)
    assert not hung
    assert returncode == 1, output
    assert 'differs from the expected table' in output
    assert cells[0] in output


@pytest.mark.real_tools
def test_example_ignores_filter_options_that_would_change_the_result(tmp_path):
    returncode, output, hung = run_cli(
        ['--example', '-o', tmp_path / 'out', '--preset', 'its1', '--min-percent-identity', '99',
         '-T', '2', '-t', '1'], timeout=180)
    assert returncode == 0, output
    assert 'matches the expected table' in output


def test_missing_expected_table_stops_the_run_before_it_starts(tmp_path, dataset):
    returncode, output, hung = run_cli(
        ['-i', dataset['input_dir'], '-d', dataset['db_path'], '-o', tmp_path / 'out',
         '--expected-table', tmp_path / 'absent.tsv'], timeout=60)
    assert returncode == 1
    assert 'absent.tsv' in output
    assert 'STEP 1' not in output
    assert not (tmp_path / 'out').exists()


def test_sample_named_like_a_taxonomy_column_is_compared(tmp_path):
    header = HEADER.replace('S1\tS2', 'Species\tS2')
    first = write_table(tmp_path / 'a.tsv', [row('A|1', 5, 0)], header=header)
    second = write_table(tmp_path / 'b.tsv', [row('A|1', 8, 0)], header=header)
    differences = OtuUtils.compare_otu_tables(first, second)
    assert len(differences) == 1 and 'Species' in differences[0]


@pytest.mark.real_tools
def test_table_left_by_an_earlier_run_is_not_compared(tmp_path, dataset):
    out_dir = tmp_path / 'out'
    arguments = ['-i', dataset['input_dir'], '-d', dataset['db_path'], '-o', out_dir, '-T', '2', '-t', '1']
    returncode, output, hung = run_cli(arguments, timeout=120)
    assert returncode == 0, output
    earlier = tmp_path / 'earlier.tsv'
    shutil.copy(out_dir / 'OTU_Table_out_reference.tsv', earlier)
    returncode, output, hung = run_cli(
        arguments + ['--min-percent-identity', '100', '--expected-table', earlier], timeout=120)
    assert returncode == 1, output
    assert 'matches the expected table' not in output
    assert not (out_dir / 'OTU_Table_out_reference.tsv').exists()


@pytest.mark.real_tools
def test_example_names_the_options_it_ignores(tmp_path):
    returncode, output, hung = run_cli(
        ['--example', '-o', tmp_path / 'out', '--preset', 'its1', '-T', '2', '-t', '1'], timeout=180)
    assert returncode == 0, output
    assert 'ignored' in output and '--preset' in output


def test_output_folder_that_cannot_be_created_is_reported(tmp_path, dataset):
    blocker = tmp_path / 'blocker'
    blocker.write_text('a file where a folder is needed')
    returncode, output, hung = run_cli(
        ['-i', dataset['input_dir'], '-d', dataset['db_path'], '-o', blocker / 'out'], timeout=60)
    assert returncode == 1
    assert 'Traceback' not in output
    assert 'output folder' in output
