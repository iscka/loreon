"""Default settings must let a first run produce results, and explain what the length filter did."""
import sys

import pytest

import OtuUtils
import metaGenomics_new
from helpers import run_cli


def parse(monkeypatch, *arguments):
    monkeypatch.setattr(sys, 'argv', ['metaGenomics_new.py', '-i', 'in', '-o', 'out', '-d', 'db.fasta',
                                      *arguments])
    return metaGenomics_new.parse_arguments()


def test_quality_filters_are_disabled_by_default(monkeypatch):
    options = parse(monkeypatch)
    assert options.min_percent_identity == 0.0
    assert options.min_ref_coverage == 0.0


def test_no_length_filter_by_default(monkeypatch):
    options = parse(monkeypatch)
    assert options.min_len is None
    assert options.max_len is None
    assert metaGenomics_new.resolve_length_window(options.min_len, options.max_len, options.preset) == (0, 0)


def test_presets_cover_the_markers_of_the_article():
    assert OtuUtils.LENGTH_PRESETS['its-lsu'] == (900, 1800)
    assert OtuUtils.LENGTH_PRESETS['rrna-operon'] == (3000, 6500)
    assert OtuUtils.LENGTH_PRESETS['its-full'] == (400, 900)
    assert set(OtuUtils.LENGTH_PRESETS) == set(OtuUtils.LENGTH_PRESET_LABELS)


def test_preset_sets_the_length_window(monkeypatch):
    options = parse(monkeypatch, '--preset', 'its-lsu')
    assert metaGenomics_new.resolve_length_window(options.min_len, options.max_len, options.preset) == (900, 1800)


def test_explicit_bounds_take_precedence_over_the_preset():
    assert metaGenomics_new.resolve_length_window(500, None, 'its-full') == (500, 900)
    assert metaGenomics_new.resolve_length_window(None, 1200, 'its-full') == (400, 1200)
    assert metaGenomics_new.resolve_length_window(450, 1000, 'its-full') == (450, 1000)


def test_length_summary_reports_percentiles():
    summary = OtuUtils.summarize_lengths(list(range(1, 101)))
    assert summary['count'] == 100
    assert summary['min'] == 1
    assert summary['max'] == 100
    assert summary['median'] in (50, 51)
    assert summary['p05'] in (5, 6)
    assert summary['p95'] in (95, 96)


def test_length_summary_honours_weights():
    summary = OtuUtils.summarize_lengths([300, 1300], weights=[1, 99])
    assert summary['median'] == 1300
    assert summary['count'] == 100


def test_advice_warns_when_the_filter_discards_most_reads():
    lengths = [1300] * 95 + [250] * 5
    text = '\n'.join(OtuUtils.length_filter_advice(lengths, None, 100, 5, 200, 300))
    assert 'WARNING' in text
    assert '95 of 100' in text
    assert '1300' in text
    assert 'its-lsu' in text


def test_advice_is_silent_when_the_filter_keeps_most_reads():
    lengths = [1300] * 95 + [250] * 5
    assert OtuUtils.length_filter_advice(lengths, None, 100, 95, 900, 1800) == []


def test_advice_notes_when_no_length_filter_is_active():
    text = '\n'.join(OtuUtils.length_filter_advice([1300] * 100, None, 100, 100, 0, 0))
    assert 'WARNING' not in text
    assert 'no read-length filter' in text.lower()
    assert '1300' in text


def test_report_shows_distribution_and_preset_shares():
    lengths = [1300] * 90 + [4500] * 10
    report = OtuUtils.format_length_report(lengths)
    assert 'median' in report.lower()
    assert 'its-lsu' in report
    assert '90.0' in report
    assert 'rrna-operon' in report
    assert '10.0' in report


def test_reads_can_be_sampled_from_an_input_folder(tmp_path):
    for sample, length in (('barcode01', 1200), ('barcode02', 600)):
        folder = tmp_path / sample
        folder.mkdir()
        with open(folder / 'reads.fastq', 'w') as handle:
            for i in range(50):
                handle.write(f"@r{i}\n{'A' * length}\n+\n{'I' * length}\n")
    lengths = OtuUtils.sample_read_lengths(tmp_path, reads_per_file=20)
    assert sorted(set(lengths)) == [600, 1200]
    assert len(lengths) == 40


@pytest.mark.real_tools
def test_first_run_with_the_three_paths_only_gives_a_non_empty_table(tmp_path, dataset):
    returncode, output, hung = run_cli(
        ['-i', dataset['input_dir'], '-o', tmp_path / 'out', '-d', dataset['db_path']], timeout=120)
    assert not hung
    assert returncode == 0, output
    table = tmp_path / 'out' / 'OTU_Table_out_reference.tsv'
    rows = table.read_text().splitlines()
    assert len(rows) == 4
    assert 'no read-length filter' in output.lower()


@pytest.mark.real_tools
def test_run_explains_when_the_length_filter_discards_every_read(tmp_path, dataset):
    returncode, output, hung = run_cli(
        ['-i', dataset['input_dir'], '-o', tmp_path / 'out', '-d', dataset['db_path'],
         '-min', '200', '-max', '300'], timeout=120)
    assert not hung
    assert returncode == 1
    assert 'discarded 96 of 96 reads' in output
    assert 'its-lsu' in output
