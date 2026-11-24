"""Option "Retain unmapped reads": the reads that are not counted are saved in separate files."""
import gzip
import io
import json
import subprocess
import sys
from html.parser import HTMLParser
from pathlib import Path

import pytest

import OtuUtils
from helpers import REPO_ROOT, make_tool, run_cli, run_guarded, write_fastq

HEADER = (
    b"@HD\tVN:1.6\tSO:coordinate\n"
    b"@SQ\tSN:refA|x\tLN:1500\n"
    b"@SQ\tSN:refB|y\tLN:900\n"
)


def sam_line(name, flag, reference, cigar, edit_distance, sequence, quality=None):
    quality = 'I' * len(sequence) if quality is None else quality
    tags = '' if edit_distance is None else f"\tNM:i:{edit_distance}\tms:i:100"
    return f"{name}\t{flag}\t{reference}\t1\t60\t{cigar}\t*\t0\t0\t{sequence}\t{quality}{tags}\n".encode()


def retained(records, **thresholds):
    """Count SAM records and return the result with the text of the two files of retained reads."""
    unmapped, below = io.BytesIO(), io.BytesIO()
    result = OtuUtils.tabulate_sam_stream(io.BytesIO(HEADER + b''.join(records)), unmapped=unmapped,
                                          below_threshold=below, **thresholds)
    return result, unmapped.getvalue().decode(), below.getvalue().decode()


def names_in(fastq_text):
    return [line[1:] for line in fastq_text.splitlines()[0::4]]


def test_reads_that_do_not_align_are_written_as_fastq_records():
    _, unmapped, below = retained([
        sam_line('counted', 0, 'refA|x', '8M', 0, 'ACGTACGT'),
        sam_line('lost', 4, '*', '*', None, 'AACCGGTT', 'ABCDEFGH'),
    ])
    assert unmapped == '@lost\nAACCGGTT\n+\nABCDEFGH\n'
    assert below == ''


def test_reads_below_the_thresholds_go_to_their_own_file():
    _, unmapped, below = retained([
        sam_line('counted', 0, 'refA|x', '100M', 2, 'A' * 100),
        sam_line('low_identity', 0, 'refA|x', '100M', 12, 'C' * 100),
        sam_line('half_aligned', 0, 'refB|y', '50M50S', 0, 'G' * 100),
        sam_line('lost', 4, '*', '*', None, 'T' * 100),
    ], min_percent_identity=95.0, min_ref_coverage=90.0)
    assert names_in(below) == ['low_identity', 'half_aligned']
    assert names_in(unmapped) == ['lost']


def test_read_aligned_on_the_reverse_strand_is_written_as_it_was_sequenced():
    _, _, below = retained([
        sam_line('reverse', 16, 'refA|x', '8M', 4, 'AACCGTTN', '12345678'),
    ], min_percent_identity=95.0)
    assert below == '@reverse\nNAACGGTT\n+\n87654321\n'


def test_secondary_and_supplementary_records_are_not_written():
    _, unmapped, below = retained([
        sam_line('read', 0, 'refA|x', '100M', 30, 'A' * 100),
        sam_line('read', 256, 'refB|y', '100M', 30, 'A' * 100),
        sam_line('read', 2048, 'refB|y', '50M50H', 30, 'A' * 50),
    ], min_percent_identity=95.0)
    assert unmapped == ''
    assert names_in(below) == ['read']


def test_counts_are_the_same_with_and_without_the_files():
    records = [
        sam_line('a', 0, 'refA|x', '100M', 2, 'A' * 100),
        sam_line('b', 16, 'refB|y', '100M', 12, 'C' * 100),
        sam_line('c', 4, '*', '*', None, 'T' * 100),
    ]
    with_files, _, _ = retained(records, min_percent_identity=95.0)
    without_files = OtuUtils.tabulate_sam_stream(io.BytesIO(HEADER + b''.join(records)),
                                                 min_percent_identity=95.0)
    assert with_files == without_files


def test_read_stored_without_qualities_gets_a_constant_quality():
    _, unmapped, _ = retained([sam_line('bare', 4, '*', '*', None, 'ACGT', '*')])
    assert unmapped == '@bare\nACGT\n+\nIIII\n'


FAILING_VIEW = """
    if [ "$1" = "view" ]; then
        for argument in "$@"; do
            if [ "$argument" = "-b" ]; then cat > /dev/null; exit 0; fi
        done
        printf '@HD\\tVN:1.6\\tSO:coordinate\\n@SQ\\tSN:refA\\tLN:1500\\n'
        printf 'lost\\t4\\t*\\t0\\t0\\t*\\t*\\t0\\t0\\tACGT\\tIIII\\n'
        printf 'low\\t0\\trefA\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\tNM:i:4\\n'
        exit 3
    fi
    exit 0
"""


def fastq_records(path):
    """Return {name: (sequence, quality)} of a FASTQ file, plain or compressed."""
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as handle:
        lines = handle.read().splitlines()
    return {lines[i][1:].split()[0]: (lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)}


def map_reads(fastq, database, out):
    (out / 'mapping').mkdir(parents=True)
    bam = OtuUtils.mapping_improved(str(fastq), str(database), 1, str(out / 'mapping'))
    assert bam is not None
    return bam


@pytest.fixture
def mapped_sample(tmp_path, dataset):
    """Map the first sample of the test dataset; return its reads, its alignment file and the output folder."""
    fastq = dataset['input_dir'] / 'barcode01' / 'reads.fastq'
    out = tmp_path / 'out'
    return fastq, map_reads(fastq, dataset['db_path'], out), out


@pytest.mark.real_tools
def test_reads_that_do_not_align_are_saved_one_file_per_sample(mapped_sample):
    fastq, bam, out = mapped_sample
    stats = OtuUtils.tabeling_improved(bam, threads=1, retain_unmapped=True)
    saved = fastq_records(out / 'unmapped_seq' / 'reads.fastq')
    original = fastq_records(fastq)
    assert sorted(saved) == [f'barcode01_random_{i}' for i in range(4)]
    assert all(saved[name] == original[name] for name in saved)
    assert not (out / 'below_threshold_seq').exists()
    assert stats['unmapped_saved'] == 4
    assert stats['unmapped_file'] == 'unmapped_seq/reads.fastq'
    assert stats['below_threshold_saved'] == 0
    assert stats['below_threshold_file'] is None


@pytest.mark.real_tools
def test_reads_below_the_thresholds_are_saved_in_a_second_folder(mapped_sample):
    fastq, bam, out = mapped_sample
    stats = OtuUtils.tabeling_improved(bam, threads=1, min_percent_identity=85.0, min_ref_coverage=50.0,
                                       retain_unmapped=True)
    original = fastq_records(fastq)
    unmapped = fastq_records(out / 'unmapped_seq' / 'reads.fastq')
    below = fastq_records(out / 'below_threshold_seq' / 'reads.fastq')
    assert (len(unmapped), len(below)) == (4, 8)
    assert not set(unmapped) & set(below)
    assert all(original[name] == record for name, record in {**unmapped, **below}.items())
    assert stats['mapped_filtered'] + len(unmapped) + len(below) == stats['total_reads'] == 48
    assert stats['below_threshold_saved'] == 8
    assert stats['below_threshold_file'] == 'below_threshold_seq/reads.fastq'


@pytest.mark.real_tools
@pytest.mark.parametrize('min_identity,min_coverage', [(95.0, 0.0), (0.0, 90.0)])
def test_one_threshold_is_enough_for_the_second_file(mapped_sample, min_identity, min_coverage):
    _, bam, out = mapped_sample
    stats = OtuUtils.tabeling_improved(bam, threads=1, min_percent_identity=min_identity,
                                       min_ref_coverage=min_coverage, retain_unmapped=True)
    below = fastq_records(out / 'below_threshold_seq' / 'reads.fastq')
    assert len(below) == stats['mapped_unfiltered'] - stats['mapped_filtered'] > 0
    assert stats['below_threshold_saved'] == len(below)


@pytest.mark.real_tools
def test_reads_sequenced_on_the_reverse_strand_are_saved_as_in_the_input(tmp_path, dataset):
    forward = fastq_records(dataset['input_dir'] / 'barcode01' / 'reads.fastq')
    complement = str.maketrans('ACGT', 'TGCA')
    fastq = tmp_path / 'reversed.fastq'
    write_fastq(fastq, [(name, sequence.translate(complement)[::-1]) for name, (sequence, _) in forward.items()])
    out = tmp_path / 'out'
    bam = map_reads(fastq, dataset['db_path'], out)
    on_reverse_strand = subprocess.run(['samtools', 'view', '-c', '-f', '16', '-F', '0x904', bam],
                                       check=True, capture_output=True, text=True).stdout
    assert int(on_reverse_strand) == 44
    OtuUtils.tabeling_improved(bam, threads=1, min_percent_identity=95.0, min_ref_coverage=90.0,
                               retain_unmapped=True)
    original = fastq_records(fastq)
    below = fastq_records(out / 'below_threshold_seq' / 'reversed.fastq')
    assert len(below) == 20
    assert all(original[name] == record for name, record in below.items())


@pytest.mark.real_tools
def test_nothing_is_saved_unless_requested(mapped_sample):
    _, bam, out = mapped_sample
    stats = OtuUtils.tabeling_improved(bam, threads=1, min_percent_identity=85.0, min_ref_coverage=50.0)
    assert not (out / 'unmapped_seq').exists()
    assert not (out / 'below_threshold_seq').exists()
    assert 'unmapped_saved' not in stats


@pytest.mark.real_tools
def test_files_are_compressed_on_request(mapped_sample):
    _, bam, out = mapped_sample
    stats = OtuUtils.tabeling_improved(bam, threads=1, min_percent_identity=85.0, min_ref_coverage=50.0,
                                       retain_unmapped=True, gzip_output=True)
    assert stats['unmapped_file'] == 'unmapped_seq/reads.fastq.gz'
    assert stats['below_threshold_file'] == 'below_threshold_seq/reads.fastq.gz'
    assert len(fastq_records(out / 'unmapped_seq' / 'reads.fastq.gz')) == 4
    assert len(fastq_records(out / 'below_threshold_seq' / 'reads.fastq.gz')) == 8


@pytest.mark.real_tools
def test_file_without_reads_is_not_left(mapped_sample):
    _, bam, out = mapped_sample
    stats = OtuUtils.tabeling_improved(bam, threads=1, min_percent_identity=1.0, retain_unmapped=True)
    assert stats['below_threshold_saved'] == 0
    assert stats['below_threshold_file'] is None
    assert not (out / 'below_threshold_seq' / 'reads.fastq').exists()
    assert stats['unmapped_saved'] == 4


@pytest.mark.real_tools
def test_files_remain_when_the_alignment_files_are_deleted(mapped_sample):
    _, bam, out = mapped_sample
    OtuUtils.tabeling_improved(bam, threads=1, retain_unmapped=True, delete_bam=True)
    assert not Path(bam).exists()
    assert len(fastq_records(out / 'unmapped_seq' / 'reads.fastq')) == 4


@pytest.mark.real_tools
def test_files_left_by_an_earlier_run_of_the_same_sample_are_replaced(mapped_sample):
    _, bam, out = mapped_sample
    OtuUtils.tabeling_improved(bam, threads=1, min_percent_identity=85.0, min_ref_coverage=50.0,
                               retain_unmapped=True, gzip_output=True)
    assert (out / 'below_threshold_seq' / 'reads.fastq.gz').exists()
    OtuUtils.tabeling_improved(bam, threads=1, retain_unmapped=True)
    assert sorted(path.name for path in (out / 'unmapped_seq').iterdir()) == ['reads.fastq']
    assert list((out / 'below_threshold_seq').iterdir()) == []


@pytest.mark.fake_tools
def test_partial_files_are_removed_when_the_sample_fails(tmp_path, fake_bin):
    make_tool(fake_bin, 'samtools', FAILING_VIEW)
    mapping_dir = tmp_path / 'out' / 'mapping'
    mapping_dir.mkdir(parents=True)
    bam = mapping_dir / 'sample_sorted_full.bam'
    bam.write_text('stand-in alignment file')
    stats = OtuUtils.tabeling_improved(str(bam), threads=1, min_percent_identity=95.0, retain_unmapped=True)
    assert stats is None
    assert not (tmp_path / 'out' / 'unmapped_seq' / 'sample.fastq').exists()
    assert not (tmp_path / 'out' / 'below_threshold_seq' / 'sample.fastq').exists()


def run_pipeline(dataset, out, *options):
    returncode, output, hung = run_cli(
        ['-i', dataset['input_dir'], '-o', out, '-d', dataset['db_path'], '-T', '2', '-t', '1', *options],
        timeout=180)
    assert returncode == 0 and not hung, output
    return json.loads((out / 'mapping_stats.json').read_text(encoding='utf-8')), output


@pytest.mark.real_tools
def test_run_saves_the_reads_that_are_not_counted(tmp_path, dataset):
    out = tmp_path / 'out'
    stats, output = run_pipeline(dataset, out, '--min-percent-identity', '85', '--min-ref-coverage', '50',
                                 '--retain-unmapped')
    assert stats['parameters']['retain_unmapped'] is True
    assert sorted(sample['sample'] for sample in stats['samples']) == ['barcode01', 'barcode02']
    for sample in stats['samples']:
        name = sample['sample']
        assert sample['unmapped_file'] == f'unmapped_seq/{name}.fastq'
        assert sample['below_threshold_file'] == f'below_threshold_seq/{name}.fastq'
        unmapped = fastq_records(out / sample['unmapped_file'])
        below = fastq_records(out / sample['below_threshold_file'])
        assert sorted(unmapped) == [f'{name}_random_{i}' for i in range(4)]
        assert len(unmapped) == sample['unmapped_saved']
        assert len(below) == sample['below_threshold_saved'] > 0
        assert not set(unmapped) & set(below)
        assert sample['mapped_filtered'] + len(unmapped) + len(below) == sample['total_reads'] == 48
    assert stats['summary']['unmapped_saved_all'] == 8
    assert stats['summary']['below_threshold_saved_all'] == sum(
        sample['below_threshold_saved'] for sample in stats['samples'])
    assert 'unmapped_seq' in output
    assert 'below_threshold_seq' in output


@pytest.mark.real_tools
def test_run_writes_nothing_more_without_the_option(tmp_path, dataset):
    out = tmp_path / 'out'
    stats, output = run_pipeline(dataset, out, '--min-percent-identity', '85')
    assert stats['parameters']['retain_unmapped'] is False
    assert 'unmapped_file' not in stats['samples'][0]
    assert 'unmapped_saved_all' not in stats['summary']
    assert not (out / 'unmapped_seq').exists()
    assert not (out / 'below_threshold_seq').exists()
    assert 'unmapped_seq' not in output


@pytest.mark.real_tools
def test_run_compresses_the_files_together_with_the_other_read_files(tmp_path, dataset):
    out = tmp_path / 'out'
    stats, _ = run_pipeline(dataset, out, '--retain-unmapped', '--gzip-intermediate')
    for sample in stats['samples']:
        assert sample['unmapped_file'] == f"unmapped_seq/{sample['sample']}.fastq.gz"
        assert len(fastq_records(out / sample['unmapped_file'])) == 4
    assert not (out / 'below_threshold_seq').exists()


def retained_stats(min_identity):
    """Mapping statistics of a run with the option, as the pipeline writes them."""
    filtering = min_identity > 0
    first = {'sample': 'barcode01', 'total_reads': 500, 'mapped_unfiltered': 490,
             'mapped_filtered': 450 if filtering else 490, 'mapping_rate_pct': 98.0,
             'quality_retention_pct': 91.84 if filtering else 100.0,
             'unmapped_saved': 10, 'unmapped_file': 'unmapped_seq/barcode01.fastq',
             'below_threshold_saved': 40 if filtering else 0,
             'below_threshold_file': 'below_threshold_seq/barcode01.fastq' if filtering else None}
    second = {'sample': 'barcode02', 'total_reads': 500, 'mapped_unfiltered': 500, 'mapped_filtered': 500,
              'mapping_rate_pct': 100.0, 'quality_retention_pct': 100.0,
              'unmapped_saved': 0, 'unmapped_file': None,
              'below_threshold_saved': 0, 'below_threshold_file': None}
    counted = first['mapped_filtered'] + second['mapped_filtered']
    return {
        'samples': [first, second],
        'parameters': {'min_percent_identity': min_identity, 'min_ref_coverage': 0.0, 'retain_unmapped': True},
        'summary': {'total_reads_all': 1000, 'total_mapped_all': 990, 'total_after_quality_filter': counted,
                    'overall_mapping_rate_pct': 99.0,
                    'overall_quality_retention_pct': round(counted / 990 * 100, 2),
                    'unmapped_saved_all': 10, 'below_threshold_saved_all': first['below_threshold_saved']},
    }


def report_text(folder, mapping_stats):
    """Run the report generator with the given mapping statistics and return the text a reader sees."""
    for module in ('pandas', 'plotly', 'jinja2'):
        pytest.importorskip(module)
    from test_report_database_note import FILTER_REPORT, OTU_TABLE, ReportText
    folder.mkdir()
    (folder / 'report_filtering_0_0.tsv').write_text(FILTER_REPORT, encoding='utf-8')
    (folder / 'OTU_Table.tsv').write_text(OTU_TABLE, encoding='utf-8')
    (folder / 'mapping_stats.json').write_text(json.dumps(mapping_stats), encoding='utf-8')
    returncode, output, hung = run_guarded(
        [sys.executable, str(REPO_ROOT / 'report_generator.py'),
         '-f', str(folder / 'report_filtering_0_0.tsv'), '-otu', str(folder / 'OTU_Table.xlsx'),
         '-o', str(folder / 'report.html'), '-ms', str(folder / 'mapping_stats.json')], timeout=120)
    assert returncode == 0 and not hung, output
    report = ReportText()
    report.feed((folder / 'report.html').read_text(encoding='utf-8'))
    return report.visible


def test_report_gives_the_path_and_the_number_of_reads_of_each_file(tmp_path):
    from test_report_database_note import DATABASE_NOTE
    text = report_text(tmp_path / 'run', retained_stats(90.0))
    assert 'Reads Not Counted, Saved to File' in text
    assert 'barcode01 10 unmapped_seq/barcode01.fastq 40 below_threshold_seq/barcode01.fastq' in text
    assert 'barcode02 0 no file 0 no file' in text
    assert text.endswith(DATABASE_NOTE)


def test_report_lists_only_the_unmapped_reads_when_no_threshold_is_set(tmp_path):
    text = report_text(tmp_path / 'run', retained_stats(0.0))
    assert 'barcode01 10 unmapped_seq/barcode01.fastq barcode02 0 no file' in text
    assert 'below_threshold_seq' not in text


def test_report_has_no_such_section_without_the_option(tmp_path):
    from test_report_database_note import MAPPING_STATS
    text = report_text(tmp_path / 'run', MAPPING_STATS)
    assert 'Saved to File' not in text
    assert 'unmapped_seq' not in text


class LineWrapCheck(HTMLParser):
    """Tell whether a text of the report sits in a block whose lines are not wrapped."""

    def __init__(self, text):
        super().__init__()
        self.text = text
        self.blocks = []
        self.kept_on_one_line = None

    def handle_starttag(self, tag, attrs):
        if tag == 'div':
            self.blocks.append((dict(attrs).get('class') or '').split())

    def handle_endtag(self, tag):
        if tag == 'div':
            self.blocks.pop()

    def handle_data(self, data):
        if self.text in data:
            self.kept_on_one_line = any('table-wrapper' in classes for classes in self.blocks)


def test_description_of_the_files_wraps_inside_the_page(tmp_path):
    report_text(tmp_path / 'run', retained_stats(90.0))
    check = LineWrapCheck('Files written by the option')
    check.feed((tmp_path / 'run' / 'report.html').read_text(encoding='utf-8'))
    assert check.kept_on_one_line is False
