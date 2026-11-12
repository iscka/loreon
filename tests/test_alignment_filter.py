"""Identity and read-coverage filter computed from CIGAR and NM, without samtools expressions."""
import io
import random

import pytest

import OtuUtils
from helpers import expected_counts, make_tool, read_counts

HEADER = (
    b"@HD\tVN:1.6\tSO:coordinate\n"
    b"@SQ\tSN:refA|x\tLN:1500\n"
    b"@SQ\tSN:refB|y\tLN:900\n"
    b"@SQ\tSN:refC|z\tLN:1200\n"
    b"@PG\tID:minimap2\tPN:minimap2\n"
)


def sam_record(name, flag, reference, cigar, edit_distance, sequence_length):
    sequence = 'A' * sequence_length
    tags = '' if edit_distance is None else f"\tNM:i:{edit_distance}\tms:i:100"
    return (
        f"{name}\t{flag}\t{reference}\t1\t60\t{cigar}\t*\t0\t0\t{sequence}\t{'I' * sequence_length}{tags}\n"
    ).encode()


def test_aligned_length_excludes_soft_clips():
    assert OtuUtils.aligned_query_length(b'10S80M10S', 100) == 80


def test_aligned_length_of_fully_aligned_read_is_read_length():
    assert OtuUtils.aligned_query_length(b'100M', 100) == 100


def test_aligned_length_ignores_hard_clips():
    assert OtuUtils.aligned_query_length(b'5H10S80M2I3D5S7H', 97) == 82


def test_aligned_length_with_long_cigar_and_single_clip():
    cigar = b'120S' + b'30M1I29M2D' * 40
    assert OtuUtils.aligned_query_length(cigar, 120 + 60 * 40) == 60 * 40


def test_aligned_length_agrees_with_a_full_cigar_parse():
    generator = random.Random(11)
    for _ in range(2000):
        operations = []
        if generator.random() < 0.3:
            operations.append((generator.randint(1, 500), 'H'))
        if generator.random() < 0.6:
            operations.append((generator.randint(1, 3000), 'S'))
        operations.append((generator.randint(1, 400), 'M'))
        for _ in range(generator.randint(0, 60)):
            operations.append((generator.randint(1, 400), generator.choice('MMMMIDN=X')))
        operations.append((generator.randint(1, 400), 'M'))
        if generator.random() < 0.6:
            operations.append((generator.randint(1, 3000), 'S'))
        if generator.random() < 0.3:
            operations.append((generator.randint(1, 500), 'H'))
        cigar = ''.join(f"{count}{operation}" for count, operation in operations).encode()
        query_length = sum(count for count, operation in operations if operation in 'MIS=X')
        soft_clipped = sum(count for count, operation in operations if operation == 'S')
        assert OtuUtils.aligned_query_length(cigar, query_length) == query_length - soft_clipped, cigar


def test_identity_threshold_is_inclusive():
    assert OtuUtils.passes_quality_filter(80, 100, 4, 95.0, 0.0)
    assert not OtuUtils.passes_quality_filter(80, 100, 5, 95.0, 0.0)


def test_coverage_threshold_uses_aligned_fraction_of_read():
    assert OtuUtils.passes_quality_filter(90, 100, 0, 0.0, 90.0)
    assert not OtuUtils.passes_quality_filter(80, 100, 0, 0.0, 90.0)


def test_both_thresholds_must_hold():
    assert OtuUtils.passes_quality_filter(95, 100, 2, 95.0, 90.0)
    assert not OtuUtils.passes_quality_filter(95, 100, 10, 95.0, 90.0)
    assert not OtuUtils.passes_quality_filter(60, 100, 0, 95.0, 90.0)


def test_missing_edit_distance_fails_only_when_identity_is_checked():
    assert not OtuUtils.passes_quality_filter(100, 100, None, 95.0, 0.0)
    assert OtuUtils.passes_quality_filter(100, 100, None, 0.0, 90.0)


def test_stream_counts_primary_alignments_per_reference_in_header_order():
    lines = io.BytesIO(
        HEADER
        + sam_record('r1', 0, 'refC|z', '100M', 1, 100)
        + sam_record('r2', 16, 'refA|x', '100M', 2, 100)
        + sam_record('r3', 0, 'refC|z', '100M', 0, 100)
        + sam_record('r4', 4, '*', '*', None, 100)
    )
    result = OtuUtils.tabulate_sam_stream(lines)
    assert result['references'] == [('refA|x', 1500), ('refB|y', 900), ('refC|z', 1200)]
    assert result['counts'] == {'refA|x': 1, 'refC|z': 2}
    assert result['total_reads'] == 4
    assert result['mapped_reads'] == 3
    assert result['passed_reads'] == 3


def test_stream_ignores_secondary_and_supplementary_records():
    lines = io.BytesIO(
        HEADER
        + sam_record('r1', 0, 'refA|x', '100M', 0, 100)
        + sam_record('r1', 256, 'refB|y', '100M', 0, 100)
        + sam_record('r1', 2048, 'refC|z', '50M50H', 0, 50)
    )
    result = OtuUtils.tabulate_sam_stream(lines)
    assert result['counts'] == {'refA|x': 1}
    assert result['total_reads'] == 1


def test_stream_applies_identity_and_coverage_thresholds():
    lines = io.BytesIO(
        HEADER
        + sam_record('keep', 0, 'refA|x', '5S95M', 2, 100)
        + sam_record('low_identity', 0, 'refA|x', '100M', 12, 100)
        + sam_record('half_aligned', 0, 'refB|y', '50M50S', 0, 100)
    )
    result = OtuUtils.tabulate_sam_stream(lines, min_percent_identity=95.0, min_ref_coverage=90.0)
    assert result['counts'] == {'refA|x': 1}
    assert result['mapped_reads'] == 3
    assert result['passed_reads'] == 1


def test_stream_forwards_header_and_passing_records_only():
    keep = sam_record('keep', 0, 'refA|x', '100M', 0, 100)
    drop = sam_record('drop', 0, 'refA|x', '100M', 30, 100)
    unmapped = sam_record('unmapped', 4, '*', '*', None, 100)
    sink = io.BytesIO()
    OtuUtils.tabulate_sam_stream(io.BytesIO(HEADER + keep + drop + unmapped),
                                 min_percent_identity=90.0, sink=sink)
    assert sink.getvalue() == HEADER + keep


@pytest.mark.real_tools
@pytest.mark.parametrize('min_identity,min_coverage', [(0.0, 0.0), (85.0, 50.0), (95.0, 90.0)])
def test_counts_match_independent_computation(tmp_path, dataset, min_identity, min_coverage):
    fastq = dataset['input_dir'] / 'barcode01' / 'reads.fastq'
    mapping_dir = tmp_path / 'out' / 'mapping'
    mapping_dir.mkdir(parents=True)
    bam = OtuUtils.mapping_improved(str(fastq), str(dataset['db_path']), 1, str(mapping_dir))
    assert bam is not None
    stats = OtuUtils.tabeling_improved(bam, threads=1, min_percent_identity=min_identity,
                                       min_ref_coverage=min_coverage)
    expected = expected_counts(bam, min_identity, min_coverage)
    observed = read_counts(tmp_path / 'out' / 'tabeling' / 'results' / 'reads.tmp.txt')
    assert observed == expected
    assert stats['mapped_filtered'] == sum(expected.values())
    assert stats['total_reads'] == 48
    assert stats['mapped_unfiltered'] == sum(expected_counts(bam, 0.0, 0.0).values())


@pytest.mark.real_tools
def test_thresholds_remove_divergent_and_partially_aligned_reads(tmp_path, dataset):
    fastq = dataset['input_dir'] / 'barcode01' / 'reads.fastq'
    mapping_dir = tmp_path / 'out' / 'mapping'
    mapping_dir.mkdir(parents=True)
    bam = OtuUtils.mapping_improved(str(fastq), str(dataset['db_path']), 1, str(mapping_dir))
    kept = {}
    for pair in ((0.0, 0.0), (85.0, 50.0), (95.0, 90.0)):
        stats = OtuUtils.tabeling_improved(bam, threads=1, min_percent_identity=pair[0],
                                           min_ref_coverage=pair[1])
        kept[pair] = stats['mapped_filtered']
    assert kept[(0.0, 0.0)] == 44
    assert kept[(85.0, 50.0)] == 36
    assert kept[(95.0, 90.0)] == 24


@pytest.mark.real_tools
def test_filtered_alignment_file_holds_exactly_the_counted_reads(tmp_path, dataset):
    import subprocess
    fastq = dataset['input_dir'] / 'barcode01' / 'reads.fastq'
    mapping_dir = tmp_path / 'out' / 'mapping'
    mapping_dir.mkdir(parents=True)
    bam = OtuUtils.mapping_improved(str(fastq), str(dataset['db_path']), 1, str(mapping_dir))
    stats = OtuUtils.tabeling_improved(bam, threads=1, min_percent_identity=95.0, min_ref_coverage=90.0)
    filtered = tmp_path / 'out' / 'tabeling' / 'reads_sorted_filtered.bam'
    idxstats = subprocess.run(['samtools', 'idxstats', str(filtered)],
                              check=True, capture_output=True, text=True).stdout
    from_index = {l.split('\t')[0]: int(l.split('\t')[2]) for l in idxstats.splitlines()
                  if l.split('\t')[0] != '*' and int(l.split('\t')[2]) > 0}
    assert from_index == read_counts(tmp_path / 'out' / 'tabeling' / 'results' / 'reads.tmp.txt')
    assert sum(from_index.values()) == stats['mapped_filtered']


@pytest.mark.real_tools
@pytest.mark.fake_tools
def test_quality_filter_works_when_samtools_rejects_filter_expressions(tmp_path, dataset, monkeypatch):
    import os
    import shutil
    real_samtools = shutil.which('samtools')
    fastq = dataset['input_dir'] / 'barcode01' / 'reads.fastq'
    mapping_dir = tmp_path / 'out' / 'mapping'
    mapping_dir.mkdir(parents=True)
    bam = OtuUtils.mapping_improved(str(fastq), str(dataset['db_path']), 1, str(mapping_dir))
    expected = expected_counts(bam, 95.0, 90.0)
    wrapper_dir = tmp_path / 'wrapper'
    make_tool(wrapper_dir, 'samtools', f"""
        for argument in "$@"; do
            if [ "$argument" = "-e" ]; then
                echo "[E::bam_sym_lookup] Unknown variable sclen" >&2
                exit 1
            fi
        done
        exec "{real_samtools}" "$@"
    """)
    monkeypatch.setenv('PATH', f"{wrapper_dir}:{os.environ['PATH']}")
    stats = OtuUtils.tabeling_improved(bam, threads=1, min_percent_identity=95.0, min_ref_coverage=90.0)
    assert stats is not None
    assert read_counts(tmp_path / 'out' / 'tabeling' / 'results' / 'reads.tmp.txt') == expected
