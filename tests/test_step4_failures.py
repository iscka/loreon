"""Step 4 must stop with an explicit error instead of waiting forever or reporting success."""
import os
import time

import pytest

import OtuUtils
import metaGenomics_new
from helpers import SYSTEM_PATH, build_dataset, make_tool, run_cli, run_python

pytestmark = pytest.mark.fake_tools

VERSION_ANSWERS = """
    if [ "$1" = "--version" ]; then
        echo "{banner}"
        exit 0
    fi
"""

FLAGSTAT_ANSWER = """
    if [ "$1" = "flagstat" ]; then
        echo "100 + 0 in total (QC-passed reads + QC-failed reads)"
        echo "100 + 0 primary"
        echo "0 + 0 secondary"
        echo "0 + 0 supplementary"
        echo "0 + 0 duplicates"
        echo "0 + 0 primary duplicates"
        echo "90 + 0 mapped (90.00% : N/A)"
        echo "90 + 0 primary mapped (90.00% : N/A)"
        exit 0
    fi
"""

SORT_WRITES_OUTPUT = """
    if [ "$1" = "sort" ]; then
        output=""
        previous=""
        for argument in "$@"; do
            if [ "$previous" = "-o" ]; then output="$argument"; fi
            previous="$argument"
        done
        cat > /dev/null
        if [ -n "$output" ]; then echo "stand-in alignment file" > "$output"; fi
        exit 0
    fi
"""


def samtools_script(view_behaviour):
    return (VERSION_ANSWERS.format(banner='samtools 1.21') + FLAGSTAT_ANSWER + SORT_WRITES_OUTPUT
            + f"""
    if [ "$1" = "view" ]; then
        {view_behaviour}
    fi
    exit 0
""")


FLOOD_AND_FAIL = 'yes "[E::bam_sym_lookup] Unknown variable sclen" | head -n 20000 >&2; exit 1'
FAIL_QUIETLY = 'echo "[E::hts_open_format] Failed to open file: disk quota exceeded" >&2; exit 3'
MINIMAP2_PRINTS_SAM = VERSION_ANSWERS.format(banner='2.28-r1209') + """
    printf '@HD\\tVN:1.6\\tSO:unsorted\\n'
    exit 0
"""


def stand_in_alignment(tmp_path):
    mapping_dir = tmp_path / 'out' / 'mapping'
    mapping_dir.mkdir(parents=True)
    bam = mapping_dir / 'sample_sorted_full.bam'
    bam.write_text('stand-in alignment file')
    return bam


def test_tabulation_stops_when_samtools_floods_its_error_stream(tmp_path, fake_bin):
    make_tool(fake_bin, 'samtools', samtools_script(FLOOD_AND_FAIL))
    bam = stand_in_alignment(tmp_path)
    returncode, output, hung = run_python(f"""
        import OtuUtils
        result = OtuUtils.tabeling_improved({str(bam)!r}, threads=1,
                                            min_percent_identity=95.0, min_ref_coverage=90.0)
        print('RESULT', result)
    """, timeout=20)
    assert not hung, 'tabulation never returned'
    assert returncode == 0
    assert 'RESULT None' in output
    assert 'Unknown variable sclen' in output


def test_mapping_stops_when_minimap2_floods_its_error_stream(tmp_path, fake_bin):
    make_tool(fake_bin, 'minimap2', VERSION_ANSWERS.format(banner='2.28-r1209') + """
        yes "[WARNING] wrong FASTQ record" | head -n 20000 >&2
        exit 1
    """)
    make_tool(fake_bin, 'samtools', samtools_script('exit 0'))
    data = build_dataset(tmp_path / 'data')
    mapping_dir = tmp_path / 'out' / 'mapping'
    mapping_dir.mkdir(parents=True)
    fastq = data['input_dir'] / 'barcode01' / 'reads.fastq'
    returncode, output, hung = run_python(f"""
        import OtuUtils
        result = OtuUtils.mapping_improved({str(fastq)!r}, {str(data['db_path'])!r}, 1, {str(mapping_dir)!r})
        print('RESULT', result)
    """, timeout=20)
    assert not hung, 'mapping never returned'
    assert returncode == 0
    assert 'RESULT None' in output
    assert 'wrong FASTQ record' in output


def test_tabulation_error_reports_tool_exit_status_and_last_error_lines(tmp_path, fake_bin, capsys):
    make_tool(fake_bin, 'samtools', samtools_script(FAIL_QUIETLY))
    bam = stand_in_alignment(tmp_path)
    assert OtuUtils.tabeling_improved(str(bam), threads=1) is None
    output = capsys.readouterr().out
    assert 'samtools view' in output
    assert 'exit status 3' in output
    assert 'disk quota exceeded' in output


def test_failed_tabulation_is_not_reported_as_success(tmp_path, fake_bin):
    make_tool(fake_bin, 'minimap2', MINIMAP2_PRINTS_SAM)
    make_tool(fake_bin, 'samtools', samtools_script(FAIL_QUIETLY))
    data = build_dataset(tmp_path / 'data')
    mapping_dir = tmp_path / 'out' / 'mapping'
    mapping_dir.mkdir(parents=True)
    fastq = data['input_dir'] / 'barcode01' / 'reads.fastq'
    result = metaGenomics_new.run_mapping_for_one_sample(
        fastq, db_path=str(data['db_path']), threads=1, output_dir=str(mapping_dir),
        kmer_size=15, window_size=10, debug=False,
        min_percent_identity=0.0, min_ref_coverage=0.0, delete_bam=False,
    )
    assert result[1] != 'Success'
    assert result[2] is None


def test_pipeline_exits_with_error_when_every_sample_fails(tmp_path, fake_bin):
    make_tool(fake_bin, 'minimap2', MINIMAP2_PRINTS_SAM)
    make_tool(fake_bin, 'samtools', samtools_script(FAIL_QUIETLY))
    data = build_dataset(tmp_path / 'data')
    env = dict(os.environ, PATH=f"{fake_bin}:{SYSTEM_PATH}")
    returncode, output, hung = run_cli(
        ['-i', data['input_dir'], '-o', tmp_path / 'out', '-d', data['db_path'], '-f', 'unite',
         '-min', '0', '-max', '0', '--min-percent-identity', '0', '--min-ref-coverage', '0',
         '-T', '2', '-t', '1'],
        timeout=60, env=env)
    assert not hung
    assert returncode == 1
    assert 'disk quota exceeded' in output
    assert 'Pipeline Complete' not in output


def test_missing_tool_is_reported_as_a_failure_not_as_an_exit(tmp_path, monkeypatch, capsys):
    empty_bin = tmp_path / 'empty'
    empty_bin.mkdir()
    monkeypatch.setenv('PATH', str(empty_bin))
    data = build_dataset(tmp_path / 'data')
    mapping_dir = tmp_path / 'out' / 'mapping'
    mapping_dir.mkdir(parents=True)
    fastq = data['input_dir'] / 'barcode01' / 'reads.fastq'
    assert OtuUtils.mapping_improved(str(fastq), str(data['db_path']), 1, str(mapping_dir)) is None
    assert 'not found' in capsys.readouterr().out


def test_sample_time_limit_stops_a_stuck_mapping(tmp_path, fake_bin):
    pid_file = tmp_path / 'minimap2.pid'
    make_tool(fake_bin, 'minimap2', VERSION_ANSWERS.format(banner='2.28-r1209') + f"""
        echo $$ > "{pid_file}"
        exec sleep 300
    """)
    make_tool(fake_bin, 'samtools', samtools_script('exit 0'))
    data = build_dataset(tmp_path / 'data')
    mapping_dir = tmp_path / 'out' / 'mapping'
    mapping_dir.mkdir(parents=True)
    fastq = data['input_dir'] / 'barcode01' / 'reads.fastq'
    started = time.monotonic()
    returncode, output, hung = run_python(f"""
        import OtuUtils
        result = OtuUtils.mapping_improved({str(fastq)!r}, {str(data['db_path'])!r}, 1,
                                           {str(mapping_dir)!r}, timeout_s=2)
        print('RESULT', result)
    """, timeout=30)
    elapsed = time.monotonic() - started
    assert not hung
    assert 'RESULT None' in output
    assert 'time limit' in output
    assert elapsed < 15
    stuck_pid = int(pid_file.read_text().strip())
    with pytest.raises(ProcessLookupError):
        os.kill(stuck_pid, 0)


def test_pipeline_prints_progress_while_samples_are_running(tmp_path, fake_bin):
    make_tool(fake_bin, 'minimap2', VERSION_ANSWERS.format(banner='2.28-r1209') + """
        sleep 4
        printf '@HD\\tVN:1.6\\tSO:unsorted\\n'
        exit 0
    """)
    make_tool(fake_bin, 'samtools', samtools_script(FAIL_QUIETLY))
    data = build_dataset(tmp_path / 'data')
    env = dict(os.environ, PATH=f"{fake_bin}:{SYSTEM_PATH}")
    returncode, output, hung = run_cli(
        ['-i', data['input_dir'], '-o', tmp_path / 'out', '-d', data['db_path'], '-f', 'unite',
         '-min', '0', '-max', '0', '--min-percent-identity', '0', '--min-ref-coverage', '0',
         '-T', '2', '-t', '1', '--progress-interval', '1'],
        timeout=60, env=env)
    assert not hung
    assert 'Step 4 progress: 0/2 samples completed' in output


@pytest.mark.real_tools
def test_failed_sample_is_left_out_and_reported_with_its_own_exit_status(tmp_path, monkeypatch):
    import shutil
    real_minimap2 = shutil.which('minimap2', path=os.environ['PATH'])
    data = build_dataset(tmp_path / 'data')
    wrapper_dir = tmp_path / 'wrapper'
    make_tool(wrapper_dir, 'minimap2', f"""
        for argument in "$@"; do
            case "$argument" in
                *barcode02*) echo "[ERROR] failed to open the query file" >&2; exit 1;;
            esac
        done
        exec "{real_minimap2}" "$@"
    """)
    out_dir = tmp_path / 'out'
    stale = out_dir / 'tabeling' / 'results' / 'barcode02.tmp.txt'
    stale.parent.mkdir(parents=True)
    stale.write_text("Alpha_one|AC000001|SH0000001.10FU|refs|x\t1500\t999\t0\n*\t0\t0\t0\n")
    env = dict(os.environ, PATH=f"{wrapper_dir}:{os.environ['PATH']}")
    returncode, output, hung = run_cli(
        ['-i', data['input_dir'], '-o', out_dir, '-d', data['db_path'], '-T', '2', '-t', '1'],
        timeout=120, env=env)
    assert not hung
    assert returncode == 3, output
    assert '1 of 2 samples failed' in output
    assert 'failed to open the query file' in output
    header = (out_dir / 'OTU_Table_out_reference.tsv').read_text().splitlines()[0].split('\t')
    assert 'barcode01' in header
    assert 'barcode02' not in header
    assert not stale.exists()


def test_error_report_names_the_tool_that_failed_first(tmp_path, fake_bin, capsys):
    make_tool(fake_bin, 'samtools', """
        for argument in "$@"; do
            if [ "$argument" = "-b" ]; then
                echo "[E::bgzf_flush] File write failed: No space left on device" >&2
                exit 1
            fi
        done
        exec yes "@CO	filler line that keeps the reader busy"
    """)
    bam = stand_in_alignment(tmp_path)
    assert OtuUtils.tabeling_improved(str(bam), threads=1) is None
    output = capsys.readouterr().out
    assert 'No space left on device' in output
    assert 'exit status 1' in output
    assert 'signal 13' not in output


def test_stopped_mapping_leaves_no_partial_files(tmp_path, fake_bin):
    make_tool(fake_bin, 'minimap2', VERSION_ANSWERS.format(banner='2.28-r1209') + """
        exec sleep 300
    """)
    make_tool(fake_bin, 'samtools', VERSION_ANSWERS.format(banner='samtools 1.21') + """
        if [ "$1" = "sort" ]; then
            output=""
            previous=""
            for argument in "$@"; do
                if [ "$previous" = "-o" ]; then output="$argument"; fi
                previous="$argument"
            done
            echo partial > "$output"
            echo partial > "$output.tmp.0000.bam"
            exec cat > /dev/null
        fi
        exit 0
    """)
    data = build_dataset(tmp_path / 'data')
    mapping_dir = tmp_path / 'out' / 'mapping'
    mapping_dir.mkdir(parents=True)
    fastq = data['input_dir'] / 'barcode01' / 'reads.fastq'
    result = OtuUtils.mapping_improved(str(fastq), str(data['db_path']), 1, str(mapping_dir),
                                       timeout_s=2)
    assert result is None
    assert list(mapping_dir.iterdir()) == []
