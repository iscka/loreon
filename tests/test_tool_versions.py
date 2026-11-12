"""Start-up check of the external tools."""
import pytest

import metaGenomics_new
from helpers import make_tool

pytestmark = pytest.mark.fake_tools


def answer_version(banner):
    return f"""
        echo "{banner}"
        exit 0
    """


def test_version_is_read_from_tool_output():
    assert metaGenomics_new.parse_tool_version('2.28-r1209') == (2, 28)
    assert metaGenomics_new.parse_tool_version('samtools 1.21\nUsing htslib 1.21') == (1, 21)
    assert metaGenomics_new.parse_tool_version('samtools 1.13') == (1, 13)
    assert metaGenomics_new.parse_tool_version('samtools 1.9') == (1, 9)
    assert metaGenomics_new.parse_tool_version('no version here') is None


def test_supported_versions_are_accepted_and_reported(fake_bin, capsys):
    make_tool(fake_bin, 'minimap2', answer_version('2.28-r1209'))
    make_tool(fake_bin, 'samtools', answer_version('samtools 1.21'))
    versions = metaGenomics_new.check_dependencies()
    output = capsys.readouterr().out
    assert versions == {'minimap2': '2.28-r1209', 'samtools': '1.21'}
    assert 'minimap2 2.28-r1209' in output
    assert 'samtools 1.21' in output


def test_oldest_supported_versions_are_accepted(fake_bin):
    make_tool(fake_bin, 'minimap2', answer_version('2.24-r1122'))
    make_tool(fake_bin, 'samtools', answer_version('samtools 1.13'))
    assert metaGenomics_new.check_dependencies() == {'minimap2': '2.24-r1122', 'samtools': '1.13'}


def test_outdated_samtools_stops_the_run_with_an_explicit_message(fake_bin, capsys):
    make_tool(fake_bin, 'minimap2', answer_version('2.28-r1209'))
    make_tool(fake_bin, 'samtools', answer_version('samtools 1.10'))
    with pytest.raises(SystemExit) as stop:
        metaGenomics_new.check_dependencies()
    output = capsys.readouterr().out
    assert stop.value.code == 1
    assert 'samtools 1.10' in output
    assert '1.13' in output


def test_outdated_minimap2_stops_the_run_with_an_explicit_message(fake_bin, capsys):
    make_tool(fake_bin, 'minimap2', answer_version('2.17-r941'))
    make_tool(fake_bin, 'samtools', answer_version('samtools 1.21'))
    with pytest.raises(SystemExit):
        metaGenomics_new.check_dependencies()
    output = capsys.readouterr().out
    assert 'minimap2 2.17' in output
    assert '2.24' in output


def test_missing_tool_stops_the_run(fake_bin, capsys):
    make_tool(fake_bin, 'minimap2', answer_version('2.28-r1209'))
    with pytest.raises(SystemExit):
        metaGenomics_new.check_dependencies()
    assert 'samtools' in capsys.readouterr().out


def test_version_output_with_bytes_that_are_not_utf8_is_accepted(fake_bin):
    make_tool(fake_bin, 'minimap2', answer_version('2.24-r1122'))
    make_tool(fake_bin, 'samtools', """
        echo "samtools 1.13"
        printf 'Samtools compilation details:\\n    CFLAGS: \\253-O2\\273\\n'
        exit 0
    """)
    assert metaGenomics_new.check_dependencies() == {'minimap2': '2.24-r1122', 'samtools': '1.13'}
