"""Analysis parameters shown in the header of the HTML report."""
import pytest

pytest.importorskip('plotly')
pytest.importorskip('jinja2')

import jinja2

import report_generator
from helpers import SOURCES_DIR


def render_header(params):
    environment = jinja2.Environment(loader=jinja2.FileSystemLoader(str(SOURCES_DIR)))
    source = (SOURCES_DIR / 'template.html').read_text(encoding='utf-8')
    start = source.index('{# --- Length Filter --- #}')
    end = source.index('{# --- Quality Filter --- #}')
    return environment.from_string(source[start:end]).render(analysis_params=params)


def test_default_run_is_shown_without_length_filter():
    params = report_generator.extract_filter_params('report_filtering_0_0.tsv', None, None, None)
    assert params['len_filter_active'] is False
    assert 'Disabled' in render_header(params)


def test_filter_switched_off_in_the_interface_is_shown_as_disabled():
    params = report_generator.extract_filter_params('report_filtering_0_999999.tsv', 0, 999999, None)
    assert params['len_filter_active'] is False


def test_window_with_both_limits_is_shown_as_a_range():
    params = report_generator.extract_filter_params('report_filtering_900_1800.tsv', 900, 1800, None)
    assert params['len_filter_active'] is True
    header = render_header(params)
    assert '900' in header and '1800' in header


def test_window_without_upper_limit_is_not_shown_as_ending_at_zero():
    params = report_generator.extract_filter_params('report_filtering_400_0.tsv', 400, 0, None)
    assert params['len_filter_active'] is True
    assert params['max_len'] == 0
    header = render_header(params)
    assert '400' in header
    assert '<strong>0</strong>' not in header
    assert 'no upper limit' in header
