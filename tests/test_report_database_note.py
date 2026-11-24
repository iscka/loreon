"""Note on the reference database printed at the bottom of the HTML report."""
import json
import sys
from html.parser import HTMLParser

import pytest

pytest.importorskip('pandas')
pytest.importorskip('plotly')
pytest.importorskip('jinja2')

from helpers import SOURCES_DIR, run_guarded

DATABASE_NOTE = (
    'These results depend on the reference database. LOREON assigns each read to the reference '
    'of its best alignment and does not fall back to a higher taxonomic rank: taxa absent from '
    'the database cannot be identified, and their reads are attributed to the closest reference '
    'available. A high mapping rate is therefore not evidence that the database is complete. '
    'Reads are left unassigned when they do not align, or when their alignment falls below the '
    'minimum identity and coverage thresholds, if set.'
)

FILTER_REPORT = (
    'barcode\tmean length\tmedian length\tbad count\tgood count\ttotal count\n'
    'barcode01\t1100.0\t1050.0\t0\t500\t500\n'
    'barcode02\t1200.0\t1150.0\t0\t500\t500\n'
)

OTU_TABLE = (
    'OTU_ID\tTaxonomy\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tDetails'
    '\tbarcode01\tbarcode02\n'
    'AC000001|SH0000001.10FU\tFungi;Ascomycota;Classa;Ordera;Familya;Alpha;Alpha_one'
    '\tFungi\tAscomycota\tClassa\tOrdera\tFamilya\tAlpha\tAlpha one\t\t300\t200\n'
    'AC000002|SH0000002.10FU\tFungi;Ascomycota;Classb;Orderb;Familyb;Beta;Beta_two'
    '\tFungi\tAscomycota\tClassb\tOrderb\tFamilyb\tBeta\tBeta two\t\t150\t250\n'
    'AC000003|SH0000003.10FU\tFungi;Basidiomycota;Classc;Orderc;Familyc;Gamma;Gamma_three'
    '\tFungi\tBasidiomycota\tClassc\tOrderc\tFamilyc\tGamma\tGamma three\t\t40\t45\n'
)

MAPPING_STATS = {
    'samples': [
        {'sample': 'barcode01', 'total_reads': 500, 'mapped_unfiltered': 490, 'mapped_filtered': 490,
         'mapping_rate_pct': 98.0, 'quality_retention_pct': 100.0},
        {'sample': 'barcode02', 'total_reads': 500, 'mapped_unfiltered': 495, 'mapped_filtered': 495,
         'mapping_rate_pct': 99.0, 'quality_retention_pct': 100.0},
    ],
    'parameters': {'min_percent_identity': 0.0, 'min_ref_coverage': 0.0},
    'summary': {
        'total_reads_all': 1000, 'total_mapped_all': 985, 'total_after_quality_filter': 985,
        'overall_mapping_rate_pct': 98.5, 'overall_quality_retention_pct': 100.0,
    },
}


class ReportText(HTMLParser):
    """Collect the text a reader sees and the attributes of the element that holds the note."""

    def __init__(self):
        super().__init__()
        self.hidden = 0
        self.last_attributes = {}
        self.pieces = []
        self.note_attributes = None

    def handle_starttag(self, tag, attrs):
        if tag in ('script', 'style', 'title'):
            self.hidden += 1
        self.last_attributes = dict(attrs)

    def handle_endtag(self, tag):
        if tag in ('script', 'style', 'title'):
            self.hidden -= 1

    def handle_data(self, data):
        if self.hidden:
            return
        self.pieces.append(data)
        if DATABASE_NOTE[:40] in ' '.join(data.split()):
            self.note_attributes = self.last_attributes

    @property
    def visible(self):
        return ' '.join(' '.join(self.pieces).split())


def read_report(folder, with_results):
    """Run the report generator on a small set of result files and parse the page it writes."""
    folder.mkdir()
    (folder / 'report_filtering_0_0.tsv').write_text(FILTER_REPORT, encoding='utf-8')
    command = [sys.executable, str(SOURCES_DIR / 'report_generator.py'),
               '-f', str(folder / 'report_filtering_0_0.tsv'),
               '-otu', str(folder / 'OTU_Table.xlsx'),
               '-o', str(folder / 'report.html')]
    if with_results:
        (folder / 'OTU_Table.tsv').write_text(OTU_TABLE, encoding='utf-8')
        (folder / 'mapping_stats.json').write_text(json.dumps(MAPPING_STATS), encoding='utf-8')
        command += ['-ms', str(folder / 'mapping_stats.json')]
    returncode, output, hung = run_guarded(command, timeout=120)
    assert returncode == 0 and not hung, output
    report = ReportText()
    report.feed((folder / 'report.html').read_text(encoding='utf-8'))
    return report


@pytest.fixture(scope='module')
def report(tmp_path_factory):
    return read_report(tmp_path_factory.mktemp('report') / 'run', with_results=True)


@pytest.fixture(scope='module')
def report_without_results(tmp_path_factory):
    return read_report(tmp_path_factory.mktemp('report') / 'run', with_results=False)


def test_report_ends_with_the_note_on_the_reference_database(report):
    assert 'Full OTU Table' in report.visible
    assert report.visible.count(DATABASE_NOTE) == 1
    assert report.visible.endswith(DATABASE_NOTE)


def test_note_closes_the_report_also_when_no_otu_table_was_produced(report_without_results):
    assert 'No OTU data available' in report_without_results.visible
    assert report_without_results.visible.count(DATABASE_NOTE) == 1
    assert report_without_results.visible.endswith(DATABASE_NOTE)


def test_note_is_printed_in_bold_at_nine_points(report):
    assert report.note_attributes is not None
    declarations = dict(
        (name.strip(), value.strip())
        for name, value in (item.split(':', 1)
                            for item in report.note_attributes.get('style', '').split(';')
                            if ':' in item)
    )
    assert declarations.get('font-size') == '9pt'
    assert declarations.get('font-weight') == 'bold'
