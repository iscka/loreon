"""Regression tests for OtuUtils quality-filter and flagstat parsing.

Covers two fixes:
  * `_parse_flagstat` now also returns the `supplementary` count (Feature C
    Tier 1).
  * the alignment quality filter expression must be a VALID htslib filter
    (the old `alen` token was rejected by samtools >=1.16, silently emptying
    the OTU table for the default 95/90 thresholds).
"""

import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

import OtuUtils


FLAGSTAT_SAMPLE = """\
1200 + 0 in total (QC-passed reads + QC-failed reads)
1000 + 0 primary
150 + 0 secondary
50 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1100 + 0 mapped (91.67% : N/A)
950 + 0 primary mapped (95.00% : N/A)
"""


def test_parse_flagstat_captures_supplementary():
    stats = OtuUtils._parse_flagstat(FLAGSTAT_SAMPLE)
    assert stats['total'] == 1000          # primary
    assert stats['mapped'] == 950          # primary mapped
    assert stats['supplementary'] == 50


SAMTOOLS = shutil.which('samtools')


@pytest.mark.skipif(SAMTOOLS is None, reason="samtools not available")
def test_quality_filter_expression_is_valid(tmp_path):
    """The exact -e expression built by tabeling_improved must be accepted by
    samtools (guards against the `alen` regression)."""
    sam = tmp_path / "in.sam"
    sam.write_text(textwrap.dedent("""\
        @HD\tVN:1.6\tSO:coordinate
        @SQ\tSN:refA\tLN:300
        q1\t0\trefA\t1\t60\t300M\t*\t0\t0\t*\t*\tNM:i:3\tAS:i:600
        q2\t0\trefA\t1\t60\t10S280M10S\t*\t0\t0\t*\t*\tNM:i:2\tAS:i:560
    """))
    bam = tmp_path / "in.bam"
    subprocess.run([SAMTOOLS, "sort", "-o", str(bam), str(sam)], check=True,
                   capture_output=True)

    # Rebuild the same expression tabeling_improved uses (identity + coverage).
    expr = ("((qlen-sclen)-[NM])*100>=95.0*(qlen-sclen) && "
            "(qlen-sclen)*100>=90.0*qlen")
    r = subprocess.run([SAMTOOLS, "view", "-c", "-F", "0x904", "-e", expr, str(bam)],
                       capture_output=True, text=True)
    assert r.returncode == 0, f"filter expression rejected: {r.stderr}"
    # q1 (100% id, full cover) passes; q2 (280/300 aligned = 93% cover) passes too.
    assert int(r.stdout.strip()) == 2


@pytest.mark.skipif(SAMTOOLS is None, reason="samtools not available")
def test_tabeling_improved_end_to_end_with_default_filter(tmp_path):
    """tabeling_improved with the default 95/90 filter must produce counts
    (previously failed outright due to the invalid `alen` token)."""
    mapping_dir = tmp_path / "mapping"
    mapping_dir.mkdir()
    sam = mapping_dir / "s1.sam"
    sam.write_text(textwrap.dedent("""\
        @HD\tVN:1.6\tSO:coordinate
        @SQ\tSN:refA\tLN:300
        @SQ\tSN:refB\tLN:300
        r1\t0\trefA\t1\t60\t300M\t*\t0\t0\t*\t*\tNM:i:1\tAS:i:598
        r2\t0\trefA\t1\t60\t300M\t*\t0\t0\t*\t*\tNM:i:2\tAS:i:596
        r3\t0\trefB\t1\t60\t300M\t*\t0\t0\t*\t*\tNM:i:0\tAS:i:600
    """))
    bam = mapping_dir / "s1_sorted_full.bam"
    subprocess.run([SAMTOOLS, "sort", "-o", str(bam), str(sam)], check=True,
                   capture_output=True)

    stats = OtuUtils.tabeling_improved(str(bam), threads=1,
                                       min_percent_identity=95.0,
                                       min_ref_coverage=90.0)
    assert stats is not None
    assert stats['mapped_filtered'] == 3
    results = tmp_path / "tabeling" / "results" / "s1.tmp.txt"
    assert results.exists()
    counts = {}
    for line in results.read_text().splitlines():
        parts = line.split('\t')
        if parts[0] not in ('*',) and len(parts) >= 3:
            counts[parts[0]] = int(parts[2])
    assert counts.get('refA') == 2
    assert counts.get('refB') == 1
