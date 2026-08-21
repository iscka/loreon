"""Tests for TaxQuality.py — taxonomic consensus, chimera tiers and EM.

Design note: the earlier version of this suite used hand-written AS values with
gaps of exactly 6 (the softmax temperature), i.e. the one regime where the
weighting looks well behaved.  It passed 18/18 while every defect found in the
audit was live.  Tests here therefore either (a) assert an invariant that must
hold for ANY input, or (b) build the fixture by running the real
minimap2 + samtools commands the pipeline itself uses.
"""

import random
import shutil
import subprocess
import textwrap
from pathlib import Path

import numpy as np
import pytest

import TaxQuality as TQ

SAMTOOLS = shutil.which('samtools')
MINIMAP2 = shutil.which('minimap2')
needs_tools = pytest.mark.skipif(
    SAMTOOLS is None or MINIMAP2 is None, reason="samtools/minimap2 not available")


# ---------------------------------------------------------------------------
# Identifier round-tripping against the DuckDB reformat expression
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("db_format,flat,raw,expected", [
    ('unite', False, "Tuber_melanosporum|MH855774|SH1533864.08FU|reps|k__Fungi",
     "MH855774|SH1533864.08FU"),
    ('unite', True, "Tuber_melanosporum|MH855774|SH1533864.08FU|reps|k__Fungi",
     "SH1533864.08FU"),
    ('eukariome', False, "ACC123;k__x;p__y", "ACC123"),
    ('cbs', False, "ACC123|extra|more", "ACC123"),
    ('none', False, "ACC123|extra", "ACC123"),
    ('silva', False, "SILVA_ACC.1", "SILVA_ACC.1"),
])
def test_reformat_ref_id(db_format, flat, raw, expected):
    assert TQ.reformat_ref_id(raw, db_format, flat) == expected


@pytest.mark.parametrize("db_format,flat", [
    ('unite', False), ('unite', True), ('eukariome', False),
    ('cbs', False), ('none', False), ('silva', False),
])
def test_reformat_matches_duckdb(db_format, flat):
    """reformat_ref_id must agree with ResultsReader's SQL for every format.

    A mismatch silently breaks the quality-column join (this is exactly how
    --flat-sh produced all-NaN quality columns)."""
    duckdb = pytest.importorskip("duckdb")
    import ResultsReader as RR
    raws = {
        'unite': ["Tuber_melanosporum|MH855774|SH1533864.08FU|reps|k__Fungi;p__Asco",
                  "Xx|AC1|SH9.01FU|reps|k__Fungi", "no_pipes_at_all"],
        'eukariome': ["ACC123;k__x;p__y", "ACC123"],
        'cbs': ["ACC123|extra|more", "ACC123"],
        'none': ["ACC123|extra", "ACC123"],
        'silva': ["SILVA_ACC.1", "ACC 123"],
    }[db_format]
    expr = RR._build_reformat_sql(db_format, flat_sh=flat)
    con = duckdb.connect()
    for raw in raws:
        esc = raw.replace("'", "''")
        sql = f"SELECT {expr} AS r FROM (SELECT '{esc}' AS otu_id)"
        assert con.execute(sql).fetchone()[0] == TQ.reformat_ref_id(raw, db_format, flat), raw
    con.close()


# ---------------------------------------------------------------------------
# Lineage splitting — positional, gaps preserved
# ---------------------------------------------------------------------------

def test_split_lineage_is_positional():
    """An empty intermediate rank must NOT shift deeper ranks up."""
    lin = TQ._split_lineage("Fungi;Ascomycota;;;Tuberaceae;Tuber;Tuber_melanosporum")
    assert lin == ('Fungi', 'Ascomycota', '', '', 'Tuberaceae', 'Tuber',
                   'Tuber melanosporum')
    assert lin[TQ.GENUS_DEPTH - 1] == 'Tuber'      # genus stays findable


def test_split_lineage_pads_and_underscores():
    assert TQ._split_lineage("Fungi;Ascomycota;Dothideomycetes") == \
        ('Fungi', 'Ascomycota', 'Dothideomycetes', '', '', '', '')
    assert TQ._split_lineage("A;B;C;D;E;F;Genus_species")[-1] == 'Genus species'


def test_split_lineage_matches_resultsreader_duckdb():
    """The module and the OTU table must agree on where each rank lives."""
    duckdb = pytest.importorskip("duckdb")
    import ResultsReader as RR
    tax = "Fungi;Ascomycota;;;Tuberaceae;Tuber;Tuber_melanosporum"
    con = duckdb.connect()
    row = con.execute(
        f"SELECT {RR._build_taxonomy_split_sql()} FROM (SELECT '{tax}' AS Taxonomy)"
    ).df().iloc[0]
    con.close()
    lin = TQ._split_lineage(tax)
    assert lin[4] == row['Family'] == 'Tuberaceae'
    assert lin[5] == row['Genus'] == 'Tuber'


# ---------------------------------------------------------------------------
# CIGAR / identity / coverage
# ---------------------------------------------------------------------------

def test_query_length_and_clips():
    assert TQ.query_length_from_cigar("10S100M5I20D5S") == 120
    assert TQ.query_length_from_cigar("10H100M10H") == 120
    assert TQ.soft_clip_length("10S100M5S") == 15
    assert TQ.aligned_query_length("10S100M5I5S") == 105


def test_samtools_identity_matches_htslib_expression():
    """Our gate must reproduce ((qlen-sclen)-NM)/(qlen-sclen) exactly."""
    cigar, nm = "10S280M10S", 14
    aln = 280
    assert TQ.samtools_identity({'NM': str(nm)}, cigar) == pytest.approx((aln - nm) / aln)
    assert TQ.samtools_identity({}, cigar) is None


def test_query_coverage():
    assert TQ.query_coverage("10S280M10S") == pytest.approx(280 / 300)
    assert TQ.query_coverage("300M") == 1.0


# ---------------------------------------------------------------------------
# Consensus
# ---------------------------------------------------------------------------

LIN = {
    'spA': ('Fungi', 'Asco', 'Pezizo', 'Pezizales', 'Tuberaceae', 'Tuber', 'Tuber melanosporum'),
    'spB': ('Fungi', 'Asco', 'Pezizo', 'Pezizales', 'Tuberaceae', 'Tuber', 'Tuber aestivum'),
    'spC': ('Fungi', 'Asco', 'Pezizo', 'Pezizales', 'Tuberaceae', 'Choiromyces', 'Choiromyces venosus'),
    'spZ': ('Bacteria', 'Proteo', 'Gamma', 'Entero', 'Enterobacteriaceae', 'Escherichia', 'Escherichia coli'),
    # same species as spA but the database is missing the Class field
    'spA_gap': ('Fungi', 'Asco', '', 'Pezizales', 'Tuberaceae', 'Tuber', 'Tuber melanosporum'),
    # sibling with no species epithet in the DB
    'spNoSp': ('Fungi', 'Asco', 'Pezizo', 'Pezizales', 'Tuberaceae', 'Tuber', ''),
}


def test_consensus_single_candidate():
    t, d, s, c, m = TQ.weighted_consensus(['spA'], np.array([1.0]), LIN)
    assert (t, d) == ('Tuber melanosporum', TQ.SPECIES_DEPTH)
    assert s == pytest.approx(1.0)


def test_consensus_same_genus_collapses_and_names_the_species_that_broke_it():
    """Falling back to genus must still say WHICH species stole the call."""
    t, d, s, c, m = TQ.weighted_consensus(['spA', 'spB'], np.array([0.55, 0.45]), LIN)
    assert (t, d) == ('Tuber', TQ.GENUS_DEPTH)
    assert c == 'Tuber aestivum'          # not '' and not a genus
    assert m == pytest.approx(0.45)


def test_consensus_climbs_to_family():
    t, d, _, _, _ = TQ.weighted_consensus(['spA', 'spC'], np.array([0.5, 0.5]), LIN)
    assert (t, d) == ('Tuberaceae', 5)


def test_consensus_missing_intermediate_rank_does_not_demote():
    """Two entries for the same species, one missing Class, must stay Species."""
    t, d, _, _, _ = TQ.weighted_consensus(
        ['spA', 'spA_gap'], np.array([0.5, 0.5]), LIN)
    assert (t, d) == ('Tuber melanosporum', TQ.SPECIES_DEPTH)


def test_consensus_unannotated_rank_is_not_counted_as_disagreement():
    """A sibling with no species name in the DB is missing data, not conflict."""
    t, d, s, _, _ = TQ.weighted_consensus(
        ['spA', 'spNoSp'], np.array([0.85, 0.15]), LIN)
    assert (t, d) == ('Tuber melanosporum', TQ.SPECIES_DEPTH)
    assert s == pytest.approx(1.0)        # 100% of the ANNOTATED mass


def test_consensus_confounder_is_never_the_consensus_itself():
    t, d, _, c, _ = TQ.weighted_consensus(
        ['spA', 'spA_gap'], np.array([0.5, 0.5]), LIN)
    assert c != t


# ---------------------------------------------------------------------------
# Resolution profile
# ---------------------------------------------------------------------------

def test_resolution_profile_partitions_all_reads():
    """pct_* must account for 100% of reads at every depth."""
    acc = TQ.OtuQualityAccumulator()
    for depth in (7, 7, 6, 5, 3, 0):
        acc.add('otu', depth, 1.0, 60, '', 0.0)
    m = acc.result()['otu']
    total = (m['pct_species'] + m['pct_genus'] +
             m['pct_above_genus'] + m['pct_unresolved'])
    assert total == pytest.approx(100.0)


def test_resolution_rank_is_modal_not_mean():
    """A bimodal OTU must not be reported at a rank no read supports."""
    acc = TQ.OtuQualityAccumulator()
    for _ in range(50):
        acc.add('otu', 7, 1.0, 60, '', 0.0)   # Species
    for _ in range(49):
        acc.add('otu', 1, 1.0, 60, '', 0.0)   # Kingdom
    m = acc.result()['otu']
    assert m['resolution_rank'] == 'Species'   # never 'Order'
    assert m['pct_species'] == pytest.approx(50.51, abs=0.02)


def test_confounder_counts_confounded_reads_only():
    acc = TQ.OtuQualityAccumulator()
    for _ in range(999):
        acc.add('otu', 7, 1.0, 60, '', 0.0)
    acc.add('otu', 7, 0.6, 60, 'Fusarium', 0.4)
    m = acc.result()['otu']
    assert m['n_reads'] == 1000
    assert m['confounded_reads'] == 1        # not 1000


# ---------------------------------------------------------------------------
# EM
# ---------------------------------------------------------------------------

def _acc(reads):
    """reads = list of (primary_raw, [(ref, score), ...])"""
    a = TQ.EMAccumulator()
    for primary, cands in reads:
        rc = TQ.ReadCandidates('q')
        rc.primary_raw = primary
        rc.primary_ref = primary
        for ref, sc in cands:
            rc.refs_raw.append(ref); rc.refs_canon.append(ref); rc.scores.append(sc)
        a.add(rc)
    return a


def test_em_conserves_read_mass():
    """Invariant: expected counts must sum to the number of reads."""
    reads = ([('A', [('A', 600.0)])] * 10 +
             [('B', [('B', 600.0)])] * 10 +
             [('A', [('A', 600.0), ('B', 600.0)])] * 20)
    counts, info = TQ.run_em(_acc(reads), tol=1e-12)
    assert sum(counts.values()) == pytest.approx(40.0)
    assert info['counts_total'] == pytest.approx(40.0)
    assert info['reads_unassignable'] == 0


def test_em_assigns_mass_to_a_reference_that_never_wins_a_primary():
    """The regression that made the EM inert.

    Two references are indistinguishable (identical AS) and minimap2 always
    tie-breaks to A.  With a primary warm start theta_B starts at 0 and can
    never recover; with the uniform default it must."""
    reads = [('A', [('A', 600.0), ('B', 600.0)])] * 100
    counts, _ = TQ.run_em(_acc(reads), tol=1e-12)
    assert counts['B'] > 0
    assert counts['A'] == pytest.approx(50.0, abs=1e-6)
    assert counts['B'] == pytest.approx(50.0, abs=1e-6)
    # the old behaviour, kept only for comparison, still collapses
    legacy, _ = TQ.run_em(_acc(reads), tol=1e-12, init='primary')
    assert legacy['B'] == pytest.approx(0.0)


def test_em_anchors_pull_ambiguous_mass():
    reads = ([('A', [('A', 600.0)])] * 30 +
             [('B', [('B', 600.0)])] * 2 +
             [('A', [('A', 600.0), ('B', 600.0)])] * 10)
    counts, _ = TQ.run_em(_acc(reads), max_iter=500, tol=1e-12)
    assert counts['A'] > counts['B']
    assert counts['A'] + counts['B'] == pytest.approx(42.0)


def test_em_convergence_criterion_is_scale_free():
    """Same generative model at two sample sizes must converge alike."""
    def build(n):
        return _acc([('A', [('A', 600.0), ('B', 594.0)])] * n +
                    [('A', [('A', 600.0)])] * n)
    _, small = TQ.run_em(build(50), max_iter=500)
    _, big = TQ.run_em(build(5000), max_iter=500)
    assert small['converged'] and big['converged']
    assert abs(small['iterations'] - big['iterations']) <= 5


def test_em_loglik_describes_the_returned_counts():
    """final_loglik must be the likelihood of the parameters actually returned."""
    a = _acc([('A', [('A', 600.0), ('B', 594.0)])] * 20 + [('A', [('A', 600.0)])] * 20)
    counts, info = TQ.run_em(a, max_iter=3)
    row_ref, row_w, row_ptr, _ = a.arrays()
    theta = np.array([counts[r] for r in sorted(a.ref_index, key=a.ref_index.get)])
    theta = theta / theta.sum()
    denom = np.add.reduceat(theta[row_ref] * row_w, row_ptr[:-1])
    assert info['final_loglik'] == pytest.approx(float(np.sum(np.log(denom))), abs=1e-9)


def test_em_empty():
    counts, info = TQ.run_em(TQ.EMAccumulator())
    assert counts == {} and info['n_reads'] == 0


def test_write_em_tsv_preserves_fractional_mass(tmp_path):
    """A long tail of sub-1.0 assignments must not be silently discarded."""
    counts = {f'r{i}': 0.45 for i in range(100)}
    counts['big'] = 55.0
    p, n, kept, dropped = TQ.write_em_tsv(counts, tmp_path / 'x.em.tsv')
    assert n == 101
    assert kept == pytest.approx(100.0)
    assert dropped == pytest.approx(0.0)
    body = [l for l in p.read_text().splitlines() if not l.startswith('*')]
    assert any(l.split('\t')[2].startswith('0.45') for l in body)


# ---------------------------------------------------------------------------
# Chimera tiers
# ---------------------------------------------------------------------------

def _rc(primary, supp=()):
    rc = TQ.ReadCandidates('q')
    rc.primary_ref = primary
    rc.primary_raw = primary
    rc.refs_raw = [primary]; rc.refs_canon = [primary]; rc.scores = [600.0]
    rc.supp_refs = set(supp)
    return rc


def test_chimera_cross_taxon_only():
    acc = TQ.ChimeraAccumulator(LIN, refined=True)
    acc.add(_rc('spA'))                 # clean
    acc.add(_rc('spA', ['spZ']))        # cross-genus -> chimera
    acc.add(_rc('spA', ['spB']))        # same genus -> not a chimera
    r = acc.result()
    assert r['n_mapped_reads'] == 3
    assert r['split_reads'] == 2
    assert r['n_candidate_chimeras'] == 1


def test_chimera_works_when_lineage_has_gaps():
    """Regression: a left-shifted lineage blanked Genus and silently scored 0."""
    lin = {'p': TQ._split_lineage("Fungi;Ascomycota;;;Tuberaceae;Tuber;Tuber_mel"),
           's': TQ._split_lineage("Bacteria;Proteo;;;Enterobacteriaceae;Escherichia;E_coli")}
    acc = TQ.ChimeraAccumulator(lin, refined=True)
    acc.add(_rc('p', ['s']))
    assert acc.result()['n_candidate_chimeras'] == 1


def test_chimera_reports_when_genus_is_unknown():
    lin = {'p': ('Fungi',) + ('',) * 6, 's': ('Bacteria',) + ('',) * 6}
    acc = TQ.ChimeraAccumulator(lin, refined=True)
    acc.add(_rc('p', ['s']))
    r = acc.result()
    assert r['n_candidate_chimeras'] == 0
    assert r['split_reads_genus_unknown'] == 1   # blind spot is surfaced


# ---------------------------------------------------------------------------
# Threshold handling
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("pct,expected", [
    (0, 0.0), (95.0, 0.95), (1, 0.01), (0.5, 0.005), (100, 1.0)])
def test_as_fraction_is_monotone(pct, expected):
    """`1` means 1%, not 100% (the old heuristic inverted small values)."""
    assert TQ._as_fraction(pct) == pytest.approx(expected)


# ---------------------------------------------------------------------------
# End-to-end against REAL minimap2 + samtools output
# ---------------------------------------------------------------------------

def _build_real_bam(tmp_path, n_reads=120, seed=5):
    """Two byte-identical accessions + one divergent, with ONT-like reads."""
    rng = random.Random(seed)
    rs = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))

    def mut(s, k):
        s = list(s)
        for p in rng.sample(range(len(s)), k):
            s[p] = rng.choice([b for b in 'ACGT' if b != s[p]])
        return ''.join(s)

    core = rs(650)
    far = mut(core, 40)
    db = tmp_path / 'db.fasta'
    # plain names: with db_format='silva' the RNAME is used verbatim, so it
    # matches the taxonomy-map key below
    db.write_text(f">identA\n{core}\n>identB\n{core}\n>far\n{far}\n")
    fq = tmp_path / 'r.fastq'
    with open(fq, 'w') as fh:
        for i in range(n_reads):
            r = mut(core, rng.randint(5, 20))
            fh.write(f"@read{i}\n{r}\n+\n{'I' * len(r)}\n")
    tax = tmp_path / 'tax.tsv'
    tax.write_text(
        "OTU_ID\tTaxonomy\n"
        "identA\tFungi;Asco;Pezizo;Pezizales;Tuberaceae;Tuber;Tuber_alpha\n"
        "identB\tFungi;Asco;Pezizo;Pezizales;Tuberaceae;Tuber;Tuber_alpha\n"
        "far\tFungi;Asco;Pezizo;Pezizales;Tuberaceae;Tuber;Tuber_beta\n")
    bam = tmp_path / 'aln.bam'
    mm = subprocess.Popen([MINIMAP2, '-ax', 'map-ont', '-N', '20', '-p', '0.5',
                           '-k', '15', '-w', '10', str(db), str(fq)],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess.run([SAMTOOLS, 'sort', '-o', str(bam), '-'], stdin=mm.stdout,
                   check=True, capture_output=True)
    mm.wait()
    return bam, tax


@needs_tools
def test_real_bam_streaming_groups_reads(tmp_path):
    """iter_read_groups must deliver every record of a read together even
    though the BAM is coordinate-sorted."""
    bam, _ = _build_real_bam(tmp_path)
    seen = {}
    for qname, recs in TQ.iter_read_groups(bam, samtools=SAMTOOLS):
        assert qname not in seen, "read group split across yields"
        seen[qname] = len(recs)
    assert len(seen) == 120
    assert max(seen.values()) >= 2, "expected secondary alignments in the fixture"


@needs_tools
def test_real_bam_em_splits_identical_references(tmp_path):
    """On real alignments, two identical accessions must share the mass."""
    bam, tax = _build_real_bam(tmp_path)
    res = TQ.run_sample(bam, tax, db_format='silva', min_percent_identity=0,
                        min_ref_coverage=0, em_out=tmp_path / 'e.tsv',
                        quality_out=tmp_path / 'q.tsv', samtools=SAMTOOLS)
    counts = {}
    for line in (tmp_path / 'e.tsv').read_text().splitlines():
        f = line.split('\t')
        if f[0] != '*':
            counts[f[0]] = float(f[2])
    assert counts['identA'] == pytest.approx(counts['identB'], rel=0.05)
    assert res['em']['counts_total'] == pytest.approx(res['em']['n_reads'], rel=1e-6)
    assert res['em']['reads_unassignable'] == 0


@needs_tools
def test_real_bam_quality_gate_matches_samtools(tmp_path):
    """The reads kept by the read-level pass must equal those kept by the
    htslib expression the main table uses."""
    bam, tax = _build_real_bam(tmp_path)
    ident, cov = 90.0, 80.0
    expr = (f"(({{}}qlen-sclen)-[NM])*100>={ident}*(qlen-sclen) && "
            f"(qlen-sclen)*100>={cov}*qlen").replace('{}', '')
    out = subprocess.run([SAMTOOLS, 'view', '-c', '-F', '0x904', '-e', expr, str(bam)],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    expected = int(out.stdout.strip())
    res = TQ.run_sample(bam, tax, db_format='silva', min_percent_identity=ident,
                        min_ref_coverage=cov, do_em=True, do_quality=True,
                        em_out=tmp_path / 'e.tsv', samtools=SAMTOOLS)
    assert res['em']['n_reads'] == expected
    assert res['taxonomic_quality']['n_reads'] == expected
