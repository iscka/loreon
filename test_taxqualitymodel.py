"""Tests for TaxQualityModel.py — the unified generative model.

Every test asserts a property the specification demands, so a regression in the
model shows up as a violated statistical guarantee rather than a changed number.
"""

import math
import shutil
import subprocess

import numpy as np
import pytest

import TaxQuality as TQ
import TaxQualityModel as TQM

SAMTOOLS = shutil.which('samtools')
MINIMAP2 = shutil.which('minimap2')
needs_tools = pytest.mark.skipif(
    SAMTOOLS is None or MINIMAP2 is None, reason="samtools/minimap2 not available")


# ---------------------------------------------------------------------------
# Sufficient statistics (§2)
# ---------------------------------------------------------------------------

def test_extract_stats_decomposes_nm():
    """NM = S + I + D, and M columns split into matches and substitutions."""
    st = TQM.extract_stats("10S100M5I3D20S", {'NM': '12'}, 0, 'ref')
    assert st.n_I == 5 and st.n_D == 3
    assert st.n_S == 12 - 5 - 3          # 4 substitutions
    assert st.n_M == 100 - 4             # remaining M columns are matches
    assert st.q_start == 10
    assert st.q_end == 10 + 100 + 5      # M and I consume the query
    assert st.qlen == 10 + 105 + 20


def test_extract_stats_reverse_strand_query_coordinates():
    """On the reverse strand the CIGAR is in reference order; query coordinates
    must be flipped so disjointness is tested in read space."""
    fwd = TQM.extract_stats("10S100M40S", {'NM': '0'}, 0, 'r')
    rev = TQM.extract_stats("10S100M40S", {'NM': '0'}, 0x10, 'r')
    assert fwd.q_start == 10
    assert rev.q_start == 40             # trailing clip becomes the leading one
    assert fwd.qlen == rev.qlen == 150


def test_extract_stats_requires_nm():
    assert TQM.extract_stats("100M", {}, 0, 'r') is None


def test_disjointness():
    a = TQM.AlignmentStats(0, 0, 0, 0, 0, 500, 1000, 'a')
    b = TQM.AlignmentStats(0, 0, 0, 0, 520, 1000, 1000, 'b')
    c = TQM.AlignmentStats(0, 0, 0, 0, 100, 600, 1000, 'c')
    assert TQM._intervals_disjoint(a, b)        #真 split
    assert not TQM._intervals_disjoint(a, c)    # heavy overlap -> not a split


# ---------------------------------------------------------------------------
# Model construction
# ---------------------------------------------------------------------------

def _mk(ref, n_M=900, n_S=30, n_I=35, n_D=35, q0=0, q1=None, qlen=1000):
    q1 = q1 if q1 is not None else qlen
    return TQM.AlignmentStats(n_M, n_S, n_I, n_D, q0, q1, qlen, ref)


def _simple_hs(reads):
    hs = TQM.HypothesisSet()
    for i, (singles, splits) in enumerate(reads):
        hs.add_read(f"q{i}", singles, splits)
    return hs.finalise()


def test_hypothesis_set_layout():
    hs = _simple_hs([([_mk('A'), _mk('B')], []),
                     ([_mk('A')], [(_mk('A', q0=0, q1=500), _mk('B', q0=520, q1=1000))])])
    assert hs.n_reads == 2
    assert hs.n_refs == 2
    assert hs.n_hyp == 2 + 2             # 2 singles + (1 single + 1 chimera)
    assert hs.is_split.tolist() == [False, True]
    assert hs.is_chim.sum() == 1


def test_clip_makes_hypotheses_comparable():
    """A hypothesis covering less of the read must be charged for the rest,
    otherwise a chimera always wins by explaining more columns."""
    hs = _simple_hs([([_mk('A', q0=0, q1=500, qlen=1000)], [])])
    # clip column count is the fifth entry
    assert hs.counts[0, 4] == 500


# ---------------------------------------------------------------------------
# EM guarantees (§7)
# ---------------------------------------------------------------------------

def test_loglik_is_monotone_non_decreasing():
    """Jensen's inequality: EM can never decrease the observed log-likelihood."""
    rng = np.random.default_rng(0)
    reads = []
    for _ in range(200):
        a = _mk('A', n_S=int(rng.integers(20, 40)))
        b = _mk('B', n_S=int(rng.integers(20, 60)))
        reads.append(([a, b], []))
    m = TQM.GenerativeModel(_simple_hs(reads)).fit(max_iter=50, tol=1e-12)
    h = np.array(m.history)
    assert np.all(np.diff(h) >= -1e-6), "log-likelihood decreased"
    assert m.info['monotone']


def test_theta_is_a_distribution_and_counts_conserve_reads():
    reads = [([_mk('A')], []) for _ in range(30)] + \
            [([_mk('A'), _mk('B')], []) for _ in range(20)]
    m = TQM.GenerativeModel(_simple_hs(reads)).fit(tol=1e-12)
    assert m.theta.sum() == pytest.approx(1.0)
    assert m.expected_counts().sum() == pytest.approx(50.0)


def test_error_model_is_estimated_from_the_data():
    """e must move away from its initial value toward the observed operation mix,
    and stay a distribution.  This is the indel-aware part of §4."""
    reads = [([_mk('A', n_M=800, n_S=20, n_I=100, n_D=80)], []) for _ in range(100)]
    m = TQM.GenerativeModel(_simple_hs(reads)).fit(tol=1e-12)
    e = m.e
    assert e.sum() == pytest.approx(1.0)
    # insertions and deletions must dominate substitutions, as in the data
    assert e[2] > e[1] and e[3] > e[1]
    assert e[2] == pytest.approx(100 / 1000, abs=0.02)


def test_error_model_can_be_frozen():
    reads = [([_mk('A', n_I=100)], []) for _ in range(50)]
    m = TQM.GenerativeModel(_simple_hs(reads), estimate_error=False).fit(tol=1e-12)
    assert np.allclose(m.e, TQM.E_INIT)


def test_uniform_init_reaches_indistinguishable_reference():
    """The regression that made the heuristic EM inert: a reference that never
    wins a primary alignment must still be able to receive mass."""
    reads = [([_mk('A'), _mk('B')], []) for _ in range(100)]   # identical stats
    m = TQM.GenerativeModel(_simple_hs(reads), init='uniform').fit(tol=1e-12)
    j = m.hs.ref_index
    assert m.theta[j['B']] == pytest.approx(0.5, abs=1e-6)
    # the specification's prescribed warm start collapses instead
    m2 = TQM.GenerativeModel(_simple_hs(reads), init='primary').fit(tol=1e-12)
    assert m2.theta[j['B']] == pytest.approx(0.0, abs=1e-9)


def test_dirichlet_prior_keeps_rare_references_alive():
    reads = [([_mk('A')], []) for _ in range(200)] + [([_mk('B')], [])]
    m = TQM.GenerativeModel(_simple_hs(reads), prior_alpha=1.0).fit(tol=1e-12)
    assert m.theta.min() > 0


# ---------------------------------------------------------------------------
# Chimera component (§5)
# ---------------------------------------------------------------------------

def _chimeric_read():
    """A realistic chimera: aligning the whole read to one reference is poor
    (the second half does not belong to it), while the two half-length segments
    each align well.  Per-segment error counts must be HALF-read sized — giving
    each segment a full read's worth of indels would make the chimeric
    hypothesis imply twice the errors and it would be correctly rejected."""
    poor_single = _mk('A', n_M=600, n_S=250, n_I=75, n_D=75)
    seg_a = _mk('A', n_M=460, n_S=15, n_I=12, n_D=13, q0=0, q1=500)
    seg_b = _mk('B', n_M=460, n_S=15, n_I=12, n_D=13, q0=500, q1=1000)
    return ([poor_single], [(seg_a, seg_b)])


def test_pi_is_estimated_not_counted():
    """pi must be a fitted parameter, and a read without a split can never be
    chimeric."""
    clean = [([_mk('A')], []) for _ in range(90)]
    split = [_chimeric_read() for _ in range(10)]
    m = TQM.GenerativeModel(_simple_hs(clean + split)).fit(tol=1e-12)
    assert 0.0 < m.pi <= 1.0
    rho = m.chimera_posterior()
    assert np.allclose(rho[:90], 0.0)          # non-split reads: rho must be 0
    assert rho[90:].sum() > 0


def test_chimera_contributes_to_both_parents():
    split = [_chimeric_read() for _ in range(50)]
    m = TQM.GenerativeModel(_simple_hs(split)).fit(tol=1e-12)
    n = m.expected_counts()
    clean = m.clean_counts()
    assert n.sum() > clean.sum()               # chimeras add to two parents
    assert m.hs.n_refs == 2 and (n > 0).all()


def test_single_source_wins_when_it_explains_the_whole_read():
    """Model selection per read (§5): if one reference explains everything with a
    good alignment, the chimeric hypothesis must not be preferred."""
    reads = [([_mk('A', n_M=980, n_S=5, n_I=5, n_D=10)],
              [(_mk('A', n_M=400, n_S=40, q0=0, q1=500),
                _mk('B', n_M=400, n_S=40, q0=500, q1=1000))]) for _ in range(50)]
    m = TQM.GenerativeModel(_simple_hs(reads)).fit(tol=1e-12)
    assert m.chimera_posterior().mean() < 0.5


# ---------------------------------------------------------------------------
# Resolution and uncertainty (§8-§9)
# ---------------------------------------------------------------------------

LIN = {
    'A': ('Fungi', 'Asco', 'Sacch', 'Sacchales', 'Sacchaceae', 'Tuber', 'Tuber alpha'),
    'B': ('Fungi', 'Asco', 'Sacch', 'Sacchales', 'Sacchaceae', 'Tuber', 'Tuber beta'),
}


def test_posterior_resolution_uses_gamma():
    reads = [([_mk('A'), _mk('B')], []) for _ in range(40)]
    m = TQM.GenerativeModel(_simple_hs(reads)).fit(tol=1e-12)
    depths, supports, taxa, _ = TQM.posterior_resolution(m, LIN, tau=0.9)
    # two indistinguishable species of one genus -> consensus at genus
    assert set(depths) == {TQ.GENUS_DEPTH}
    assert taxa[0] == 'Tuber'
    assert supports[0] == pytest.approx(1.0)


def test_bootstrap_intervals_bracket_the_point_estimate():
    reads = [([_mk('A')], []) for _ in range(60)] + \
            [([_mk('B')], []) for _ in range(40)]
    hs = _simple_hs(reads)
    m = TQM.GenerativeModel(hs).fit(tol=1e-12)
    lo, hi, draws = TQM.bootstrap(hs, B=60, seed=1, max_iter=50)
    assert draws.shape == (60, hs.n_refs)
    assert np.all(lo <= m.theta + 1e-9) and np.all(hi >= m.theta - 1e-9)
    assert np.all(hi >= lo)


def test_calibration_coverage_increases_with_nominal_level():
    reads = [([_mk('A')], []) for _ in range(70)] + \
            [([_mk('B')], []) for _ in range(30)]
    hs = _simple_hs(reads)
    _, _, draws = TQM.bootstrap(hs, B=100, seed=2, max_iter=50)
    truth = np.zeros(hs.n_refs)
    truth[hs.ref_index['A']] = 0.7
    truth[hs.ref_index['B']] = 0.3
    cov = TQM.calibration_curve(draws, truth)
    levels = sorted(cov)
    vals = [cov[l] for l in levels]
    assert all(b >= a - 1e-9 for a, b in zip(vals, vals[1:]))   # monotone
    assert cov[0.95] >= cov[0.5]


def test_adaptive_units_reports_the_deepest_supported_node():
    lo = np.array([0.4, 0.0])                 # only A's interval excludes zero
    units = TQM.adaptive_units(lo, ['A', 'B'], LIN)
    assert any(u[0] == TQ.SPECIES_DEPTH for u in units.values())


def test_shannon_ci():
    draws = np.array([[0.5, 0.5], [0.6, 0.4], [0.45, 0.55]])
    s = TQM._shannon_ci(draws)
    assert s['ci_low'] <= s['mean'] <= s['ci_high']
    assert s['mean'] == pytest.approx(0.68, abs=0.02)


def test_empty_input_is_safe():
    m = TQM.GenerativeModel(TQM.HypothesisSet().finalise()).fit()
    assert m.info['n_reads'] == 0


# ---------------------------------------------------------------------------
# End-to-end on a real BAM
# ---------------------------------------------------------------------------

@needs_tools
def test_end_to_end_real_bam(tmp_path):
    """Full fit on real minimap2 output, including the two identical references
    that the heuristic module could only resolve after its init was fixed."""
    import random
    rng = random.Random(11)
    rs = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))

    def mut(s, k):
        s = list(s)
        for p in rng.sample(range(len(s)), k):
            s[p] = rng.choice([b for b in 'ACGT' if b != s[p]])
        return ''.join(s)

    core = rs(900)
    db = tmp_path / 'db.fasta'
    db.write_text(f">identA\n{core}\n>identB\n{core}\n>far\n{mut(core,60)}\n")
    fq = tmp_path / 'r.fastq'
    with open(fq, 'w') as fh:
        for i in range(150):
            r = mut(core, rng.randint(10, 30))
            fh.write(f"@read{i}\n{r}\n+\n{'I'*len(r)}\n")
    tax = tmp_path / 'tax.tsv'
    tax.write_text("OTU_ID\tTaxonomy\n"
                   "identA\tFungi;Asco;S;S;S;Tuber;Tuber alpha\n"
                   "identB\tFungi;Asco;S;S;S;Tuber;Tuber alpha\n"
                   "far\tFungi;Asco;S;S;S;Tuber;Tuber beta\n")
    bam = tmp_path / 'a.bam'
    mm = subprocess.Popen([MINIMAP2, '-ax', 'map-ont', '-N', '20', '-p', '0.8',
                           '-k', '15', '-w', '10', str(db), str(fq)],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess.run([SAMTOOLS, 'sort', '-o', str(bam), '-'], stdin=mm.stdout,
                   check=True, capture_output=True)
    mm.wait()

    res = TQM.run_sample(bam, tax, db_format='silva', em_out=tmp_path / 'e.tsv',
                         quality_out=tmp_path / 'q.tsv', ci_out=tmp_path / 'c.tsv',
                         bootstrap_B=25, samtools=SAMTOOLS)
    assert res['model']['converged']
    assert res['model']['monotone']
    # error model learned from the data must remain a distribution
    assert sum(res['model']['error_model'].values()) == pytest.approx(1.0)
    counts = {}
    for line in (tmp_path / 'e.tsv').read_text().splitlines():
        f = line.split('\t')
        if f[0] != '*':
            counts[f[0]] = float(f[2])
    assert counts['identA'] == pytest.approx(counts['identB'], rel=0.15)
    assert 'uncertainty' in res and (tmp_path / 'c.tsv').exists()
    sh = res['uncertainty']['shannon']
    assert sh['ci_low'] <= sh['mean'] + 1e-9 <= sh['ci_high'] + 1e-9
