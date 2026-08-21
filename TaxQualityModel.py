#!/usr/bin/env python3
"""TaxQualityModel.py — Unified generative model for long-read amplicon
quantification (LOREON / TaxQuality).

This module implements the model formalised in
``LOREON_TaxQuality_specifica_formale.tex``.  It is a SEPARATE, OPT-IN
alternative to the heuristic pipeline in ``TaxQuality.py``: nothing here runs
unless the user explicitly asks for it (``--em-model generative``).  The
heuristic module is untouched and remains the default.

What the specification requires, and where it lives here
--------------------------------------------------------
  §4  indel-aware error model, e = (e_M, e_S, e_I, e_D) estimated from the data
        -> ``AlignmentStats``, ``HypothesisSet.log_likelihood``, ``_m_step_error``
  §5  chimera component inside the likelihood, pi estimated (not counted)
        -> ``HypothesisSet.build`` (disjointness enforced), ``_m_step_pi``
  §6  joint mixture likelihood, eq. (5)
        -> ``GenerativeModel.observed_loglik``
  §7  EM: E-step eq. (6)-(7), M-step eq. (9)-(11)
        -> ``GenerativeModel.fit``
  §7.3 Dirichlet MAP
        -> ``prior_alpha``
  §8  posterior taxonomic resolution d_r, adaptive units
        -> ``posterior_resolution``, ``adaptive_units``
  §9  bootstrap credible intervals and calibration
        -> ``bootstrap``, ``calibration_curve``

Difference from the heuristic module, in one line: there the per-candidate
weights are a softmax over minimap2's AS with a fixed temperature; here they are
true posterior responsibilities of a model whose error parameters are estimated
from the data.

Deliberate deviation from the specification
-------------------------------------------
§7.4 prescribes warm-starting theta from normalised primary counts.  That
initialisation has an absorbing-zero flaw: a reference that never wins a primary
alignment gets theta_j = 0, and since gamma is proportional to theta it can never
recover, no matter how strong its alignment evidence.  This is exactly the
multi-mapping case the model exists to resolve.  The default here is uniform
(as kallisto/RSEM/Bracken); ``init='primary'`` reproduces the specification's
behaviour for comparison.

Extension beyond the specification
-----------------------------------
Soft/hard-clipped query bases are charged a fixed penalty (the alternative
mentioned in the §4 note).  Without it a chimeric hypothesis, which explains more
of the read, would always beat a single-source hypothesis simply by covering more
columns; with it every hypothesis is scored over the FULL read length, so the
per-read model selection in §5 is a fair comparison.  The penalty is fixed, not
fitted: when it was allowed to float as a fifth category of the estimated
simplex it absorbed 71% of the mass on real ITS data and inverted the incentives.
"""

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

import TaxQuality as TQ

# Estimated error model: a categorical over ALIGNED columns only (specification
# eq. 1).  Clipped query bases are handled by a separate FIXED penalty and are
# deliberately NOT part of this simplex.
OPS = ('M', 'S', 'I', 'D')
N_OPS = len(OPS)
N_COUNT_COLS = N_OPS + 1          # the fifth column holds clipped bases

# Literature starting point for ONT R10 (indel-dominated); re-estimated by EM.
E_INIT = np.array([0.935, 0.020, 0.022, 0.023])

# Probability assigned to a query base that the hypothesis does not align at all.
# It must be SMALL and FIXED.  Letting it float inside the simplex lets it absorb
# most of the mass on real data (measured: e_C = 0.71, e_M = 0.27), which inverts
# the incentives — matching a base becomes more expensive than leaving it
# unexplained — and collapses the chimera component.  A clipped base is about as
# unlikely as a substituted one, which is what this default encodes.
CLIP_PENALTY = 0.01

_LOG13 = math.log(1.0 / 3.0)     # choice of substituted base
_LOG14 = math.log(1.0 / 4.0)     # choice of inserted base

RANK_NAMES = TQ.RANK_NAMES
N_RANKS = TQ.N_RANKS


# ---------------------------------------------------------------------------
# Sufficient statistics (§2)
# ---------------------------------------------------------------------------

class AlignmentStats:
    """Per-alignment sufficient statistics and query interval."""
    __slots__ = ('n_M', 'n_S', 'n_I', 'n_D', 'q_start', 'q_end', 'qlen', 'ref')

    def __init__(self, n_M, n_S, n_I, n_D, q_start, q_end, qlen, ref):
        self.n_M, self.n_S, self.n_I, self.n_D = n_M, n_S, n_I, n_D
        self.q_start, self.q_end, self.qlen, self.ref = q_start, q_end, qlen, ref

    @property
    def counts(self):
        """(M,S,I,D) — clip is added later, relative to the whole read."""
        return (self.n_M, self.n_S, self.n_I, self.n_D)


def extract_stats(cigar, tags, flag, ref):
    """Sufficient statistics from a CIGAR + NM tag.

    NM = n_S + n_I + n_D, so substitutions are recovered as NM - I - D and
    matches as (aligned M columns) - n_S.  Query coordinates are returned in
    ORIGINAL read orientation, so that disjointness between a primary and its
    supplementary segment can be tested (specification §5).
    """
    if not cigar or cigar == '*':
        return None
    n_I = n_D = m_cols = 0
    lead_clip = 0
    trail_clip = 0
    seen_aligned = False
    for length, op in TQ._parse_cigar(cigar):
        if op in 'SH':
            if seen_aligned:
                trail_clip += length
            else:
                lead_clip += length
        elif op == 'I':
            n_I += length; seen_aligned = True
        elif op == 'D' or op == 'N':
            n_D += length; seen_aligned = True
        elif op in 'M=X':
            m_cols += length; seen_aligned = True
    try:
        nm = int(tags['NM'])
    except (KeyError, ValueError):
        return None
    n_S = max(0, nm - n_I - n_D)
    n_M = max(0, m_cols - n_S)

    aligned_q = m_cols + n_I
    qlen = lead_clip + aligned_q + trail_clip
    if flag & 0x10:            # reverse strand: CIGAR is in reference order
        q_start = trail_clip
    else:
        q_start = lead_clip
    return AlignmentStats(n_M, n_S, n_I, n_D, q_start, q_start + aligned_q, qlen, ref)


def _intervals_disjoint(a, b, max_overlap_frac=0.2):
    """Two segments explain different parts of the read (specification §5).

    The heuristic module omits this test, which is why it reported a 72%
    inter-taxon chimera rate on real ITS data: roughly half of the multi-segment
    reads had OVERLAPPING alignments, which are repeats or alternative partial
    alignments, not split structure.
    """
    inter = max(0, min(a.q_end, b.q_end) - max(a.q_start, b.q_start))
    shorter = min(a.q_end - a.q_start, b.q_end - b.q_start)
    if shorter <= 0:
        return False
    return inter < max_overlap_frac * shorter


# ---------------------------------------------------------------------------
# Hypothesis set (§6)
# ---------------------------------------------------------------------------

class HypothesisSet:
    """Flattened per-read hypotheses: singles (Cand_r) + chimeras (Split_r).

    Layout is CSR-like so both EM steps are single vectorised sweeps:
      counts  (H, 5)  operation counts per hypothesis, clip included
      const   (H,)    base-choice terms, constant across iterations
      ref1/2  (H,)    reference indices (ref2 = -1 for singles)
      is_chim (H,)    hypothesis type
      row_ptr (R+1,)  hypothesis span of each read
      is_split(R,)    whether the read has a usable split (indicator in eq. 5)
    """

    def __init__(self):
        self.ref_index = {}
        self.ref_names = []
        self._counts = []
        self._const = []
        self._ref1 = []
        self._ref2 = []
        self._is_chim = []
        self.row_ptr = [0]
        self.is_split = []
        self.primary_ref = []
        self.qnames = []
        self.n_dropped_overlap = 0

    def _ref_id(self, name):
        j = self.ref_index.get(name)
        if j is None:
            j = len(self.ref_names)
            self.ref_index[name] = j
            self.ref_names.append(name)
        return j

    def add_read(self, qname, singles, splits):
        """singles: list[AlignmentStats]; splits: list[(AlignmentStats, AlignmentStats)]"""
        if not singles:
            return False
        qlen = max(s.qlen for s in singles)
        n0 = len(self._counts)

        for st in singles:
            clip = max(0, qlen - (st.q_end - st.q_start))
            self._counts.append((st.n_M, st.n_S, st.n_I, st.n_D, clip))
            self._const.append(st.n_S * _LOG13 + st.n_I * _LOG14)
            self._ref1.append(self._ref_id(st.ref))
            self._ref2.append(-1)
            self._is_chim.append(False)

        for a, b in splits:
            cov = (a.q_end - a.q_start) + (b.q_end - b.q_start)
            clip = max(0, qlen - cov)
            self._counts.append((a.n_M + b.n_M, a.n_S + b.n_S,
                                 a.n_I + b.n_I, a.n_D + b.n_D, clip))
            self._const.append((a.n_S + b.n_S) * _LOG13 + (a.n_I + b.n_I) * _LOG14)
            self._ref1.append(self._ref_id(a.ref))
            self._ref2.append(self._ref_id(b.ref))
            self._is_chim.append(True)

        if len(self._counts) == n0:
            return False
        self.row_ptr.append(len(self._counts))
        self.is_split.append(bool(splits))
        self.primary_ref.append(self._ref_id(singles[0].ref))
        self.qnames.append(qname)
        return True

    def finalise(self):
        self.counts = np.asarray(self._counts, dtype=np.float64) if self._counts \
            else np.zeros((0, N_COUNT_COLS))
        self.const = np.asarray(self._const, dtype=np.float64)
        self.ref1 = np.asarray(self._ref1, dtype=np.int64)
        self.ref2 = np.asarray(self._ref2, dtype=np.int64)
        self.is_chim = np.asarray(self._is_chim, dtype=bool)
        self.row_ptr = np.asarray(self.row_ptr, dtype=np.int64)
        self.is_split = np.asarray(self.is_split, dtype=bool)
        self.primary_ref = np.asarray(self.primary_ref, dtype=np.int64)
        self.seg_len = np.diff(self.row_ptr)
        self._counts = self._const = self._ref1 = self._ref2 = self._is_chim = None
        return self

    @property
    def n_reads(self):
        return len(self.row_ptr) - 1

    @property
    def n_refs(self):
        return len(self.ref_names)

    @property
    def n_hyp(self):
        return self.counts.shape[0]

    def log_likelihood(self, log_e, clip_penalty=CLIP_PENALTY):
        """log L(r | h) for every hypothesis — specification eq. (1), plus the
        soft-clip penalty term of the §4 note.

        ``log_e`` covers the four aligned operations; the clip column is scored
        with a fixed penalty so that every hypothesis is evaluated over the full
        read and single-source vs chimeric hypotheses stay comparable."""
        full = np.concatenate([log_e, [math.log(clip_penalty)]])
        return self.counts @ full + self.const


def build_hypotheses(bam_path, lineages, db_format='unite', flat_sh=False,
                     threads=1, min_percent_identity=0.0, min_ref_coverage=0.0,
                     max_overlap_frac=0.2, samtools='samtools', debug=False):
    """One streaming pass over the BAM producing the hypothesis set.

    Reuses the name-collated reader from the heuristic module so the two share
    exactly one BAM-parsing implementation.
    """
    hs = HypothesisSet()
    min_id = TQ._as_fraction(min_percent_identity)
    min_cov = TQ._as_fraction(min_ref_coverage)
    n_seen = n_gated = 0

    for qname, records in TQ.iter_read_groups(bam_path, threads=threads,
                                              samtools=samtools):
        n_seen += 1
        primary = None
        singles = {}
        supps = []
        gate_ok = True
        for flag, rname, mapq, cigar, tags in records:
            # Store the RAW reference name: write_em_tsv output is consumed by
            # makeOtu_duckdb, which applies the format-specific reformat itself.
            # Writing already-canonical ids would make DuckDB canonicalise twice
            # and the taxonomy join would match nothing.
            st = extract_stats(cigar, tags, flag, rname)
            if st is None:
                continue
            if flag & 0x800:
                supps.append(st)
                continue
            ident = TQ.samtools_identity(tags, cigar)
            cov = TQ.query_coverage(cigar)
            if not (flag & 0x100):
                primary = st
                ok_id = (min_id <= 0) or (ident is None) or (ident >= min_id)
                gate_ok = ok_id and ((min_cov <= 0) or (cov >= min_cov))
            else:
                if min_id > 0 and ident is not None and ident < min_id:
                    continue
                if min_cov > 0 and cov < min_cov:
                    continue
            prev = singles.get(st.ref)
            if prev is None or (st.n_M - st.n_S) > (prev.n_M - prev.n_S):
                singles[st.ref] = st

        if primary is None or not singles:
            continue
        if not gate_ok:
            n_gated += 1
            continue

        ordered = [singles.pop(primary.ref)] if primary.ref in singles else []
        ordered += list(singles.values())
        if not ordered:
            ordered = [primary]

        # Chimera hypotheses: primary paired with each supplementary segment that
        # explains a DISJOINT part of the read (specification §5).
        splits = []
        for sp in supps:
            if sp.ref == primary.ref:
                continue
            if _intervals_disjoint(primary, sp, max_overlap_frac):
                splits.append((primary, sp))
            else:
                hs.n_dropped_overlap += 1
        hs.add_read(qname, ordered, splits)

    hs.finalise()
    # canonical id per reference, for taxonomy lookups only
    hs.canonical = [TQ.reformat_ref_id(n, db_format, flat_sh) for n in hs.ref_names]
    hs.n_reads_seen = n_seen
    hs.n_reads_gated = n_gated
    if debug:
        print(f"[Model] {hs.n_reads:,} reads, {hs.n_refs:,} refs, {hs.n_hyp:,} hypotheses, "
              f"{int(hs.is_split.sum()):,} with a usable split "
              f"({hs.n_dropped_overlap:,} overlapping segments rejected)", flush=True)
    return hs


# ---------------------------------------------------------------------------
# The model (§7)
# ---------------------------------------------------------------------------

class GenerativeModel:
    """Joint EM over abundances theta, error model e, and chimera rate pi."""

    def __init__(self, hs, init='uniform', prior_alpha=0.0, estimate_error=True,
                 fixed_pi=None, clip_penalty=CLIP_PENALTY):
        self.hs = hs
        self.init = init
        self.prior_alpha = float(prior_alpha)
        self.estimate_error = estimate_error
        self.fixed_pi = fixed_pi
        self.clip_penalty = float(clip_penalty)
        J = hs.n_refs
        if init == 'primary' and hs.n_reads:
            counts = np.bincount(hs.primary_ref, minlength=J).astype(float)
            self.theta = counts / counts.sum() if counts.sum() > 0 else np.full(J, 1.0 / J)
        else:
            self.theta = np.full(J, 1.0 / J) if J else np.zeros(0)
        self.e = E_INIT.copy()
        n_split = int(hs.is_split.sum())
        self.pi = (0.5 if fixed_pi is None else float(fixed_pi)) if n_split else 0.0
        self.history = []

    # -- E-step (eq. 6-7) ---------------------------------------------------
    def _log_joint(self):
        hs = self.hs
        logL = hs.log_likelihood(np.log(np.maximum(self.e, 1e-300)), self.clip_penalty)
        lt = np.log(np.maximum(self.theta, 1e-300))
        out = np.empty(hs.n_hyp)
        single = ~hs.is_chim
        out[single] = lt[hs.ref1[single]] + logL[single]
        out[hs.is_chim] = (lt[hs.ref1[hs.is_chim]] + lt[hs.ref2[hs.is_chim]]
                           + logL[hs.is_chim])
        # mixture weights: (1-pi) for single hypotheses of a split read, pi for
        # chimeric ones.  A read without a split cannot be chimeric.
        pi = min(max(self.pi, 1e-12), 1 - 1e-12)
        read_of = np.repeat(np.arange(hs.n_reads), hs.seg_len)
        split_read = hs.is_split[read_of]
        out[single & split_read] += math.log1p(-pi)
        out[hs.is_chim] += math.log(pi)
        return out

    def _responsibilities(self):
        hs = self.hs
        lj = self._log_joint()
        mx = np.maximum.reduceat(lj, hs.row_ptr[:-1])
        shifted = np.exp(lj - np.repeat(mx, hs.seg_len))
        denom = np.add.reduceat(shifted, hs.row_ptr[:-1])
        ll_per_read = np.log(denom) + mx
        gamma = shifted / np.repeat(denom, hs.seg_len)
        return gamma, float(ll_per_read.sum())

    # -- M-step (eq. 9-11) --------------------------------------------------
    def _m_step_theta(self, gamma):
        hs = self.hs
        J = hs.n_refs
        n = np.bincount(hs.ref1, weights=gamma, minlength=J)
        if hs.is_chim.any():
            # a chimera contributes to BOTH parents (eq. 9)
            n += np.bincount(hs.ref2[hs.is_chim], weights=gamma[hs.is_chim], minlength=J)
        if self.prior_alpha > 0:
            a = self.prior_alpha
            self.theta = (n + a) / (n.sum() + J * a)
        else:
            self.theta = n / n.sum() if n.sum() > 0 else np.full(J, 1.0 / J)
        return n

    def _m_step_pi(self, gamma):
        hs = self.hs
        if self.fixed_pi is not None or not hs.is_split.any():
            return
        rho = np.zeros(hs.n_reads)
        if hs.is_chim.any():
            read_of = np.repeat(np.arange(hs.n_reads), hs.seg_len)
            np.add.at(rho, read_of[hs.is_chim], gamma[hs.is_chim])
        self.pi = float(np.clip(rho[hs.is_split].mean(), 1e-6, 1 - 1e-6))

    def _m_step_error(self, gamma):
        """Weighted MLE of the categorical over ALIGNED columns — eq. (11).

        The clip column is excluded from both numerator and denominator: it is a
        fixed penalty, not a fitted parameter."""
        if not self.estimate_error:
            return
        num = self.hs.counts[:, :N_OPS].T @ gamma
        tot = num.sum()
        if tot > 0:
            self.e = np.maximum(num / tot, 1e-12)
            self.e /= self.e.sum()

    # -- driver -------------------------------------------------------------
    def fit(self, max_iter=200, tol=1e-6, verbose=False):
        """Run EM until the per-read log-likelihood change falls below tol."""
        if self.hs.n_reads == 0:
            self.info = {'n_reads': 0, 'iterations': 0, 'converged': True,
                         'loglik': 0.0, 'monotone': True}
            return self
        prev = -np.inf
        monotone = True
        it = 0
        for it in range(1, max_iter + 1):
            gamma, ll = self._responsibilities()
            if np.isfinite(prev) and ll < prev - 1e-6:
                monotone = False          # EM guarantees this cannot happen (§7.4)
            self._m_step_theta(gamma)
            self._m_step_pi(gamma)
            self._m_step_error(gamma)
            delta = (ll - prev) / self.hs.n_reads if np.isfinite(prev) else np.inf
            self.history.append(ll)
            prev = ll
            if verbose:
                print(f"    iter {it:>3}  logL={ll:,.2f}  d/read={delta:.3e}  "
                      f"pi={self.pi:.4f}  e={np.round(self.e,4)}", flush=True)
            if abs(delta) < tol:
                break
        gamma, ll = self._responsibilities()
        self.gamma = gamma
        self.info = {
            'n_reads': int(self.hs.n_reads),
            'n_refs': int(self.hs.n_refs),
            'n_hypotheses': int(self.hs.n_hyp),
            'iterations': int(it),
            'converged': bool(abs(delta) < tol) if np.isfinite(delta) else False,
            'loglik': float(ll),
            'monotone': bool(monotone),
            'pi': float(self.pi),
            'error_model': {o: float(v) for o, v in zip(OPS, self.e)},
            'init': self.init,
            'reads_with_split': int(self.hs.is_split.sum()),
            'overlapping_segments_rejected': int(self.hs.n_dropped_overlap),
        }
        return self

    def expected_counts(self):
        """n_j — expected number of source molecules per reference (eq. 9)."""
        hs = self.hs
        n = np.bincount(hs.ref1, weights=self.gamma, minlength=hs.n_refs)
        if hs.is_chim.any():
            n += np.bincount(hs.ref2[hs.is_chim], weights=self.gamma[hs.is_chim],
                             minlength=hs.n_refs)
        return n

    def clean_counts(self):
        """Abundance excluding chimeric contributions (option noted in §7.2)."""
        hs = self.hs
        single = ~hs.is_chim
        return np.bincount(hs.ref1[single], weights=self.gamma[single],
                           minlength=hs.n_refs)

    def chimera_posterior(self):
        """rho_r — posterior probability that each read is chimeric."""
        hs = self.hs
        rho = np.zeros(hs.n_reads)
        if hs.is_chim.any():
            read_of = np.repeat(np.arange(hs.n_reads), hs.seg_len)
            np.add.at(rho, read_of[hs.is_chim], self.gamma[hs.is_chim])
        return rho


# ---------------------------------------------------------------------------
# Posterior taxonomic resolution and adaptive units (§8)
# ---------------------------------------------------------------------------

def posterior_resolution(model, lineages, tau=0.9):
    """Per-read consensus depth from the POSTERIOR gamma (not a softmax).

    Returns (depths, supports, primary_refs, taxa) aligned with the reads."""
    hs = model.hs
    single = ~hs.is_chim
    depths = np.zeros(hs.n_reads, dtype=int)
    supports = np.zeros(hs.n_reads)
    taxa = [''] * hs.n_reads
    confounders = [''] * hs.n_reads

    canon = getattr(hs, 'canonical', hs.ref_names)
    lin_cache = {}
    def lin_of(j):
        v = lin_cache.get(j)
        if v is None:
            v = lineages.get(canon[j], ('',) * N_RANKS)
            lin_cache[j] = v
        return v

    for r in range(hs.n_reads):
        s, e = hs.row_ptr[r], hs.row_ptr[r + 1]
        sel = np.arange(s, e)[single[s:e]]
        if sel.size == 0:
            continue
        w = model.gamma[sel]
        tot = w.sum()
        if tot <= 0:
            continue
        w = w / tot                       # normalise over the non-chimeric part
        ids = [canon[j] for j in hs.ref1[sel]]
        lins = {canon[j]: lin_of(j) for j in hs.ref1[sel]}
        taxon, depth, support, conf, _ = TQ.weighted_consensus(ids, w, lins, tau=tau)
        depths[r], supports[r], taxa[r], confounders[r] = depth, support, taxon, conf
    return depths, supports, taxa, confounders


def adaptive_units(ci_low, ref_names, lineages, ranks=None):
    """Deepest taxonomic node whose credible interval excludes zero (§8).

    Turns the identifiability limit into a declared property: when two
    references cannot be told apart by the marker their individual abundances
    are not identifiable but their SUM is, so the taxon is reported at the
    finest rank that is actually supported."""
    ranks = ranks or list(range(N_RANKS, 0, -1))
    node_low = defaultdict(float)
    for j, name in enumerate(ref_names):
        lin = lineages.get(name, ('',) * N_RANKS)
        for d in ranks:
            label = lin[d - 1]
            if label:
                node_low[(d, tuple(p for p in lin[:d] if p))] += ci_low[j]
    units = {}
    for (d, path), low in node_low.items():
        if low > 0:
            key = path
            if key not in units or d > units[key][0]:
                units[key] = (d, low)
    return units


# ---------------------------------------------------------------------------
# Uncertainty: bootstrap and calibration (§9)
# ---------------------------------------------------------------------------

def bootstrap(hs, B=200, alpha=0.05, seed=0, max_iter=100, tol=1e-5,
              prior_alpha=0.0, init='uniform', estimate_error=True, verbose=False):
    """Non-parametric bootstrap over reads -> credible intervals on theta.

    Resampling is done by reweighting reads (multinomial counts), which is
    equivalent to resampling with replacement but avoids copying the hypothesis
    table B times."""
    rng = np.random.default_rng(seed)
    J = hs.n_refs
    R = hs.n_reads
    if R == 0 or J == 0:
        return np.zeros(J), np.zeros(J), np.zeros((0, J))
    draws = np.zeros((B, J))
    base = GenerativeModel(hs, init=init, prior_alpha=prior_alpha,
                           estimate_error=estimate_error).fit(max_iter=max_iter, tol=tol)
    for b in range(B):
        w_read = rng.multinomial(R, np.full(R, 1.0 / R)).astype(float)
        m = GenerativeModel(hs, init=init, prior_alpha=prior_alpha,
                            estimate_error=estimate_error)
        m.theta = base.theta.copy(); m.e = base.e.copy(); m.pi = base.pi   # warm start
        w_hyp = np.repeat(w_read, hs.seg_len)
        prev = -np.inf
        for _ in range(max_iter):
            gamma, _ = m._responsibilities()
            g = gamma * w_hyp
            n = np.bincount(hs.ref1, weights=g, minlength=J)
            if hs.is_chim.any():
                n += np.bincount(hs.ref2[hs.is_chim], weights=g[hs.is_chim], minlength=J)
            new = (n + prior_alpha) / (n.sum() + J * prior_alpha) if prior_alpha > 0 \
                else (n / n.sum() if n.sum() > 0 else m.theta)
            shift = np.abs(new - m.theta).max()
            m.theta = new
            if shift < 1e-7:
                break
        draws[b] = m.theta
        if verbose and (b + 1) % 25 == 0:
            print(f"    bootstrap {b+1}/{B}", flush=True)
    lo = np.percentile(draws, 100 * alpha / 2, axis=0)
    hi = np.percentile(draws, 100 * (1 - alpha / 2), axis=0)
    return lo, hi, draws


def calibration_curve(draws, truth, levels=(0.5, 0.8, 0.9, 0.95, 0.99)):
    """Empirical coverage of the credible intervals against a known truth (§9).

    ``truth`` is a vector aligned with the reference index.  Returns
    {nominal: empirical_coverage}; a well-calibrated method sits on the
    diagonal."""
    out = {}
    truth = np.asarray(truth, dtype=float)
    for lv in levels:
        a = 1 - lv
        lo = np.percentile(draws, 100 * a / 2, axis=0)
        hi = np.percentile(draws, 100 * (1 - a / 2), axis=0)
        inside = (truth >= lo) & (truth <= hi)
        out[lv] = float(inside.mean())
    return out


# ---------------------------------------------------------------------------
# Per-sample orchestrator
# ---------------------------------------------------------------------------

def run_sample(bam_path, tax_map_path, db_format='unite', flat_sh=False, threads=1,
               min_percent_identity=0.0, min_ref_coverage=0.0, tau=0.9,
               max_iter=200, tol=1e-6, prior_alpha=0.0, init='uniform',
               estimate_error=True, chimera=True, bootstrap_B=0, alpha=0.05, seed=0,
               abundance='all', em_out=None, quality_out=None, ci_out=None,
               samtools='samtools', debug=False):
    """Fit the generative model to one sample and write LOREON-compatible output."""
    bam_path = Path(bam_path)
    lineages = TQ.load_taxonomy_lineages(tax_map_path, db_format, flat_sh)

    hs = build_hypotheses(bam_path, lineages, db_format=db_format, flat_sh=flat_sh,
                          threads=threads, min_percent_identity=min_percent_identity,
                          min_ref_coverage=min_ref_coverage, samtools=samtools,
                          debug=debug)
    # chimera=False pins pi to 0, reducing the model to a pure abundance EM with
    # the indel-aware error model.  Worth having as an ablation: on ITS data the
    # chimera component drives pi towards 1 (most reads split because no single
    # UNITE reference explains them end-to-end, which is a database gap rather
    # than a PCR chimera) and, since eq. (9) credits BOTH parents of a chimera,
    # that inflates false positives.
    model = GenerativeModel(hs, init=init, prior_alpha=prior_alpha,
                            estimate_error=estimate_error,
                            fixed_pi=None if chimera else 0.0).fit(
        max_iter=max_iter, tol=tol, verbose=debug)

    counts = model.clean_counts() if abundance == 'clean' else model.expected_counts()
    result = {'sample_bam': bam_path.name, 'model': model.info, 'paths': {},
              'reads_seen': getattr(hs, 'n_reads_seen', hs.n_reads),
              'reads_excluded_by_quality_filter': getattr(hs, 'n_reads_gated', 0)}

    if em_out is not None:
        by_ref = {hs.ref_names[j]: float(counts[j]) for j in range(hs.n_refs)}
        _, n_written, kept, dropped = TQ.write_em_tsv(by_ref, em_out)
        result['paths']['em'] = str(em_out)
        result['model'].update({'refs_written': n_written, 'mass_written': kept,
                                'mass_dropped': dropped})

    # Feature A from true posteriors
    depths, supports, taxa, confs = posterior_resolution(model, lineages, tau=tau)
    acc = TQ.OtuQualityAccumulator()
    canon = getattr(hs, 'canonical', hs.ref_names)
    for r in range(hs.n_reads):
        otu = canon[hs.primary_ref[r]]
        acc.add(otu, int(depths[r]), float(supports[r]), 0,
                confs[r], 1.0 - float(supports[r]))
    otu_quality = acc.result()
    if quality_out is not None:
        TQ.write_quality_tsv(otu_quality, quality_out)
        result['paths']['quality'] = str(quality_out)
    hist = Counter(int(d) for d in depths)
    result['taxonomic_quality'] = TQ.quality_stats_block(
        hist, float(supports.sum()), hs.n_reads, otu_quality)

    # Feature C from the posterior, not from counting supplementary records
    rho = model.chimera_posterior()
    n_split = int(hs.is_split.sum())
    result['noise'] = {
        'n_mapped_reads': int(hs.n_reads),
        'split_reads': n_split,
        'split_read_pct': round(n_split / hs.n_reads * 100, 3) if hs.n_reads else 0.0,
        'chimera_rate_posterior': round(float(model.pi), 5),
        'expected_chimeric_reads': round(float(rho.sum()), 1),
        'chimera_rate_refined_pct': round(float(rho.sum()) / hs.n_reads * 100, 3)
        if hs.n_reads else 0.0,
        'n_candidate_chimeras': int((rho > 0.5).sum()),
        'overlapping_segments_rejected': int(hs.n_dropped_overlap),
        'refined_available': True,
    }

    if bootstrap_B:
        if debug:
            print(f"[Model] bootstrap B={bootstrap_B} ...", flush=True)
        lo, hi, draws = bootstrap(hs, B=bootstrap_B, alpha=alpha, seed=seed,
                                  prior_alpha=prior_alpha, init=init,
                                  estimate_error=estimate_error, verbose=debug)
        N = hs.n_reads
        if ci_out is not None:
            with open(ci_out, 'w', encoding='utf-8') as fh:
                fh.write("otu_id\ttheta\tci_low\tci_high\tcount\tcount_low\tcount_high\n")
                for j, name in enumerate(hs.ref_names):
                    if counts[j] <= 0.01 and hi[j] <= 0:
                        continue
                    fh.write(f"{name}\t{model.theta[j]:.8f}\t{lo[j]:.8f}\t{hi[j]:.8f}\t"
                             f"{counts[j]:.4f}\t{lo[j]*N:.4f}\t{hi[j]*N:.4f}\n")
            result['paths']['ci'] = str(ci_out)
        units = adaptive_units(lo, getattr(hs, 'canonical', hs.ref_names), lineages)
        result['uncertainty'] = {
            'bootstrap_B': int(bootstrap_B),
            'alpha': alpha,
            'refs_with_ci_excluding_zero': int((lo > 0).sum()),
            'adaptive_units': len(units),
            'shannon': _shannon_ci(draws, alpha),
        }
    return result


def _shannon_ci(draws, alpha=0.05):
    """Shannon index with a credible interval, from the bootstrap draws (§9)."""
    with np.errstate(divide='ignore', invalid='ignore'):
        H = np.array([-np.nansum(np.where(t > 0, t * np.log(t), 0.0)) for t in draws])
    return {'mean': float(H.mean()),
            'ci_low': float(np.percentile(H, 100 * alpha / 2)),
            'ci_high': float(np.percentile(H, 100 * (1 - alpha / 2)))}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Unified generative model (LOREON TaxQuality specification).")
    p.add_argument('--bam', required=True)
    p.add_argument('--tax-map', required=True)
    p.add_argument('--format', default='unite',
                   choices=['unite', 'silva', 'eukariome', 'cbs', 'none'])
    p.add_argument('--flat-sh', action='store_true')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--min-percent-identity', type=float, default=0.0)
    p.add_argument('--min-ref-coverage', type=float, default=0.0)
    p.add_argument('--consensus-tau', type=float, default=0.9)
    p.add_argument('--max-iter', type=int, default=200)
    p.add_argument('--tol', type=float, default=1e-6)
    p.add_argument('--prior-alpha', type=float, default=0.0)
    p.add_argument('--init', choices=['uniform', 'primary'], default='uniform')
    p.add_argument('--no-estimate-error', action='store_true',
                   help="keep the error model fixed at its initial value")
    p.add_argument('--no-chimera', action='store_true',
                   help="disable the chimera component (pi fixed at 0)")
    p.add_argument('--abundance', choices=['all', 'clean'], default='all',
                   help="'clean' excludes chimeric contributions from theta")
    p.add_argument('--bootstrap', type=int, default=0, metavar='B')
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--em-out', default=None)
    p.add_argument('--quality-out', default=None)
    p.add_argument('--ci-out', default=None)
    p.add_argument('--stats-out', default=None)
    p.add_argument('--debug', action='store_true')
    return p.parse_args(argv)


def main(argv=None):
    a = _parse_args(argv)
    res = run_sample(
        bam_path=a.bam, tax_map_path=a.tax_map, db_format=a.format, flat_sh=a.flat_sh,
        threads=a.threads, min_percent_identity=a.min_percent_identity,
        min_ref_coverage=a.min_ref_coverage, tau=a.consensus_tau,
        max_iter=a.max_iter, tol=a.tol, prior_alpha=a.prior_alpha, init=a.init,
        estimate_error=not a.no_estimate_error, chimera=not a.no_chimera,
        bootstrap_B=a.bootstrap,
        alpha=a.alpha, seed=a.seed, abundance=a.abundance,
        em_out=a.em_out, quality_out=a.quality_out, ci_out=a.ci_out, debug=a.debug)
    if a.stats_out:
        Path(a.stats_out).write_text(json.dumps(res, indent=2), encoding='utf-8')
    else:
        json.dump(res, sys.stdout, indent=2); print()
    return 0


if __name__ == '__main__':
    sys.exit(main())
