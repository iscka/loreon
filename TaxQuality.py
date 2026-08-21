#!/usr/bin/env python3
"""TaxQuality.py — Taxonomic quality, noise/chimera rate and fractional (EM)
assignment for the LOREON long-read metabarcoding pipeline.

All three features consume the *same* datum: a single read-level pass over the
full per-sample BAM (`mapping/<sample>_sorted_full.bam`).  The pass derives:

  * Feature A — a weighted *adaptive taxonomic consensus* per read (LCA-style),
    aggregated into a per-OTU resolution profile (how deep each OTU is
    credible, and with whom it is confused).
  * Feature C — a split/chimera rate refined by comparing the taxon of the
    primary alignment with the taxa of its supplementary segments.
  * Feature B — fractional abundance via Expectation-Maximisation (kallisto /
    RSEM / Bracken style) that redistributes multi-mapping reads by maximum
    likelihood.

Design constraints honoured here:
  * No pysam (PyInstaller-hostile).  The BAM is read by streaming the textual
    output of `samtools view`, name-collated so records of one read arrive
    together and nothing needs to be held for the whole sample.
  * Everything is additive and opt-in.
  * Only numpy is used for the numerical core.

Read-set consistency
--------------------
The main OTU table counts reads surviving `samtools view -F 0x904` plus the
identity/coverage expression built in ``OtuUtils.tabeling_improved``.  This
module applies the *same* gates to primary alignments so that the EM table and
the quality columns describe the same read set as the primary table.  The gate
formula deliberately mirrors the htslib expression, including its treatment of
deletions, rather than using a more "correct" identity — consistency with the
table it sits beside matters more than which identity definition is nicer.
"""

import argparse
import gzip
import json
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Taxonomy ranks (7 core ranks, matching ResultsReader.TAXONOMY_RANKS[:7])
# ---------------------------------------------------------------------------
RANK_NAMES = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
N_RANKS = len(RANK_NAMES)                 # 7
SPECIES_DEPTH = N_RANKS                    # 7
GENUS_DEPTH = N_RANKS - 1                  # 6

# Softmax temperature, in raw minimap2 alignment-score (AS) units.
#
# CALIBRATION NOTE.  map-ont scoring is match=+2 / mismatch=-4, so one extra
# mismatch moves AS by ~6 points.  T=6 therefore makes a one-mismatch-worse
# candidate e^-1 as likely as the best one.  A likelihood-calibrated value
# (Karlin-Altschul lambda for this scoring scheme) would be nearer T~2, i.e.
# T=6 is deliberately CONSERVATIVE: it reports more residual ambiguity than
# the alignment scores strictly support.  Exposed as --score-temperature.
#
# AS is intentionally NOT normalised by read length.  Within a read the softmax
# is invariant to a constant divisor, and *across* reads the length-dependence
# is correct behaviour, not a bug: a 1% divergence over 900 bp really is ~6x
# more mismatches of evidence than over 150 bp, and the estimate should be
# correspondingly more decisive.  Long reads cannot dominate the abundance
# estimate regardless, because every read contributes exactly 1.0 of
# responsibility mass to the EM.
SCORE_TEMPERATURE = 6.0

# Candidates whose normalised within-read weight falls below this are dropped
# as numerically negligible; the primary is always retained.
MIN_CANDIDATE_WEIGHT = 1e-4

# A read counts as "confounded" at its consensus rank when the runner-up group
# carries at least this much weight mass.
CONFOUNDER_MIN_MASS = 0.01

_CIGAR_OPS = set('MIDNSHP=X')


# ---------------------------------------------------------------------------
# Reference-id → canonical OTU_ID reformatting
# ---------------------------------------------------------------------------

def reformat_ref_id(raw: str, db_format: str, flat_sh: bool = False) -> str:
    """Map a raw reference name (BAM RNAME) to the canonical OTU_ID used as the
    taxonomy-map key and as the OTU table index.

    Mirrors ``OtuUtils._reformat_single_file`` / ``ResultsReader._build_reformat_sql``
    exactly — including the UNITE Flat-SH collapse — so ids produced here join
    1:1 against the taxonomy map and the OTU table.
    """
    if not raw or raw == '*':
        return raw
    if db_format == 'unite':
        # Replicate DuckDB's split_part() semantics EXACTLY, including for
        # malformed headers with fewer than 5 fields: split_part returns '' when
        # the index is out of range, so a header with no '|' becomes '|'.
        # Returning the raw id instead (the defensive choice) would key the
        # quality row differently from the OTU table row and silently drop that
        # OTU's quality data in the join.  The table's index is the source of
        # truth, so consistency wins over defensiveness here.
        parts = raw.split('|')
        get = lambda i: parts[i] if i < len(parts) else ''
        # flat_sh collapses accessions sharing a Species Hypothesis
        return get(2) if flat_sh else f"{get(1)}|{get(2)}"
    if db_format == 'eukariome':
        return raw.split(';', 1)[0]
    if db_format in ('cbs', 'none'):
        return raw.split('|', 1)[0]
    # silva and anything else: unchanged
    return raw


def load_taxonomy_lineages(tax_map_path, db_format: str = 'unite',
                           flat_sh: bool = False) -> dict:
    """Load ``<db>_taxonomy_map.tsv`` into ``{canonical_OTU_ID: lineage}``.

    The map is keyed by the non-flat canonical id (``acc|SH`` for UNITE).  When
    ``flat_sh`` is active the OTU table is indexed by the SH alone, so the keys
    are re-projected to match.
    """
    lineages = {}
    with open(Path(tax_map_path), 'r', encoding='utf-8') as fh:
        fh.readline()  # header: OTU_ID\tTaxonomy
        for line in fh:
            line = line.rstrip('\n\r')
            if not line:
                continue
            otu_id, _, taxonomy = line.partition('\t')
            if flat_sh and db_format == 'unite':
                # cached key is "acc|SH" -> project to "SH"
                parts = otu_id.split('|')
                if len(parts) >= 2:
                    otu_id = parts[1]
            lineages[otu_id] = _split_lineage(taxonomy)
    return lineages


def _split_lineage(taxonomy: str):
    """Split a semicolon-separated taxonomy string into a 7-tuple of ranks.

    POSITIONAL: an empty field keeps its slot.  Dropping empty fields (the
    previous behaviour) left-shifted every rank below a gap, so e.g.
    ``Fungi;Ascomycota;;;Tuberaceae;Tuber;Tuber_melanosporum`` filed the family
    as a Class and blanked Genus/Species.  Only trailing ranks are padded.
    """
    if not taxonomy:
        return ('',) * N_RANKS
    parts = [p.strip().replace('_', ' ') for p in taxonomy.split(';')]
    if len(parts) >= N_RANKS:
        return tuple(parts[:N_RANKS])
    return tuple(parts) + ('',) * (N_RANKS - len(parts))


# ---------------------------------------------------------------------------
# CIGAR / tag helpers
# ---------------------------------------------------------------------------

def _parse_cigar(cigar: str):
    """Yield (length, op) pairs. Hand-rolled: faster than re for hot loops."""
    if not cigar or cigar == '*':
        return
    num = 0
    for ch in cigar:
        if ch.isdigit():
            num = num * 10 + (ord(ch) - 48)
        elif ch in _CIGAR_OPS:
            yield num, ch
            num = 0
        else:
            return


def query_length_from_cigar(cigar: str) -> int:
    """Full read length implied by a CIGAR (query-consuming ops incl. clips).

    M/I/S/=/X consume the query; H is part of the original read too, so it is
    counted — this recovers the true read length from any single record."""
    return sum(n for n, op in _parse_cigar(cigar) if op in 'MIS=XH')


def soft_clip_length(cigar: str) -> int:
    """Soft-clipped query bases (htslib `sclen`)."""
    return sum(n for n, op in _parse_cigar(cigar) if op == 'S')


def aligned_query_length(cigar: str) -> int:
    """Query bases inside the alignment (M/I/=/X — excludes clips)."""
    return sum(n for n, op in _parse_cigar(cigar) if op in 'MI=X')


def parse_tags(fields):
    """Extract the optional tags we need (TAG:TYPE:VALUE)."""
    tags = {}
    for f in fields:
        if f[:5] in ('AS:i:', 'NM:i:', 'de:f:', 'SA:Z:'):
            tags[f[:2]] = f[5:]
    return tags


def samtools_identity(tags: dict, cigar: str):
    """Identity replicating the htslib filter used by ``tabeling_improved``:

        ((qlen - sclen) - NM) / (qlen - sclen)

    Returns a fraction in [0,1], or None when NM or the CIGAR is unavailable.
    This intentionally matches the main pipeline's definition (deletions inflate
    NM without enlarging the denominator) so that the same reads pass both
    gates.  ``gap_compressed_identity`` is the biologically nicer measure and is
    used for reporting, never for gating.
    """
    if not cigar or cigar == '*' or 'NM' not in tags:
        return None
    aln = query_length_from_cigar(cigar) - soft_clip_length(cigar)
    if aln <= 0:
        return None
    try:
        nm = int(tags['NM'])
    except ValueError:
        return None
    return max(0.0, (aln - nm) / aln)


def query_coverage(cigar: str) -> float:
    """Fraction of the read aligned, matching htslib ``(qlen-sclen)/qlen``."""
    qlen = query_length_from_cigar(cigar)
    if qlen <= 0:
        return 0.0
    return (qlen - soft_clip_length(cigar)) / qlen


def gap_compressed_identity(tags: dict, cigar: str):
    """minimap2 `de` tag (gap-compressed divergence) → identity, else None.

    Used for reporting only, never for gating (see samtools_identity)."""
    if 'de' in tags:
        try:
            return max(0.0, 1.0 - float(tags['de']))
        except ValueError:
            return None
    return None


def supplementary_targets_from_sa(sa_value: str):
    """Yield reference names from an SA:Z tag (rname,pos,strand,CIGAR,mapQ,NM;...)."""
    for seg in sa_value.split(';'):
        if not seg:
            continue
        rname = seg.split(',', 1)[0]
        if rname:
            yield rname


# ---------------------------------------------------------------------------
# Streaming, name-collated BAM traversal
# ---------------------------------------------------------------------------

def iter_read_groups(bam_path, threads=1, samtools='samtools'):
    """Yield ``(qname, [record, ...])`` with all records of a read together.

    The BAM is coordinate-sorted, so records of one read are scattered.  Rather
    than holding the whole sample in RAM, the stream is passed through
    ``samtools collate`` (a fast name-grouping shuffle, no full sort) so that
    consecutive records share a qname and the caller can process one read at a
    time.

    Each record is ``(flag, rname, mapq, cigar, tags)``.  Unmapped records are
    skipped.  stderr is redirected to a temp file and only read on failure —
    reading it after the stdout loop (the previous behaviour) deadlocks as soon
    as samtools writes more than one pipe buffer of warnings.
    """
    bam_path = str(bam_path)
    with tempfile.TemporaryFile() as errfile:
        collate = subprocess.Popen(
            [samtools, 'collate', f'-@{threads}', '-u', '-O', bam_path],
            stdout=subprocess.PIPE, stderr=errfile)
        view = subprocess.Popen(
            [samtools, 'view', f'-@{threads}', '-'],
            stdin=collate.stdout, stdout=subprocess.PIPE, stderr=errfile,
            text=True, encoding='utf-8', errors='replace', bufsize=1024 * 1024)
        collate.stdout.close()   # let collate see SIGPIPE if view dies

        cur_q = None
        batch = []
        try:
            for line in view.stdout:
                if not line or line[0] == '@':
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 11:
                    continue
                flag = int(fields[1])
                if flag & 0x4:
                    continue
                rname = fields[2]
                if rname == '*':
                    continue
                rec = (flag, rname, int(fields[4]), fields[5], parse_tags(fields[11:]))
                qname = fields[0]
                if qname != cur_q:
                    if cur_q is not None and batch:
                        yield cur_q, batch
                    cur_q, batch = qname, [rec]
                else:
                    batch.append(rec)
            if cur_q is not None and batch:
                yield cur_q, batch
        finally:
            if view.stdout:
                view.stdout.close()
            view.wait()
            collate.wait()
            if view.returncode not in (0, None) or collate.returncode not in (0, None):
                errfile.seek(0)
                err = errfile.read().decode('utf-8', errors='replace')
                raise RuntimeError(
                    f"samtools collate|view failed on {bam_path} "
                    f"(collate={collate.returncode}, view={view.returncode}): {err.strip()[:2000]}")


# ---------------------------------------------------------------------------
# Per-read candidate assembly
# ---------------------------------------------------------------------------

class ReadCandidates:
    """Compact per-read summary produced while streaming."""
    __slots__ = ('qname', 'primary_raw', 'primary_ref', 'primary_mapq',
                 'refs_raw', 'refs_canon', 'scores', 'supp_refs', 'passes_gate')

    def __init__(self, qname):
        self.qname = qname
        self.primary_raw = None
        self.primary_ref = None
        self.primary_mapq = 0
        self.refs_raw = []
        self.refs_canon = []
        self.scores = []
        self.supp_refs = set()
        self.passes_gate = False


def build_read(qname, records, db_format, flat_sh, min_identity, min_coverage):
    """Turn one read's records into a ReadCandidates, applying the quality gates.

    ``passes_gate`` mirrors the main pipeline: the PRIMARY alignment must clear
    both the identity and the coverage thresholds, exactly as
    ``samtools view -F 0x904 -e ...`` does for the primary count table.
    Secondary candidates are additionally gated so that noise does not enter the
    consensus or the EM.  Candidates on the same reference are merged, keeping
    the best score, so a repeat-induced double alignment cannot inflate that
    reference's weight.
    """
    rc = ReadCandidates(qname)
    best_by_ref = {}

    for flag, rname, mapq, cigar, tags in records:
        canonical = reformat_ref_id(rname, db_format, flat_sh)

        if flag & 0x800:                       # supplementary: split segment
            rc.supp_refs.add(canonical)
            if 'SA' in tags:
                for sa in supplementary_targets_from_sa(tags['SA']):
                    rc.supp_refs.add(reformat_ref_id(sa, db_format, flat_sh))
            continue

        try:
            as_score = float(tags['AS']) if 'AS' in tags else None
        except ValueError:
            as_score = None
        ident = samtools_identity(tags, cigar)
        cov = query_coverage(cigar)

        if not (flag & 0x100):                 # primary
            rc.primary_raw = rname
            rc.primary_ref = canonical
            rc.primary_mapq = mapq
            ok_id = (min_identity <= 0.0) or (ident is None) or (ident >= min_identity)
            ok_cov = (min_coverage <= 0.0) or (cov >= min_coverage)
            rc.passes_gate = ok_id and ok_cov
            if 'SA' in tags:
                for sa in supplementary_targets_from_sa(tags['SA']):
                    rc.supp_refs.add(reformat_ref_id(sa, db_format, flat_sh))
        else:                                  # secondary candidate
            if min_identity > 0.0 and ident is not None and ident < min_identity:
                continue
            if min_coverage > 0.0 and cov < min_coverage:
                continue

        if as_score is None:
            # No AS: cannot place this candidate on the score scale.  Keep the
            # primary (it defines the read's OTU) but never let a scoreless
            # secondary enter the softmax with a fabricated score of 0.
            if flag & 0x100:
                continue
            as_score = 0.0

        prev = best_by_ref.get(rname)
        if prev is None or as_score > prev[0]:
            best_by_ref[rname] = (as_score, canonical)

    for raw, (score, canon) in best_by_ref.items():
        rc.refs_raw.append(raw)
        rc.refs_canon.append(canon)
        rc.scores.append(score)
    return rc


# ---------------------------------------------------------------------------
# Weights (softmax over AS) — shared by consensus (Feature A) and EM (Feature B)
# ---------------------------------------------------------------------------

def candidate_weights(scores, temperature=SCORE_TEMPERATURE):
    """Normalised softmax weights over alignment scores."""
    s = np.asarray(scores, dtype=float)
    if s.size == 0:
        return s
    s = s - s.max()
    w = np.exp(s / temperature)
    total = w.sum()
    if total <= 0:
        return np.full(s.shape, 1.0 / s.size)
    return w / total


def prune_candidates(refs_raw, refs_canon, scores, temperature=SCORE_TEMPERATURE,
                     min_weight=MIN_CANDIDATE_WEIGHT):
    """Drop numerically negligible candidates, always keeping the best one.

    This is the tail filter that MIN_CANDIDATE_WEIGHT documents; previously the
    constant was declared but never applied."""
    if len(scores) <= 1:
        return refs_raw, refs_canon, scores
    w = candidate_weights(scores, temperature)
    keep = w >= min_weight
    keep[int(np.argmax(w))] = True
    if keep.all():
        return refs_raw, refs_canon, scores
    return ([r for r, k in zip(refs_raw, keep) if k],
            [c for c, k in zip(refs_canon, keep) if k],
            [s for s, k in zip(scores, keep) if k])


# ---------------------------------------------------------------------------
# Feature A — weighted adaptive taxonomic consensus
# ---------------------------------------------------------------------------

def weighted_consensus(cand_ids, weights, lineages, tau=0.9):
    """Adaptive weighted-LCA consensus for one read.

    Returns ``(taxon, depth, support, confounder, confounded_mass)``:
      * taxon            — consensus taxon label at the chosen rank
      * depth            — rank depth (Species=7 … Kingdom=1; 0 if unresolved)
      * support          — share of *annotated* mass backing the consensus
      * confounder       — the competing taxon AT THE LEVEL WHERE THE CONSENSUS
                           BROKE (one rank deeper), i.e. who stole the finer call
      * confounded_mass  — that competitor's mass share

    Two corrections over a naive implementation:

    1. ``support`` is measured against the mass that is *annotated at that rank*,
       not against 1.0.  Otherwise a candidate whose database entry simply lacks
       a species name drags the read's consensus up a rank, reporting reference
       incompleteness as ambiguity of the sample.
    2. Grouping uses the lineage path with empty ranks treated as wildcards, so
       two references that agree wherever both are annotated are not split apart
       by a missing intermediate rank.
    """
    n = len(cand_ids)
    if n == 0:
        return ('Unclassified', 0, 0.0, '', 0.0)

    lins = [lineages.get(cid, ('',) * N_RANKS) for cid in cand_ids]
    w = np.asarray(weights, dtype=float)

    deeper_runner = ('', 0.0)     # runner-up seen one level below the accepted one

    for depth in range(N_RANKS, 0, -1):
        idx = depth - 1
        members = [j for j in range(n) if lins[j][idx] != '']
        if not members:
            continue
        annotated = float(w[members].sum())
        if annotated <= 0:
            continue

        # Key on the ancestor ranks that are annotated for EVERY candidate at
        # this level.  A rank missing from even one database entry carries no
        # discriminating information, so including it would split two records
        # of the same taxon apart purely because one lineage has a gap.
        # (Filtering empties out of each key individually does NOT work: it
        # left-shifts the tuple and makes the keys unequal again.)
        usable = [i for i in range(depth)
                  if all(lins[j][i] != '' for j in members)]
        groups = defaultdict(float)
        labels = {}
        for j in members:
            key = tuple(lins[j][i] for i in usable)
            groups[key] += w[j]
            labels[key] = lins[j][idx]

        ranked = sorted(groups.items(), key=lambda kv: kv[1], reverse=True)
        top_key, top_mass = ranked[0]
        support = top_mass / annotated
        if support >= tau:
            if deeper_runner[0]:
                conf, conf_mass = deeper_runner
            elif len(ranked) > 1:
                conf, conf_mass = labels[ranked[1][0]], ranked[1][1] / annotated
            else:
                conf, conf_mass = '', 0.0
            if conf == labels[top_key]:
                conf, conf_mass = '', 0.0
            return (labels[top_key], depth, float(support), conf, float(conf_mass))

        # consensus failed here — remember who the competition was, then go up
        if len(ranked) > 1:
            deeper_runner = (labels[ranked[1][0]], ranked[1][1] / annotated)
        else:
            deeper_runner = (labels[top_key], 0.0)

    return ('Unclassified', 0, 0.0, deeper_runner[0], deeper_runner[1])


def _depth_to_rank(depth):
    d = int(depth)
    if 1 <= d <= N_RANKS:
        return RANK_NAMES[d - 1]
    return 'Unclassified'


class OtuQualityAccumulator:
    """Streaming accumulator for the per-OTU resolution profile."""

    def __init__(self):
        self.depth_counts = defaultdict(Counter)     # otu -> {depth: n}
        self.support_sum = defaultdict(float)
        self.mapq_sum = defaultdict(float)
        self.n = defaultdict(int)
        self.confounders = defaultdict(Counter)      # otu -> {taxon: confounded reads}
        self.conf_mass = defaultdict(lambda: defaultdict(float))

    def add(self, otu, depth, support, mapq, confounder, confounded_mass):
        self.depth_counts[otu][depth] += 1
        self.support_sum[otu] += support
        self.mapq_sum[otu] += mapq
        self.n[otu] += 1
        if confounder and confounded_mass >= CONFOUNDER_MIN_MASS:
            self.confounders[otu][confounder] += 1
            self.conf_mass[otu][confounder] += confounded_mass

    def result(self):
        out = {}
        for otu, n in self.n.items():
            dc = self.depth_counts[otu]
            pct = lambda pred: round(
                sum(c for d, c in dc.items() if pred(d)) / n * 100.0, 2)
            modal_depth = max(dc.items(), key=lambda kv: (kv[1], kv[0]))[0]
            conf = self.confounders[otu]
            top_conf, top_conf_reads = ('', 0)
            if conf:
                top_conf, top_conf_reads = conf.most_common(1)[0]
            out[otu] = {
                'n_reads': n,
                # These four partition the reads exactly; they sum to 100.
                'pct_species': pct(lambda d: d >= SPECIES_DEPTH),
                'pct_genus': pct(lambda d: d == GENUS_DEPTH),
                'pct_above_genus': pct(lambda d: 1 <= d < GENUS_DEPTH),
                'pct_unresolved': pct(lambda d: d <= 0),
                # Modal rank, not a median: averaging an ordinal rank scale
                # invents ranks that no read supports (a 50/50 species/kingdom
                # OTU used to be reported as "Order").
                'resolution_rank': _depth_to_rank(modal_depth),
                'resolution_depth_mode': int(modal_depth),
                'mean_support': round(self.support_sum[otu] / n, 4),
                'mean_mapq': round(self.mapq_sum[otu] / n, 2),
                'top_confounder': top_conf,
                'confounded_reads': int(top_conf_reads),
                'confounded_mass': round(self.conf_mass[otu].get(top_conf, 0.0), 3),
            }
        return out


QUALITY_COLUMNS = ['n_reads', 'pct_species', 'pct_genus', 'pct_above_genus',
                   'pct_unresolved', 'resolution_rank', 'resolution_depth_mode',
                   'mean_support', 'mean_mapq', 'top_confounder',
                   'confounded_reads', 'confounded_mass']


def write_quality_tsv(otu_quality, out_path):
    """Write the per-OTU quality table joined later by makeOtu_duckdb."""
    out_path = Path(out_path)
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write('\t'.join(['otu_id'] + QUALITY_COLUMNS) + '\n')
        for otu, m in sorted(otu_quality.items()):
            fh.write(otu + '\t' + '\t'.join(str(m[c]) for c in QUALITY_COLUMNS) + '\n')
    return out_path


def quality_stats_block(depth_hist, support_sum, n_reads, otu_quality):
    """Build the ``taxonomic_quality`` block for mapping_stats.json."""
    edges = Counter()
    for otu, m in otu_quality.items():
        if m['top_confounder']:
            # weight by the reads actually confounded, not by the OTU's total
            edges[(otu, m['top_confounder'])] += m['confounded_reads']
    top_confusion = [{'otu': a, 'confounder': b, 'confounded_reads': int(c)}
                     for (a, b), c in edges.most_common(15) if c > 0]
    pct_species = (round(depth_hist.get(SPECIES_DEPTH, 0) / n_reads * 100.0, 2)
                   if n_reads else 0.0)
    return {
        'n_reads': int(n_reads),
        'pct_reads_species': pct_species,
        'mean_support': round(support_sum / n_reads, 4) if n_reads else 0.0,
        'depth_histogram': {_depth_to_rank(d): int(depth_hist.get(d, 0))
                            for d in range(N_RANKS, 0, -1)},
        'unresolved_reads': int(depth_hist.get(0, 0)),
        'top_confusion': top_confusion,
    }


# ---------------------------------------------------------------------------
# Feature C — split / cross-taxon chimera statistics
# ---------------------------------------------------------------------------

class ChimeraAccumulator:
    def __init__(self, lineages, refined=True, keep_qnames=False, max_qnames=1000):
        self.lineages = lineages
        self.refined = refined
        self.n_mapped = 0
        self.split_reads = 0
        self.inter_taxon = 0
        self.no_genus = 0
        self.qnames = [] if keep_qnames else None
        self.max_qnames = max_qnames

    def add(self, rc):
        if rc.primary_ref is None:
            return
        self.n_mapped += 1
        if not rc.supp_refs:
            return
        self.split_reads += 1
        if not self.refined:
            return
        pg = self.lineages.get(rc.primary_ref, ('',) * N_RANKS)[GENUS_DEPTH - 1]
        saw_comparable = False
        for supp in rc.supp_refs:
            if supp == rc.primary_ref:
                continue
            sg = self.lineages.get(supp, ('',) * N_RANKS)[GENUS_DEPTH - 1]
            if pg and sg:
                saw_comparable = True
                if pg != sg:
                    self.inter_taxon += 1
                    if self.qnames is not None and len(self.qnames) < self.max_qnames:
                        self.qnames.append(rc.qname)
                    return
        if not saw_comparable:
            self.no_genus += 1

    def result(self):
        n = self.n_mapped
        r = {
            'n_mapped_reads': n,
            'split_reads': self.split_reads,
            'split_read_pct': round(self.split_reads / n * 100, 3) if n else 0.0,
            'n_candidate_chimeras': self.inter_taxon,
            'chimera_rate_refined_pct': round(self.inter_taxon / n * 100, 3) if n else 0.0,
            # split reads whose genus could not be compared (missing annotation)
            'split_reads_genus_unknown': self.no_genus,
            'refined_available': self.refined,
        }
        return r


# ---------------------------------------------------------------------------
# Feature B — Expectation-Maximisation fractional assignment
# ---------------------------------------------------------------------------

class EMAccumulator:
    """Collects the compact sparse structure the EM needs, while streaming.

    Stores int32 reference indices and float32 weights rather than per-read
    Python objects: ~8 bytes per candidate instead of several hundred."""

    def __init__(self, temperature=SCORE_TEMPERATURE):
        self.temperature = temperature
        self.ref_index = {}
        self.row_ref = []
        self.row_w = []
        self.row_ptr = [0]
        self.primary_counts = []
        self.n_reads = 0

    def add(self, rc):
        if not rc.refs_raw:
            return
        w = candidate_weights(rc.scores, self.temperature)
        for raw, wj in zip(rc.refs_raw, w):
            j = self.ref_index.get(raw)
            if j is None:
                j = len(self.ref_index)
                self.ref_index[raw] = j
                self.primary_counts.append(0.0)
            self.row_ref.append(j)
            self.row_w.append(wj)
        self.row_ptr.append(len(self.row_ref))
        if rc.primary_raw is not None:
            j = self.ref_index.get(rc.primary_raw)
            if j is not None:
                self.primary_counts[j] += 1.0
        self.n_reads += 1

    def arrays(self):
        return (np.asarray(self.row_ref, dtype=np.int64),
                np.asarray(self.row_w, dtype=np.float64),
                np.asarray(self.row_ptr, dtype=np.int64),
                np.asarray(self.primary_counts, dtype=np.float64))


def run_em(acc, max_iter=200, tol=1e-6, prior=0.0, init='uniform'):
    """EM over the accumulated candidates → per-reference expected counts.

    Model: theta_j = latent relative abundance of reference j (sum = 1).
      E-step: gamma_rj = theta_j L_rj / sum_k theta_k L_rk
      M-step: theta_j = (sum_r gamma_rj + prior) / (N + J*prior)

    ``init`` defaults to **uniform** (as kallisto/RSEM/Bracken do).  Warm-starting
    from primary counts — the previous behaviour — creates an absorbing zero
    state: any reference that never wins a primary alignment has theta_j = 0, and
    since gamma is proportional to theta it can never receive mass however strong
    its alignment evidence.  That is precisely the multi-mapping case the EM
    exists to resolve, so the warm start defeated the feature.  ``init='primary'``
    is kept for comparison but is not recommended.

    Convergence is tested on the *per-read* mean log-likelihood change, so the
    criterion does not tighten as the sample grows.

    Returns ``(counts_by_raw_ref, info)`` keyed by raw BAM reference name, so
    the result can be written in .tmp.txt schema and reuse makeOtu_duckdb.
    """
    row_ref, row_w, row_ptr, primary_counts = acc.arrays()
    J = len(acc.ref_index)
    n_reads = acc.n_reads
    if J == 0 or n_reads == 0:
        return {}, {'n_reads': 0, 'n_refs': 0, 'iterations': 0, 'converged': True,
                    'final_loglik': 0.0, 'delta_per_read': 0.0,
                    'counts_total': 0.0, 'reads_unassignable': 0}

    if init == 'primary' and primary_counts.sum() > 0:
        theta = primary_counts / primary_counts.sum()
    else:
        theta = np.full(J, 1.0 / J)

    seg_len = np.diff(row_ptr)
    prev_ll = -np.inf
    iterations = 0
    converged = False
    delta_per_read = 0.0
    counts = np.zeros(J, dtype=float)
    unassignable = 0

    for iterations in range(1, max_iter + 1):
        contrib = theta[row_ref] * row_w
        denom = np.add.reduceat(contrib, row_ptr[:-1])
        # A read whose candidates all have theta == 0 has zero likelihood.  With
        # a uniform init this cannot happen, but a user-supplied init or an
        # extreme prior could produce it; count such reads instead of silently
        # scoring them as log(1) = 0, which used to hide the failure inside the
        # very number offered as the sanity check.
        bad = denom <= 0
        unassignable = int(bad.sum())
        denom_safe = np.where(bad, 1.0, denom)
        gamma = contrib / np.repeat(denom_safe, seg_len)
        counts = np.bincount(row_ref, weights=gamma, minlength=J)
        ll = float(np.sum(np.log(denom_safe[~bad]))) if not bad.all() else float('-inf')

        assigned = float(counts.sum())
        if prior > 0.0:
            theta = (counts + prior) / (assigned + J * prior)
        elif assigned > 0:
            theta = counts / assigned
        delta_per_read = (ll - prev_ll) / n_reads if np.isfinite(prev_ll) else np.inf
        prev_ll = ll
        if abs(delta_per_read) < tol:
            converged = True
            break

    # Report the log-likelihood OF THE RETURNED parameters: recompute once more
    # after the final M-step so final_loglik and counts describe the same theta.
    contrib = theta[row_ref] * row_w
    denom = np.add.reduceat(contrib, row_ptr[:-1])
    bad = denom <= 0
    final_ll = float(np.sum(np.log(np.where(bad, 1.0, denom)[~bad]))) if not bad.all() else float('-inf')

    counts_by_ref = {raw: float(counts[idx]) for raw, idx in acc.ref_index.items()}
    info = {
        'n_reads': int(n_reads),
        'n_refs': int(J),
        'iterations': int(iterations),
        'converged': bool(converged),
        'final_loglik': final_ll,
        'delta_per_read': float(delta_per_read) if np.isfinite(delta_per_read) else None,
        'counts_total': round(float(counts.sum()), 4),
        'reads_unassignable': unassignable,
        'init': init,
    }
    return counts_by_ref, info


def write_em_tsv(counts_by_ref, out_path, min_count=0.01):
    """Write EM counts in the .tmp.txt schema (``id  length  mapped  unmapped``).

    Counts are written as FRACTIONAL values.  An earlier version rounded to int
    and dropped anything below 0.5, which (a) contradicted the whole point of a
    fractional-assignment table and (b) silently destroyed read mass — a long
    tail of sub-1.0 assignments could lose a large share of the sample with no
    warning.  ``makeOtu_duckdb`` reads this column as DOUBLE for the EM path.

    ``min_count`` defaults to 0.01 reads: the EM leaves a long tail of
    numerically-negligible responsibilities (on a real ITS sample, 327 rows were
    written for 73 references carrying any real mass, the rest being values that
    round to zero in the output).  Dropping below a hundredth of a read removed
    254 such rows while retaining 100% of the assigned mass.  The dropped mass is
    reported so the loss is never silent.

    Returns (path, n_written, mass_written, mass_dropped).
    """
    out_path = Path(out_path)
    written = 0
    kept = 0.0
    dropped = 0.0
    with open(out_path, 'w', encoding='utf-8') as fh:
        for raw, count in sorted(counts_by_ref.items(), key=lambda kv: -kv[1]):
            if count <= min_count:
                dropped += count
                continue
            fh.write(f"{raw}\t0\t{count:.6f}\t0\n")
            kept += count
            written += 1
        fh.write("*\t0\t0\t0\n")
    return out_path, written, round(kept, 4), round(dropped, 4)


# ---------------------------------------------------------------------------
# Per-sample orchestrator
# ---------------------------------------------------------------------------

def run_sample(bam_path, tax_map_path, db_format='unite', flat_sh=False, threads=1,
               min_percent_identity=95.0, min_ref_coverage=90.0, tau=0.9,
               temperature=SCORE_TEMPERATURE, do_quality=True, do_em=True,
               chimera_refined=True, em_max_iter=200, em_tol=1e-6, em_prior=0.0,
               em_init='uniform', quality_out=None, em_out=None, readcand_out=None,
               samtools='samtools', debug=False):
    """Single streaming pass over one sample's full BAM producing all outputs."""
    bam_path = Path(bam_path)
    if not bam_path.exists():
        raise FileNotFoundError(f"Full BAM not found for tax-quality pass: {bam_path}")

    min_identity = _as_fraction(min_percent_identity)
    min_coverage = _as_fraction(min_ref_coverage)
    lineages = load_taxonomy_lineages(tax_map_path, db_format, flat_sh)

    if debug:
        print(f"[TaxQuality:{bam_path.name}] streaming read-level pass "
              f"(identity>={min_identity:.3f}, coverage>={min_coverage:.3f}, "
              f"tau={tau}, T={temperature})...", flush=True)

    qual = OtuQualityAccumulator() if do_quality else None
    chim = ChimeraAccumulator(lineages, refined=chimera_refined)
    em = EMAccumulator(temperature) if do_em else None
    depth_hist = Counter()
    support_sum = 0.0
    n_consensus = 0
    n_seen = 0
    n_gated_out = 0

    cand_fh = None
    if readcand_out is not None:
        opener = gzip.open if str(readcand_out).endswith('.gz') else open
        cand_fh = opener(readcand_out, 'wt', encoding='utf-8')
        cand_fh.write('qname\tis_primary\tref_raw\tref_canonical\tAS\tweight\tlineage\n')

    try:
        for qname, records in iter_read_groups(bam_path, threads=threads,
                                               samtools=samtools):
            rc = build_read(qname, records, db_format, flat_sh,
                            min_identity, min_coverage)
            n_seen += 1
            # Chimera statistics are a property of the RAW library, so they are
            # computed over every mapped read regardless of the quality gate.
            chim.add(rc)

            if rc.primary_ref is None or not rc.refs_raw:
                continue
            if not rc.passes_gate:
                n_gated_out += 1
                continue

            rc.refs_raw, rc.refs_canon, rc.scores = prune_candidates(
                rc.refs_raw, rc.refs_canon, rc.scores, temperature)
            w = candidate_weights(rc.scores, temperature)

            if qual is not None:
                taxon, depth, support, conf, conf_mass = weighted_consensus(
                    rc.refs_canon, w, lineages, tau=tau)
                qual.add(rc.primary_ref, depth, support, rc.primary_mapq,
                         conf, conf_mass)
                depth_hist[depth] += 1
                support_sum += support
                n_consensus += 1

            if em is not None:
                em.add(rc)

            if cand_fh is not None:
                for i, (raw, cid) in enumerate(zip(rc.refs_raw, rc.refs_canon)):
                    lin = ';'.join(lineages.get(cid, ('',) * N_RANKS))
                    cand_fh.write(f"{qname}\t{1 if raw == rc.primary_raw else 0}\t"
                                  f"{raw}\t{cid}\t{rc.scores[i]:.1f}\t{w[i]:.4f}\t{lin}\n")
    finally:
        if cand_fh is not None:
            cand_fh.close()

    result = {
        'sample_bam': bam_path.name,
        'paths': {},
        'reads_seen': n_seen,
        'reads_excluded_by_quality_filter': n_gated_out,
        'filters_applied': {
            'min_percent_identity': min_percent_identity,
            'min_ref_coverage': min_ref_coverage,
        },
    }

    if do_quality:
        otu_quality = qual.result()
        if quality_out is not None:
            write_quality_tsv(otu_quality, quality_out)
            result['paths']['quality'] = str(quality_out)
        result['taxonomic_quality'] = quality_stats_block(
            depth_hist, support_sum, n_consensus, otu_quality)

    result['noise'] = chim.result()

    if do_em:
        counts_by_ref, info = run_em(em, max_iter=em_max_iter, tol=em_tol,
                                     prior=em_prior, init=em_init)
        if em_out is not None:
            _, n_written, kept, dropped = write_em_tsv(counts_by_ref, em_out)
            info['refs_written'] = n_written
            info['mass_written'] = kept
            info['mass_dropped'] = dropped
            result['paths']['em'] = str(em_out)
        result['em'] = info
        if debug:
            print(f"[TaxQuality:{bam_path.name}] EM: {info['iterations']} iters, "
                  f"converged={info['converged']}, logL={info['final_loglik']:.2f}, "
                  f"mass={info['counts_total']}/{info['n_reads']}", flush=True)

    return result


def _as_fraction(pct):
    """Convert a percentage threshold to a fraction.

    Values are percentages by contract (95 == 95%).  The previous heuristic
    ('>1 means percent') inverted the meaning of small values: --min-percent-identity 1
    became a 100% gate.  Only a value of exactly 0 disables the gate."""
    if pct is None:
        return 0.0
    if pct <= 0:
        return 0.0
    return pct / 100.0


# ---------------------------------------------------------------------------
# Stand-alone CLI (one sample)
# ---------------------------------------------------------------------------

def _parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Taxonomic-quality / chimera / EM pass over one sample BAM.")
    p.add_argument('--bam', required=True, help='Path to <sample>_sorted_full.bam')
    p.add_argument('--tax-map', required=True, help='Taxonomy map TSV')
    p.add_argument('--format', default='unite',
                   choices=['unite', 'silva', 'eukariome', 'cbs', 'none'])
    p.add_argument('--flat-sh', action='store_true')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--min-percent-identity', type=float, default=95.0)
    p.add_argument('--min-ref-coverage', type=float, default=90.0)
    p.add_argument('--consensus-tau', type=float, default=0.9)
    p.add_argument('--score-temperature', type=float, default=SCORE_TEMPERATURE)
    p.add_argument('--no-quality', action='store_true')
    p.add_argument('--no-em', action='store_true')
    p.add_argument('--chimera-refined', action='store_true', default=True)
    p.add_argument('--em-max-iter', type=int, default=200)
    p.add_argument('--em-tol', type=float, default=1e-6)
    p.add_argument('--em-prior', type=float, default=0.0)
    p.add_argument('--em-init', choices=['uniform', 'primary'], default='uniform')
    p.add_argument('--quality-out', default=None)
    p.add_argument('--em-out', default=None)
    p.add_argument('--readcand-out', default=None)
    p.add_argument('--stats-out', default=None)
    p.add_argument('--debug', action='store_true')
    return p.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)
    result = run_sample(
        bam_path=args.bam, tax_map_path=args.tax_map, db_format=args.format,
        flat_sh=args.flat_sh, threads=args.threads,
        min_percent_identity=args.min_percent_identity,
        min_ref_coverage=args.min_ref_coverage,
        tau=args.consensus_tau, temperature=args.score_temperature,
        do_quality=not args.no_quality, do_em=not args.no_em,
        chimera_refined=args.chimera_refined, em_max_iter=args.em_max_iter,
        em_tol=args.em_tol, em_prior=args.em_prior, em_init=args.em_init,
        quality_out=args.quality_out, em_out=args.em_out,
        readcand_out=args.readcand_out, debug=args.debug,
    )
    if args.stats_out:
        with open(args.stats_out, 'w', encoding='utf-8') as fh:
            json.dump(result, fh, indent=2)
    else:
        json.dump(result, sys.stdout, indent=2)
        print()
    return 0


if __name__ == '__main__':
    sys.exit(main())
