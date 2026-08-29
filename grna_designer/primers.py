"""Flanking validation primer design and QC for CRISPR guide candidates."""
from __future__ import annotations

from typing import Dict, Tuple

from .cloning import gc_percent, nearest_neighbor_tm
from .models import Guide
from .sequtil import rev_comp

DEFAULT_AMPLICON_MIN = 400
DEFAULT_AMPLICON_MAX = 600
TARGET_TM_MIN = 58.0
TARGET_TM_MAX = 62.0
TARGET_LENGTH_MIN = 18
TARGET_LENGTH_MAX = 24
TARGET_GC_MIN = 40.0
TARGET_GC_MAX = 60.0


def _primer_quality(seq: str) -> Dict[str, object]:
    seq = seq.upper()
    if not seq:
        return {"length": 0, "gc": 0.0, "tm": 0.0, "passes": False}
    gc = gc_percent(seq)
    tm = nearest_neighbor_tm(seq)
    length = len(seq)
    passes = (
        TARGET_LENGTH_MIN <= length <= TARGET_LENGTH_MAX
        and TARGET_GC_MIN <= gc <= TARGET_GC_MAX
        and TARGET_TM_MIN <= tm <= TARGET_TM_MAX
    )
    return {"length": length, "gc": gc, "tm": tm, "passes": passes}


def _select_candidate_primer(seq: str, center: int, direction: str) -> str:
    seq = seq.upper()
    n = len(seq)
    if not seq:
        return ""
    center = max(0, min(center, n - 1))
    best: Tuple[float, float, str] | None = None

    search_start = max(0, center - 35)
    search_end = min(n, center + 35)
    for start in range(search_start, search_end):
        for length in range(TARGET_LENGTH_MAX, TARGET_LENGTH_MIN - 1, -1):
            end = start + length
            if end > n:
                continue
            candidate = seq[start:end]
            if direction == "reverse":
                candidate = rev_comp(candidate)
            tm = nearest_neighbor_tm(candidate)
            gc = gc_percent(candidate)
            score = abs(tm - 60.0) + abs(gc - 50.0) * 0.12
            if TARGET_TM_MIN <= tm <= TARGET_TM_MAX and TARGET_GC_MIN <= gc <= TARGET_GC_MAX:
                return candidate
            if best is None or score < best[0]:
                best = (score, abs(tm - 60.0), candidate)

    if best is not None:
        return best[2]
    return seq[max(0, center - TARGET_LENGTH_MIN + 1):center + 1][:TARGET_LENGTH_MIN]


def primer_pair_quality(fwd_primer: str, rev_primer: str) -> Dict[str, bool]:
    """Flag risky primer-primer interactions based on 3' complementarity and dimer checks."""
    fw = str(fwd_primer).upper()
    rv = str(rev_primer).upper()
    rc_rv = rev_comp(rv)
    rc_fw = rev_comp(fw)

    fwd_3 = fw[-6:]
    rev_3 = rv[-6:]
    fwd_self = fwd_3 in rc_fw[:-1] or fw[-8:] in rc_fw
    rev_self = rev_3 in rc_rv[:-1] or rv[-8:] in rc_rv

    dimer = False
    for primer, other in ((fw, rv), (rv, fw)):
        for offset in range(1, min(len(primer), len(other)) - 3):
            end = primer[-offset - 3:-offset] if offset else primer[-3:]
            if len(end) < 3:
                continue
            if end in rev_comp(other[-len(end):]):
                dimer = True
                break
        if dimer:
            break

    return {
        "primer_dimer_risk": bool(dimer),
        "self_complementarity_risk": bool(fwd_self or rev_self),
    }


def _amplicon_bounds(seq: str, cut_site: int, target_size: int = 500) -> Tuple[int, int]:
    seq_length = len(seq)
    start = max(0, cut_site - (target_size // 2))
    end = min(seq_length, cut_site + (target_size // 2))

    if end - start < DEFAULT_AMPLICON_MIN:
        padding = (DEFAULT_AMPLICON_MIN - (end - start)) // 2
        start = max(0, start - padding)
        end = min(seq_length, end + padding)
    if end - start > DEFAULT_AMPLICON_MAX:
        overshoot = end - start - DEFAULT_AMPLICON_MAX
        start += overshoot // 2
        end -= overshoot - overshoot // 2
    return start, end


def design_validation_primers(seq: str, guide: Guide) -> Dict[str, object]:
    """Design flanking genomic PCR validation primers around the cut site."""
    seq = str(seq).upper()
    if not seq:
        raise ValueError(
            "Sequence cannot be empty for validation primer design.")

    cut_site = int(getattr(guide, "pos", 0))
    if hasattr(guide, "guide") and getattr(guide, "guide"):
        cut_site = max(
            0, min(len(seq) - 1, int(guide.pos) + len(guide.guide) // 2))

    amplicon_start, amplicon_end = _amplicon_bounds(
        seq, cut_site, target_size=500)
    amplicon_length = max(0, amplicon_end - amplicon_start)
    if amplicon_length < DEFAULT_AMPLICON_MIN:
        amplicon_start = max(
            0, min(len(seq) - DEFAULT_AMPLICON_MIN, amplicon_start))
        amplicon_end = min(len(seq), amplicon_start + DEFAULT_AMPLICON_MIN)
        amplicon_length = amplicon_end - amplicon_start
    if amplicon_length > DEFAULT_AMPLICON_MAX:
        amplicon_end = min(len(seq), amplicon_start + DEFAULT_AMPLICON_MAX)
        amplicon_length = amplicon_end - amplicon_start

    fwd_primer = _select_candidate_primer(seq, cut_site, "forward")
    rev_primer = _select_candidate_primer(seq, cut_site, "reverse")
    pair_quality = primer_pair_quality(fwd_primer, rev_primer)

    if len(fwd_primer) < TARGET_LENGTH_MIN:
        fwd_primer = (fwd_primer + seq[:TARGET_LENGTH_MIN])[:TARGET_LENGTH_MIN]
    if len(rev_primer) < TARGET_LENGTH_MIN:
        rev_primer = rev_comp(
            (rev_primer + seq[:TARGET_LENGTH_MIN])[:TARGET_LENGTH_MIN])
    if len(fwd_primer) > TARGET_LENGTH_MAX:
        fwd_primer = fwd_primer[:TARGET_LENGTH_MAX]
    if len(rev_primer) > TARGET_LENGTH_MAX:
        rev_primer = rev_primer[:TARGET_LENGTH_MAX]

    fwd_quality = _primer_quality(fwd_primer)
    rev_quality = _primer_quality(rev_primer)

    if not fwd_quality["passes"]:
        fwd_primer = _select_candidate_primer(seq, cut_site, "forward")
    if not rev_quality["passes"]:
        rev_primer = _select_candidate_primer(seq, cut_site, "reverse")

    fwd_quality = _primer_quality(fwd_primer)
    rev_quality = _primer_quality(rev_primer)
    return {
        "fwd_primer": fwd_primer,
        "rev_primer": rev_primer,
        "tm_fwd": float(fwd_quality["tm"]),
        "tm_rev": float(rev_quality["tm"]),
        "gc_fwd": float(fwd_quality["gc"]),
        "gc_rev": float(rev_quality["gc"]),
        "amplicon_start": amplicon_start,
        "amplicon_end": amplicon_end,
        "amplicon_length": amplicon_length,
        "primer_dimer_risk": bool(pair_quality["primer_dimer_risk"]),
        "self_complementarity_risk": bool(pair_quality["self_complementarity_risk"]),
    }


def attach_validation_primers(guide: Guide, sequence: str) -> Guide:
    """Attach validation primer design to a guide candidate object."""
    design = design_validation_primers(sequence, guide)
    guide.flanking_fwd_primer = design["fwd_primer"]
    guide.flanking_rev_primer = design["rev_primer"]
    guide.tm_fwd = design["tm_fwd"]
    guide.tm_rev = design["tm_rev"]
    guide.amplicon_length = design["amplicon_length"]
    guide.amplicon_start = design["amplicon_start"]
    guide.amplicon_end = design["amplicon_end"]
    return guide
