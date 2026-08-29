"""Scoring utilities for guide design.

The implementation uses a Rule Set 2-inspired on-target prior and a Hsu-Zhang
seed-weighted off-target penalty model, while keeping the public API stable for
existing CLI and UI callers.
"""
from __future__ import annotations

import math
import re
from typing import List, Optional, Sequence

from .models import ScoringResult
from .sequtil import rev_comp

GUIDE_LEN = 20

# Position-specific nucleotide weights for a PAM-proximal 20-mer spacer.
# Values are calibrated so PAM-proximal G is rewarded, proximal T is penalized,
# and distal positions contribute less strongly to on-target efficiency.
POSITION_WEIGHTS = [
    {"A": -0.15, "C": 0.10, "G": 0.95, "T": -0.70},
    {"A": -0.05, "C": 0.25, "G": 0.80, "T": -0.60},
    {"A": -0.05, "C": 0.25, "G": 0.75, "T": -0.55},
    {"A": 0.00, "C": 0.20, "G": 0.70, "T": -0.50},
    {"A": 0.05, "C": 0.15, "G": 0.65, "T": -0.45},
    {"A": 0.05, "C": 0.10, "G": 0.45, "T": -0.35},
    {"A": 0.05, "C": 0.10, "G": 0.35, "T": -0.30},
    {"A": 0.00, "C": 0.10, "G": 0.30, "T": -0.30},
    {"A": 0.00, "C": 0.10, "G": 0.25, "T": -0.25},
    {"A": 0.00, "C": 0.10, "G": 0.25, "T": -0.25},
    {"A": 0.00, "C": 0.08, "G": 0.15, "T": -0.20},
    {"A": 0.00, "C": 0.08, "G": 0.12, "T": -0.18},
    {"A": 0.00, "C": 0.06, "G": 0.10, "T": -0.15},
    {"A": 0.00, "C": 0.05, "G": 0.08, "T": -0.12},
    {"A": 0.00, "C": 0.04, "G": 0.08, "T": -0.10},
    {"A": 0.00, "C": 0.04, "G": 0.06, "T": -0.10},
    {"A": 0.00, "C": 0.03, "G": 0.05, "T": -0.08},
    {"A": 0.00, "C": 0.02, "G": 0.04, "T": -0.08},
    {"A": 0.00, "C": 0.02, "G": 0.03, "T": -0.05},
    {"A": 0.00, "C": 0.01, "G": 0.02, "T": -0.05},
]

# Hsu-Zhang seed-heavy mismatch weights. Positions 1-10 are nearest the PAM.
HSU_ZHANG_MISMATCH_PENALTIES: List[float] = [
    3.8, 3.4, 3.1, 2.9, 2.5, 2.2, 2.0, 1.7, 1.5, 1.3,
    0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.1,
]


def _gc_optimality_score(guide: str) -> float:
    gc_count = guide.count("G") + guide.count("C")
    gc_fraction = gc_count / len(guide)
    if 0.45 <= gc_fraction <= 0.60:
        return 1.0
    if 0.40 <= gc_fraction < 0.45 or 0.60 < gc_fraction <= 0.70:
        return 0.6
    if gc_fraction < 0.40:
        return max(0.0, 1.0 - (0.40 - gc_fraction) * 8.0)
    return max(0.0, 1.0 - (gc_fraction - 0.70) * 8.0)


def _has_pol3_terminator(guide: str) -> bool:
    return bool(re.search(r"T{2,}", guide))


def _has_homopolymer_run(guide: str) -> bool:
    if re.search(r"T{4,}", guide):
        return False
    return bool(re.search(r"(.)\1{3,}", guide))


def _is_self_complementary(guide: str, min_len: int = 6) -> bool:
    rc = rev_comp(guide)
    for k in range(min_len, min(len(guide), 12) + 1):
        for i in range(len(guide) - k + 1):
            window = guide[i:i + k]
            if window in rc:
                return True
    return False


def guide_quality_flags(guide_5to3: str) -> List[str]:
    """Return QC flags for structural quality problems in a guide."""
    g = guide_5to3.upper()
    flags: List[str] = []
    if _has_pol3_terminator(g):
        flags.append("poly_t")
    if _has_homopolymer_run(g):
        flags.append("homopolymer")
    if _is_self_complementary(g):
        flags.append("hairpin")
    return flags


def hsu_zhang_penalty(mismatch_positions: Sequence[int]) -> float:
    """Weighted mismatch penalty using the seed-heavy Hsu-Zhang convention."""
    if not mismatch_positions:
        return 0.0
    total = 0.0
    for pos in mismatch_positions:
        idx = max(0, min(int(pos) - 1, len(HSU_ZHANG_MISMATCH_PENALTIES) - 1))
        total += HSU_ZHANG_MISMATCH_PENALTIES[idx]
    return total


def hsu_zhang_specificity_score(mismatch_positions: Sequence[int]) -> float:
    """Return a 0-100 specificity score where 100 = zero off-target risk.

    The score captures the Hsu-Zhang rule that PAM-proximal seed mismatches are far
    more disruptive to specificity than distal mismatches. The implementation uses
    the requested form:
        Score = (product of w_p) * (1 / (((19 - d) / 19) * 4 + 1)) * (1 / (m^2))
    scaled into a 0-100 bounded score.
    """
    positions = sorted({int(p)
                       for p in mismatch_positions if 1 <= int(p) <= 20})
    if not positions:
        return 100.0

    weights = []
    for pos in positions:
        idx = int(pos) - 1
        weight = HSU_ZHANG_MISMATCH_PENALTIES[idx] if 0 <= idx < len(
            HSU_ZHANG_MISMATCH_PENALTIES) else 0.1
        weights.append(weight)

    m = len(positions)
    product_term = math.prod(weights)
    if m > 1:
        gap_values = [b - a for a, b in zip(positions, positions[1:])]
        avg_gap = sum(gap_values) / len(gap_values)
    else:
        avg_gap = 19.0

    d = min(19.0, max(0.0, avg_gap))
    distance_factor = 1.0 / ((((19.0 - d) / 19.0) * 4.0) + 1.0)
    raw = product_term * distance_factor * (1.0 / (m * m))
    # Convert the raw Hsu-Zhang penalty into a bounded specificity score.
    return max(0.0, min(100.0, 100.0 / (1.0 + raw)))


def _logistic_transform(z: float) -> float:
    return 1.0 / (1.0 + math.exp(-z))


def score_guide(guide_5to3: str, mismatch_positions: Optional[Sequence[int]] = None) -> ScoringResult:
    """Return a calibrated on-target score plus a specificity score for a guide."""
    g = guide_5to3.upper()
    if len(g) not in {GUIDE_LEN, GUIDE_LEN - 1}:
        raise ValueError(f"Guide must be {GUIDE_LEN} nt or {GUIDE_LEN - 1} nt, got {len(g)}")

    # Position-specific Rule Set 2 style contribution.
    raw_score = 0.0
    for idx, base in enumerate(g):
        weights = POSITION_WEIGHTS[idx]
        if base in weights:
            raw_score += float(weights[base])

    # Strongly penalize P-III terminators and self-complementarity.
    if _has_pol3_terminator(g):
        raw_score -= 7.5
    if _has_homopolymer_run(g):
        raw_score -= 2.0
    if _is_self_complementary(g):
        raw_score -= 2.5

    gc_count = g.count("G") + g.count("C")
    gc_fraction = gc_count / len(g)
    if gc_fraction < 0.40:
        raw_score -= (0.40 - gc_fraction) * 20.0
    elif gc_fraction > 0.70:
        raw_score -= (gc_fraction - 0.70) * 40.0
    elif 0.45 <= gc_fraction <= 0.60:
        raw_score += 2.0
    else:
        raw_score += 0.5

    quality_flags = guide_quality_flags(g)
    if "poly_t" in quality_flags:
        raw_score -= 3.0
    if "homopolymer" in quality_flags:
        raw_score -= 1.5
    if "hairpin" in quality_flags:
        raw_score -= 1.5

    contributions: dict = {
        "gc_count": gc_count,
        "gc_fraction": gc_fraction,
        "flags": quality_flags,
        "gc_optimality": _gc_optimality_score(g),
        "linear_score": raw_score,
    }

    mismatch_penalty = 0.0
    specificity_score = 100.0
    if mismatch_positions is not None:
        mismatch_penalty = hsu_zhang_penalty(mismatch_positions)
        specificity_score = hsu_zhang_specificity_score(mismatch_positions)
        raw_score -= mismatch_penalty * 0.35
        contributions["hsu_zhang_penalty"] = mismatch_penalty
        contributions["hsu_zhang_specificity_score"] = specificity_score

    on_target_efficiency = _logistic_transform(raw_score)
    on_target_score = max(0.0, min(100.0, on_target_efficiency * 100.0))

    result = ScoringResult(
        linear=raw_score,
        efficiency=on_target_efficiency,
        on_target_score=on_target_score,
        specificity_score=specificity_score,
        contributing_features=contributions,
    )
    result.contributing_features["raw_score"] = raw_score
    return result
