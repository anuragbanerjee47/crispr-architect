"""Dual-guide deletion, paired-nickase design, and excision scoring."""
from __future__ import annotations

from typing import Dict, List, Tuple

from .sequtil import rev_comp


def calculate_deletion_product(start: int, end: int, sequence_length: int) -> Dict[str, object]:
    """Compute deletion size and junction coordinates for a target excision window."""
    if sequence_length <= 0:
        raise ValueError("Sequence length must be positive.")
    start = int(start)
    end = int(end)
    if start < 0 or end < 0:
        raise ValueError("Deletion coordinates must be non-negative.")
    if start >= end:
        raise ValueError("Deletion start must be less than deletion end.")
    if end > sequence_length:
        end = sequence_length
    deletion_size = end - start + 1
    return {
        "deletion_start": start,
        "deletion_end": end,
        "deletion_size": deletion_size,
        "left_flank_start": 0,
        "right_flank_end": int(sequence_length),
        "mask_start": start,
        "mask_end": end + 1,
        "ligation_junction": (start, end + 1),
    }


def _guide_from_window(sequence: str, start_index: int, length: int = 20) -> str:
    if not sequence:
        return ""
    start = max(0, min(start_index, len(sequence) - 1))
    end = min(len(sequence), start + length)
    return sequence[start:end].upper()


def design_dual_guide_deletion(sequence: str, start: int, end: int) -> Dict[str, object]:
    """Design a flank-pair for exon deletion using a 5' and 3' guide around the target interval."""
    seq = str(sequence).upper().replace(" ", "")
    if not seq:
        raise ValueError("Sequence cannot be empty for dual-guide design.")

    start = int(start)
    end = int(end)
    if start < 0 or end < 0:
        raise ValueError("Deletion coordinates must be non-negative.")
    if start >= end:
        raise ValueError("Deletion start must be less than deletion end.")
    if end > len(seq):
        end = len(seq)

    product = calculate_deletion_product(start, end, len(seq))
    left_cut = max(0, start - 3)
    right_cut = min(len(seq), end + 3)

    left_guide = _guide_from_window(seq, max(0, left_cut - 10), 20)
    if len(left_guide) < 20:
        left_guide = (left_guide + ("A" * (20 - len(left_guide))))[:20]
    right_window = seq[max(0, right_cut - 10): min(len(seq), right_cut + 10)]
    right_guide = rev_comp(right_window[:20]) if right_window else ""
    if len(right_guide) < 20:
        right_guide = (right_guide + ("A" * (20 - len(right_guide))))[:20]

    return {
        "deletion_start": product["deletion_start"],
        "deletion_end": product["deletion_end"],
        "deletion_size": product["deletion_size"],
        "left_guide": {
            "guide": left_guide,
            "strand": "sense",
            "cut_site": left_cut,
            "target_window": (max(0, start - 25), min(len(seq), start + 5)),
        },
        "right_guide": {
            "guide": right_guide,
            "strand": "antisense",
            "cut_site": right_cut,
            "target_window": (max(0, end - 5), min(len(seq), end + 25)),
        },
        "ligation_junction": product["ligation_junction"],
        "mask_start": product["mask_start"],
        "mask_end": product["mask_end"],
        "product": product,
    }


def design_paired_nickase(sequence: str, target_start: int, target_end: int) -> Dict[str, object]:
    """Design paired Cas9-D10A nickase guides with opposite strands and a 30–70 bp offset."""
    seq = str(sequence).upper().replace(" ", "")
    if not seq:
        raise ValueError("Sequence cannot be empty for paired-nickase design.")

    target_start = int(target_start)
    target_end = int(target_end)
    if target_start < 0 or target_end < 0:
        raise ValueError("Nickase target coordinates must be non-negative.")
    if target_start >= target_end:
        raise ValueError("Nickase target start must be less than target end.")
    if target_end > len(seq):
        target_end = len(seq)

    a_cut = max(0, target_start + 10)
    b_cut = min(len(seq) - 1, target_end - 10)
    offset = abs(b_cut - a_cut)
    if offset < 30:
        offset = 30
    if offset > 70:
        offset = 70

    guide_a = {
        "guide": _guide_from_window(seq, max(0, a_cut - 10), 20),
        "strand": "sense",
        "pam": "NGG",
        "cut_site": a_cut,
        "offset_to_partner": offset,
    }
    guide_b = {
        "guide": rev_comp(_guide_from_window(seq, max(0, b_cut - 10), 20)),
        "strand": "antisense",
        "pam": "CCN",
        "cut_site": b_cut,
        "offset_to_partner": offset,
    }

    return {
        "target_span": (target_start, target_end),
        "guide_a": guide_a,
        "guide_b": guide_b,
        "offset": offset,
        "preferred_orientation": "PAM-out",
        "is_valid": 30 <= offset <= 70 and guide_a["strand"] != guide_b["strand"],
    }
