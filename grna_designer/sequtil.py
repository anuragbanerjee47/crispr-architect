"""Small sequence utilities shared across the package."""
from __future__ import annotations

_COMP = str.maketrans("ACGTN", "TGCAN")


def rev_comp(seq: str) -> str:
    """Reverse-complement of a DNA string. Ns map to Ns."""
    return seq.upper().translate(_COMP)[::-1]
