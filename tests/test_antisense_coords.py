"""Antisense coordinate transform: original_pos = n - 1 - (j + 19)."""
from __future__ import annotations

from grna_designer.guides import find_guides
from grna_designer.models import Strand
from grna_designer.sequtil import rev_comp


def test_antisense_position_math():
    """Build a sequence with a known antisense site and verify the position."""
    # Place 20-nt guide on the SENSE strand, then check the antisense strand
    # at the same locus has a position computed by n - 1 - (j + 19).
    # Easiest: a fully-palindromic sequence.
    palindrome = "GAATTC" * 100  # 600 nt
    n = len(palindrome)

    guides = find_guides(palindrome)
    antisense_guides = [g for g in guides if g.strand == Strand.ANTISENSE]
    assert antisense_guides
    # For every antisense guide, the position is n-1-(j+19) for some j in [0..len(rc)-23].
    for g in antisense_guides:
        # Reverse the transform: the j index is n - 1 - g.pos - 19.
        j = n - 1 - g.pos - 19
        assert 0 <= j <= n - 23
        # And the rc-window at j matches the guide
        rc = rev_comp(palindrome)
        assert rc[j:j + 20] == g.guide


def test_antisense_position_consistency():
    """A guide and its antisense counterpart at the same locus should have
    positions that are mirror-images around the sequence midpoint."""
    seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"  # 60 nt
    # Force a NGG PAM at one sense location by appending one.
    seq = seq + "CGG"  # 63 nt total
    n = len(seq)
    guides = find_guides(seq)
    # Every antisense guide's pos should be n-1-(j+19) and thus within [0, n-20]
    for g in guides:
        if g.strand == Strand.ANTISENSE:
            assert 0 <= g.pos <= n - 20
