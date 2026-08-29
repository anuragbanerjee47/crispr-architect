"""PAM scanning and protospacer detection."""
from __future__ import annotations

from grna_designer.guides import find_guides
from grna_designer.models import Strand


def test_planted_ngg_at_position_601_sense():
    """The BRCA1 row-2 guide GAACTCTGAGGACAAAGCAG should be found at 1-based
    position 601 (0-based 600) on the sense strand."""
    seq = (Path := __import__("pathlib").Path(
        __file__).parent / "fixtures" / "nm_007294_first5kb.fasta").read_text()
    # Strip header
    seq = "".join(line for line in seq.splitlines() if not line.startswith(">"))
    guides = find_guides(seq)
    sense_at_600 = [g for g in guides if g.pos == 600 and g.strand == Strand.SENSE]
    assert sense_at_600, "expected a sense guide at 0-based 600"
    g = sense_at_600[0]
    assert g.guide == "GAACTCTGAGGACAAAGCAG"
    assert g.pam == "CGG"


def test_planted_ngg_at_position_200_sense():
    """The second planted guide AATCACCCCTCAAGGAACCA at 1-based 200."""
    seq = (Path := __import__("pathlib").Path(
        __file__).parent / "fixtures" / "nm_007294_first5kb.fasta").read_text()
    seq = "".join(line for line in seq.splitlines() if not line.startswith(">"))
    guides = find_guides(seq)
    sense_at_199 = [g for g in guides if g.pos == 199 and g.strand == Strand.SENSE]
    assert sense_at_199
    g = sense_at_199[0]
    assert g.guide == "AATCACCCCTCAAGGAACCA"
    assert g.pam == "GGG"


def test_non_ngg_pam_is_rejected():
    """The planted AAT PAM at 1-based 400 must NOT produce a guide."""
    seq = (Path := __import__("pathlib").Path(
        __file__).parent / "fixtures" / "nm_007294_first5kb.fasta").read_text()
    seq = "".join(line for line in seq.splitlines() if not line.startswith(">"))
    guides = find_guides(seq)
    # No guide should be at 0-based 399 (the AAT PAM position)
    assert not any(g.pos == 399 and g.strand == Strand.SENSE for g in guides)


def test_boundary_overhang_skipped():
    """find_guides must not produce entries that overhang the sequence end."""
    seq = "ACGTACGTACGTACGTACGTACGTACGT"  # 28 nt; needs 23 to fit one guide
    guides = find_guides(seq)
    for g in guides:
        assert 0 <= g.pos <= len(seq) - 23, f"overhang at {g.pos}"
