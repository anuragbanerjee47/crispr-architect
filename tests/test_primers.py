from grna_designer.models import Guide, Strand
from grna_designer.primers import design_validation_primers, primer_pair_quality


def test_validation_primers_have_valid_length_and_tm_ranges():
    seq = "A" * 120 + "GCGATGCTAACGTTGACGAA" * 12 + "T" * 120
    guide = Guide(
        guide="GCGATGCTAACGTTGACGAA",
        pam="CGG",
        strand=Strand.SENSE,
        pos=120,
    )
    design = design_validation_primers(seq, guide)

    assert 18 <= len(design["fwd_primer"]) <= 24
    assert 18 <= len(design["rev_primer"]) <= 24
    assert 58.0 <= design["tm_fwd"] <= 62.0
    assert 58.0 <= design["tm_rev"] <= 62.0
    assert 40.0 <= design["gc_fwd"] <= 60.0
    assert 40.0 <= design["gc_rev"] <= 60.0
    assert design["amplicon_length"] >= 400
    assert design["amplicon_length"] <= 600


def test_amplicon_alignment_tracks_the_cut_site():
    seq = "A" * 400 + "GCGATGCTAACGTTGACGAA" + "T" * 400
    guide = Guide(
        guide="GCGATGCTAACGTTGACGAA",
        pam="CGG",
        strand=Strand.SENSE,
        pos=400,
    )
    design = design_validation_primers(seq, guide)

    cut_site = guide.pos + len(guide.guide) - 3
    assert design["amplicon_start"] <= cut_site <= design["amplicon_end"]
    assert design["amplicon_length"] <= 600
    assert design["amplicon_length"] >= 400


def test_boundary_cut_sites_still_return_valid_primers():
    seq = "ATGC" * 60
    guide = Guide(
        guide="ATGCGATGCGATGCGATGCG",
        pam="CGG",
        strand=Strand.SENSE,
        pos=2,
    )
    design = design_validation_primers(seq, guide)

    assert 0 <= design["amplicon_start"] < len(seq)
    assert design["amplicon_end"] <= len(seq)
    assert design["amplicon_length"] > 0
    assert len(design["fwd_primer"]) >= 18
    assert len(design["rev_primer"]) >= 18


def test_primer_pair_is_not_self_complementary_or_dimeric():
    fwd = "GCGATGCTAACGTTGACGAA"
    rev = "TTCGTCAACGTTAGCATCGC"
    quality = primer_pair_quality(fwd, rev)

    assert quality["primer_dimer_risk"] is False
    assert quality["self_complementarity_risk"] is False
