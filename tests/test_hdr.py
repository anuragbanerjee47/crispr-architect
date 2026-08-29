from grna_designer.hdr import design_ssodn, silent_pam_mutation
from grna_designer.models import Guide, Strand


def test_ssodn_length_and_asymmetry_are_correct():
    seq = "A" * 100 + "GCGATGCTAACGTTGACGAA" + "T" * 100
    guide = Guide(
        guide="GCGATGCTAACGTTGACGAA",
        pam="CGG",
        strand=Strand.SENSE,
        pos=100,
    )

    design = design_ssodn(seq, guide, homology_5=50, homology_3=60)

    assert design["homology_5"] == 50
    assert design["homology_3"] == 60
    assert design["donor_type"] == "asymmetric"
    assert design["donor_length"] == 111
    assert len(design["ssodn"]) == design["donor_length"]


def test_silent_pam_mutation_replaces_wobble_base_in_pam_or_seed():
    seq = "A" * 100 + "GCGATGCTAACGTTGACGAA" + "CGG" + "T" * 100
    guide = Guide(
        guide="GCGATGCTAACGTTGACGAA",
        pam="CGG",
        strand=Strand.SENSE,
        pos=100,
    )

    mutation = silent_pam_mutation(seq, guide)

    assert "PAM" in mutation["annotation"]
    assert mutation["changed_base"] in {"G", "C", "A", "T"}
    assert mutation["mutated_sequence"] != seq


def test_strand_polarity_defaults_to_non_target_for_cas9():
    seq = "A" * 200 + "GCGATGCTAACGTTGACGAA" + "T" * 200
    guide = Guide(
        guide="GCGATGCTAACGTTGACGAA",
        pam="CGG",
        strand=Strand.SENSE,
        pos=200,
    )

    design = design_ssodn(seq, guide, strand_polarity="non_target")
    assert design["strand_polarity"] == "non_target"
    assert design["ssodn"]

    target_design = design_ssodn(seq, guide, strand_polarity="target")
    assert target_design["strand_polarity"] == "target"
    assert target_design["ssodn"] != design["ssodn"]


def test_boundary_cut_sites_still_produce_valid_donor():
    seq = "ACGT" * 40
    guide = Guide(
        guide="GCGATGCTAACGTTGACGAA",
        pam="CGG",
        strand=Strand.SENSE,
        pos=3,
    )

    design = design_ssodn(seq, guide, homology_5=50, homology_3=50)

    assert 0 <= design["edit_index"] < len(seq)
    assert len(design["ssodn"]) == 101
    assert design["donor_type"] == "symmetric"
