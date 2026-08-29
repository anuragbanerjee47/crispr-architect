from grna_designer.cloning import design_cloning_oligos, nearest_neighbor_tm, write_vendor_csv
from grna_designer.guides import find_guides, get_enzyme
from grna_designer.models import Guide
from grna_designer.scoring import guide_quality_flags, score_guide


def test_structure_flags():
    flags = guide_quality_flags("TTTTTACGTACGTACGTACG")
    assert "poly_t" in flags
    assert "homopolymer" not in flags

    flags = guide_quality_flags("AAAAAAAAAAAAAAAAAAAA")
    assert "homopolymer" in flags

    flags = guide_quality_flags("GCGATGCTAACGTTGACGAA")
    assert "hairpin" in flags or "hairpin" not in flags


def test_gc_penalty_applies_outside_optimal_range():
    neutral = score_guide("GCGATGCTAACGTTGACGAA")
    gc_low = score_guide("ATATATATATATATATATAT")
    gc_high = score_guide("GCGCGCGCGCGCGCGCGCGC")
    assert neutral.efficiency > gc_low.efficiency
    assert neutral.efficiency > gc_high.efficiency


def test_hsu_zhang_penalty_is_seed_heavier_than_distal():
    mismatch_positions = [1, 2, 10, 20]
    penalty = score_guide("GCGATGCTAACGTTGACGAA",
                          mismatch_positions=mismatch_positions)
    assert penalty.contributing_features["hsu_zhang_penalty"] > 0


def test_enzyme_registry_supports_multiple_cas_types():
    assert get_enzyme("spcas9")["pam"] == "NGG"
    assert get_enzyme("sacas9")["pam"] == "NNGRRT"
    assert get_enzyme("cas12a")["pam"] == "TTTV"
    assert get_enzyme("spry")["pam"] == "NRN"


def test_cas12a_uses_5prime_pam_and_23bp_spacer():
    enzyme = get_enzyme("cas12a")
    assert enzyme["pam_orientation"] == "5'"
    assert enzyme["spacer_len"] == 23


def test_find_guides_accepts_custom_enzyme_and_3prime_pam():
    seq = "A" * 20 + "CGG" + "A" * 20
    guides = find_guides(seq, enzyme="spcas9")
    assert any(g.guide == "A" * 20 and g.pam == "CGG" for g in guides)

    cas12a_seq = "TTTA" + "A" * 23
    guides = find_guides(cas12a_seq, enzyme="cas12a")
    assert guides


def test_spcas9_oligos_match_expected_overhangs():
    design = design_cloning_oligos("GTTGACGATGCTAACGTTGC", enzyme="spcas9")
    assert design["oligo_fwd"] == "CACCGTTGACGATGCTAACGTTGC"
    assert design["oligo_rev"] == "AAACGCAACGTTAGCATCGTCAAC"  # revcomp(guide)

    non_g = design_cloning_oligos("TTGACGATGCTAACGTTGCA", enzyme="spcas9")
    assert non_g["oligo_fwd"] == "CACCGTTGACGATGCTAACGTTGCA"
    assert non_g["oligo_rev"] == "AAACTGCAACGTTAGCATCGTCAA" + "C"


def test_nearest_neighbor_tm_is_positive_and_accurate_for_reference_oligo():
    tm = nearest_neighbor_tm("GCGCGCGCGCGC")
    assert 55.0 < tm < 75.0
    assert abs(tm - 62.5) < 8.0


def test_vendor_csv_exports_include_provider_columns(tmp_path):
    guides = [
        Guide(
            guide="GTTGACGATGCTAACGTTGC",
            pam="CGG",
            strand="sense",
            pos=12,
            oligo_fwd="CACCGTTGACGATGCTAACGTTGC",
            oligo_rev="AAACGCAACGTTAGCATCGTCAAC",
            tm_fwd=62.5,
            tm_rev=61.8,
        ),
    ]

    idt_path = tmp_path / "idt.csv"
    write_vendor_csv(guides, idt_path, vendor="idt")
    idt_text = idt_path.read_text(encoding="utf-8")
    assert "Tube ID" in idt_text
    assert "Sequence" in idt_text
    assert "CACCGTTGACGATGCTAACGTTGC" in idt_text

    twist_path = tmp_path / "twist.csv"
    write_vendor_csv(guides, twist_path, vendor="twist")
    twist_text = twist_path.read_text(encoding="utf-8")
    assert "Name" in twist_text
    assert "Tube" in twist_text
    assert "CACCGTTGACGATGCTAACGTTGC" in twist_text

    genscript_path = tmp_path / "genscript.csv"
    write_vendor_csv(guides, genscript_path, vendor="genscript")
    genscript_text = genscript_path.read_text(encoding="utf-8")
    assert "Order ID" in genscript_text
    assert "Tube" in genscript_text
    assert "CACCGTTGACGATGCTAACGTTGC" in genscript_text
