import pytest

from grna_designer.dual_guide import (
    calculate_deletion_product,
    design_dual_guide_deletion,
    design_paired_nickase,
)


def test_dual_guide_deletion_coordinates_and_size():
    seq = "A" * 100 + "GCGATGCTAACGTTGACGAA" + "T" * 100
    design = design_dual_guide_deletion(seq, start=120, end=140)

    assert design["deletion_start"] == 120
    assert design["deletion_end"] == 140
    assert design["deletion_size"] == 21
    assert design["left_guide"]
    assert design["right_guide"]
    assert design["ligation_junction"] == (120, 141)


def test_paired_nickase_uses_opposite_strands_and_valid_offset():
    seq = "A" * 200 + "GCGATGCTAACGTTGACGAA" * 3 + "T" * 200
    design = design_paired_nickase(seq, target_start=240, target_end=280)

    assert design["guide_a"]["strand"] != design["guide_b"]["strand"]
    assert 30 <= design["offset"] <= 70
    assert design["target_span"] == (240, 280)


def test_empty_and_edge_cases_are_handled_cleanly():
    with pytest.raises(ValueError):
        design_dual_guide_deletion("", start=0, end=0)

    edge = design_dual_guide_deletion("ACGTACGT", start=0, end=3)
    assert edge["deletion_size"] == 4
    assert edge["mask_start"] == 0
    assert edge["mask_end"] == 4

    with pytest.raises(ValueError):
        design_paired_nickase("", target_start=0, target_end=0)


def test_deletion_product_calculation_matches_span():
    product = calculate_deletion_product(100, 130, 500)
    assert product["deletion_size"] == 31
    assert product["left_flank_start"] == 0
    assert product["right_flank_end"] == 500
