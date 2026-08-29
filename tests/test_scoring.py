from grna_designer.scoring import (
    guide_quality_flags,
    hsu_zhang_specificity_score,
    score_guide,
)


def test_poly_t_guide_is_penalized_on_target_score():
    good = score_guide("GCGATGCTAACGTTGACGAA")
    bad = score_guide("GCGATGCTAACGTTTTTTTT")
    assert good.on_target_score > bad.on_target_score


def test_seed_mismatches_reduce_specificity_more_than_distal_mismatches():
    seed_risk = hsu_zhang_specificity_score([1, 2, 10, 15])
    distal_risk = hsu_zhang_specificity_score([11, 12, 16, 19])
    assert seed_risk < distal_risk


def test_guide_quality_flags_catch_pol3_terminators():
    flags = guide_quality_flags("GTTTTGCGATGCTAACGTTGC")
    assert "poly_t" in flags
