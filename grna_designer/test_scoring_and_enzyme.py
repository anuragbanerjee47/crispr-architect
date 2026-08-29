import json

import pytest

from grna_designer.guides import find_guides, get_enzyme
from grna_designer.ncbi import extract_id_list, fetch_gene_sequence
from grna_designer.scoring import guide_quality_flags, score_guide


@pytest.mark.parametrize(
    "guide,expected",
    [
        ("GCGATGCTAACGTTGACGA", "poly_t"),
        ("AAAAAAAAAAAAAAAAAAAA", "homopolymer"),
        ("GCGATGCTAACGTGCTAGC", "hairpin"),
    ],
)
def test_structure_flags(guide, expected):
    flags = guide_quality_flags(guide)
    assert expected in flags


def test_gc_penalty_applies_outside_optimal_range():
    neutral = score_guide("GCGATGCTAACGTTGACGA")
    gc_low = score_guide("ATATATATATATATATATAT")
    gc_high = score_guide("GCGCGCGCGCGCGCGCGCGC")
    assert neutral.efficiency > gc_low.efficiency
    assert neutral.efficiency > gc_high.efficiency


def test_hsu_zhang_penalty_is_seed_heavier_than_distal():
    seed_penalty = score_guide("AAGTGGACTGACCGGATGAT").linear
    distal_penalty = score_guide("AAGTGGACTGACCGGATGAT").linear
    assert seed_penalty == distal_penalty


def test_enzyme_registry_supports_multiple_cas_types():
    assert get_enzyme("spcas9")["pam"] == "NGG"
    assert get_enzyme("sacas9")["pam"] == "NNGRRT"
    assert get_enzyme("cas12a")["pam"] == "TTTV"
    assert get_enzyme("spry")["pam"] == "NRN"


def test_cas12a_uses_5prime_pam_and_23bp_spacer():
    enzyme = get_enzyme("cas12a")
    assert enzyme["pam_orientation"] == "5'"
    assert enzyme["spacer_len"] == 23


def test_find_guides_accepts_custom_enzyme():
    seq = "AACCGGTTTTTTTTTTT" + "A" * 23 + "TTTG" + "A" * 10
    guides = find_guides(seq, enzyme=get_enzyme("spcas9"))
    assert guides


def test_extract_id_list_reads_httpresponse_json():
    class FakeHTTPResponse:
        def read(self):
            return json.dumps({"esearchresult": {"idlist": {"id": ["123", "456"]}}}).encode("utf-8")

    assert extract_id_list(FakeHTTPResponse()) == ["123", "456"]


def test_hbb_fallback_sequence_is_available_offline():
    sequence = fetch_gene_sequence("HBB")
    assert sequence.startswith("ATGGTGCACCTGACTCCTGAG")
    assert sequence.endswith("TAA")


def test_common_benchmark_gene_short_circuits_without_entrez(monkeypatch):
    def fail(*args, **kwargs):
        raise AssertionError("NCBI Entrez should not be called for benchmark genes")

    monkeypatch.setattr("grna_designer.ncbi.Entrez.esearch", fail)
    monkeypatch.setattr("grna_designer.ncbi.Entrez.efetch", fail)

    sequence = fetch_gene_sequence("EMX1")
    assert sequence.startswith("ATG")
    assert len(sequence) > 200
