import pytest

from grna_designer.security import (
    MAX_SEQUENCE_LENGTH,
    sanitize_accession,
    sanitize_batch_text,
    sanitize_gene_symbol,
    sanitize_sequence,
)


def test_sanitize_sequence_accepts_standard_iupac_bases():
    seq = sanitize_sequence(">seq\nACGTN\nACGT")
    assert seq == "ACGTNACGT"


def test_sanitize_sequence_rejects_invalid_base_characters():
    with pytest.raises(ValueError):
        sanitize_sequence("ACGTX")


def test_sanitize_sequence_rejects_malicious_payloads():
    with pytest.raises(ValueError):
        sanitize_sequence("ACGT<script>alert(1)</script>")


def test_sanitize_accession_accepts_standard_ncbi_patterns():
    assert sanitize_accession("NM_007294") == "NM_007294"
    assert sanitize_accession("NC_000001") == "NC_000001"


def test_sanitize_accession_rejects_malformed_or_injected_values():
    with pytest.raises(ValueError):
        sanitize_accession("; rm -rf /")
    with pytest.raises(ValueError):
        sanitize_accession("invalid")


def test_sanitize_gene_symbol_rejects_injection_payloads():
    with pytest.raises(ValueError):
        sanitize_gene_symbol("BRCA1; DROP TABLE genes")


def test_sanitize_batch_text_rejects_oversized_payloads():
    with pytest.raises(ValueError):
        sanitize_batch_text("A" * (MAX_SEQUENCE_LENGTH + 1))


def test_sanitize_sequence_enforces_max_size_limit():
    with pytest.raises(ValueError):
        sanitize_sequence("A" * (MAX_SEQUENCE_LENGTH + 1))
