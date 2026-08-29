"""Security and input sanitization helpers for gRNA designer.

The project keeps all validation logic in memory and never writes user-supplied
sequence data to disk. This module is deliberately strict: only canonical IUPAC
DNA/RNA symbols are accepted for nucleotide inputs, and all web inputs are
bounded to prevent untrusted payloads from exhausting memory.
"""
from __future__ import annotations

import re
from typing import List

IUPAC_BASES = set("ACGTNUR")
MAX_SEQUENCE_LENGTH = 5_000_000
MAX_BATCH_BYTES = 5_000_000

ACCESSION_RE = re.compile(
    r"^(?:[A-Z]{1,2}\d{5,}|(?:NM|NC|XM|NR|NP|NG|NT|NW|NZ|AP|AC|AF|AY|EU|U[0-9]|M[0-9])_[A-Z0-9.]+)$",
    re.IGNORECASE,
)
GENE_SYMBOL_RE = re.compile(r"^[A-Za-z0-9._-]+(?:\s+[A-Za-z0-9._-]+)*$")
BLOCKED_TOKENS = (
    "--",
    ";",
    "||",
    "&&",
    "<script",
    "DROP TABLE",
    "DELETE FROM",
    "INSERT INTO",
    "SELECT *",
    "javascript:",
    "data:",
)


def _contains_blocked_payload(value: str) -> bool:
    candidate = value.upper()
    for token in BLOCKED_TOKENS:
        if token in candidate:
            return True
    return False


def _strip_fasta_headers(raw: str) -> str:
    """Extract nucleotide sequence text from FASTA or raw sequence input."""
    if raw is None:
        raise ValueError("Input is empty.")
    text = str(raw).strip()
    if not text:
        raise ValueError("Input is empty.")
    if len(text.encode("utf-8")) > MAX_BATCH_BYTES:
        raise ValueError(
            f"Input exceeds the safe memory limit of {MAX_BATCH_BYTES} bytes.")
    if _contains_blocked_payload(text):
        raise ValueError("Input contains blocked payload content.")

    sequence_chunks: List[str] = []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith(">"):
            continue
        sequence_chunks.append(s)

    if not sequence_chunks:
        raise ValueError("No DNA/RNA sequence content was found.")

    sequence = "".join(sequence_chunks)
    return re.sub(r"\s+", "", sequence)


def sanitize_sequence(raw: str, *, field_name: str = "sequence", max_length: int = MAX_SEQUENCE_LENGTH) -> str:
    """Validate and normalize a DNA/RNA sequence. Only IUPAC symbols are accepted."""
    raw_text = _strip_fasta_headers(raw)
    raw_text = raw_text.upper().replace("U", "T")

    if len(raw_text) > max_length:
        raise ValueError(
            f"{field_name} exceeds the maximum allowed length of {max_length} bases.")
    if not raw_text:
        raise ValueError(f"{field_name} is empty.")
    if not re.fullmatch(r"[ACGTN]+", raw_text):
        raise ValueError(
            f"{field_name} contains unsupported characters. Only A, C, G, T, N are allowed.")
    return raw_text


def validate_sequence(raw: str, *, field_name: str = "sequence", max_length: int = MAX_SEQUENCE_LENGTH) -> str:
    """Compatibility wrapper kept for callers expecting a validation function name."""
    return sanitize_sequence(raw, field_name=field_name, max_length=max_length)


def sanitize_accession(raw: str) -> str:
    """Sanitize and validate an NCBI accession ID."""
    text = str(raw or "").strip()
    if not text:
        raise ValueError("An NCBI accession is required.")
    if len(text) > 32:
        raise ValueError("Accession ID is too long.")
    if _contains_blocked_payload(text):
        raise ValueError("Accession contains invalid or malicious content.")
    if not ACCESSION_RE.fullmatch(text):
        raise ValueError(
            "Invalid NCBI accession format. Use patterns such as NM_123456 or NC_000001.")
    return text.upper()


def sanitize_gene_symbol(raw: str) -> str:
    """Validate user-facing gene symbols or search terms."""
    text = str(raw or "").strip()
    if not text:
        raise ValueError("A gene symbol or search term is required.")
    if _contains_blocked_payload(text):
        raise ValueError("Gene symbol contains invalid or malicious content.")
    if len(text) > 128:
        raise ValueError("Gene symbol exceeds the maximum supported length.")
    if not GENE_SYMBOL_RE.fullmatch(text):
        raise ValueError(
            "Gene symbol contains unsupported characters. Use letters, numbers, underscores, hyphens, dots, or spaces only.")
    return text.strip()


def sanitize_batch_text(raw: str, *, max_entries: int = 200) -> List[str]:
    """Validate a batch of accessions or sequences kept in memory only."""
    text = str(raw or "").strip()
    if not text:
        raise ValueError("Batch input is empty.")
    if len(text.encode("utf-8")) > MAX_BATCH_BYTES:
        raise ValueError(
            f"Batch input exceeds the safe memory limit of {MAX_BATCH_BYTES} bytes.")

    entries: List[str] = []
    for line in text.splitlines():
        value = line.strip()
        if not value:
            continue
        if value.startswith(">"):
            if len("".join(entries)) + len(value) > MAX_BATCH_BYTES:
                raise ValueError(
                    f"Batch input exceeds the safe memory limit of {MAX_BATCH_BYTES} bytes.")
            entries.append(value)
            continue
        if " " in value and not re.fullmatch(r"[A-Za-z0-9._-]+(?:\s+[A-Za-z0-9._-]+)*", value):
            raise ValueError("Batch input contains unsupported characters.")
        if not re.fullmatch(r"[A-Za-z0-9_\-.]+", value):
            raise ValueError("Batch input contains unsupported characters.")
        if len(value) > MAX_SEQUENCE_LENGTH:
            raise ValueError(
                f"Batch entry exceeds the maximum sequence length of {MAX_SEQUENCE_LENGTH} bases.")
        if len("".join(entries)) + len(value) > MAX_BATCH_BYTES:
            raise ValueError(
                f"Batch input exceeds the safe memory limit of {MAX_BATCH_BYTES} bytes.")
        entries.append(value)
    if not entries:
        raise ValueError("No valid batch entries were found.")
    if len(entries) > max_entries:
        raise ValueError(
            f"Batch input exceeds the maximum of {max_entries} entries.")
    return entries
