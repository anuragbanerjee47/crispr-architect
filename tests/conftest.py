"""Test fixtures and offline mocks.

This file is auto-loaded by pytest. It must not hit the network.
"""
from __future__ import annotations

import io
import os
from pathlib import Path
from typing import Callable

import pytest

# Ensure the package is importable when running pytest from the project root.
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in os.sys.path:
    os.sys.path.insert(0, str(_ROOT))

FIXTURES = Path(__file__).resolve().parent / "fixtures"


@pytest.fixture(autouse=True)
def _set_entrez_email():
    """Set a placeholder NCBI email so any code path that touches Entrez
    doesn't fail the policy check."""
    from Bio import Entrez
    Entrez.email = "pytest@example.com"
    yield


@pytest.fixture
def fixture_fasta_path() -> Path:
    return FIXTURES / "nm_007294_first5kb.fasta"


@pytest.fixture
def fixture_fasta_text() -> str:
    return (FIXTURES / "nm_007294_first5kb.fasta").read_text(encoding="utf-8")


@pytest.fixture
def fixture_feature_table_text() -> str:
    return (FIXTURES / "nm_007294_feature_table.txt").read_text(encoding="utf-8")


@pytest.fixture
def fixture_cas_offinder_output() -> str:
    return (FIXTURES / "cas_offinder_sample_output.tsv").read_text(encoding="utf-8")


@pytest.fixture
def mock_entrez(monkeypatch, fixture_fasta_text, fixture_feature_table_text) -> None:
    """Patch `EntrezClient._efetch` to return fixture text instead of hitting NCBI."""
    from grna_designer import ncbi

    def fake_efetch(self, **kwargs):
        db = kwargs.get("db")
        rettype = kwargs.get("rettype")
        if db == "nuccore" and rettype == "fasta":
            return fixture_fasta_text
        if db == "nuccore" and rettype == "ft":
            return fixture_feature_table_text
        raise AssertionError(f"Unexpected efetch: {db} {rettype}")

    monkeypatch.setattr(ncbi.EntrezClient, "_efetch", fake_efetch)


@pytest.fixture
def mock_cas_offinder(monkeypatch, tmp_path, fixture_cas_offinder_output) -> Callable:
    """Patch CasOffinderSearch._run_one to copy the fixture output file in
    place of running a real binary."""
    from grna_designer import offtarget

    def fake_run_one(self, guide_23nt, in_tsv, out_tsv):
        out_tsv.write_text(fixture_cas_offinder_output, encoding="utf-8")

    monkeypatch.setattr(offtarget.CasOffinderSearch, "_run_one", fake_run_one)
    return fake_run_one
