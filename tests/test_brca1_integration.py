"""End-to-end: with Entrez mocked, run the full pipeline and assert a known
BRCA1 guide lands in the top-N of the combined-score ranking.

This is the canary test. If this passes, the whole pipeline is wired correctly.
"""
from __future__ import annotations

import io
import sys
from pathlib import Path

import pytest

from grna_designer.cli import main as cli_main
from grna_designer.models import Strand


def test_brca1_known_guide_appears_in_top_30(mock_entrez, capsys):
    """Run the full CLI on NM_007294 (with mocks) and verify a known
    guide from brca1_known_guides.csv appears in the top 30."""
    known = Path(__file__).parent / "fixtures" / "brca1_known_guides.csv"
    known_guides = []
    for line in known.read_text().splitlines()[1:]:
        guide_seq, band, src = line.split(",")
        known_guides.append(guide_seq)

    # We do not run the off-target search (faster; not relevant to this test).
    # Score-sort is still applied.
    rc = cli_main([
        "NM_007294",
        "--window", "1000",
        "--top", "30",
        "--no-offtarget",
    ])
    assert rc == 0, f"CLI exited with {rc}"

    captured = capsys.readouterr().out
    # At least one known guide should appear in the printed table.
    assert any(g in captured for g in known_guides), (
        f"None of the known BRCA1 guides appeared in the top 30. "
        f"Output:\n{captured}"
    )


def test_brca1_known_guide_appears_in_csv(mock_entrez, tmp_path, capsys):
    """Same as above, but check the CSV output."""
    known = Path(__file__).parent / "fixtures" / "brca1_known_guides.csv"
    known_guides = [line.split(",")[0] for line in known.read_text().splitlines()[1:]]

    out_csv = tmp_path / "guides.csv"
    rc = cli_main([
        "NM_007294",
        "--window", "1000",
        "--top", "30",
        "--no-offtarget",
        "--csv", str(out_csv),
    ])
    assert rc == 0
    csv_text = out_csv.read_text()
    assert any(g in csv_text for g in known_guides), csv_text
