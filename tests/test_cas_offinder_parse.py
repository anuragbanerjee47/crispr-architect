"""Cas-OFFinder output parsing."""
from __future__ import annotations

from pathlib import Path

from grna_designer.offtarget import _parse_cas_offinder_output


def test_parse_cas_offinder_sample_output():
    text = (Path(__file__).parent / "fixtures" / "cas_offinder_sample_output.tsv").read_text()
    tmp = Path(__file__).parent / "fixtures" / "_tmp_parsed.tsv"
    tmp.write_text(text, encoding="utf-8")
    try:
        rows = _parse_cas_offinder_output(tmp)
    finally:
        tmp.unlink(missing_ok=True)
    # 12 rows in the fixture
    assert len(rows) == 12
    # First row is the on-target at chr1:601, + strand, 0 mm
    assert rows[0] == (601, "+", 0, 0)
    # Second row is the on-target antisense, 0 mm
    assert rows[1] == (601, "-", 0, 0)
    # Some 3-mm row exists
    threes = [r for r in rows if r[2] == 3]
    assert len(threes) >= 1
