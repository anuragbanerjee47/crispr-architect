"""Exon feature-table parsing and guide filtering."""
from __future__ import annotations

from pathlib import Path

from grna_designer.exons import parse_feature_table, filter_guides_to_exons
from grna_designer.models import ExonInterval, Guide, Strand


def test_parse_feature_table_extracts_cds_intervals():
    text = (Path(__file__).parent / "fixtures" / "nm_007294_feature_table.txt").read_text()
    intervals = parse_feature_table(text)
    cds = [e for e in intervals if e.feature_type == "CDS"]
    assert cds == [ExonInterval(150, 700, "CDS")]
    utr5 = [e for e in intervals if e.feature_type == "5UTR"]
    assert utr5 == [ExonInterval(1, 100, "5UTR")]


def test_filter_guides_to_cds_drops_utr_guides():
    """A guide at 1-based position 80 (in 5'UTR) should be dropped under the
    default (no --include-utr) policy."""
    guide_utr = Guide(guide="A" * 20, pam="AGG",
                      strand=Strand.SENSE, pos=79)  # 0-based 79 -> 1-based 80
    guide_cds = Guide(guide="C" * 20, pam="CGG",
                      strand=Strand.SENSE, pos=200)  # 0-based 200 -> 1-based 201, in CDS
    intervals = [
        ExonInterval(1, 100, "5UTR"),
        ExonInterval(150, 700, "CDS"),
    ]
    kept = filter_guides_to_exons(
        [guide_utr, guide_cds], intervals, include_utr=False
    )
    assert kept == [guide_cds]


def test_filter_guides_with_include_utr_keeps_utr_guides():
    guide_utr = Guide(guide="A" * 20, pam="AGG",
                      strand=Strand.SENSE, pos=79)
    guide_cds = Guide(guide="C" * 20, pam="CGG",
                      strand=Strand.SENSE, pos=200)
    intervals = [
        ExonInterval(1, 100, "5UTR"),
        ExonInterval(150, 700, "CDS"),
    ]
    kept = filter_guides_to_exons(
        [guide_utr, guide_cds], intervals, include_utr=True
    )
    assert guide_utr in kept and guide_cds in kept


def test_filter_guides_no_intervals_is_pass_through():
    guides = [Guide(guide="A" * 20, pam="AGG", strand=Strand.SENSE, pos=0)]
    assert filter_guides_to_exons(guides, [], include_utr=False) == guides
