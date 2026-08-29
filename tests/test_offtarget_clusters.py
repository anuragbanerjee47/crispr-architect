"""Off-target cluster merging logic.

Sense and antisense hits at the same genomic locus must collapse into one
cluster; on-target cluster is subtracted.
"""
from __future__ import annotations

from grna_designer.offtarget import merge_clusters, compute_specificity, InTranscriptOffTarget
from grna_designer.models import Guide, OffTargetCounts, Strand


def test_sense_and_antisense_at_same_locus_collapse_to_one():
    """Two intervals overlapping at the same position but on different strands
    must merge into a single cluster."""
    sense = [(100, 120)]
    anti = [(80, 100)]  # overlaps the sense interval at position 100
    # On-target is at position 50, far from the cluster.
    assert merge_clusters(sense, anti, on_target_start=50) == 1


def test_disjoint_clusters_counted_individually():
    sense = [(100, 120), (500, 520), (800, 820)]
    anti = []
    # On-target at 100 -- one cluster contains the on-target, so result is 2.
    assert merge_clusters(sense, anti, on_target_start=110) == 2
    # On-target not inside any cluster -- all 3 remain.
    assert merge_clusters(sense, anti, on_target_start=99999) == 3


def test_empty_input_returns_zero():
    assert merge_clusters([], [], on_target_start=0) == 0


def test_touching_intervals_merge():
    """Two intervals that touch at a single base (e <= s) merge into one cluster."""
    sense = [(100, 120), (120, 140)]
    assert merge_clusters(sense, [], on_target_start=99999) == 1


def test_in_transcript_search_on_synthetic_sequence():
    """Plant a guide with a known off-target and verify InTranscriptOffTarget
    finds it at 0 mismatches."""
    # 60-nt sequence; guide at position 0 (1-based 1) with NGG PAM
    seq = "GAACTCTGAGGACAAAGCAG" + "CGG" + "A" * 37
    # Plant an exact-match off-target at position 30 (1-based 31) on sense.
    # We replace positions 30..50 with the same guide.
    seq_list = list(seq)
    for i, base in enumerate("GAACTCTGAGGACAAAGCAG"):
        seq_list[30 + i] = base
    seq = "".join(seq_list)

    guide = Guide(guide="GAACTCTGAGGACAAAGCAG", pam="CGG",
                  strand=Strand.SENSE, pos=0)
    engine = InTranscriptOffTarget(seq, mismatches=3)
    counts = engine.search(guide)
    # Expect at least 1 zero-mm cluster (the planted off-target).
    assert counts.zero >= 1


def test_compute_specificity_pure_function():
    """compute_specificity is independent of the search backend."""
    c0 = compute_specificity(OffTargetCounts(zero=0))
    c1 = compute_specificity(OffTargetCounts(zero=1))
    assert c0 > c1
    assert 0.0 <= c1 <= 1.0
    # Pure monotonicity: more off-targets -> lower specificity.
    assert compute_specificity(OffTargetCounts(zero=10)) < c1
