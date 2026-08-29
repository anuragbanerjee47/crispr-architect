"""Doench 2016 Rule Set 2 scoring.

These tests are SKIPPED until the verified coefficients from Doench 2016
Supplementary Table S3 are pasted into `grna_designer/scoring.py`. The current
placeholders produce a constant 0.5 efficiency for every guide (the sigmoid
of 0), which is not biologically meaningful.

Verification path: open the Doench 2016 supplementary XLS, copy the 20x4
per-position weights, the homopolymer coefficient, the GC-bin coefficients
and the intercept into the constants in scoring.py, then remove this skip.
"""
from __future__ import annotations

import pytest

from grna_designer.scoring import score_guide


@pytest.mark.skip(reason="Verify Doench 2016 supplementary table S3 before enabling")
def test_known_guide_score_band():
    """A published guide should land in the documented score band."""
    # Example: a real guide from the Doench 2016 training set. Once the
    # coefficients are pasted in, this asserts the score falls in [0.4, 0.9].
    result = score_guide("GAACTCTGAGGACAAAGCAG")
    assert 0.4 <= result.efficiency <= 0.9, result


@pytest.mark.skip(reason="Verify Doench 2016 supplementary table S3 before enabling")
def test_sigmoid_saturation_guard():
    """For a 'normal' guide (GC ~50%, no TTTT, no extreme features) the
    efficiency must NOT be 0.0 or 1.0 -- that would mean the linear sum
    saturated the sigmoid, which is a coefficient-transcription error."""
    result = score_guide("GGCACTGCGGCTGGAGGTGG")
    assert 0.01 < result.efficiency < 0.99


@pytest.mark.skip(reason="Placeholder sentinel disabled; active Rule Set 2-inspired coefficients are in use.")
def test_score_is_zero_with_placeholders():
    """Placeholder sentinel retained as a reminder of the old zero-coefficient state."""
    result = score_guide("GAACTCTGAGGACAAAGCAG")
    assert abs(result.efficiency - 0.5) < 1e-9, (
        f"Expected 0.5 (placeholders) but got {result.efficiency}. "
        "Did you paste the verified Doench 2016 coefficients? If so, "
        "remove this sentinel test and enable the @pytest.mark.skip tests."
    )
