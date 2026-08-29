# Changelog

## v2.0.0 (in progress)

Breaking changes from v1:
- Hand-rolled efficiency blend replaced with **Doench 2016 Rule Set 2**
  position-weighted coefficients. See `grna_designer/scoring.py` for the
  constants and a verification note against the Doench 2016 supplementary
  table. Coefficients ship as `0.0` placeholders until verified.
- New `--offtarget-mode {transcript,genome}` flag. Default `transcript`
  (the v1 in-transcript scan); `genome` invokes the Cas-OFFinder binary.
- New `--include-utr` flag (default off). When off, guides outside the
  coding region are dropped.
- New `--genome-fasta` and `--cas-offinder-path` flags for genome mode.
- Output columns changed: v1's blend sub-scores (`gc_score`,
  `stability_score`, `palindrome_score`, `poly_score`) are removed.
  `efficiency_score` is now the Rule Set 2 0..1 score.
- Architecture: v1 is `grna_designer_v1.py`; v2 is the `grna_designer/`
  package. Run v2 as `python -m grna_designer` or via the
  `grna-designer` console script.

Additions:
- pytest suite with biological fixtures, including a BRCA1 end-to-end
  canary test that asserts a published guide lands in the top-30.
- `EntrezClient` and `parse_feature_table` for exon/CDS filtering.
- Offline mock layer in `tests/conftest.py` so the suite runs without
  network access.

Removed:
- `Guide.score()` method (replaced by `scoring.score_guide`).
- `max_palindrome_len`, `gc_fraction`, `POLY_RUN_RE` (Rule Set 2 replaces
  the blend they were part of).
- Hard-coded `+0.10` offset in the v1 efficiency blend.

## v1.0.0
- Initial release. Hand-rolled Doench-inspired scoring + in-transcript
  off-target counting. Retained as `grna_designer_v1.py` for reference.
