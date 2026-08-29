"""NCBI feature-table parsing and exon/CDS filtering of guides."""
from __future__ import annotations

import re
from typing import List

from .models import ExonInterval, Guide

# A single line of an NCBI feature-table interval. e.g. "     5UTR     1..226"
# or "     exon     1..226" or "            226..427" (continuation).
_HEADER_KEYS = {
    "5UTR", "3UTR", "UTR", "exon", "CDS", "gene", "mRNA", "misc_feature",
    "precursor_RNA", "primary_transcript", "stem_loop", "regulatory",
}


def parse_feature_table(text: str) -> List[ExonInterval]:
    """Parse the NCBI `ft` rettype into a list of ExonInterval.

    The format is plain text. Each feature starts with a header line whose
    first whitespace-separated token is a known feature key. Subsequent lines
    that begin with whitespace + an interval belong to the same feature.
    Coordinates are 1-based inclusive.
    """
    intervals: List[ExonInterval] = []
    current_key: str | None = None

    # Interval pattern: 123..456  (NCBI uses '<' and '>' to indicate partials)
    int_re = re.compile(r"<?(\d+)\.\.>?(\d+)")
    comp_re = re.compile(r"complement\(<?(\d+)\.\.>?(\d+)\)")

    for raw in text.splitlines():
        line = raw.rstrip()
        if not line or line.startswith(">"):
            # New record header like ">Feature NM_007294.4"
            current_key = None
            continue
        # Strip leading whitespace; split on whitespace.
        stripped = line.lstrip()
        if not stripped:
            continue
        tokens = stripped.split()
        first = tokens[0]

        if first in _HEADER_KEYS:
            current_key = first
            # Header line may itself carry the first interval.
            rest = stripped[len(first):].strip()
        else:
            # Continuation interval under the current feature key.
            rest = stripped

        if current_key is None:
            continue

        # Find intervals in the rest of the line.
        for m in int_re.finditer(rest):
            start, end = int(m.group(1)), int(m.group(2))
            intervals.append(ExonInterval(start=start, end=end,
                                          feature_type=current_key))
        # Ignore complement() for v2 (we only deal with mRNA-coordinate
        # space; antisense coordinates are already handled in guides.py).
        _ = comp_re  # explicit no-op for clarity

    return intervals


def filter_guides_to_exons(
    guides: List[Guide],
    intervals: List[ExonInterval],
    *,
    include_utr: bool = False,
) -> List[Guide]:
    """Drop guides outside exons, or outside the CDS when include_utr=False.

    Behavior:
    - `include_utr=False` (default): keep guides whose 1-based position falls
      inside a CDS interval.
    - `include_utr=True`: keep guides inside any exon, regardless of feature
      type. Guides still outside all exon intervals are dropped.
    """
    if not intervals:
        return list(guides)

    if include_utr:
        keep_intervals = [e for e in intervals if e.feature_type in {"exon", "CDS", "5UTR", "3UTR"}]
    else:
        keep_intervals = [e for e in intervals if e.feature_type == "CDS"]
        if not keep_intervals:
            # Fallback: no CDS annotations in this record; fall back to exon.
            keep_intervals = [e for e in intervals if e.feature_type == "exon"]

    def _keeps(g: Guide) -> bool:
        pos1 = g.pos + 1
        return any(e.contains(pos1) for e in keep_intervals)

    return [g for g in guides if _keeps(g)]
