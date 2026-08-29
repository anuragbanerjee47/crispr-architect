"""grna_designer v2.

Public surface:
    Guide, score_guide, find_guides, get_enzyme, OffTargetCounts, OffTargetMode, Strand.
"""
from .cloning import annotate_guide, design_cloning_oligos, nearest_neighbor_tm, write_vendor_csv
from .dual_guide import calculate_deletion_product, design_dual_guide_deletion, design_paired_nickase
from .guides import find_guides, get_enzyme
from .hdr import attach_hdr_donor, design_ssodn, silent_pam_mutation
from .models import Guide, OffTargetCounts, OffTargetMode, ScoringResult, Strand
from .primers import attach_validation_primers, design_validation_primers, primer_pair_quality
from .scoring import (
    guide_quality_flags,
    hsu_zhang_penalty,
    hsu_zhang_specificity_score,
    score_guide,
)

__version__ = "2.0.0"

__all__ = [
    "Guide",
    "OffTargetCounts",
    "OffTargetMode",
    "ScoringResult",
    "Strand",
    "annotate_guide",
    "attach_hdr_donor",
    "attach_validation_primers",
    "calculate_deletion_product",
    "design_cloning_oligos",
    "design_dual_guide_deletion",
    "design_paired_nickase",
    "design_ssodn",
    "design_validation_primers",
    "find_guides",
    "get_enzyme",
    "guide_quality_flags",
    "hsu_zhang_penalty",
    "hsu_zhang_specificity_score",
    "nearest_neighbor_tm",
    "primer_pair_quality",
    "score_guide",
    "silent_pam_mutation",
    "write_vendor_csv",
    "__version__",
]
