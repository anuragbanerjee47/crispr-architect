"""Shared dataclasses and enums for gRNA designer v2."""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional


class Strand(str, Enum):
    SENSE = "sense"
    ANTISENSE = "antisense"


class OffTargetMode(str, Enum):
    TRANSCRIPT = "transcript"
    GENOME = "genome"


@dataclass
class OffTargetCounts:
    """Number of off-target loci at 0/1/2/3 mismatches (post-cluster-merge)."""
    zero: int = 0
    one: int = 0
    two: int = 0
    three: int = 0

    def total(self) -> int:
        return self.zero + self.one + self.two + self.three


@dataclass
class ScoringResult:
    """Output of Rule Set 2 scoring for a single guide."""
    linear: float
    efficiency: float
    on_target_score: float = 0.0
    specificity_score: float = 100.0
    contributing_features: dict = field(default_factory=dict)


@dataclass
class ExonInterval:
    """1-based inclusive mRNA-coordinate interval, in the same frame as the
    sequence returned by `Entrez.efetch(db="nuccore", rettype="fasta")`."""
    start: int
    end: int
    feature_type: str  # "exon", "CDS", "5UTR", "3UTR", "misc_feature", ...

    def contains(self, pos_1based: int) -> bool:
        return self.start <= pos_1based <= self.end


@dataclass
class Guide:
    """A candidate SpCas9 protospacer + PAM."""
    guide: str                # 20-nt, 5'->3' as written on the oligo
    pam: str                  # 3-nt NGG
    strand: Strand
    pos: int                  # 0-based start in the original fetched sequence
    # Populated by scoring
    efficiency: float = 0.0
    # Populated by off-target search
    offtarget: OffTargetCounts = field(default_factory=OffTargetCounts)
    specificity: float = 1.0
    # Populated by exon filter
    in_exon: bool = True
    # Human-readable flags from QC passes
    issues: str = ""
    # Wet-lab cloning oligos and Tm values
    oligo_fwd: str = ""
    oligo_rev: str = ""
    tm_fwd: float = 0.0
    tm_rev: float = 0.0
    flanking_fwd_primer: str = ""
    flanking_rev_primer: str = ""
    amplicon_length: int = 0
    amplicon_start: int = 0
    amplicon_end: int = 0
    # HDR donor / ssODN builder metadata
    ssodn_sequence: str = ""
    donor_type: str = ""
    homology_5: int = 0
    homology_3: int = 0
    edit_index: int = 0
    strand_polarity: str = "non_target"
    donor_notes: str = ""

    def __post_init__(self) -> None:
        if isinstance(self.strand, str):
            self.strand = Strand(self.strand.lower())

    @property
    def strand_name(self) -> str:
        return self.strand.value if isinstance(self.strand, Strand) else str(self.strand).lower()

    def combined_score(self, *, efficiency_weight: float = 0.6) -> float:
        return round(
            efficiency_weight * self.efficiency
            + (1.0 - efficiency_weight) * self.specificity,
            3,
        )
