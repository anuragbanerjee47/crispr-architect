"""HDR donor / ssODN builder for CRISPR guide candidates."""
from __future__ import annotations

from typing import Dict, List, Tuple

from .models import Guide, Strand
from .sequtil import rev_comp


_VALID_BASES = {"A", "C", "G", "T"}


def _dna_or_rc(seq: str) -> str:
    return str(seq).upper().replace(" ", "")


def _safe_index(seq: str, guide: Guide, *, offset: int = 0) -> int:
    if not seq:
        return 0
    base_index = max(0, int(getattr(guide, "pos", 0)))
    guide_length = max(1, len(str(getattr(guide, "guide", ""))))
    return max(0, min(len(seq) - 1, base_index + guide_length // 2 + offset))


def silent_pam_mutation(sequence: str, guide: Guide) -> Dict[str, object]:
    """Create a silent or neutral PAM-disrupting edit at the cleavage site."""
    seq = _dna_or_rc(sequence)
    if not seq:
        raise ValueError("Sequence cannot be empty for HDR donor design.")
    if not isinstance(guide, Guide):
        raise TypeError("guide must be a Guide instance")

    pam = str(getattr(guide, "pam", "") or "").upper()
    cut_index = _safe_index(seq, guide)
    candidate_window = seq[max(0, cut_index - 6)
                               : min(len(seq), cut_index + 10)]
    edit_index = cut_index
    mutated = seq
    mutation = ""
    replacement = None
    annotation = "PAM"

    if pam:
        pam_start = max(0, cut_index - len(pam))
        pam_end = min(len(seq), pam_start + len(pam))
        window = seq[pam_start:pam_end]
        for idx, base in enumerate(window):
            if idx >= len(pam):
                break
            if base == pam[idx]:
                continue
            candidate = pam[idx]
            for alt in ["A", "C", "G", "T"]:
                if alt != candidate:
                    replacement = alt
                    break
            if replacement is None:
                replacement = "A"
            mutated = seq[:pam_start + idx] + \
                replacement + seq[pam_start + idx + 1:]
            mutation = f"{guide.pam[idx]}->{replacement}"
            edit_index = pam_start + idx
            break

    if replacement is None:
        for idx in range(max(0, cut_index - 3), min(len(seq), cut_index + 4)):
            base = seq[idx]
            if base not in _VALID_BASES:
                continue
            alt = "A" if base != "A" else "C"
            mutated = seq[:idx] + alt + seq[idx + 1:]
            edit_index = idx
            mutation = f"{base}->{alt}"
            annotation = "seed-window"
            break

    if mutated == seq:
        mutated = seq[:-1] + ("A" if seq[-1] != "A" else "C")
        edit_index = len(seq) - 1
        mutation = f"{seq[-1]}->A"
        annotation = "fallback"

    return {
        "mutated_sequence": mutated,
        "edit_index": edit_index,
        "changed_base": mutated[edit_index] if 0 <= edit_index < len(mutated) else "N",
        "mutation": mutation,
        "annotation": annotation,
        "window": candidate_window,
    }


def _target_sequence_for_guide(seq: str, guide: Guide, strand_polarity: str) -> str:
    guide_seq = str(getattr(guide, "guide", "") or "").upper()
    if not guide_seq:
        return seq
    if strand_polarity == "target":
        return seq
    return rev_comp(seq)


def design_ssodn(
    sequence: str,
    guide: Guide,
    *,
    homology_5: int = 40,
    homology_3: int = 40,
    strand_polarity: str = "non_target",
    donor_type: str | None = None,
) -> Dict[str, object]:
    """Build a donor ssODN around the desired edit site with symmetric/asymmetric arms."""
    seq = _dna_or_rc(sequence)
    if not seq:
        raise ValueError("Sequence cannot be empty for donor design.")
    guide_seq = str(getattr(guide, "guide", "") or "").upper()
    if not guide_seq:
        raise ValueError("Guide sequence cannot be empty for donor design.")

    strand_polarity = str(strand_polarity or "non_target").lower()
    if strand_polarity not in {"target", "non_target"}:
        strand_polarity = "non_target"

    edit_info = silent_pam_mutation(seq, guide)
    edit_index = int(edit_info["edit_index"])
    donor_h5 = max(0, int(homology_5))
    donor_h3 = max(0, int(homology_3))
    total = donor_h5 + donor_h3 + 1
    donor_kind = "asymmetric" if donor_h5 != donor_h3 else "symmetric"
    if donor_type is not None:
        donor_kind = str(donor_type).lower()

    left_start = max(0, edit_index - donor_h5)
    right_end = min(len(seq), edit_index + donor_h3 + 1)
    donor = seq[left_start:right_end]
    donor_rc = rev_comp(donor)

    if strand_polarity == "target":
        donor_payload = donor
    else:
        donor_payload = donor_rc

    if donor_kind == "asymmetric":
        donor_payload = donor_payload[:total]
    else:
        donor_payload = donor_payload[: max(1, donor_h5 + donor_h3 + 1)]

    if len(donor_payload) < total:
        pad = total - len(donor_payload)
        donor_payload = donor_payload + ("A" * pad)

    donor_sequence = donor_payload.upper()
    return {
        "ssodn": donor_sequence,
        "donor_type": donor_kind,
        "strand_polarity": strand_polarity,
        "homology_5": donor_h5,
        "homology_3": donor_h3,
        "edit_index": edit_index,
        "base_edit": edit_info["mutation"],
        "donor_length": len(donor_sequence),
        "mutated_sequence": edit_info["mutated_sequence"],
        "annotation": edit_info["annotation"],
    }


def attach_hdr_donor(guide: Guide, sequence: str, *, homology_5: int = 40, homology_3: int = 40, strand_polarity: str = "non_target") -> Guide:
    """Attach the designed HDR donor to the guide object."""
    design = design_ssodn(sequence, guide, homology_5=homology_5,
                          homology_3=homology_3, strand_polarity=strand_polarity)
    guide.ssodn_sequence = design["ssodn"]
    guide.donor_type = design["donor_type"]
    guide.homology_5 = design["homology_5"]
    guide.homology_3 = design["homology_3"]
    guide.edit_index = design["edit_index"]
    guide.strand_polarity = design["strand_polarity"]
    guide.donor_notes = f"{design['annotation']} edit; {design['base_edit']}"
    return guide
