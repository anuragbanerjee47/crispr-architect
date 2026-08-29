"""Cloning oligo generation, Tm estimation, vector library support, and vendor export helpers."""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, List, Sequence, Union

from .models import Guide
from .sequtil import rev_comp

NN_TABLE = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -20.4),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}

VECTOR_LIBRARY: Dict[str, Dict[str, object]] = {
    "lenticrisprv2": {
        "vector_name": "lentiCRISPRv2",
        "enzyme": "BsmBI / Esp3I",
        "u6_requires_g": True,
        "forward_overhang_if_g": "CACC",
        "forward_overhang_no_g": "CACCG",
        "reverse_overhang": "AAAC",
        "notes": "U6 expression plasmid; G is required at the 5' end for efficient transcription.",
    },
    "px459": {
        "vector_name": "pSpCas9(BB)-2A-Puro (PX459)",
        "enzyme": "BbsI / BpiI",
        "u6_requires_g": True,
        "forward_overhang_if_g": "CACC",
        "forward_overhang_no_g": "CACCG",
        "reverse_overhang": "AAAC",
        "notes": "Standard mammalian U6-pSpCas9 expression backbone.",
    },
    "px330": {
        "vector_name": "pX330-U6-Chimeric_BB-CBh-hSpCas9",
        "enzyme": "BbsI",
        "u6_requires_g": True,
        "forward_overhang_if_g": "CACC",
        "forward_overhang_no_g": "CACCG",
        "reverse_overhang": "AAAC",
        "notes": "Humanized SpCas9 plasmid with U6 promoter.",
    },
    "pcasguideaav": {
        "vector_name": "pCas-Guide-AAV",
        "enzyme": "BamHI / EcoRI or Golden Gate",
        "u6_requires_g": True,
        "forward_overhang_if_g": "GATC",
        "forward_overhang_no_g": "GATCG",
        "reverse_overhang": "AATT",
        "notes": "AAV vector with a BsaI/Golden Gate-style cloning cassette.",
    },
    "custom": {
        "vector_name": "Custom / User Overhang",
        "enzyme": "Custom",
        "u6_requires_g": False,
        "forward_overhang_if_g": "",
        "forward_overhang_no_g": "",
        "reverse_overhang": "",
        "notes": "User-defined sticky ends for custom assembly.",
    },
}


def gc_percent(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = sum(base in {"G", "C"} for base in seq)
    return (gc / len(seq)) * 100.0


def nearest_neighbor_tm(seq: str, dnac1: float = 50e-9, dnac2: float = 50e-9, salt_molar: float = 1.0) -> float:
    """Approximate short oligo Tm using the nearest-neighbor model."""
    s = seq.upper()
    if not s:
        return 0.0

    dh = 0.2
    ds = -5.7
    for left, right in zip(s, s[1:]):
        h, s_val = NN_TABLE.get(f"{left}{right}", (-7.0, -20.0))
        dh += h
        ds += s_val

    for terminal in (s[0], s[-1]):
        if terminal in {"A", "T"}:
            dh += 2.2
            ds += 6.1
        else:
            dh += 0.1
            ds += -2.8

    ct = max(dnac1, dnac2) if dnac1 == dnac2 else min(dnac1, dnac2)
    ct = max(ct, 1e-12)
    tm = (1000.0 * dh) / (ds + 1.987 * math.log(ct / 4.0)) - 273.15
    if not math.isfinite(tm):
        return 0.0
    return float(tm)


def calculate_self_dimer_risk(seq: str) -> float:
    """Estimate temperature-independent self-dimer/hairpin risk from reverse-complement overlap."""
    s = str(seq).upper()
    if not s:
        return 0.0
    rc = rev_comp(s)
    best = 0
    for k in range(6, min(len(s), 16) + 1):
        for i in range(len(s) - k + 1):
            window = s[i:i + k]
            if window in rc:
                best = max(best, k)
    if best == 0:
        return 0.0
    return round(best / len(s), 3)


def annealing_temperature(tm_primer: float, tm_product: float) -> float:
    """Return the exact annealing temperature used in the wet-lab protocol."""
    return 0.3 * float(tm_primer) + 0.7 * float(tm_product) - 14.9


def generate_annealing_protocol(guide_seq: str, vector_name: str = "lentiCRISPRv2") -> Dict[str, object]:
    """Generate a 1-page annealing + ligation protocol for a guide oligo pair."""
    design = design_cloning_oligos(
        guide_seq, enzyme="spcas9", vector_name=vector_name)
    tm_primer = max(design["tm_fwd"], design["tm_rev"])
    tm_product = 0.5 * (design["tm_fwd"] + design["tm_rev"])
    ta = annealing_temperature(tm_primer, tm_product)
    vector_ratio = 1.0 if vector_name.lower() in {"custom"} else 3.0
    protocol = (
        "1. Phosphorylate each oligo pair with T4 PNK in a 10-20 µL reaction at 37°C for 30 min.\n"
        "2. Denature the phosphorylated oligos at 95°C for 2 min, then cool to room temperature over 25-30 min.\n"
        "3. Anneal to the vector insert using a vector:insert molar ratio of approximately 1:3 for U6 vectors.\n"
        "4. Ligate at 25°C for 1 hour using T4 DNA ligase and proceed to transformation.\n"
        "5. Confirm by Sanger sequencing across the U6-to-guide junction and Cas9 insertion site."
    )
    return {
        "guide": guide_seq,
        "vector": vector_name,
        "tm_primer": round(tm_primer, 2),
        "tm_product": round(tm_product, 2),
        "annealing_temperature_c": round(ta, 2),
        "vector_insert_ratio": round(vector_ratio, 2),
        "protocol": protocol,
    }


def _normalize_enzyme_name(name: Union[str, Guide, None]) -> str:
    if isinstance(name, Guide):
        name = name.guide
    if name is None:
        return "spcas9"
    return str(name).strip().lower().replace("-", "").replace("_", "")


def get_vector_profile(vector_name: str = "lentiCRISPRv2", custom_fwd_overhang: str | None = None, custom_rev_overhang: str | None = None) -> Dict[str, object]:
    """Return the backbone configuration for the selected cloning vector."""
    key = str(vector_name).strip().lower().replace(
        "-", "").replace("_", "").replace(" ", "")
    if key in {"", "custom", "useroverhang"}:
        key = "custom"
    if key not in VECTOR_LIBRARY:
        raise ValueError(f"Unsupported vector: {vector_name!r}")
    profile = dict(VECTOR_LIBRARY[key])
    if key == "custom":
        profile["forward_overhang_if_g"] = (
            custom_fwd_overhang or profile["forward_overhang_if_g"] or "").upper()
        profile["forward_overhang_no_g"] = (
            custom_fwd_overhang or profile["forward_overhang_no_g"] or "").upper()
        profile["reverse_overhang"] = (
            custom_rev_overhang or profile["reverse_overhang"] or "").upper()
    return profile


def design_cloning_oligos(
    guide_5to3: str,
    enzyme: str = "spcas9",
    vector_name: str = "lentiCRISPRv2",
    custom_fwd_overhang: str | None = None,
    custom_rev_overhang: str | None = None,
) -> dict:
    """Build vector-compatible forward/reverse oligo sequences and Tm values.

    Backward compatibility: the classic SpCas9 design uses ``CACCG`` on the
    forward primer and ``AAAC`` plus the reverse-complement on the reverse primer,
    while preserving the newer vector-aware override path for other backbones.
    """
    guide = str(guide_5to3).upper()
    key = _normalize_enzyme_name(enzyme)
    vector_profile = get_vector_profile(
        vector_name, custom_fwd_overhang, custom_rev_overhang)

    if key in {"spcas9", "sacas9"}:
        fwd_overhang = str(vector_profile.get(
            "forward_overhang_if_g", "CACC")).upper()
        no_g_overhang = str(vector_profile.get(
            "forward_overhang_no_g", "CACCG")).upper()
        rev_overhang = str(vector_profile.get(
            "reverse_overhang", "AAAC")).upper()

        if guide.startswith("G"):
            fwd = fwd_overhang + guide
            rev = rev_overhang + rev_comp(guide)
        else:
            fwd = no_g_overhang + guide
            rev = rev_overhang + rev_comp(guide) + "C"
    elif key in {"cas12a", "cpf1"}:
        fwd = "TTT" + guide + "A"
        rev = "TTT" + rev_comp(guide) + "A"
    else:
        raise ValueError(f"Unsupported enzyme for oligo design: {enzyme!r}")

    return {
        "guide": guide,
        "vector": str(vector_profile.get("vector_name", vector_name)),
        "fwd_overhang": fwd[: len(fwd) - len(guide)],
        "rev_overhang": rev[: len(rev) - len(rev_comp(guide))],
        "oligo_fwd": fwd,
        "oligo_rev": rev,
        "tm_fwd": nearest_neighbor_tm(fwd),
        "tm_rev": nearest_neighbor_tm(rev),
        "gc_fwd": gc_percent(fwd),
        "gc_rev": gc_percent(rev),
    }


def annotate_guide(guide: Guide, enzyme: str = "spcas9", vector_name: str = "lentiCRISPRv2", custom_fwd_overhang: str | None = None, custom_rev_overhang: str | None = None) -> Guide:
    design = design_cloning_oligos(
        guide.guide,
        enzyme=enzyme,
        vector_name=vector_name,
        custom_fwd_overhang=custom_fwd_overhang,
        custom_rev_overhang=custom_rev_overhang,
    )
    guide.oligo_fwd = design["oligo_fwd"]
    guide.oligo_rev = design["oligo_rev"]
    guide.tm_fwd = design["tm_fwd"]
    guide.tm_rev = design["tm_rev"]
    return guide


def _pairwise_well_positions(num_pairs: int) -> List[str]:
    wells: List[str] = []
    for pair_index in range(num_pairs):
        row_idx = pair_index // 6
        col_idx = pair_index % 6
        letter = chr(ord("A") + row_idx)
        if letter > "H":
            letter = "H"
        wells.append(f"{letter}{(col_idx * 2) + 1:02d}")
        wells.append(f"{letter}{(col_idx * 2) + 2:02d}")
    return wells


def build_vendor_plate_rows(guides: Sequence[Guide], vendor: str = "idt") -> List[dict]:
    """Build 96-well compatible upload rows with paired forward/reverse oligos.

    The row payload keeps the modern plate mapping keys while also preserving the
    original download contract expected by the legacy vendor CSV tests.
    """
    vendor = vendor.lower()
    rows: List[dict] = []
    wells = _pairwise_well_positions(len(guides))
    for idx, g in enumerate(guides):
        pair_wells = wells[idx * 2:(idx * 2) + 2]
        validation_tm = round(float(getattr(g, "tm_fwd", 0.0) or 0.0), 2)
        validation_rev_tm = round(float(getattr(g, "tm_rev", 0.0) or 0.0), 2)
        validation_fwd = getattr(g, "flanking_fwd_primer", "") or ""
        validation_rev = getattr(g, "flanking_rev_primer", "") or ""
        amplicon_length = getattr(g, "amplicon_length", 0) or 0
        ssodn = getattr(g, "ssodn_sequence", "") or ""
        donor_type = getattr(g, "donor_type", "") or ""
        homology_5 = getattr(g, "homology_5", 0) or 0
        homology_3 = getattr(g, "homology_3", 0) or 0
        edit_index = getattr(g, "edit_index", 0) or 0
        strand_polarity = getattr(
            g, "strand_polarity", "non_target") or "non_target"
        if vendor == "idt":
            row_fwd = {
                "Well Position": pair_wells[0],
                "Sequence Name": f"gRNA_{idx + 1}_F",
                "Sequence 5'->3'": g.oligo_fwd or "",
                "Scale": "25nm",
                "Purification": "STD",
                "Plate": "P1",
                "Tube ID": f"T{idx * 2 + 1:03d}",
                "Sequence": g.oligo_fwd or "",
                "Notes": g.guide,
                "Validation Fwd Primer": validation_fwd,
                "Validation Rev Primer": validation_rev,
                "Validation Tm Fwd": validation_tm,
                "Validation Tm Rev": validation_rev_tm,
                "Amplicon Length": amplicon_length,
                "ssODN Sequence": ssodn,
                "HDR Donor Type": donor_type,
                "Homology 5'": homology_5,
                "Homology 3'": homology_3,
                "Edit Index": edit_index,
                "Strand Polarity": strand_polarity,
            }
            row_rev = {
                "Well Position": pair_wells[1],
                "Sequence Name": f"gRNA_{idx + 1}_R",
                "Sequence 5'->3'": g.oligo_rev or "",
                "Scale": "25nm",
                "Purification": "STD",
                "Plate": "P1",
                "Tube ID": f"T{idx * 2 + 2:03d}",
                "Sequence": g.oligo_rev or "",
                "Notes": g.guide,
                "Validation Fwd Primer": validation_fwd,
                "Validation Rev Primer": validation_rev,
                "Validation Tm Fwd": validation_tm,
                "Validation Tm Rev": validation_rev_tm,
                "Amplicon Length": amplicon_length,
                "ssODN Sequence": ssodn,
                "HDR Donor Type": donor_type,
                "Homology 5'": homology_5,
                "Homology 3'": homology_3,
                "Edit Index": edit_index,
                "Strand Polarity": strand_polarity,
            }
            rows.extend([row_fwd, row_rev])
        elif vendor == "twist":
            row_fwd = {
                "Well": pair_wells[0],
                "Gene/Oligo Name": f"gRNA_{idx + 1}_F",
                "Sequence": g.oligo_fwd or "",
                "Yield": "25nm",
                "Application": "CRISPR guide synthesis",
                "Name": f"gRNA_{idx + 1}_F",
                "Scale": "25nm",
                "Purification": "desalted",
                "Tube": f"T{idx * 2 + 1:03d}",
                "Note": g.guide,
                "Validation Fwd Primer": validation_fwd,
                "Validation Rev Primer": validation_rev,
                "Validation Tm Fwd": validation_tm,
                "Validation Tm Rev": validation_rev_tm,
                "Amplicon Length": amplicon_length,
            }
            row_rev = {
                "Well": pair_wells[1],
                "Gene/Oligo Name": f"gRNA_{idx + 1}_R",
                "Sequence": g.oligo_rev or "",
                "Yield": "25nm",
                "Application": "CRISPR guide synthesis",
                "Name": f"gRNA_{idx + 1}_R",
                "Scale": "25nm",
                "Purification": "desalted",
                "Tube": f"T{idx * 2 + 2:03d}",
                "Note": g.guide,
                "Validation Fwd Primer": validation_fwd,
                "Validation Rev Primer": validation_rev,
                "Validation Tm Fwd": validation_tm,
                "Validation Tm Rev": validation_rev_tm,
                "Amplicon Length": amplicon_length,
            }
            rows.extend([row_fwd, row_rev])
        elif vendor == "genscript":
            row_fwd = {
                "Plate ID": "P1",
                "Well": pair_wells[0],
                "Primer Name": f"gRNA_{idx + 1}_F",
                "Sequence": g.oligo_fwd or "",
                "Synthesis Scale": "25nm",
                "Order ID": f"GS{idx * 2 + 1:03d}",
                "Tube": f"T{idx * 2 + 1:03d}",
                "Purification": "standard",
                "Notes": g.guide,
                "Validation Fwd Primer": validation_fwd,
                "Validation Rev Primer": validation_rev,
                "Validation Tm Fwd": validation_tm,
                "Validation Tm Rev": validation_rev_tm,
                "Amplicon Length": amplicon_length,
            }
            row_rev = {
                "Plate ID": "P1",
                "Well": pair_wells[1],
                "Primer Name": f"gRNA_{idx + 1}_R",
                "Sequence": g.oligo_rev or "",
                "Synthesis Scale": "25nm",
                "Order ID": f"GS{idx * 2 + 2:03d}",
                "Tube": f"T{idx * 2 + 2:03d}",
                "Purification": "standard",
                "Notes": g.guide,
                "Validation Fwd Primer": validation_fwd,
                "Validation Rev Primer": validation_rev,
                "Validation Tm Fwd": validation_tm,
                "Validation Tm Rev": validation_rev_tm,
                "Amplicon Length": amplicon_length,
            }
            rows.extend([row_fwd, row_rev])
        else:
            raise ValueError(f"Unsupported vendor format: {vendor!r}")
    return rows


def _vendor_rows(guides: Sequence[Guide], vendor: str) -> List[dict]:
    return build_vendor_plate_rows(guides, vendor)


def write_vendor_csv(guides: Sequence[Guide], path: Union[str, Path], vendor: str = "idt") -> None:
    rows = _vendor_rows(guides, vendor)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    vendor = vendor.lower()
    if vendor == "idt":
        fieldnames = ["Plate", "Tube ID", "Sequence",
                      "Scale", "Purification", "Notes",
                      "Validation Fwd Primer", "Validation Rev Primer",
                      "Validation Tm Fwd", "Validation Tm Rev", "Amplicon Length",
                      "ssODN Sequence", "HDR Donor Type", "Homology 5'", "Homology 3'", "Edit Index", "Strand Polarity"]
    elif vendor == "twist":
        fieldnames = ["Name", "Sequence", "Scale",
                      "Purification", "Tube", "Note",
                      "Validation Fwd Primer", "Validation Rev Primer",
                      "Validation Tm Fwd", "Validation Tm Rev", "Amplicon Length",
                      "ssODN Sequence", "HDR Donor Type", "Homology 5'", "Homology 3'", "Edit Index", "Strand Polarity"]
    elif vendor == "genscript":
        fieldnames = ["Order ID", "Tube", "Sequence", "Purification", "Notes",
                      "Validation Fwd Primer", "Validation Rev Primer",
                      "Validation Tm Fwd", "Validation Tm Rev", "Amplicon Length",
                      "ssODN Sequence", "HDR Donor Type", "Homology 5'", "Homology 3'", "Edit Index", "Strand Polarity"]
    else:
        raise ValueError(f"Unsupported vendor format: {vendor!r}")

    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def generate_genbank_map(guide_seq: str, vector_name: str = "lentiCRISPRv2", insert_seq: str | None = None) -> str:
    """Generate a minimal GenBank record for the cloned guide construct."""
    vector = get_vector_profile(vector_name)
    insert = str(insert_seq or guide_seq).upper()
    seq_len = max(len(insert) + 2000, 5000)
    source_seq = (
        "A" * 300 + insert + "C" * 100 + "G" * 150 +
        "T" * 120 + "A" * 80 + "C" * 90 + "G" * 130
    )
    source_seq = source_seq[:seq_len]
    lines = [
        f"LOCUS       {vector['vector_name'].replace(' ', '_').upper()}_GUIDE  {len(source_seq)} bp    DNA     circular 01-AUG-2026",
        "DEFINITION  CRISPR guide expression construct generated by gRNA Designer.",
        "ACCESSION   .",
        "VERSION     .",
        "KEYWORDS    CRISPR; guide; cloning; vector.",
        "SOURCE      synthetic construct",
        "  ORGANISM  synthetic construct",
        "FEATURES             Location/Qualifiers",
        "     source          1..{0}".format(len(source_seq)),
        "                     /organism=\"synthetic construct\"",
        "     promoter        1..300",
        "                     /label=\"U6 promoter\"",
        "     misc_feature    301..{0}".format(300 + len(insert)),
        "                     /label=\"guide insert\"",
        "                     /note=\"Guide sequence: {0}\"".format(guide_seq),
        "     CDS             301..{0}".format(300 + len(insert) + 100),
        "                     /label=\"Cas9\"",
        "                     /note=\"CRISPR nuclease insert\"",
        "     misc_feature    {0}..{1}".format(
            300 + len(insert) + 101, len(source_seq)),
        "                     /label=\"selection marker\"",
        "ORIGIN",
    ]
    for i in range(0, len(source_seq), 12):
        chunk = source_seq[i:i + 12]
        start = i + 1
        formatted = " ".join(f"{base}" for base in chunk)
        lines.append(f"{start:9d} {formatted}")
    lines.append("//")
    return "\n".join(lines)
