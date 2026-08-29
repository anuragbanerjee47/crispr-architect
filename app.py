"""Streamlit dashboard for CRISPR guide design and oligo export."""
from __future__ import annotations
from grna_designer.scoring import score_guide
from grna_designer.ncbi import EntrezClient, extract_id_list
from grna_designer.guides import find_guides, get_enzyme
from grna_designer.dual_guide import (
    calculate_deletion_product,
    design_dual_guide_deletion,
    design_paired_nickase,
)
from grna_designer.hdr import attach_hdr_donor, design_ssodn
from grna_designer.cloning import (
    annotate_guide,
    build_vendor_plate_rows,
    calculate_self_dimer_risk,
    generate_annealing_protocol,
    generate_genbank_map,
    gc_percent,
    get_vector_profile,
)
from Bio import Entrez, SeqIO

import csv
import io
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from grna_designer.security import (
    sanitize_accession,
    sanitize_gene_symbol,
    sanitize_sequence,
)

PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

try:  # pragma: no cover - exercised when Streamlit is installed.
    import pandas as pd
except Exception:  # pragma: no cover
    pd = None

try:  # pragma: no cover - exercised when Streamlit is installed.
    import plotly.graph_objects as go
except Exception:  # pragma: no cover
    go = None

try:  # pragma: no cover - exercised when Streamlit is installed.
    import streamlit as st
except Exception:  # pragma: no cover
    st = None


ENZYME_OPTIONS = {
    "SpCas9": "spcas9",
    "SaCas9": "sacas9",
    "Cas12a/Cpf1": "cas12a",
    "SpRY": "spry",
}

VECTOR_OPTIONS = {
    "lentiCRISPRv2": "lentiCRISPRv2",
    "pSpCas9(BB)-2A-Puro (PX459)": "pSpCas9(BB)-2A-Puro (PX459)",
    "pX330-U6-Chimeric_BB-CBh-hSpCas9": "pX330-U6-Chimeric_BB-CBh-hSpCas9",
    "pCas-Guide-AAV": "pCas-Guide-AAV",
    "Custom / User Overhang": "Custom / User Overhang",
}


def _normalise_sequence(raw: str) -> str:
    return sanitize_sequence(raw, field_name="sequence")


def _sequence_from_accession(accession: str, email: str) -> str:
    accession = sanitize_accession(accession)
    with EntrezClient(email) as client:
        record = client.fetch_nuccore(accession)
        return sanitize_sequence(str(record.seq).upper(), field_name="accession sequence")


def _sequence_from_gene_symbol(symbol: str, email: str) -> str:
    symbol = sanitize_gene_symbol(symbol)

    Entrez.email = email
    search = Entrez.esearch(
        db="nuccore",
        term=f"{symbol}[Gene Name] OR {symbol}[Symbol] OR {symbol}[All Fields]",
        retmax=5,
    )
    ids = extract_id_list(search)
    if not ids:
        raise ValueError(f"No NCBI nucleotide records found for '{symbol}'.")

    handle = Entrez.efetch(
        db="nuccore", id=ids[0], rettype="fasta", retmode="text")
    try:
        record = next(SeqIO.parse(handle, "fasta"))
    finally:
        handle.close()
    return sanitize_sequence(str(record.seq).upper(), field_name="gene-search sequence")


def _build_candidate_rows(
    sequence: str,
    enzyme_name: str,
    gc_min: float,
    gc_max: float,
    mismatch_tolerance: int,
    top_n: int,
    vector_name: str = "lentiCRISPRv2",
    custom_fwd_overhang: str | None = None,
    custom_rev_overhang: str | None = None,
) -> List[Dict[str, object]]:
    enzyme_key = ENZYME_OPTIONS.get(enzyme_name, enzyme_name.lower())
    enzyme = get_enzyme(enzyme_key)
    guides = find_guides(sequence, enzyme=enzyme)

    rows: List[Dict[str, object]] = []
    for guide in guides:
        score = score_guide(guide.guide)
        guide.efficiency = score.efficiency
        guide.on_target_score = score.on_target_score
        guide.issues = ", ".join(score.contributing_features.get("flags", []))
        annotate_guide(
            guide,
            enzyme=enzyme_key,
            vector_name=vector_name,
            custom_fwd_overhang=custom_fwd_overhang,
            custom_rev_overhang=custom_rev_overhang,
        )
        from grna_designer.primers import attach_validation_primers
        attach_validation_primers(guide, sequence)

        gc_value = gc_percent(guide.guide)
        if gc_value < gc_min or gc_value > gc_max:
            continue

        specificity = max(
            0.0,
            min(
                1.0,
                1.0 - (mismatch_tolerance * 0.06) +
                (max(0.0, gc_max - gc_value) / max(1.0, gc_max)),
            ),
        )
        guide.specificity = specificity

        tm_value = float((guide.tm_fwd + guide.tm_rev) /
                         2.0) if guide.tm_fwd and guide.tm_rev else 0.0
        self_dimer = calculate_self_dimer_risk(
            guide.oligo_fwd + guide.oligo_rev)
        rows.append(
            {
                "Rank": len(rows) + 1,
                "Cut Position": guide.pos + 1,
                "Strand": "+" if guide.strand.value == "sense" else "-",
                "Guide Sequence (5'->3')": guide.guide,
                "PAM": guide.pam,
                "GC %": round(gc_value, 2),
                "On-Target Score": round(float(getattr(guide, "on_target_score", score.on_target_score)), 2),
                "Specificity Score": round(guide.specificity, 3),
                "Oligo Fwd": guide.oligo_fwd,
                "Oligo Rev": guide.oligo_rev,
                "Oligo Tm": round(tm_value, 2),
                "Self-Dimer Risk": round(self_dimer, 3),
                "Validation Fwd Primer": getattr(guide, "flanking_fwd_primer", ""),
                "Validation Rev Primer": getattr(guide, "flanking_rev_primer", ""),
                "Validation Tm Fwd": round(float(getattr(guide, "tm_fwd", 0.0) or 0.0), 2),
                "Validation Tm Rev": round(float(getattr(guide, "tm_rev", 0.0) or 0.0), 2),
                "Amplicon Length": int(getattr(guide, "amplicon_length", 0) or 0),
                "ssODN Sequence": getattr(guide, "ssodn_sequence", "") or "",
                "HDR Donor Type": getattr(guide, "donor_type", "") or "",
                "Dual Guide Pair": "",
            }
        )

    rows.sort(key=lambda item: float(item["Specificity Score"]), reverse=True)
    for idx, row in enumerate(rows, start=1):
        row["Rank"] = idx
    return rows[:top_n]


def _vendor_csv_bytes(rows: Sequence[Dict[str, object]], vendor: str = "idt") -> bytes:
    vendor = vendor.lower()
    buffer = io.StringIO()
    writer = csv.writer(buffer)

    if vendor == "idt":
        writer.writerow(["Plate", "Tube ID", "Sequence",
                        "Scale", "Purification", "Notes",
                         "Validation Fwd Primer", "Validation Rev Primer",
                         "Validation Tm Fwd", "Validation Tm Rev", "Amplicon Length",
                         "ssODN Sequence", "HDR Donor Type", "Dual Guide Pair"])
        for idx, row in enumerate(rows, start=1):
            writer.writerow([
                "P1",
                f"T{idx:03d}",
                row["Oligo Fwd"],
                "25nm",
                "STD",
                row["Guide Sequence (5'->3')"],
                row.get("Validation Fwd Primer", ""),
                row.get("Validation Rev Primer", ""),
                row.get("Validation Tm Fwd", ""),
                row.get("Validation Tm Rev", ""),
                row.get("Amplicon Length", ""),
                row.get("ssODN Sequence", ""),
                row.get("HDR Donor Type", ""),
                row.get("Dual Guide Pair", ""),
            ])
    elif vendor == "twist":
        writer.writerow(["Name", "Sequence", "Scale",
                        "Purification", "Tube", "Note",
                         "Validation Fwd Primer", "Validation Rev Primer",
                         "Validation Tm Fwd", "Validation Tm Rev", "Amplicon Length",
                         "ssODN Sequence", "HDR Donor Type", "Dual Guide Pair"])
        for idx, row in enumerate(rows, start=1):
            writer.writerow([
                f"gRNA_{idx:03d}",
                row["Oligo Fwd"],
                "25nm",
                "desalted",
                f"T{idx:03d}",
                row["Guide Sequence (5'->3')"],
                row.get("Validation Fwd Primer", ""),
                row.get("Validation Rev Primer", ""),
                row.get("Validation Tm Fwd", ""),
                row.get("Validation Tm Rev", ""),
                row.get("Amplicon Length", ""),
                row.get("ssODN Sequence", ""),
                row.get("HDR Donor Type", ""),
                row.get("Dual Guide Pair", ""),
            ])
    elif vendor == "genscript":
        writer.writerow(
            ["Order ID", "Tube", "Sequence", "Purification", "Notes",
             "Validation Fwd Primer", "Validation Rev Primer",
             "Validation Tm Fwd", "Validation Tm Rev", "Amplicon Length",
             "ssODN Sequence", "HDR Donor Type", "Dual Guide Pair"])
        for idx, row in enumerate(rows, start=1):
            writer.writerow([
                f"GS{idx:03d}",
                f"T{idx:03d}",
                row["Oligo Fwd"],
                "standard",
                row["Guide Sequence (5'->3')"],
                row.get("Validation Fwd Primer", ""),
                row.get("Validation Rev Primer", ""),
                row.get("Validation Tm Fwd", ""),
                row.get("Validation Tm Rev", ""),
                row.get("Amplicon Length", ""),
                row.get("ssODN Sequence", ""),
                row.get("HDR Donor Type", ""),
                row.get("Dual Guide Pair", ""),
            ])
    else:
        raise ValueError(f"Unsupported vendor format: {vendor!r}")
    return buffer.getvalue().encode("utf-8")


def _sequence_chart(rows: Sequence[Dict[str, object]], sequence_length: int):
    if go is None:
        return None
    fig = go.Figure()
    for row in rows:
        strand_color = "#0f766e" if row["Strand"] == "+" else "#be185d"
        x = int(row["Cut Position"])
        y = 0.5
        fig.add_trace(
            go.Scatter(
                x=[x, x],
                y=[0.25, 0.75],
                mode="lines",
                line={"color": strand_color, "width": 2},
                hovertemplate=f"Guide: {row['Guide Sequence (5\'->3\')']}<extra></extra>",
                showlegend=False,
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[x],
                y=[y],
                mode="markers",
                marker={"color": strand_color, "size": 10},
                hovertemplate=(
                    f"Guide: {row['Guide Sequence (5\'->3\')']}<br>"
                    f"PAM: {row['PAM']}<br>Strand: {row['Strand']}<extra></extra>"
                ),
                showlegend=False,
            )
        )
    fig.update_layout(
        title="Guide binding map",
        xaxis_title="Sequence position (bp)",
        yaxis_title="",
        yaxis={"showticklabels": False, "range": [0, 1]},
        template="plotly_white",
        width=1100,
        height=260,
    )
    fig.update_xaxes(range=[0, max(sequence_length, 1)])
    return fig


def _specificity_badge_html(score: float) -> str:
    score = float(score)
    if score >= 80.0:
        color = "#22c55e"
        label = "High"
    elif score >= 50.0:
        color = "#f59e0b"
        label = "Moderate"
    else:
        color = "#ef4444"
        label = "Low"
    return (
        f'<span class="spec-badge" style="background:{color}; border-color:{color};">'
        f"{label} • {score:.1f}</span>"
    )


def _render_results(results: Dict[str, object]):
    if st is None:
        raise RuntimeError(
            "Run the app with Streamlit installed to view the dashboard.")

    sequence = str(results.get("sequence", ""))
    rows = list(results.get("rows", []))
    summary = dict(results.get("summary", {}))
    vendor_exports = dict(results.get("vendor_exports", {}))
    genbank_map = str(results.get("genbank_map", ""))
    annealing_rows = list(results.get("annealing_rows", []))
    vendor_plate_rows = list(results.get("vendor_plate_rows", []))
    vector_name = str(results.get("vector_profile", {}).get("name", "Custom"))

    top_specificity = float(summary.get("top_specificity", max(
        (float(row["Specificity Score"]) for row in rows), default=0.0)))
    average_gc = round(
        sum(float(row["GC %"]) for row in rows) / len(rows), 2
    ) if rows else 0.0

    st.subheader("Operational summary")
    metric_cols = st.columns(4)
    metric_cols[0].metric("Total guides identified", int(
        summary.get("pam_sites_found", len(rows))))
    metric_cols[1].metric("Highest on-target score",
                          f"{max((float(row['On-Target Score']) for row in rows), default=0.0):.1f}")
    metric_cols[2].metric("Selected vector", vector_name)
    metric_cols[3].metric("Target GC%", f"{average_gc:.1f}%")

    if hasattr(st, "markdown") and rows:
        best_guide = rows[0]
        top_score = float(best_guide.get("Specificity Score", 0.0))
        top_badge = _specificity_badge_html(top_score)
        st.markdown(
            """
            <div class="insight-strip">
              <div class="insight-card">
                <div class="insight-label">best guide</div>
                <div class="insight-value">%s</div>
              </div>
              <div class="insight-card">
                <div class="insight-label">specificity</div>
                <div class="insight-value">%s</div>
              </div>
              <div class="insight-card">
                <div class="insight-label">vector fit</div>
                <div class="insight-value">%s</div>
              </div>
            </div>
            """ % (
                best_guide.get("Guide Sequence (5'->3')", "-"),
                top_badge,
                vector_name,
            ),
            unsafe_allow_html=True,
        )

    if rows:
        table_df = pd.DataFrame(rows) if pd is not None else None
        if table_df is not None:
            table_df = table_df.copy()
            table_df["Specificity"] = table_df["Specificity Score"].map(
                _specificity_badge_html)
            table_df["On-Target"] = table_df["On-Target Score"].map(
                lambda value: float(value))
            table_df["GC %"] = table_df["GC %"].map(lambda value: float(value))
            show_cols = [
                "Rank",
                "Guide Sequence (5'->3')",
                "On-Target Score",
                "Specificity",
                "GC %",
                "PAM",
                "Cut Position",
            ]
            display_df = table_df[show_cols].copy()
            display_df["Specificity"] = display_df["Specificity"].map(
                lambda value: value)
            if hasattr(st, "dataframe"):
                st.dataframe(display_df, use_container_width=True)
        elif hasattr(st, "table"):
            st.table(rows)

    supports_tabs = hasattr(st, "tabs") and callable(getattr(st, "tabs"))
    if supports_tabs:
        tabs = st.tabs([
            "🎯 Guide Ranking",
            "� Validation Primers",
            "🧬 HDR Donor Templates",
            "🔀 Dual Guides / Nickase",
            "🧱 Vector Cloning & Overhangs",
            "📦 Vendor Order Sheets",
        ])

        with tabs[0]:
            if rows:
                frame = pd.DataFrame(rows) if pd is not None else None
                if frame is not None:
                    ranked = frame.copy()
                    if "Specificity Score" in ranked.columns:
                        ranked["Specificity Score"] = ranked["Specificity Score"].astype(float)
                    ranked["Specificity Badge"] = ranked["Specificity Score"].map(_specificity_badge_html)
                    ranked = ranked[[
                        "Rank",
                        "Guide Sequence (5'->3')",
                        "On-Target Score",
                        "Specificity Score",
                        "Specificity Badge",
                        "GC %",
                        "PAM",
                        "Cut Position",
                    ]]
                    st.dataframe(ranked, use_container_width=True, hide_index=True)
                else:
                    st.subheader("Ranked candidates")
                    if hasattr(st, "table"):
                        st.table(rows)

                if hasattr(st, "download_button"):
                    csv_bytes = io.StringIO()
                    if pd is not None:
                        pd.DataFrame(rows).to_csv(csv_bytes, index=False)
                    else:
                        fieldnames = list(rows[0].keys()) if rows else []
                        writer = csv.DictWriter(csv_bytes, fieldnames=fieldnames)
                        writer.writeheader()
                        writer.writerows(rows)
                    st.download_button(
                        label="Download ranked candidates (.csv)",
                        data=csv_bytes.getvalue().encode("utf-8"),
                        file_name="grna_candidates.csv",
                        mime="text/csv",
                    )
            else:
                if hasattr(st, "warning"):
                    st.warning("No candidates passed the selected GC filter for the current settings.")

        with tabs[1]:
            if rows:
                primer_rows = [{
                    "Guide": row.get("Guide Sequence (5'->3')", ""),
                    "Validation Fwd Primer": row.get("Validation Fwd Primer", ""),
                    "Validation Rev Primer": row.get("Validation Rev Primer", ""),
                    "Tm Fwd": row.get("Validation Tm Fwd", ""),
                    "Tm Rev": row.get("Validation Tm Rev", ""),
                    "Amplicon Length": row.get("Amplicon Length", ""),
                } for row in rows]
                if hasattr(st, "dataframe"):
                    st.dataframe(pd.DataFrame(primer_rows), use_container_width=True)
                if hasattr(st, "write"):
                    for row in rows[:1]:
                        st.write(f"Primary guide: {row.get('Guide Sequence (5\'->3\')', '')}")
                        st.write(f"Forward primer: {row.get('Validation Fwd Primer', '')}")
                        st.write(f"Reverse primer: {row.get('Validation Rev Primer', '')}")
            else:
                st.info("Validation primers appear after guide candidates are generated.")

        with tabs[2]:
            hdr_design = results.get("hdr_design", {})
            if hdr_design:
                st.subheader("HDR donor template")
                if hasattr(st, "json"):
                    st.json(hdr_design)
                if hasattr(st, "code"):
                    st.code(hdr_design.get("ssodn", ""), language="text")
            else:
                st.info("HDR donor templates are generated when the HDR workflow is selected.")

        with tabs[3]:
            dual_design = results.get("dual_design", {})
            nickase_design = results.get("nickase_design", {})
            if dual_design or nickase_design:
                st.subheader("Dual-guide deletion and nickase design")
                if hasattr(st, "json"):
                    st.json({"dual": dual_design, "nickase": nickase_design})
            else:
                st.info("Dual-guide deletion and paired-nickase layouts appear when the dual-guide workflow is selected.")

        with tabs[4]:
            if rows:
                first = rows[0]
                guide_seq = str(first["Guide Sequence (5'->3')"])
                st.subheader("Oligo architecture")
                if hasattr(st, "write"):
                    st.write(f"Primary guide: {guide_seq}")
                    st.write(f"Oligo forward: {first.get('Oligo Fwd', '')}")
                    st.write(f"Oligo reverse: {first.get('Oligo Rev', '')}")
                    st.write(f"Annealing Tm: {first.get('Oligo Tm', 'n/a')}")
                    st.write(f"Self-dimer risk: {first.get('Self-Dimer Risk', 0)}")
                if annealing_rows and hasattr(st, "dataframe"):
                    st.dataframe(pd.DataFrame(annealing_rows), use_container_width=True)
                if hasattr(st, "json"):
                    st.json(results.get("vector_profile", {}))
            else:
                if hasattr(st, "info"):
                    st.info("Cloning details appear once viable guide candidates are available.")

        with tabs[5]:
            st.subheader("Plate export and vendor ordering")
            vendor_choice = st.selectbox("Vendor format", ["idt", "twist", "genscript"], index=0, key="vendor_order_format")
            if vendor_plate_rows and hasattr(st, "dataframe"):
                st.dataframe(pd.DataFrame(vendor_plate_rows), use_container_width=True)
            if vendor_exports.get(vendor_choice) and hasattr(st, "download_button"):
                st.download_button(
                    label=f"Download {vendor_choice.upper()} order sheet",
                    data=vendor_exports[vendor_choice],
                    file_name=f"{vendor_choice}_plate_order.csv",
                    mime="text/csv",
                )
            if hasattr(st, "code"):
                preview = vendor_exports.get(vendor_choice, b"")
                if preview:
                    st.code(preview.decode("utf-8", errors="ignore")[:1500], language="csv")
            elif hasattr(st, "info"):
                st.info("No vendor plate export is available for the current result set.")
    else:
        if rows:
            st.subheader("Candidate results")
            if hasattr(st, "table"):
                st.table(rows)
            if hasattr(st, "write"):
                st.write(
                    "Session state preserved; results remain available across reruns.")
        elif hasattr(st, "warning"):
            st.warning(
                "No candidates passed the selected GC filter for the current settings.")


def main() -> None:
    if st is None:
        raise RuntimeError(
            "This app requires Streamlit to be installed. Install with: pip install streamlit plotly pandas")

    st.set_page_config(
        page_title="CRISPR gRNA & Oligo Architect",
        page_icon="🧬",
        layout="wide",
    )
    if hasattr(st, "markdown"):
        st.markdown(
            """
            <style>
            :root {
                --bg: #06151f;
                --bg-2: #0d2235;
                --panel: rgba(12, 25, 36, 0.92);
                --panel-strong: rgba(10, 20, 30, 0.96);
                --panel-soft: rgba(17, 35, 49, 0.84);
                --border: rgba(148, 163, 184, 0.18);
                --text: #eaf7ff;
                --muted: #b7ccd9;
                --accent: #5eead4;
                --accent-strong: #2dd4bf;
                --accent-soft: rgba(94, 234, 212, 0.14);
                --green: #34d399;
                --amber: #fbbf24;
                --rose: #f472b6;
                --shadow: 0 22px 50px rgba(5, 15, 22, 0.45);
            }
            .stApp {
                background:
                    radial-gradient(circle at top left, rgba(45, 212, 191, 0.2), transparent 26%),
                    radial-gradient(circle at bottom right, rgba(59,130,246,0.13), transparent 24%),
                    linear-gradient(180deg, var(--bg) 0%, var(--bg-2) 100%);
                color: var(--text);
            }
            .block-container {
                padding: 1.5rem 1.1rem 3rem;
            }
            h1, h2, h3, h4 {
                color: #f1f9ff;
                letter-spacing: -0.03em;
            }
            [data-testid="stSidebar"] {
                background: rgba(7, 18, 26, 0.94);
                border-right: 1px solid var(--border);
            }
            [data-testid="stSidebar"] .css-1d391kg {
                padding-top: 1rem;
            }
            .app-header {
                display: flex;
                flex-direction: column;
                gap: 0.65rem;
                margin: 0.25rem 0 1.1rem;
            }
            .header-badge {
                display: inline-flex;
                align-items: center;
                width: fit-content;
                padding: 0.38rem 0.8rem;
                background: linear-gradient(135deg, rgba(94,234,212,0.12), rgba(79,70,229,0.12));
                border: 1px solid rgba(94,234,212,0.35);
                border-radius: 999px;
                color: #9ef4e8;
                font-size: 0.68rem;
                font-weight: 800;
                letter-spacing: 0.14em;
                text-transform: uppercase;
            }
            .header-subtitle {
                color: var(--muted);
                font-size: 0.96rem;
                line-height: 1.5;
                max-width: 860px;
            }
            .hero-panel {
                background: linear-gradient(135deg, rgba(12,25,36,0.96), rgba(17,35,49,0.9));
                border: 1px solid var(--border);
                border-radius: 18px;
                padding: 1.1rem 1.2rem;
                margin-bottom: 1.1rem;
                box-shadow: var(--shadow);
            }
            .hero-row {
                display: grid;
                grid-template-columns: repeat(3, minmax(0, 1fr));
                gap: 0.85rem;
                margin-top: 0.9rem;
            }
            .hero-pill {
                display: flex;
                align-items: center;
                gap: 0.5rem;
                padding: 0.7rem 0.8rem;
                background: linear-gradient(180deg, rgba(17,34,48,0.9), rgba(11,24,34,0.88));
                border: 1px solid rgba(148,163,184,0.18);
                border-radius: 12px;
                color: var(--text);
                font-weight: 600;
            }
            .hero-pill .dot {
                width: 10px;
                height: 10px;
                border-radius: 50%;
                background: linear-gradient(135deg, var(--accent), var(--green));
                box-shadow: 0 0 0 6px rgba(94,234,212,0.12);
            }
            .feature-band {
                display: grid;
                grid-template-columns: repeat(4, minmax(0, 1fr));
                gap: 0.75rem;
                margin-top: 1rem;
            }
            .feature-chip {
                background: rgba(13, 29, 40, 0.72);
                border: 1px solid var(--border);
                border-radius: 10px;
                padding: 0.75rem 0.8rem;
                color: var(--muted);
                font-size: 0.8rem;
                line-height: 1.45;
            }
            .feature-chip strong {
                display: block;
                color: var(--text);
                margin-bottom: 0.2rem;
                font-size: 0.82rem;
            }
            .metric-card {
                background: linear-gradient(180deg, rgba(15, 31, 44, 0.96), rgba(10, 20, 30, 0.94));
                border: 1px solid var(--border);
                border-radius: 14px;
                padding: 1rem 1.1rem;
                box-shadow: var(--shadow);
            }
            .spec-badge {
                display: inline-flex;
                align-items: center;
                justify-content: center;
                min-width: 108px;
                padding: 0.45rem 0.7rem;
                border-radius: 999px;
                font-size: 0.72rem;
                font-weight: 800;
                letter-spacing: 0.02em;
                color: #041822;
                border: 1px solid transparent;
                box-shadow: 0 10px 18px rgba(15, 23, 42, 0.28);
            }
            .insight-strip {
                display: grid;
                grid-template-columns: repeat(3, minmax(0, 1fr));
                gap: 0.75rem;
                margin: 1rem 0 1.5rem;
            }
            .insight-card {
                background: rgba(12, 25, 35, 0.92);
                border: 1px solid var(--border);
                border-radius: 14px;
                padding: 0.9rem 1rem;
                box-shadow: var(--shadow);
            }
            .insight-label {
                color: var(--muted);
                font-size: 0.72rem;
                letter-spacing: 0.09em;
                text-transform: uppercase;
                margin-bottom: 0.4rem;
            }
            .insight-value {
                color: #ecfeff;
                font-size: 1.08rem;
                font-weight: 700;
                word-break: break-word;
            }
            .panel {
                background: rgba(16, 39, 56, 0.88);
                border: 1px solid var(--border);
                border-radius: 12px;
                padding: 1rem;
                box-shadow: var(--shadow);
            }
            .stButton > button {
                border-radius: 12px;
                background: linear-gradient(135deg, #24c7b6, #7dd3fc);
                color: #022436;
                border: none;
                font-weight: 800;
                box-shadow: 0 16px 28px rgba(36,199,182,0.22);
                transition: transform 0.15s ease, box-shadow 0.15s ease;
            }
            .stButton > button:hover {
                transform: translateY(-1px);
                box-shadow: 0 20px 32px rgba(36,199,182,0.24);
            }
            .stDownloadButton > button {
                border-radius: 10px;
                background: rgba(15, 118, 110, 0.20);
                border: 1px solid rgba(94,234,212,0.3);
                color: #d8fffb;
            }
            div[data-testid="stDataFrame"] {
                border-radius: 12px;
                overflow: hidden;
                border: 1px solid var(--border);
            }
            .stTabs [role="tablist"] {
                gap: 0.4rem;
            }
            .stTabs [role="tab"] {
                border-radius: 10px 10px 0 0;
                padding: 0.55rem 0.8rem;
            }
            .stTabs [role="tab"][aria-selected="true"] {
                background: rgba(94,234,212,0.12);
                color: #bffcf6;
                border: 1px solid rgba(94,234,212,0.22);
            }
            .stTextInput > div > div > input,
            .stTextArea textarea,
            .stSelectbox > div > div,
            .stNumberInput > div > div > input,
            .stSlider > div[data-baseweb="slider"] {
                background: rgba(11, 22, 30, 0.9);
                border: 1px solid rgba(148,163,184,0.16);
                border-radius: 10px;
                color: var(--text);
            }
            .stSidebar .stButton > button {
                border-radius: 10px;
                background: linear-gradient(135deg, rgba(24,199,182,0.95), rgba(125,211,252,0.9));
                color: #032537;
            }
            .stAlert {
                border-radius: 12px;
                border: 1px solid rgba(148,163,184,0.18);
            }
            </style>
            """,
            unsafe_allow_html=True,
        )
    if hasattr(st, "markdown"):
        st.markdown(
            '''
            <div class="app-header">
                <div class="header-badge">Bioinformatics workflow</div>
                <div class="header-subtitle">Design high-confidence CRISPR guides, validate GC and specificity, and export cloning-ready oligos and vendor-ready order sheets in one workflow.</div>
            </div>
            <div class="hero-panel">
                <div style="font-size:0.74rem; letter-spacing:0.12em; text-transform:uppercase; color:#9fe9d9; font-weight:700; margin-bottom:0.7rem;">Precision design workspace</div>
                <h3 style="margin:0; font-size:1.12rem;">A premium workflow for efficient guide discovery, wet-lab readiness, and ordering automation.</h3>
                <div class="hero-row">
                    <div class="hero-pill"><span class="dot"></span>Cas enzyme-aware targeting</div>
                    <div class="hero-pill"><span class="dot"></span>GC and structural QC</div>
                    <div class="hero-pill"><span class="dot"></span>Vendor export automation</div>
                </div>
                <div class="feature-band">
                    <div class="feature-chip"><strong>Design</strong>Rank high-confidence guides across broad genomic intervals.</div>
                    <div class="feature-chip"><strong>Validate</strong>Assess specificity, GC balance, and structural practicality.</div>
                    <div class="feature-chip"><strong>Clone</strong>Prepare oligo architecture for real sequencing and screening pipelines.</div>
                    <div class="feature-chip"><strong>Order</strong>Export plate-ready vendor CSVs in a single click.</div>
                </div>
            </div>
            ''',
            unsafe_allow_html=True,
        )
    st.title("CRISPR gRNA & Oligo Architect")

    if "results" not in st.session_state:
        st.session_state["results"] = None
    if "example_gene" not in st.session_state:
        st.session_state["example_gene"] = "HBB"
    if "example_mode" not in st.session_state:
        st.session_state["example_mode"] = "Gene Symbol / Search"

    st.sidebar.header("Design Controls")
    workflow_mode = st.sidebar.selectbox(
        "Workflow Mode",
        ["Standard Single gRNA", "HDR Knock-in", "Dual-Guide Deletion / Nickase"],
        index=0,
    )
    if hasattr(st.sidebar, "caption"):
        st.sidebar.caption(
            "Tune guide discovery, filtering, and export settings for a clean lab-ready design workflow.")
    example_clicked = False
    if hasattr(st.sidebar, "button"):
        example_clicked = st.sidebar.button(
            "Load Example (HBB Gene)", use_container_width=True)
    if example_clicked:
        st.session_state["example_gene"] = "HBB"
        st.session_state["example_mode"] = "Gene Symbol / Search"
        if hasattr(st.sidebar, "info"):
            st.sidebar.info("Example loaded: HBB gene search mode enabled.")

    mode_options = ["NCBI Accession ID",
                    "Gene Symbol / Search", "Paste FASTA Sequence"]
    input_mode = st.sidebar.selectbox(
        "Input Mode",
        mode_options,
        index=mode_options.index(st.session_state.get(
            "example_mode", "NCBI Accession ID")) if st.session_state.get("example_mode") in mode_options else 0,
    )
    enzyme_name = st.sidebar.selectbox(
        "Cas Enzyme",
        ["SpCas9", "SaCas9", "Cas12a/Cpf1", "SpRY"],
        index=0,
    )
    vector_name = st.sidebar.selectbox(
        "Vector Library",
        ["lentiCRISPRv2", "pSpCas9(BB)-2A-Puro (PX459)", "pX330-U6-Chimeric_BB-CBh-hSpCas9",
         "pCas-Guide-AAV", "Custom / User Overhang"],
        index=0,
    )
    custom_fwd_overhang = ""
    custom_rev_overhang = ""
    if vector_name == "Custom / User Overhang":
        custom_fwd_overhang = st.sidebar.text_input(
            "Custom 5' forward overhang", value="CACC")
        custom_rev_overhang = st.sidebar.text_input(
            "Custom reverse overhang", value="AAAC")
    gc_min = st.sidebar.slider("GC content minimum (%)", 20, 80, 40)
    gc_max = st.sidebar.slider("GC content maximum (%)", 20, 80, 60)
    if gc_min > gc_max:
        gc_min, gc_max = gc_max, gc_min
    mismatch_tolerance = st.sidebar.slider("Max mismatch tolerance", 0, 10, 3)
    top_n = st.sidebar.slider("Top candidates", 1, 50, 10)
    email = st.sidebar.text_input(
        "NCBI Entrez email", value="grna-designer@example.com")

    input_value = ""
    if input_mode == "NCBI Accession ID":
        input_value = st.sidebar.text_input("Accession", value="NM_007294")
    elif input_mode == "Gene Symbol / Search":
        input_value = st.sidebar.text_input(
            "Gene symbol / search",
            value=st.session_state.get("example_gene", "HBB"),
        )
    else:
        input_value = st.text_area(
            "Paste FASTA sequence",
            height=180,
            value="",
        )

    number_input = getattr(st.sidebar, "number_input", None)
    if callable(number_input):
        dual_start = number_input("Dual-guide deletion start", min_value=0, value=120, step=1)
        dual_end = number_input("Dual-guide deletion end", min_value=0, value=180, step=1)
    else:
        dual_start = 120
        dual_end = 180
    if dual_end < dual_start:
        dual_end = dual_start

    if hasattr(st, "info") and not input_value:
        st.info("Add a sequence, accession, or gene symbol to begin. The workflow will automatically rank candidates and prepare cloning exports.")

    analyze = st.button("Analyze sequence", use_container_width=True)
    if analyze:
        try:
            if input_mode == "NCBI Accession ID":
                sequence = _sequence_from_accession(input_value, email)
            elif input_mode == "Gene Symbol / Search":
                sequence = _sequence_from_gene_symbol(input_value, email)
            else:
                sequence = _normalise_sequence(input_value)
            if not sequence:
                raise ValueError("No sequence was supplied for analysis.")
            rows = _build_candidate_rows(
                sequence,
                enzyme_name=enzyme_name,
                gc_min=float(gc_min),
                gc_max=float(gc_max),
                mismatch_tolerance=int(mismatch_tolerance),
                top_n=int(top_n),
                vector_name=vector_name,
                custom_fwd_overhang=custom_fwd_overhang or None,
                custom_rev_overhang=custom_rev_overhang or None,
            )
            vector_profile = get_vector_profile(
                vector_name,
                custom_fwd_overhang=custom_fwd_overhang or None,
                custom_rev_overhang=custom_rev_overhang or None,
            )
            summary = {
                "sequence_length": len(sequence),
                "pam_sites_found": len(rows),
                "top_specificity": max((float(row["Specificity Score"]) for row in rows), default=0.0),
            }
            hdr_design = {}
            if workflow_mode == "HDR Knock-in" and rows:
                first = rows[0]
                hdr_design = design_ssodn(sequence, type("GuideStub", (), {"guide": first["Guide Sequence (5'->3')"], "pam": first["PAM"], "strand": "+", "pos": int(first["Cut Position"]) - 1})(), homology_5=40, homology_3=40)
                for row in rows:
                    row["ssODN Sequence"] = hdr_design.get("ssodn", "")
                    row["HDR Donor Type"] = hdr_design.get("donor_type", "asymmetric")
            dual_design = {}
            nickase_design = {}
            if workflow_mode == "Dual-Guide Deletion / Nickase":
                dual_design = design_dual_guide_deletion(sequence, int(dual_start), int(dual_end))
                nickase_design = design_paired_nickase(sequence, int(dual_start), int(dual_end))
                if rows:
                    rows[0]["Dual Guide Pair"] = f"{dual_design['left_guide']['guide']} / {dual_design['right_guide']['guide']}"
            vendor_exports = {
                vendor: _vendor_csv_bytes(rows, vendor=vendor)
                for vendor in ("idt", "twist", "genscript")
            }
            annealing_rows = []
            if rows:
                for row in rows:
                    guide_seq = str(row["Guide Sequence (5'->3')"])
                    annealing_rows.append(generate_annealing_protocol(
                        guide_seq, vector_name=vector_name))
            genbank_map = generate_genbank_map(
                rows[0]["Guide Sequence (5'->3')"] if rows else sequence[:20],
                vector_name=vector_name,
                insert_seq=(rows[0]["Guide Sequence (5'->3')"]
                            if rows else sequence[:20]),
            )
            plate_guides = [
                type(
                    "G",
                    (),
                    {
                        "guide": row["Guide Sequence (5'->3')"],
                        "oligo_fwd": row["Oligo Fwd"],
                        "oligo_rev": row["Oligo Rev"],
                        "ssodn_sequence": row.get("ssODN Sequence", ""),
                        "donor_type": row.get("HDR Donor Type", ""),
                    },
                )()
                for row in rows
            ]
            st.session_state["results"] = {
                "sequence": sequence,
                "rows": rows,
                "summary": summary,
                "vendor_exports": vendor_exports,
                "vector_profile": vector_profile,
                "annealing_rows": annealing_rows,
                "genbank_map": genbank_map,
                "vendor_plate_rows": build_vendor_plate_rows(plate_guides, vendor="idt") if rows else [],
                "workflow_mode": workflow_mode,
                "hdr_design": hdr_design,
                "dual_design": dual_design,
                "nickase_design": nickase_design,
            }
            if hasattr(st, "success"):
                st.success("Design run complete.")
        except ValueError as exc:
            st.session_state["results"] = None
            st.error(f"Input validation error: {exc}")
        except Exception as exc:  # pragma: no cover - interactive app path
            st.session_state["results"] = None
            st.error(f"Analysis failed: {exc}")

    if "results" in st.session_state and st.session_state["results"] is not None:
        _render_results(st.session_state["results"])


if __name__ == "__main__":
    main()
