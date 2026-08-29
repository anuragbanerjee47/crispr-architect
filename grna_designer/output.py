"""Printed table and CSV output for gRNA results."""
from __future__ import annotations

import csv
from typing import List

from .models import Guide


def rating(score: float) -> str:
    if score >= 0.75:
        return "Excellent"
    if score >= 0.55:
        return "Good"
    if score >= 0.40:
        return "Marginal"
    return "Poor"


def print_table(guides: List[Guide], top_n: int) -> None:
    top = sorted(guides, key=lambda g: g.combined_score(),
                 reverse=True)[:top_n]
    if not top:
        print("No guides found.")
        return

    print(
        f"\nTop {len(top)} gRNA candidates (ranked by efficiency + specificity):\n")
    header = (
        f"{'Pos':>6}  {'Strand':<10}  {'Guide (5->3)':<22}  {'PAM':<3}  "
        f"{'OnT':>6}  {'Spec':>6}  {'Eff':>5}  {'OT0':>4}  {'OT1':>4}  {'OT2':>4}  {'OT3':>4}  "
        f"{'Oligo_Fwd':<24}  {'Tm_Fwd':>6}  {'Oligo_Rev':<24}  {'Tm_Rev':>6}  "
        f"{'Score':>5}  {'Rating':<10}  {'Flags'}"
    )
    print(header)
    print("-" * len(header))
    for g in top:
        cs = g.combined_score()
        fwd = g.oligo_fwd or ""
        rev = g.oligo_rev or ""
        on_target = getattr(g, "on_target_score", 100.0 * g.efficiency)
        specificity = getattr(g, "specificity_score", g.specificity * 100.0)
        print(
            f"{g.pos + 1:>6}  {g.strand_name:<10}  {g.guide:<22}  {g.pam:<3}  "
            f"{on_target:>6.2f}  {specificity:>6.2f}  {g.efficiency:>5.3f}  "
            f"{g.offtarget.zero:>4}  {g.offtarget.one:>4}  {g.offtarget.two:>4}  {g.offtarget.three:>4}  "
            f"{fwd:<24}  {g.tm_fwd:>6.2f}  {rev:<24}  {g.tm_rev:>6.2f}  "
            f"{cs:>5.3f}  {rating(cs):<10}  {g.issues}"
        )


def write_csv(guides: List[Guide], top_n: int, path: str) -> None:
    top = sorted(guides, key=lambda g: g.combined_score(),
                 reverse=True)[:top_n]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "position",
            "strand",
            "guide_5to3",
            "PAM",
            "on_target_score",
            "specificity_score",
            "efficiency_score",
            "combined_score",
            "offtarget_0mm",
            "offtarget_1mm",
            "offtarget_2mm",
            "offtarget_3mm",
            "Oligo_Fwd",
            "Oligo_Rev",
            "Tm_Fwd",
            "Tm_Rev",
            "rating",
            "flags",
        ])
        for g in top:
            cs = g.combined_score()
            on_target = getattr(g, "on_target_score", 100.0 * g.efficiency)
            specificity = getattr(g, "specificity_score",
                                  g.specificity * 100.0)
            w.writerow([
                g.pos + 1,
                g.strand_name,
                g.guide,
                g.pam,
                round(float(on_target), 3),
                round(float(specificity), 3),
                round(g.efficiency, 4),
                cs,
                g.offtarget.zero,
                g.offtarget.one,
                g.offtarget.two,
                g.offtarget.three,
                g.oligo_fwd,
                g.oligo_rev,
                round(g.tm_fwd, 2),
                round(g.tm_rev, 2),
                rating(cs),
                g.issues,
            ])
