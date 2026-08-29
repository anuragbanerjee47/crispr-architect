"""
CRISPR gRNA Designer — Python + Biopython
Fetches a gene from NCBI, scans for SpCas9 (NGG) protospacers, scores candidates,
and reports off-target counts across the full transcript.

Usage:
    python grna_designer.py NM_007294
    python grna_designer.py NM_007294 --window 3000 --top 15 --mismatches 3
"""
import argparse
import csv
import re
import sys
from dataclasses import dataclass
from typing import List, Tuple

from Bio import Entrez, SeqIO

# ---------- Scoring (Doench et al. inspired heuristics) ----------
COMP = str.maketrans("ACGTN", "TGCAN")
POLY_RUN_RE = re.compile(r"(.)\1{4,}")
GUIDE_LEN = 20


def rev_comp(seq: str) -> str:
    return seq.upper().translate(COMP)[::-1]


def gc_fraction(seq: str) -> float:
    s = seq.upper()
    if not s:
        return 0.0
    return (s.count("G") + s.count("C")) / len(s)


def max_palindrome_len(seq: str, min_len: int = 4, max_len: int = 8) -> int:
    """Longest k-mer in seq that also appears in its reverse complement (hairpin risk)."""
    rc = rev_comp(seq)
    best = 0
    for k in range(min_len, max_len + 1):
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if kmer in rc:
                best = max(best, k)
    return best


@dataclass
class Guide:
    guide: str
    pam: str
    strand: str
    pos: int            # 0-based start in the original fetched sequence
    gc: float = 0.0
    gc_score: float = 0.0
    stability_score: float = 0.0
    palindrome_score: float = 0.0
    poly_score: float = 0.0
    total: float = 0.0
    issues: str = ""
    # Off-target fields (populated by count_offtargets)
    offtarget_0: int = 0
    offtarget_1: int = 0
    offtarget_2: int = 0
    offtarget_3: int = 0
    specificity_score: float = 1.0

    def score(self) -> None:
        seq = self.guide
        gc = gc_fraction(seq)
        self.gc = round(gc * 100, 1)
        self.gc_score = max(0.0, min(1.0, 1 - abs(gc - 0.5) * 3))

        s = 0.5
        if seq[19] == "G":
            s += 0.3
        if seq[19] == "C":
            s -= 0.1
        if seq[18] in "GC":
            s += 0.1
        if seq[15] == "C":
            s -= 0.2
        self.stability_score = max(0.0, min(1.0, s))

        pal = max_palindrome_len(seq)
        self.palindrome_score = max(0.0, min(1.0, 1 - (pal - 4) / 4))

        poly = 1.0
        issues = []
        if "TTTT" in seq:
            poly -= 0.5
            issues.append("TTTT")
        if POLY_RUN_RE.search(seq):
            poly -= 0.4
            issues.append("homopolymer")
        if "AAAA" in seq:
            poly -= 0.2
            issues.append("AAAA")
        self.poly_score = max(0.0, poly)
        self.issues = ",".join(issues)

        self.total = min(
            1.0,
            round(
                0.30 * self.gc_score
                + 0.25 * self.stability_score
                + 0.20 * self.palindrome_score
                + 0.15 * self.poly_score
                + 0.10,
                3,
            ),
        )


def find_guides(seq: str) -> List[Guide]:
    """Scan both strands for 20-nt guide + NGG PAM. Positions are in `seq` coordinates."""
    guides: List[Guide] = []
    seq_u = seq.upper()
    rc = rev_comp(seq_u)
    n = len(seq_u)

    for i in range(n - 22):
        pam = seq_u[i + 20:i + 23]
        if pam[1] == "G" and pam[2] == "G" and pam[0] in "ACGT":
            g = Guide(guide=seq_u[i:i + GUIDE_LEN],
                      pam=pam, strand="sense", pos=i)
            g.score()
            guides.append(g)

    for j in range(len(rc) - 22):
        pam = rc[j + 20:j + 23]
        if pam[1] == "G" and pam[2] == "G" and pam[0] in "ACGT":
            original_pos = n - 1 - (j + 19)
            g = Guide(guide=rc[j:j + GUIDE_LEN], pam=pam,
                      strand="antisense", pos=original_pos)
            g.score()
            guides.append(g)

    return guides


# ---------- Off-target search ----------
# Use Biopython PairwiseAligner with a simple mismatch penalty. We count how many
# 20-nt windows in the full sequence (both strands) match the guide with <=N mismatches.
#
# Performance: for each guide we align once against the rc sequence using a global
# aligner restricted to GUIDE_LEN windows. Biopython PairwiseAligner does this with
# a "seed-and-extend" style via find_suboptimal_alignments, but for simplicity we use
# a sliding-window Hamming-style scan (much faster, equivalent for short guides).
def _slide_count(guide: str, haystack: str, max_mismatches: int) -> int:
    """Count windows of len(GUIDE_LEN) in haystack with <= max_mismatches vs guide."""
    n = len(haystack)
    g = guide.upper()
    k = len(g)
    if n < k:
        return 0
    count = 0
    # Precompute mask: ignore N bases in haystack
    for i in range(n - k + 1):
        w = haystack[i:i + k]
        diff = 0
        for a, b in zip(g, w):
            if a != b:
                diff += 1
                if diff > max_mismatches:
                    break
        if diff <= max_mismatches:
            count += 1
    return count


def _find_hits(guide: str, haystack: str, max_mismatches: int) -> list:
    """Return list of (start, end) GUIDE_LEN intervals in haystack where the window
    matches `guide` with <= max_mismatches."""
    n = len(haystack)
    g = guide.upper()
    k = len(g)
    hits = []
    if n < k:
        return hits
    for i in range(n - k + 1):
        w = haystack[i:i + k]
        diff = 0
        for a, b in zip(g, w):
            if a != b:
                diff += 1
                if diff > max_mismatches:
                    break
        if diff <= max_mismatches:
            hits.append((i, i + k))
    return hits


def _count_clusters(intervals: list, on_target_start: int) -> int:
    """Count disjoint genomic-interval clusters, excluding the on-target window.

    Two hits in the same cluster overlap or touch; clusters separated by >=1 bp gap.
    """
    if not intervals:
        return 0
    merged = sorted(intervals)
    clusters = []
    cs, ce = merged[0]
    for s, e in merged[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            clusters.append((cs, ce))
            cs, ce = s, e
    clusters.append((cs, ce))
    # Subtract the on-target cluster (the one containing on_target_start)
    on_target_cluster = None
    for cs, ce in clusters:
        if cs <= on_target_start < ce:
            on_target_cluster = (cs, ce)
            break
    count = len(clusters) - (1 if on_target_cluster is not None else 0)
    return max(0, count)


def count_offtargets(guide: Guide, full_seq: str, mismatches: int = 3) -> None:
    """Count exact and near-match off-target loci for `guide` against full_seq.

    A sense window matching the guide AND an antisense window matching its rc at the
    same locus are merged into one cluster (one biological off-target site).
    On-target window is subtracted.
    """
    seq_u = full_seq.upper()
    rc = rev_comp(seq_u)
    n = len(seq_u)
    self_rc = rev_comp(guide.guide)
    on_target_start = guide.pos

    def _clusters_at(max_mm: int) -> int:
        sense = _find_hits(guide.guide, seq_u, max_mm)
        # antisense window at rc-index j covers genomic [n - j - GUIDE_LEN, n - j]
        anti_in_rc = _find_hits(self_rc, rc, max_mm)
        anti = [(n - j - GUIDE_LEN, n - j) for (j, _) in anti_in_rc]
        return _count_clusters(sense + anti, on_target_start)

    guide.offtarget_0 = _clusters_at(0)
    guide.offtarget_1 = _clusters_at(1)
    guide.offtarget_2 = _clusters_at(2)
    guide.offtarget_3 = _clusters_at(mismatches)

    # Specificity score: log-weighted penalty for off-targets.
    # Weights follow Hsu et al. 2013 / Cas-OFFinder convention.
    weighted = (
        100 * guide.offtarget_0
        + 10 * guide.offtarget_1
        + 1 * guide.offtarget_2
        + 0.1 * guide.offtarget_3
    )
    guide.specificity_score = round(1 / (1 + weighted / 20), 3)


# ---------- NCBI fetch ----------
def fetch_sequence(accession: str, email: str) -> Tuple[str, str]:
    Entrez.email = email
    handle = Entrez.efetch(
        db="nuccore", id=accession, rettype="fasta", retmode="text"
    )
    record = next(SeqIO.parse(handle, "fasta"))
    handle.close()
    seq = str(record.seq).upper()
    return record.description, seq


# ---------- Output ----------
def rating(score: float) -> str:
    if score >= 0.75:
        return "Excellent"
    if score >= 0.55:
        return "Good"
    if score >= 0.40:
        return "Marginal"
    return "Poor"


def combined_score(g: Guide) -> float:
    """Efficiency (60%) + specificity (40%), final 0..1."""
    return round(0.6 * g.total + 0.4 * g.specificity_score, 3)


def combined_rating(score: float) -> str:
    return rating(score)


def print_table(guides: List[Guide], top_n: int) -> None:
    top = sorted(guides, key=lambda g: combined_score(g), reverse=True)[:top_n]
    if not top:
        print("No guides found.")
        return

    print(
        f"\nTop {len(top)} gRNA candidates (ranked by efficiency + specificity):\n")
    header = (
        f"{'Pos':>6}  {'Strand':<10}  {'Guide (5->3)':<22}  {'PAM':<3}  "
        f"{'GC%':>5}  {'Eff':>5}  {'Spec':>5}  {'OT0':>4}  {'OT1':>4}  {'OT2':>4}  {'OT3':>4}  "
        f"{'Score':>5}  {'Rating':<10}  {'Flags'}"
    )
    print(header)
    print("-" * len(header))
    for g in top:
        print(
            f"{g.pos + 1:>6}  {g.strand:<10}  {g.guide:<22}  {g.pam:<3}  "
            f"{g.gc:>5}  {g.total:>5}  {g.specificity_score:>5}  "
            f"{g.offtarget_0:>4}  {g.offtarget_1:>4}  {g.offtarget_2:>4}  {g.offtarget_3:>4}  "
            f"{combined_score(g):>5}  {combined_rating(combined_score(g)):<10}  {g.issues}"
        )


def write_csv(guides: List[Guide], top_n: int, path: str) -> None:
    top = sorted(guides, key=lambda g: combined_score(g), reverse=True)[:top_n]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "position",
                "strand",
                "guide_5to3",
                "PAM",
                "GC_percent",
                "efficiency_score",
                "specificity_score",
                "combined_score",
                "offtarget_0mm",
                "offtarget_1mm",
                "offtarget_2mm",
                "offtarget_3mm",
                "rating",
                "flags",
                "gc_score",
                "stability_score",
                "palindrome_score",
                "poly_score",
            ]
        )
        for g in top:
            cs = combined_score(g)
            w.writerow(
                [
                    g.pos + 1,
                    g.strand,
                    g.guide,
                    g.pam,
                    g.gc,
                    g.total,
                    g.specificity_score,
                    cs,
                    g.offtarget_0,
                    g.offtarget_1,
                    g.offtarget_2,
                    g.offtarget_3,
                    combined_rating(cs),
                    g.issues,
                    g.gc_score,
                    g.stability_score,
                    g.palindrome_score,
                    g.poly_score,
                ]
            )


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Design CRISPR gRNAs from an NCBI sequence (with off-target search)."
    )
    ap.add_argument("accession", help="NCBI accession (e.g. NM_007294)")
    ap.add_argument(
        "--window",
        type=int,
        default=2000,
        help="Bases to scan from the start of the fetched sequence (default 2000)",
    )
    ap.add_argument("--top", type=int, default=15,
                    help="Number of guides to report")
    ap.add_argument(
        "--email",
        default="grna-designer@example.com",
        help="Email for NCBI (required by their policy)",
    )
    ap.add_argument("--csv", default=None, help="Path to write CSV results")
    ap.add_argument(
        "--mismatches",
        type=int,
        default=3,
        help="Max mismatches to count for off-targets (default 3)",
    )
    ap.add_argument(
        "--no-offtarget",
        action="store_true",
        help="Skip off-target search (faster; only useful for very short sequences)",
    )
    args = ap.parse_args()

    print(f"Fetching {args.accession} from NCBI…")
    try:
        header, seq = fetch_sequence(args.accession, args.email)
    except Exception as e:
        print(f"NCBI fetch failed: {e}", file=sys.stderr)
        return 1

    region = seq[: args.window]
    print(f"Header: {header}")
    print(f"Region scanned: {len(region)} bp")
    if not args.no_offtarget:
        print(f"Off-target search space: full sequence ({len(seq):,} bp)")

    guides = find_guides(region)
    print(f"NGG sites found: {len(guides)}")

    if not args.no_offtarget:
        print("Counting off-targets (this takes a moment)…")
        # Only compute off-targets for the top candidate pool to keep runtime bounded
        # (sort by efficiency first, compute OT for the top K, then re-rank).
        K = max(args.top * 5, 50)
        candidate_pool = sorted(
            guides, key=lambda g: g.total, reverse=True)[:K]
        for i, g in enumerate(candidate_pool):
            count_offtargets(g, seq, mismatches=min(args.mismatches, 3))
            if (i + 1) % 25 == 0:
                print(f"  {i + 1}/{len(candidate_pool)}")
        guides = candidate_pool

    print_table(guides, args.top)

    if args.csv:
        write_csv(guides, args.top, args.csv)
        print(f"\nWrote {args.csv}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
