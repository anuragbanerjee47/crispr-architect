"""Command-line interface for gRNA designer v2."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any, Iterable

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

try:  # pragma: no cover - optional dependency for progress bars
    from tqdm import tqdm
except ImportError:  # pragma: no cover
    def tqdm(iterable: Iterable[Any], *args: Any, **kwargs: Any) -> Iterable[Any]:
        del args, kwargs
        return iterable

from .cloning import annotate_guide, write_vendor_csv
from .exons import filter_guides_to_exons, parse_feature_table
from .guides import find_guides, get_enzyme
from .models import Guide, OffTargetMode
from .ncbi import EntrezClient, extract_id_list
from .offtarget import CasOffinderSearch, InTranscriptOffTarget
from .output import print_table, write_csv
from .scoring import score_guide


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="grna-designer",
        description="Design CRISPR gRNAs from an NCBI sequence. v2: rule-based scoring, "
                    "enzyme-aware PAM scanning, and optional genome-wide off-target search.",
    )
    p.add_argument("accession", nargs="?",
                   help="NCBI accession (e.g. NM_007294)")
    p.add_argument("--window", type=int, default=2000,
                   help="Bases to scan from the start of the fetched sequence (default 2000)")
    p.add_argument("--top", type=int, default=15,
                   help="Number of guides to report")
    p.add_argument("--email", default="grna-designer@example.com",
                   help="Email for NCBI (required by their policy)")
    p.add_argument("--csv", default=None, help="Path to write CSV results")
    p.add_argument("--vendor", choices=["idt", "twist", "genscript"], default=None,
                   help="Export a vendor-formatted oligo CSV for synthesis ordering")
    p.add_argument("--cas", default="spcas9",
                   help="CRISPR enzyme: spcas9, sacas9, cas12a, or spry (default: spcas9)")
    p.add_argument("--mismatches", type=int, default=3,
                   help="Max mismatches for off-target counting (default 3)")
    p.add_argument("--offtarget-mode", choices=[m.value for m in OffTargetMode],
                   default=OffTargetMode.TRANSCRIPT.value,
                   help="Off-target search space: transcript (default) or genome (via Cas-OFFinder)")
    p.add_argument("--genome-fasta", default=None,
                   help="Genome FASTA path (required for --offtarget-mode genome)")
    p.add_argument("--cas-offinder-path", default=None,
                   help="Path to cas-offinder binary (or set CAS_OFFINDER_BIN)")
    p.add_argument("--include-utr", action="store_true",
                   help="Keep guides in 5' and 3' UTRs. Default: coding-region only.")
    p.add_argument("--no-offtarget", action="store_true",
                   help="Skip off-target search entirely (faster).")
    p.add_argument("--batch", "--input-list", dest="batch", default=None,
                   help="Path to a batch file of accessions, gene symbols, or FASTA files.")
    p.add_argument("--output-dir", default="batch_output",
                   help="Output directory for per-target CSVs and batch summary reports.")
    return p


def _score_all(guides: list[Guide], *, enzyme_name: str = "spcas9") -> None:
    for g in guides:
        try:
            result = score_guide(g.guide)
            g.efficiency = result.efficiency
            g.issues = ", ".join(result.contributing_features.get("flags", []))
            oligo = annotate_guide(g, enzyme=enzyme_name)
            g.oligo_fwd = oligo.oligo_fwd
            g.oligo_rev = oligo.oligo_rev
            g.tm_fwd = oligo.tm_fwd
            g.tm_rev = oligo.tm_rev
        except ValueError:
            g.efficiency = 0.0


def _run_offtarget(
    guides: list[Guide],
    full_seq: str,
    *,
    mode: OffTargetMode,
    mismatches: int,
    genome_fasta: str | None,
    cas_offinder_path: str | None,
) -> None:
    if mode == OffTargetMode.TRANSCRIPT:
        engine = InTranscriptOffTarget(full_seq, mismatches=mismatches)
    else:
        engine = CasOffinderSearch(
            binary=cas_offinder_path,
            genome_fasta=genome_fasta,
            mismatches=mismatches,
        )
    for g in guides:
        engine.search(g)


def _safe_target_name(value: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in {
                      "_", "-"} else "_" for ch in str(value).strip())
    return cleaned or "target"


def _parse_batch_targets(batch_path: str) -> list[str]:
    path = Path(batch_path)
    if not path.exists():
        raise FileNotFoundError(f"Batch input file not found: {batch_path}")

    entries: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            entry = raw.strip()
            if entry and not entry.startswith("#"):
                entries.append(entry)
    return entries


def _resolve_target_sequence(target: str, *, email: str) -> str:
    value = str(target).strip()
    if not value:
        raise ValueError("Batch target entry is empty.")

    candidate_path = Path(value)
    if candidate_path.exists() and candidate_path.suffix.lower() in {".fa", ".fasta", ".fna"}:
        record = next(SeqIO.parse(str(candidate_path), "fasta"))
        return str(record.seq).upper()

    if value.startswith(">"):
        record = next(SeqIO.parse(value.splitlines(), "fasta"))
        return str(record.seq).upper()

    if any(ch.isdigit() for ch in value):
        with EntrezClient(email) as client:
            record: SeqRecord = client.fetch_nuccore(value)
            return str(record.seq).upper()

    from Bio import Entrez
    Entrez.email = email
    search = Entrez.esearch(
        db="nuccore",
        term=f"{value}[Gene Name] OR {value}[Symbol] OR {value}[All Fields]",
        retmax=1,
    )
    ids = extract_id_list(search)
    if not ids:
        raise ValueError(
            f"No NCBI record found for gene symbol or search term: {value}")
    handle = Entrez.efetch(
        db="nuccore", id=ids[0], rettype="fasta", retmode="text")
    try:
        record = next(SeqIO.parse(handle, "fasta"))
    finally:
        handle.close()
    return str(record.seq).upper()


def _process_single_batch_target(
    target: str,
    *,
    email: str,
    enzyme_name: str,
    gc_min: float,
    gc_max: float,
    window: int,
    top_n: int,
    output_dir: Path,
) -> dict:
    seq = _resolve_target_sequence(target, email=email)
    enzyme = get_enzyme(enzyme_name)
    guides = find_guides(seq[:window], enzyme=enzyme)
    _score_all(guides, enzyme_name=enzyme_name)

    filtered: list[Guide] = []
    for g in guides:
        if not g.guide:
            continue
        gc = (sum(base in {"G", "C"}
              for base in g.guide) / len(g.guide)) * 100.0
        if gc_min <= gc <= gc_max:
            filtered.append(g)

    ranked = sorted(filtered or guides,
                    key=lambda g: g.efficiency, reverse=True)[:top_n]
    csv_path = output_dir / f"{_safe_target_name(target)}.csv"
    write_csv(ranked, len(ranked), str(csv_path))

    return {
        "target": target,
        "sequence_length": len(seq),
        "pam_sites_found": len(guides),
        "selected_candidates": len(ranked),
        "best_efficiency": round(max((g.efficiency for g in ranked), default=0.0), 4),
        "best_guide": ranked[0].guide if ranked else "",
    }


def process_batch_file(
    batch_path: str,
    *,
    email: str,
    enzyme_name: str = "spcas9",
    gc_min: float = 40.0,
    gc_max: float = 60.0,
    window: int = 2000,
    top_n: int = 10,
    output_dir: str = "batch_output",
) -> list[dict[str, Any]]:
    targets = _parse_batch_targets(batch_path)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    batch_rows: list[dict[str, Any]] = []
    for target in tqdm(targets, desc="Processing batch targets"):
        batch_rows.append(
            _process_single_batch_target(
                target,
                email=email,
                enzyme_name=enzyme_name,
                gc_min=gc_min,
                gc_max=gc_max,
                window=window,
                top_n=top_n,
                output_dir=out_dir,
            )
        )

    summary_path = out_dir / "batch_summary.csv"
    fieldnames = [
        "target",
        "sequence_length",
        "pam_sites_found",
        "selected_candidates",
        "best_efficiency",
        "best_guide",
    ]
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(batch_rows)

    print(f"Completed batch processing for {len(batch_rows)} targets.")
    print(f"Summary report written to: {summary_path}")
    return batch_rows


def _single_target_main(args: argparse.Namespace) -> int:
    if not args.accession:
        raise ValueError(
            "An accession or FASTA target is required unless --batch is used.")

    enzyme = get_enzyme(args.cas)
    print(f"Fetching {args.accession} from NCBI...")
    try:
        with EntrezClient(args.email) as client:
            record: SeqRecord = client.fetch_nuccore(args.accession)
            seq = str(record.seq).upper()
            header = record.description
            print(f"Header: {header}")
            print(
                f"Region scanned: {min(args.window, len(seq))} bp of {len(seq):,} total")
            print(
                f"Selected Cas enzyme: {enzyme['name']} ({enzyme['pam']} PAM, {enzyme['pam_orientation']} orientation)")
            if not args.no_offtarget and args.offtarget_mode == OffTargetMode.GENOME:
                print("Off-target search: genome mode")
            elif not args.no_offtarget:
                print(f"Off-target search: transcript ({len(seq):,} bp)")

            exon_intervals = []
            if not args.include_utr:
                try:
                    ft_text = client.fetch_feature_table(args.accession)
                    exon_intervals = parse_feature_table(ft_text)
                    print(f"Exon feature rows parsed: {len(exon_intervals)}")
                except (AttributeError, OSError, RuntimeError, TimeoutError, TypeError, ValueError) as exc:  # pragma: no cover - network fallback path
                    print(
                        f"Warning: could not fetch feature table ({exc}); exon filter disabled.", file=sys.stderr)
                    exon_intervals = []
    except (AttributeError, OSError, RuntimeError, TimeoutError, TypeError, ValueError) as exc:
        print(f"NCBI fetch failed: {exc}", file=sys.stderr)
        return 1

    region = seq[: args.window]
    guides = find_guides(
        region, exon_filter=exon_intervals or None, enzyme=enzyme)
    print(f"PAM sites found in window: {len(guides)}")

    if exon_intervals:
        before = len(guides)
        guides = filter_guides_to_exons(
            guides, exon_intervals, include_utr=args.include_utr)
        print(
            f"After exon filter: {len(guides)} (dropped {before - len(guides)})")

    _score_all(guides, enzyme_name=args.cas)

    if not args.no_offtarget:
        print("Counting off-targets (this takes a moment)...")
        K = max(args.top * 5, 50)
        pool = sorted(guides, key=lambda g: g.efficiency, reverse=True)[:K]
        try:
            _run_offtarget(
                pool,
                seq,
                mode=OffTargetMode(args.offtarget_mode),
                mismatches=args.mismatches,
                genome_fasta=args.genome_fasta,
                cas_offinder_path=args.cas_offinder_path,
            )
        except FileNotFoundError as exc:
            print(f"\nOff-target backend error: {exc}", file=sys.stderr)
            return 2
        except RuntimeError as exc:
            print(f"\nOff-target backend failed: {exc}", file=sys.stderr)
            return 2
        guides = pool

    print_table(guides, args.top)

    if args.csv:
        write_csv(guides, args.top, args.csv)
        print(f"\nWrote {args.csv}")

    if args.vendor:
        vendor_path = args.csv if args.csv else f"vendor_{args.vendor}.csv"
        write_vendor_csv(guides, vendor_path, vendor=args.vendor)
        print(f"Wrote vendor CSV ({args.vendor}) to {vendor_path}")

    return 0


def main(argv: list[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)
    if args.batch:
        summary = process_batch_file(
            args.batch,
            email=args.email,
            enzyme_name=args.cas,
            gc_min=40.0,
            gc_max=60.0,
            window=args.window,
            top_n=args.top,
            output_dir=args.output_dir,
        )
        print(f"Batch summary rows: {len(summary)}")
        return 0
    return _single_target_main(args)


if __name__ == "__main__":
    sys.exit(main())
