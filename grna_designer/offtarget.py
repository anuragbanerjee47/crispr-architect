"""Off-target search: transcript-mode (Hamming scan) and Cas-OFFinder subprocess."""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Protocol, Tuple

from .models import Guide, OffTargetCounts
from .sequtil import rev_comp

GUIDE_LEN = 20

# Hsu et al. 2013 / Cas-OFFinder weighting convention.
_OT_WEIGHTS = {0: 100.0, 1: 10.0, 2: 1.0, 3: 0.1}


def compute_specificity(counts: OffTargetCounts) -> float:
    """Hsu/Cas-OFFinder weighted penalty -> [0, 1] specificity score."""
    weighted = (
        _OT_WEIGHTS[0] * counts.zero
        + _OT_WEIGHTS[1] * counts.one
        + _OT_WEIGHTS[2] * counts.two
        + _OT_WEIGHTS[3] * counts.three
    )
    return round(1.0 / (1.0 + weighted / 20.0), 3)


def merge_clusters(
    sense: List[Tuple[int, int]],
    antisense: List[Tuple[int, int]],
    on_target_start: int,
) -> int:
    """Promoted from v1 `_count_clusters`.

    Merge all genomic intervals (sense + antisense at the same locus collapse
    into one cluster), then subtract the on-target cluster. Returns the count
    of off-target clusters.
    """
    intervals = list(sense) + list(antisense)
    if not intervals:
        return 0
    merged = sorted(intervals)
    clusters: List[Tuple[int, int]] = []
    cs, ce = merged[0]
    for s, e in merged[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            clusters.append((cs, ce))
            cs, ce = s, e
    clusters.append((cs, ce))

    on_target_cluster = None
    for cs, ce in clusters:
        if cs <= on_target_start < ce:
            on_target_cluster = (cs, ce)
            break
    return max(0, len(clusters) - (1 if on_target_cluster is not None else 0))


def _find_hits(guide: str, haystack: str, max_mismatches: int) -> List[Tuple[int, int]]:
    """Return GUIDE_LEN intervals in `haystack` matching `guide` with <= max_mm."""
    g = guide.upper()
    n = len(haystack)
    k = len(g)
    hits: List[Tuple[int, int]] = []
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


class OffTargetSearch(Protocol):
    """Pluggable off-target backend. Implementations populate guide.offtarget
    and guide.specificity in place."""

    def search(self, guide: Guide) -> OffTargetCounts: ...


class InTranscriptOffTarget:
    """v1-style sliding-window Hamming scan, restricted to one transcript."""

    def __init__(self, full_seq: str, *, mismatches: int = 3) -> None:
        self.full_seq = full_seq.upper()
        self.rc = rev_comp(self.full_seq)
        self.n = len(self.full_seq)
        self.mismatches = min(mismatches, 3)

    def _clusters_at(self, guide: Guide, max_mm: int) -> int:
        sense = _find_hits(guide.guide, self.full_seq, max_mm)
        anti_in_rc = _find_hits(rev_comp(guide.guide), self.rc, max_mm)
        anti = [(self.n - j - GUIDE_LEN, self.n - j) for (j, _) in anti_in_rc]
        return merge_clusters(sense, anti, guide.pos)

    def search(self, guide: Guide) -> OffTargetCounts:
        counts = OffTargetCounts(
            zero=self._clusters_at(guide, 0),
            one=self._clusters_at(guide, 1),
            two=self._clusters_at(guide, 2),
            three=self._clusters_at(guide, self.mismatches),
        )
        guide.offtarget = counts
        guide.specificity = compute_specificity(counts)
        return counts


def _resolve_binary(cli_value: Optional[str]) -> str:
    """Find the Cas-OFFinder binary. Resolution order:
    1) --cas-offinder-path CLI value (must exist)
    2) CAS_OFFINDER_BIN env var
    3) shutil.which lookup
    """
    if cli_value:
        p = Path(cli_value)
        if not p.exists():
            raise FileNotFoundError(f"--cas-offinder-path not found: {p}")
        return str(p)
    env = os.environ.get("CAS_OFFINDER_BIN")
    if env and Path(env).exists():
        return env
    found = shutil.which("cas-offinder") or shutil.which("cas-offinder.exe")
    if not found:
        raise FileNotFoundError(
            "cas-offinder not on PATH. Pass --cas-offinder-path, set "
            "CAS_OFFINDER_BIN, or install from "
            "https://github.com/snugel/cas-offinder/releases"
        )
    return found


def _resolve_genome(cli_value: Optional[str]) -> str:
    if cli_value:
        p = Path(cli_value)
        if not p.exists():
            raise FileNotFoundError(f"--genome-fasta not found: {p}")
        return str(p)
    env = os.environ.get("GRNA_GENOME_FASTA")
    if env and Path(env).exists():
        return env
    raise FileNotFoundError(
        "Genome FASTA not provided. Pass --genome-fasta PATH or set "
        "GRNA_GENOME_FASTA. (hg38 is ~3 GB; there is no default.)"
    )


def _parse_cas_offinder_output(path: Path) -> List[Tuple[int, str, int, int]]:
    """Parse a Cas-OFFinder output TSV.

    Expected columns: chromosome, position, guide_seq, target_seq, strand,
    mismatches, bulges. We return (genomic_position, strand, mismatches, bulges)
    for each row. NOTE: column order is verified against the v2.4 release; if
    Cas-OFFinder's format changes, update this parser.
    """
    rows: List[Tuple[int, str, int, int]] = []
    if not path.exists():
        return rows
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 7:
            # Fall back to whitespace split for older releases.
            parts = line.split()
        if len(parts) < 7:
            continue
        try:
            position = int(parts[1])
        except ValueError:
            continue
        strand = parts[4]
        try:
            mismatches = int(parts[5])
            bulges = int(parts[6])
        except ValueError:
            continue
        rows.append((position, strand, mismatches, bulges))
    return rows


class CasOffinderSearch:
    """Genome-wide off-target search via the Cas-OFFinder binary.

    One subprocess per guide. The binary is fast (seconds for hg38 at 3 mm),
    so per-guide invocation keeps the parsing simple. If you need to batch,
    accumulate guides into one TSV and split the output by guide prefix.
    """

    def __init__(
        self,
        *,
        binary: Optional[str] = None,
        genome_fasta: Optional[str] = None,
        mismatches: int = 3,
        work_dir: Optional[Path] = None,
    ) -> None:
        self.binary = _resolve_binary(binary)
        self.genome = _resolve_genome(genome_fasta)
        self.mismatches = min(mismatches, 3)
        self.work_dir = Path(work_dir) if work_dir else Path(tempfile.gettempdir())

    def _run_one(self, guide_23nt: str, in_tsv: Path, out_tsv: Path) -> None:
        in_tsv.write_text(
            f"{guide_23nt}\t{self.mismatches}\t0\t0\n",
            encoding="utf-8",
        )
        try:
            subprocess.run(
                [self.binary, self.genome, str(in_tsv), str(out_tsv)],
                check=True,
                capture_output=True,
            )
        except FileNotFoundError as e:
            raise FileNotFoundError(f"cas-offinder binary not executable: {self.binary}") from e
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"cas-offinder failed (exit {e.returncode}): {e.stderr.decode(errors='replace')}"
            ) from e

    def _rows_to_clusters(self, rows: List[Tuple[int, str, int, int]],
                          guide: Guide) -> OffTargetCounts:
        """Group rows by (mismatch_bucket) -> disjoint genomic clusters.

        For each mismatch bucket (0..self.mismatches), antisense rows are
        mirrored: their genomic position is converted to a 0..genome-length
        window so merge_clusters can collapse a sense+antisense site.
        """
        bucket_rows = {0: [], 1: [], 2: [], 3: []}
        for position, strand, mismatches, bulges in rows:
            # max(mismatches, bulge_size) for bucket assignment
            bucket = min(3, max(mismatches, bulges))
            if bucket not in bucket_rows:
                continue
            if strand == "-":
                # Antisense: the "start" of the window is `position - GUIDE_LEN`.
                # For cluster merging, represent the interval as
                # (position - GUIDE_LEN, position).
                start = max(0, position - GUIDE_LEN)
                bucket_rows[bucket].append((start, position, True))
            else:
                bucket_rows[bucket].append((position, position + GUIDE_LEN, False))

        on_target_start = guide.pos  # mRNA coordinate; will not match genome pos.
        # The mRNA coordinate won't match the genomic coordinate, so the
        # on-target cluster will not be subtracted in genome mode unless
        # `guide.pos` is in genomic coordinates. The caller can subtract
        # manually if it has both views.
        counts = OffTargetCounts()
        for b, items in bucket_rows.items():
            sense = [(s, e) for s, e, anti in items if not anti]
            anti = [(s, e) for s, e, anti in items if anti]
            counts_at_b = merge_clusters(sense, anti, on_target_start)
            if b == 0:
                counts.zero = counts_at_b
            elif b == 1:
                counts.one = counts_at_b
            elif b == 2:
                counts.two = counts_at_b
            elif b == 3:
                counts.three = counts_at_b
        return counts

    def search(self, guide: Guide) -> OffTargetCounts:
        guide_23nt = guide.guide.upper() + "NGG"  # place PAM; verify against release format
        with tempfile.NamedTemporaryFile(
            delete=False, suffix=".tsv", dir=str(self.work_dir)
        ) as inf:
            in_tsv = Path(inf.name)
        with tempfile.NamedTemporaryFile(
            delete=False, suffix=".tsv.out", dir=str(self.work_dir)
        ) as outf:
            out_tsv = Path(outf.name)
        try:
            self._run_one(guide_23nt, in_tsv, out_tsv)
            rows = _parse_cas_offinder_output(out_tsv)
            counts = self._rows_to_clusters(rows, guide)
        finally:
            for p in (in_tsv, out_tsv):
                try:
                    p.unlink()
                except FileNotFoundError:
                    pass
        guide.offtarget = counts
        guide.specificity = compute_specificity(counts)
        return counts
