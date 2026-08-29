# CRISPR gRNA Designer v2

Python package that fetches a gene from NCBI, scans for CRISPR guide
candidates using a configurable Cas-enzyme registry, scores candidates with
an efficiency + specificity model, filters to the coding region by default,
and (optionally) counts off-target sites genome-wide via the **Cas-OFFinder**
binary.

## Install

```
pip install -e .
```

Or run without installing:

```
python -m grna_designer NM_007294 --top 15
```

## Run

```
python -m grna_designer NM_007294 --window 2000 --top 15 --csv guides.csv
```

### Custom Cas enzyme

```bash
python -m grna_designer NM_007294 --cas cas12a --top 10
python -m grna_designer NM_007294 --cas sacas9 --top 10
```

Supported enzymes:

- `spcas9`: PAM `NGG`, 20 nt spacer, 3' PAM orientation
- `sacas9`: PAM `NNGRRT`, 21 nt spacer, 3' PAM orientation
- `cas12a`: PAM `TTTV`, 23 nt spacer, 5' PAM orientation
- `spry`: PAM `NRN`, 20 nt spacer, 3' PAM orientation

    --offtarget-mode genome \
    --genome-fasta hg38.fa \
    --cas-offinder-path cas-offinder.exe \
    --top 15
```

Download Cas-OFFinder from <https://github.com/snugel/cas-offinder/releases>
(Windows, macOS, and Linux binaries are all available).

## Cloning oligos and vendor export

Each scored guide is annotated with cloning oligos and approximate melting
values before export:

- `Oligo_Fwd`
- `Oligo_Rev`
- `Tm_Fwd`
- `Tm_Rev`

For SpCas9 / SaCas9 cloning the overhangs follow the standard lentiCRISPRv2
style design:

- guide starts with `G`: `CACC + guide`
- guide does not start with `G`: `CACCG + guide`
- reverse: `AAAC + revcomp(guide)` and, when needed, a terminal `C`

The CLI can also export provider-ready CSVs for synthesis ordering:

```bash
python -m grna_designer NM_007294 --vendor idt --csv guides.csv
python -m grna_designer NM_007294 --vendor twist --csv guides.csv
python -m grna_designer NM_007294 --vendor genscript --csv guides.csv
```

## Options

| Flag | Default | Meaning |
|------|---------|---------|
| `accession` | — | NCBI ID, e.g. `NM_007294` |
| `--window` | 2000 | bases to scan from start of fetched sequence |
| `--top` | 15 | number of guides to print |
| `--email` | grna-designer@example.com | NCBI requires an email |
| `--csv` | — | path to write CSV results |
| `--cas` | `spcas9` | enzyme registry: `spcas9`, `sacas9`, `cas12a`, `spry` |
| `--mismatches` | 3 | max mismatches for off-target counting |
| `--offtarget-mode` | transcript | `transcript` or `genome` (Cas-OFFinder) |
| `--genome-fasta` | — | genome FASTA (required for `--offtarget-mode genome`) |
| `--cas-offinder-path` | — | path to `cas-offinder` binary; or set `CAS_OFFINDER_BIN` |
| `--include-utr` | off | keep guides in 5'/3' UTRs (default: CDS only) |
| `--no-offtarget` | off | skip off-target search entirely |

## Scoring

**Efficiency** (0..1) — guide quality + activity prior:
- lightweight Doench-inspired position prior
- GC-content penalty outside the preferred 40–60% range
- structural flags for poly-T / Pol III terminators, strong homopolymers,
  and hairpin/self-complementary sequences
- logistic sigmoid over the linearized score

**Specificity** (0..1) — Hsu-Zhang seed-heavy mismatch penalty model:
- mismatch positions are weighted more heavily in the seed region (bases 1–10)
  than in distal positions (11–20)
- this scoring is designed to discourage low-specificity guides with strong seed
  mismatch penalties

**Combined** = 0.6 x efficiency + 0.4 x specificity

Ratings: >=0.75 Excellent / >=0.55 Good / >=0.40 Marginal / else Poor.

## Example accessions

| Accession | Gene |
|-----------|------|
| NM_007294 | BRCA1 |
| NM_000546 | TP53 |
| NM_005228 | EGFR |
| NM_001005484 | OR1D2 |

## Tests

```
pip install -e .[test]
python -m pytest tests/
```

The suite runs offline; Entrez and Cas-OFFinder are mocked in
`tests/conftest.py`.

## Files

- `grna_designer/` - v2 package
- `grna_designer_v1.py` - v1, retained for reference
- `tests/` - pytest suite
- `brca1_guides.csv` - v1 sample output
- `CHANGELOG.md` - what changed in v2

## Citations

- Doench et al. 2016, "Optimized sgRNA design to maximize activity and
  minimize off-target effects of CRISPR-Cas9", *Nature Biotechnology*
  34(2):184-191.
- Bae et al. 2014, "Cas-OFFinder: a fast and versatile algorithm that
  searches for potential off-target sites of Cas9 RNA-guided endonucleases",
  *Bioinformatics* 30(10):1473-1475.
