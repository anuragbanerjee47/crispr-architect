from __future__ import annotations

from pathlib import Path

from grna_designer.cli import _build_argparser, _parse_batch_targets


def test_batch_arg_parser_accepts_batch_input(tmp_path):
    parser = _build_argparser()
    args = parser.parse_args(
        ["--batch", "targets.txt", "--output-dir", "tmp_batch"])
    assert args.batch == "targets.txt"
    assert args.output_dir == "tmp_batch"


def test_parse_batch_targets_reads_entries(tmp_path):
    batch = tmp_path / "targets.txt"
    batch.write_text("# comment\nNM_007294\nBRCA1\n", encoding="utf-8")
    assert _parse_batch_targets(str(batch)) == ["NM_007294", "BRCA1"]
