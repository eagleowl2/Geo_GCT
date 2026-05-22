"""Tests for gct_pipeline.exporter."""

from __future__ import annotations

import io

import numpy as np

from gct_pipeline.exporter import export_gct
from gct_pipeline.models import SeriesMatrix


def _make_matrix(
    n_probes: int = 3,
    n_samples: int = 3,
) -> SeriesMatrix:
    probe_ids = tuple(f"probe_{i}" for i in range(n_probes))
    sample_ids = tuple(f"GSM{i:03d}" for i in range(n_samples))
    expression = np.arange(n_probes * n_samples, dtype=np.float32).reshape(
        n_probes, n_samples
    )
    return SeriesMatrix(
        metadata={"Series_title": "Test"},
        sample_ids=sample_ids,
        probe_ids=probe_ids,
        expression=expression,
    )


def _export_to_str(
    matrix: SeriesMatrix, gene_mapping: dict[str, str]
) -> str:
    buf = io.StringIO()
    export_gct(matrix, gene_mapping, buf)
    return buf.getvalue()


# ---------------------------------------------------------------------------
# Header format
# ---------------------------------------------------------------------------


def test_first_line_is_version() -> None:
    matrix = _make_matrix()
    output = _export_to_str(matrix, {})
    assert output.splitlines()[0] == "#1.2"


def test_second_line_dimensions() -> None:
    matrix = _make_matrix(n_probes=5, n_samples=4)
    output = _export_to_str(matrix, {})
    assert output.splitlines()[1] == "5\t4"


def test_column_header_format() -> None:
    matrix = _make_matrix()
    header = output_lines(_export_to_str(matrix, {}))[2]
    assert header.startswith("Name\tDescription\t")
    assert "GSM000" in header
    assert "GSM002" in header


def output_lines(text: str) -> list[str]:
    return text.splitlines()


# ---------------------------------------------------------------------------
# Row count
# ---------------------------------------------------------------------------


def test_total_line_count() -> None:
    matrix = _make_matrix(n_probes=4, n_samples=2)
    lines = output_lines(_export_to_str(matrix, {}))
    # 3 header lines + 4 data rows
    assert len(lines) == 7


# ---------------------------------------------------------------------------
# Gene symbol in output
# ---------------------------------------------------------------------------


def test_gene_symbol_written_in_description_column() -> None:
    matrix = _make_matrix(n_probes=2, n_samples=2)
    mapping = {"probe_0": "TP53", "probe_1": "BRCA1"}
    lines = output_lines(_export_to_str(matrix, mapping))
    assert lines[3].split("\t")[1] == "TP53"
    assert lines[4].split("\t")[1] == "BRCA1"


def test_unmapped_probe_gets_empty_description() -> None:
    matrix = _make_matrix(n_probes=1, n_samples=1)
    lines = output_lines(_export_to_str(matrix, {}))
    assert lines[3].split("\t")[1] == ""


# ---------------------------------------------------------------------------
# Float formatting
# ---------------------------------------------------------------------------


def test_float_values_formatted_to_4_decimals() -> None:
    matrix = _make_matrix(n_probes=1, n_samples=1)
    lines = output_lines(_export_to_str(matrix, {}))
    # expression[0,0] = 0.0 → "0.0000"
    assert lines[3].split("\t")[2] == "0.0000"


def test_non_integer_float_formatting() -> None:
    probe_ids = ("p1",)
    sample_ids = ("S1",)
    expr = np.array([[3.14159]], dtype=np.float32)
    matrix = SeriesMatrix(
        metadata={},
        sample_ids=sample_ids,
        probe_ids=probe_ids,
        expression=expr,
    )
    lines = output_lines(_export_to_str(matrix, {}))
    val = lines[3].split("\t")[2]
    assert val.startswith("3.141")


# ---------------------------------------------------------------------------
# File path output
# ---------------------------------------------------------------------------


def test_export_to_file_path(tmp_path: object) -> None:
    import os
    import tempfile
    matrix = _make_matrix()
    with tempfile.NamedTemporaryFile(mode="w", suffix=".gct", delete=False) as f:
        path = f.name
    try:
        export_gct(matrix, {}, path)
        with open(path, encoding="utf-8") as f:
            content = f.read()
        assert content.startswith("#1.2")
    finally:
        os.unlink(path)
