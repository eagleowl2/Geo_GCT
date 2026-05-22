"""Tests for gct_pipeline.parser."""

from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pytest

from gct_pipeline.models import ParseError, SeriesMatrix
from gct_pipeline.parser import parse

from .conftest import (
    EMPTY_TABLE_TEXT,
    MALFORMED_FLOAT_TEXT,
    MINI_MATRIX_TEXT,
    NO_TABLE_BEGIN_TEXT,
    NO_TABLE_END_TEXT,
)

# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------


def test_parse_stream_returns_series_matrix() -> None:
    result = parse(io.StringIO(MINI_MATRIX_TEXT))
    assert isinstance(result, SeriesMatrix)


def test_parse_metadata() -> None:
    result = parse(io.StringIO(MINI_MATRIX_TEXT))
    assert result.metadata["Series_title"] == "Test Series"
    assert result.metadata["Series_geo_accession"] == "GSE999"


def test_parse_sample_ids() -> None:
    result = parse(io.StringIO(MINI_MATRIX_TEXT))
    assert result.sample_ids == ("GSM001", "GSM002", "GSM003")


def test_parse_sample_titles() -> None:
    result = parse(io.StringIO(MINI_MATRIX_TEXT))
    assert result.sample_titles == ("Control 1", "Control 2", "Treated 1")


def test_parse_no_sample_titles_gives_empty_tuple() -> None:
    text_without_titles = MINI_MATRIX_TEXT.replace(
        "!Sample_title = Control 1\tControl 2\tTreated 1\n", ""
    )
    result = parse(io.StringIO(text_without_titles))
    assert result.sample_titles == ()


def test_parse_probe_ids() -> None:
    result = parse(io.StringIO(MINI_MATRIX_TEXT))
    assert result.probe_ids == ("1007_s_at", "1053_at", "121_at")


def test_parse_expression_shape() -> None:
    result = parse(io.StringIO(MINI_MATRIX_TEXT))
    assert result.expression.shape == (3, 3)


def test_parse_expression_dtype_float32() -> None:
    result = parse(io.StringIO(MINI_MATRIX_TEXT))
    assert result.expression.dtype == np.float32


def test_parse_expression_values() -> None:
    result = parse(io.StringIO(MINI_MATRIX_TEXT))
    np.testing.assert_allclose(result.expression[0], [5.1, 6.2, 7.3], rtol=1e-5)
    np.testing.assert_allclose(result.expression[2], [1.0, 2.0, 3.0], rtol=1e-5)


def test_parse_result_is_frozen() -> None:
    result = parse(io.StringIO(MINI_MATRIX_TEXT))
    with pytest.raises((AttributeError, TypeError)):
        result.metadata = {}  # type: ignore[misc]


# ---------------------------------------------------------------------------
# File path input
# ---------------------------------------------------------------------------


def test_parse_from_path(tmp_path: Path) -> None:
    matrix_file = tmp_path / "test.txt"
    matrix_file.write_text(MINI_MATRIX_TEXT, encoding="utf-8")
    result = parse(matrix_file)
    assert result.probe_ids == ("1007_s_at", "1053_at", "121_at")


def test_parse_from_str_path(tmp_path: Path) -> None:
    matrix_file = tmp_path / "test.txt"
    matrix_file.write_text(MINI_MATRIX_TEXT, encoding="utf-8")
    result = parse(str(matrix_file))
    assert len(result.probe_ids) == 3


# ---------------------------------------------------------------------------
# Empty expression table
# ---------------------------------------------------------------------------


def test_parse_empty_table() -> None:
    result = parse(io.StringIO(EMPTY_TABLE_TEXT))
    assert result.probe_ids == ()
    assert result.expression.shape == (0, 2)


# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------


def test_parse_missing_table_end_raises() -> None:
    with pytest.raises(ParseError, match="table_end"):
        parse(io.StringIO(NO_TABLE_END_TEXT))


def test_parse_missing_table_begin_raises() -> None:
    with pytest.raises(ParseError, match="table_begin"):
        parse(io.StringIO(NO_TABLE_BEGIN_TEXT))


def test_parse_malformed_float_raises() -> None:
    with pytest.raises(ParseError, match="float"):
        parse(io.StringIO(MALFORMED_FLOAT_TEXT))
