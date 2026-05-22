"""Shared fixtures and helpers for the gct_pipeline test suite."""

from __future__ import annotations

import io
import json
from typing import Any
from unittest.mock import MagicMock

import pytest

# ---------------------------------------------------------------------------
# Synthetic GEO Series Matrix content
# ---------------------------------------------------------------------------

MINI_MATRIX_TEXT = """\
!Series_title = Test Series
!Series_geo_accession = GSE999
!Sample_geo_accession = GSM001\tGSM002\tGSM003
!Sample_title = Control 1\tControl 2\tTreated 1
!series_matrix_table_begin
ID_REF\tGSM001\tGSM002\tGSM003
1007_s_at\t5.1\t6.2\t7.3
1053_at\t8.4\t9.5\t10.6
121_at\t1.0\t2.0\t3.0
!series_matrix_table_end
"""

EMPTY_TABLE_TEXT = """\
!Series_title = Empty
!Sample_geo_accession = GSM001\tGSM002
!series_matrix_table_begin
ID_REF\tGSM001\tGSM002
!series_matrix_table_end
"""

NO_TABLE_END_TEXT = """\
!Series_title = Broken
!Sample_geo_accession = GSM001
!series_matrix_table_begin
ID_REF\tGSM001
1007_s_at\t5.1
"""

NO_TABLE_BEGIN_TEXT = """\
!Series_title = NoTable
!Sample_geo_accession = GSM001
"""

MALFORMED_FLOAT_TEXT = """\
!Series_title = Bad
!Sample_geo_accession = GSM001\tGSM002
!series_matrix_table_begin
ID_REF\tGSM001\tGSM002
1007_s_at\t5.1\tnot_a_number
!series_matrix_table_end
"""


@pytest.fixture
def mini_matrix_text() -> str:
    return MINI_MATRIX_TEXT


@pytest.fixture
def mini_stream() -> io.StringIO:
    return io.StringIO(MINI_MATRIX_TEXT)


# ---------------------------------------------------------------------------
# MyGene.info mock helpers
# ---------------------------------------------------------------------------


def make_mygene_response(hits: list[dict[str, Any]]) -> MagicMock:
    """Build a fake urllib urlopen response returning *hits* as JSON."""
    body = json.dumps(hits).encode()
    mock_resp = MagicMock()
    mock_resp.read.return_value = body
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    return mock_resp


def standard_hits(probe_ids: list[str]) -> list[dict[str, Any]]:
    """Return canned MyGene.info hits mapping each probe to a gene symbol."""
    symbols = {
        "1007_s_at": "DDR1",
        "1053_at": "RFC2",
        "121_at": "PAX8",
    }
    return [
        {"query": p, "symbol": symbols.get(p, f"GENE_{p}")}
        for p in probe_ids
    ]
