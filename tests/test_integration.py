"""End-to-end integration test: parse → convert (mocked) → export → verify."""

from __future__ import annotations

import io
from pathlib import Path
from unittest.mock import patch

from gct_pipeline.exporter import export_gct
from gct_pipeline.parser import parse

from .conftest import MINI_MATRIX_TEXT, make_mygene_response, standard_hits


def test_full_pipeline_produces_valid_gct(tmp_path: Path) -> None:
    # 1. Parse
    matrix = parse(io.StringIO(MINI_MATRIX_TEXT))
    assert len(matrix.probe_ids) == 3

    # 2. Convert (mocked HTTP)
    probes = list(matrix.probe_ids)
    hits = standard_hits(probes)
    mock_resp = make_mygene_response(hits)
    cache_file = tmp_path / "cache.json"

    with patch("urllib.request.urlopen", return_value=mock_resp):
        from gct_pipeline.converter import convert_probes

        gene_mapping = convert_probes(probes, cache_path=cache_file)

    assert gene_mapping["1007_s_at"] == "DDR1"

    # 3. Export
    gct_buf = io.StringIO()
    export_gct(matrix, gene_mapping, gct_buf)
    gct_content = gct_buf.getvalue()

    # 4. Verify GCT structure
    lines = gct_content.splitlines()
    assert lines[0] == "#1.2"
    assert lines[1] == "3\t3"
    assert lines[2].startswith("Name\tDescription\tGSM001")

    data_lines = lines[3:]
    assert len(data_lines) == 3

    first_row = data_lines[0].split("\t")
    assert first_row[0] == "1007_s_at"
    assert first_row[1] == "DDR1"
    # Values match original (within float32 precision)
    assert abs(float(first_row[2]) - 5.1) < 0.01


def test_cache_makes_second_run_skip_http(tmp_path: Path) -> None:
    matrix = parse(io.StringIO(MINI_MATRIX_TEXT))
    probes = list(matrix.probe_ids)
    cache_file = tmp_path / "cache.json"

    hits = standard_hits(probes)
    mock_resp = make_mygene_response(hits)

    with patch("urllib.request.urlopen", return_value=mock_resp) as mock_url:
        from gct_pipeline.converter import convert_probes

        convert_probes(probes, cache_path=cache_file)
        assert mock_url.call_count == 1

        # Second call — cache is warm, no HTTP
        convert_probes(probes, cache_path=cache_file)
        assert mock_url.call_count == 1  # still 1
