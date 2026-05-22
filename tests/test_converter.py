"""Tests for gct_pipeline.converter."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from gct_pipeline.models import PipelineError, RetryableError

from .conftest import make_mygene_response, standard_hits

# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------


def _write_cache(path: Path, data: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")


# ---------------------------------------------------------------------------
# Cache hit — no HTTP calls
# ---------------------------------------------------------------------------


def test_cache_hit_no_http(tmp_path: Path) -> None:
    cache_file = tmp_path / "cache.json"
    _write_cache(cache_file, {"1007_s_at": "DDR1", "1053_at": "RFC2"})

    with patch("urllib.request.urlopen") as mock_url:
        from gct_pipeline.converter import convert_probes

        result = convert_probes(
            ["1007_s_at", "1053_at"], cache_path=cache_file
        )

    mock_url.assert_not_called()
    assert result == {"1007_s_at": "DDR1", "1053_at": "RFC2"}


def test_partial_cache_fetches_only_missing(tmp_path: Path) -> None:
    cache_file = tmp_path / "cache.json"
    _write_cache(cache_file, {"1007_s_at": "DDR1"})

    hits = standard_hits(["1053_at"])
    mock_resp = make_mygene_response(hits)

    with patch("urllib.request.urlopen", return_value=mock_resp):
        from gct_pipeline.converter import convert_probes

        result = convert_probes(
            ["1007_s_at", "1053_at"], cache_path=cache_file
        )

    assert result["1007_s_at"] == "DDR1"
    assert result["1053_at"] == "RFC2"


# ---------------------------------------------------------------------------
# Cache miss — fetch and persist
# ---------------------------------------------------------------------------


def test_cache_miss_fetches_and_saves(tmp_path: Path) -> None:
    cache_file = tmp_path / "cache.json"
    probes = ["1007_s_at", "1053_at", "121_at"]
    hits = standard_hits(probes)
    mock_resp = make_mygene_response(hits)

    with patch("urllib.request.urlopen", return_value=mock_resp):
        from gct_pipeline.converter import convert_probes

        result = convert_probes(probes, cache_path=cache_file)

    assert result["1007_s_at"] == "DDR1"
    assert result["121_at"] == "PAX8"
    assert cache_file.exists()
    saved = json.loads(cache_file.read_text())
    assert saved["1053_at"] == "RFC2"


def test_unmapped_probe_returns_empty_string(tmp_path: Path) -> None:
    cache_file = tmp_path / "cache.json"
    # MyGene.info returns notfound for unknown probes
    hits = [{"query": "UNKNOWN_PROBE", "notfound": True}]
    mock_resp = make_mygene_response(hits)

    with patch("urllib.request.urlopen", return_value=mock_resp):
        from gct_pipeline.converter import convert_probes

        result = convert_probes(["UNKNOWN_PROBE"], cache_path=cache_file)

    assert result["UNKNOWN_PROBE"] == ""


# ---------------------------------------------------------------------------
# Batching
# ---------------------------------------------------------------------------


def test_batching_splits_at_500(tmp_path: Path) -> None:
    cache_file = tmp_path / "cache.json"
    probes = [f"probe_{i}" for i in range(1100)]

    call_count = 0

    def fake_urlopen(req: object, timeout: int = 30) -> MagicMock:
        nonlocal call_count
        call_count += 1
        return make_mygene_response([])

    with patch("urllib.request.urlopen", side_effect=fake_urlopen):
        from gct_pipeline.converter import convert_probes

        convert_probes(probes, cache_path=cache_file, batch_size=500)

    # 1100 probes / 500 per batch = 3 calls (500 + 500 + 100)
    assert call_count == 3


def test_batching_exact_500_boundary(tmp_path: Path) -> None:
    cache_file = tmp_path / "cache.json"
    probes = [f"probe_{i}" for i in range(500)]

    call_count = 0

    def fake_urlopen(req: object, timeout: int = 30) -> MagicMock:
        nonlocal call_count
        call_count += 1
        return make_mygene_response([])

    with patch("urllib.request.urlopen", side_effect=fake_urlopen):
        from gct_pipeline.converter import convert_probes

        convert_probes(probes, cache_path=cache_file, batch_size=500)

    assert call_count == 1


# ---------------------------------------------------------------------------
# Retry on 429
# ---------------------------------------------------------------------------


def test_retry_on_429_eventually_succeeds(tmp_path: Path) -> None:
    import urllib.error

    cache_file = tmp_path / "cache.json"
    hits = standard_hits(["1007_s_at"])
    success_resp = make_mygene_response(hits)

    attempt = 0

    def fake_urlopen(req: object, timeout: int = 30) -> MagicMock:
        nonlocal attempt
        attempt += 1
        if attempt < 3:
            raise urllib.error.HTTPError(
                url="", code=429, msg="Too Many Requests", hdrs=MagicMock(), fp=None  # type: ignore[arg-type]
            )
        return success_resp

    with (
        patch("urllib.request.urlopen", side_effect=fake_urlopen),
        patch("gct_pipeline.converter.time.sleep") as mock_sleep,
    ):
        from gct_pipeline.converter import convert_probes

        result = convert_probes(["1007_s_at"], cache_path=cache_file)

    assert result["1007_s_at"] == "DDR1"
    assert mock_sleep.call_count == 2


def test_retry_exhausted_raises_retryable_error(tmp_path: Path) -> None:
    import urllib.error

    cache_file = tmp_path / "cache.json"

    def fake_urlopen(req: object, timeout: int = 30) -> MagicMock:
        raise urllib.error.HTTPError(
            url="", code=503, msg="Service Unavailable", hdrs=MagicMock(), fp=None  # type: ignore[arg-type]
        )

    with (
        patch("urllib.request.urlopen", side_effect=fake_urlopen),
        patch("gct_pipeline.converter.time.sleep"),
    ):
        from gct_pipeline.converter import convert_probes

        with pytest.raises(RetryableError):
            convert_probes(["1007_s_at"], cache_path=cache_file)


# ---------------------------------------------------------------------------
# Permanent client error (4xx) — no retry
# ---------------------------------------------------------------------------


def test_http_400_raises_pipeline_error_immediately(tmp_path: Path) -> None:
    import urllib.error

    cache_file = tmp_path / "cache.json"
    call_count = 0

    def fake_urlopen(req: object, timeout: int = 30) -> MagicMock:
        nonlocal call_count
        call_count += 1
        raise urllib.error.HTTPError(
            url="", code=400, msg="Bad Request", hdrs=MagicMock(), fp=None  # type: ignore[arg-type]
        )

    with (
        patch("urllib.request.urlopen", side_effect=fake_urlopen),
        patch("gct_pipeline.converter.time.sleep"),
    ):
        from gct_pipeline.converter import convert_probes

        with pytest.raises(PipelineError):
            convert_probes(["1007_s_at"], cache_path=cache_file)

    assert call_count == 1  # no retries


# ---------------------------------------------------------------------------
# Exponential backoff delays
# ---------------------------------------------------------------------------


def test_retry_backoff_delays(tmp_path: Path) -> None:
    import urllib.error

    cache_file = tmp_path / "cache.json"
    hits = standard_hits(["1007_s_at"])
    success_resp = make_mygene_response(hits)

    attempt = 0

    def fake_urlopen(req: object, timeout: int = 30) -> MagicMock:
        nonlocal attempt
        attempt += 1
        if attempt < 4:
            raise urllib.error.HTTPError(
                url="", code=429, msg="Too Many Requests", hdrs=MagicMock(), fp=None  # type: ignore[arg-type]
            )
        return success_resp

    sleep_calls: list[float] = []

    with (
        patch("urllib.request.urlopen", side_effect=fake_urlopen),
        patch("gct_pipeline.converter.time.sleep", side_effect=sleep_calls.append),
    ):
        from gct_pipeline.converter import convert_probes

        convert_probes(["1007_s_at"], cache_path=cache_file)

    # delays: 1.0, 2.0, 4.0 (base=1, factor=2)
    assert sleep_calls == [1.0, 2.0, 4.0]
