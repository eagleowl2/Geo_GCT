"""Affymetrix probe ID → gene symbol converter via MyGene.info with caching."""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Sequence
from pathlib import Path

from gct_pipeline.models import PipelineError, RetryableError

_DEFAULT_CACHE = Path.home() / ".gct_pipeline" / "probe_cache.json"
_MYGENE_URL = "https://mygene.info/v3/query"


def convert_probes(
    probe_ids: Sequence[str],
    cache_path: Path | None = None,
    batch_size: int = 500,
    platform: str = "hgu133plus2",
) -> dict[str, str]:
    """Map Affymetrix probe IDs to gene symbols using MyGene.info.

    Results are persisted in a local JSON cache so repeat runs are instant.
    Unmapped probes receive an empty string.

    Args:
        probe_ids: Probe identifiers to look up.
        cache_path: JSON file used as persistent cache. Defaults to
            ~/.gct_pipeline/probe_cache.json.
        batch_size: Number of probes per HTTP POST (MyGene.info cap is ~1000;
            500 is conservative to stay well under URL/body limits).
        platform: Affymetrix platform slug passed as ``scopes`` hint.

    Returns:
        Dict mapping every input probe_id to its gene symbol (or "").
    """
    resolved_cache = cache_path or _DEFAULT_CACHE
    cache = _load_cache(resolved_cache)

    unique_ids = list(dict.fromkeys(probe_ids))  # deduplicate, preserve order
    missing = [p for p in unique_ids if p not in cache]

    if missing:
        fetched = _fetch_all(missing, batch_size=batch_size, platform=platform)
        cache.update(fetched)
        _save_cache(resolved_cache, cache)

    return {p: cache.get(p, "") for p in probe_ids}


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------


def _load_cache(path: Path) -> dict[str, str]:
    if path.exists():
        try:
            data: dict[str, str] = json.loads(path.read_text(encoding="utf-8"))
            return data
        except (json.JSONDecodeError, OSError):
            return {}
    return {}


def _save_cache(path: Path, mapping: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(mapping, indent=2), encoding="utf-8")


# ---------------------------------------------------------------------------
# Fetching
# ---------------------------------------------------------------------------


def _fetch_all(
    probe_ids: list[str],
    batch_size: int,
    platform: str,
) -> dict[str, str]:
    result: dict[str, str] = {}
    for i in range(0, len(probe_ids), batch_size):
        batch = probe_ids[i : i + batch_size]

        def _call(b: list[str] = batch, p: str = platform) -> dict[str, str]:
            return _post_mygene(b, p)

        partial = _with_retry(_call)
        result.update(partial)
        # Ensure every probe in the batch has an entry (even if unmapped)
        for probe in batch:
            result.setdefault(probe, "")
    return result


def _post_mygene(probe_batch: list[str], platform: str) -> dict[str, str]:
    """POST a single batch to MyGene.info and return probe→symbol mapping."""
    body = urllib.parse.urlencode(
        {
            "q": ",".join(probe_batch),
            "scopes": f"reporter,{platform}",
            "fields": "symbol,name",
            "species": "human",
            "size": len(probe_batch),
        }
    ).encode()

    req = urllib.request.Request(
        _MYGENE_URL,
        data=body,
        method="POST",
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )

    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            payload = json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        if exc.code == 429 or exc.code >= 500:
            raise RetryableError(f"MyGene.info HTTP {exc.code}") from exc
        raise PipelineError(f"MyGene.info HTTP {exc.code}") from exc
    except urllib.error.URLError as exc:
        raise RetryableError(f"Network error: {exc.reason}") from exc

    mapping: dict[str, str] = {}
    for hit in payload:
        if isinstance(hit, dict) and "notfound" not in hit:
            probe = hit.get("query", "")
            symbol = hit.get("symbol", "")
            if probe:
                mapping[probe] = symbol
    return mapping


# ---------------------------------------------------------------------------
# Retry
# ---------------------------------------------------------------------------


def _with_retry(
    fn: Callable[[], dict[str, str]],
    max_attempts: int = 5,
    base: float = 1.0,
    factor: float = 2.0,
    cap: float = 30.0,
) -> dict[str, str]:
    """Call *fn* with exponential backoff on RetryableError."""
    last_exc: RetryableError | None = None
    for attempt in range(max_attempts):
        try:
            return fn()
        except RetryableError as exc:
            last_exc = exc
            if attempt == max_attempts - 1:
                break
            delay = min(base * (factor**attempt), cap)
            time.sleep(delay)

    assert last_exc is not None
    raise last_exc
