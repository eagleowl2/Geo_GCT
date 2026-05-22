# Test Suite Notes

## Overview

The test suite uses **pytest** and covers all three pipeline modules plus a full end-to-end integration path. It runs in under 2 seconds and makes **zero real network calls** — all HTTP is intercepted at the `urllib.request.urlopen` layer.

```
35 passed in 1.45s   |   94% line coverage
```

---

## Test Files

| File | Tests | Focus |
|---|---|---|
| `tests/test_parser.py` | 13 | GEO Series Matrix parsing |
| `tests/test_converter.py` | 10 | Probe conversion, caching, retry |
| `tests/test_exporter.py` | 9 | GCT 1.2 output format |
| `tests/test_integration.py` | 2 | Full parse → convert → export pipeline |
| `tests/conftest.py` | — | Shared fixtures and mock helpers |

---

## Coverage Report

```
Name                            Stmts   Miss  Cover   Missing
-------------------------------------------------------------
src/gct_pipeline/__init__.py        5      0   100%
src/gct_pipeline/converter.py      76      4    95%   65-66, 126-127
src/gct_pipeline/exporter.py       19      0   100%
src/gct_pipeline/models.py         19      1    95%   33
src/gct_pipeline/parser.py         84      8    90%   52, 65, 68, 75, 84, 92, 119, 136
-------------------------------------------------------------
TOTAL                             203     13    94%
```

Uncovered lines are defensive branches that require contrived internal state (e.g. shape mismatch in `SeriesMatrix.__post_init__`, table-header guard clauses reachable only through direct internal calls).

---

## Key Design Decisions

### No real network calls

All tests patch `urllib.request.urlopen` using `unittest.mock.patch`. The `conftest.py` helper `make_mygene_response()` builds a minimal fake response object that satisfies the context-manager protocol (`__enter__` / `__exit__`) and returns canned JSON via `.read()`. This means:

- Tests run offline.
- Tests are deterministic — no rate-limit flakiness.
- Test speed is bounded by CPU, not network latency.

### No real sleeps in retry tests

`time.sleep` inside `gct_pipeline.converter` is patched to a no-op for all retry tests:

```python
with patch("gct_pipeline.converter.time.sleep") as mock_sleep:
    ...
assert mock_sleep.call_count == 2          # two retries happened
assert sleep_calls == [1.0, 2.0, 4.0]     # backoff delays verified
```

The exact delay sequence is asserted — not just that sleep was called — so regressions in backoff arithmetic are caught immediately.

### Isolated cache directories

Every converter test receives an isolated `tmp_path` fixture directory and passes it as `cache_path`. Tests never touch `~/.gct_pipeline/probe_cache.json`, so:

- Tests are idempotent regardless of local cache state.
- Tests can run in parallel without cache collisions.

### Synthetic input via `io.StringIO`

The parser is tested entirely through in-memory streams. No fixture files are written to disk except in the two tests that explicitly exercise the `Path` / `str` input path (using pytest's `tmp_path`). This keeps the test suite self-contained and fast.

---

## Test Scenarios by Module

### Parser (`test_parser.py`)

| Test | What it checks |
|---|---|
| `test_parse_stream_returns_series_matrix` | Return type is `SeriesMatrix` |
| `test_parse_metadata` | `!Series_*` lines parsed into metadata dict |
| `test_parse_sample_ids` | Tab-separated GSM IDs extracted correctly |
| `test_parse_probe_ids` | Probe IDs match declaration order |
| `test_parse_expression_shape` | ndarray shape is `(n_probes, n_samples)` |
| `test_parse_expression_dtype_float32` | dtype is `float32`, not `float64` |
| `test_parse_expression_values` | Numeric values match source with `rtol=1e-5` |
| `test_parse_result_is_frozen` | Mutation raises `AttributeError` / `TypeError` |
| `test_parse_from_path` | `Path` object accepted as input |
| `test_parse_from_str_path` | `str` path accepted as input |
| `test_parse_empty_table` | Zero-probe table → shape `(0, n_samples)` |
| `test_parse_missing_table_end_raises` | `ParseError` on truncated file |
| `test_parse_missing_table_begin_raises` | `ParseError` when table block absent |
| `test_parse_malformed_float_raises` | `ParseError` on non-numeric cell |

### Converter (`test_converter.py`)

| Test | What it checks |
|---|---|
| `test_cache_hit_no_http` | Fully cached input never calls `urlopen` |
| `test_partial_cache_fetches_only_missing` | Only uncached probes are fetched |
| `test_cache_miss_fetches_and_saves` | Results written to disk after fetch |
| `test_unmapped_probe_returns_empty_string` | MyGene.info `notfound` → `""` |
| `test_batching_splits_at_500` | 1100 probes → exactly 3 HTTP calls |
| `test_batching_exact_500_boundary` | 500 probes → exactly 1 HTTP call |
| `test_retry_on_429_eventually_succeeds` | Two 429s then success → correct result |
| `test_retry_exhausted_raises_retryable_error` | 5 consecutive 503s → `RetryableError` |
| `test_http_400_raises_pipeline_error_immediately` | 400 → `PipelineError`, no retries |
| `test_retry_backoff_delays` | Sleep delays are exactly `[1.0, 2.0, 4.0]` |

### Exporter (`test_exporter.py`)

| Test | What it checks |
|---|---|
| `test_first_line_is_version` | Output starts with `#1.2` |
| `test_second_line_dimensions` | Second line is `{n_probes}\t{n_samples}` |
| `test_column_header_format` | Header starts with `Name\tDescription\t` |
| `test_total_line_count` | 3 header lines + 1 per probe |
| `test_gene_symbol_written_in_description_column` | Gene symbol in column 2 |
| `test_unmapped_probe_gets_empty_description` | Missing mapping → `""` not crash |
| `test_float_values_formatted_to_4_decimals` | `0.0` → `"0.0000"` |
| `test_non_integer_float_formatting` | `3.14159` starts with `"3.141"` |
| `test_export_to_file_path` | File path output roundtrip |

### Integration (`test_integration.py`)

| Test | What it checks |
|---|---|
| `test_full_pipeline_produces_valid_gct` | parse → convert → export → verify GCT content |
| `test_cache_makes_second_run_skip_http` | Second `convert_probes` call hits cache only |

---

## Running the Tests

```bash
# Standard run
pytest

# With coverage report
pytest --cov=gct_pipeline --cov-report=term-missing

# Single module
pytest tests/test_converter.py -v

# Single test
pytest tests/test_parser.py::test_parse_empty_table -v
```

---

## What Is Not Tested

- **Real MyGene.info responses** — intentional. Integration against the live API belongs in a separate, opt-in test marked with a custom pytest marker (e.g. `@pytest.mark.live`) and excluded from CI.
- **Files larger than a few KB** — the streaming parser is architecturally validated; large-file performance is better confirmed with a profiling script than a unit test.
- **Cache corruption recovery** — `_load_cache` silently returns `{}` on `json.JSONDecodeError`; the branch exists but is not exercised. A future test could write a malformed JSON file and assert the converter continues without error.
