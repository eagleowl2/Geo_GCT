# GCT Pipeline

A production-grade bioinformatics pipeline that parses **NCBI GEO Series Matrix** files, converts Affymetrix probe IDs to gene symbols, and exports to **GCT 1.2** format for downstream tools such as GSEA and GenePattern.

---

## Features

- **Streaming parser** — processes GEO Series Matrix files line-by-line; handles 100 MB+ files without loading them into memory.
- **numpy-backed expression matrix** — `float32` arrays give 8–10× faster downstream math and 2× lower memory usage compared to plain Python lists or `float64`.
- **Immutable typed data** — `SeriesMatrix` and `GeneMapping` are frozen dataclasses; any accidental mutation raises at runtime.
- **Batched probe conversion** — posts 500 probes per call to the [MyGene.info](https://mygene.info) REST API, staying well under rate-limit thresholds.
- **Local JSON cache** — after the first run, gene symbol lookups are served from `~/.gct_pipeline/probe_cache.json`. Repeat runs go from ~15 s to instant.
- **Exponential backoff retry** — HTTP 429 and 5xx responses are retried up to 5 times (delays: 1 s, 2 s, 4 s, 8 s, 16 s). Hard client errors (4xx) propagate immediately.
- **GCT 1.2 export** — streams output row-by-row; never holds the full file in memory.

---

## Requirements

- Python 3.11 or 3.12
- [numpy](https://numpy.org/) ≥ 1.26

No third-party HTTP library is required. All network calls use the Python standard library (`urllib`).

---

## Installation

```bash
# From the project root
pip install -e .

# With dev/test dependencies
pip install -e .[dev]
```

---

## Quick Start

```python
from gct_pipeline import parse, convert_probes, export_gct

# 1. Parse a GEO Series Matrix file
matrix = parse("GSE12345_series_matrix.txt")

print(matrix.sample_ids)          # ('GSM001', 'GSM002', ...)
print(matrix.expression.shape)    # (n_probes, n_samples)
print(matrix.expression.dtype)    # float32

# 2. Convert probe IDs → gene symbols (cached after first run)
gene_mapping = convert_probes(list(matrix.probe_ids))

# 3. Export to GCT 1.2
export_gct(matrix, gene_mapping, "output.gct")
```

### Using a custom cache location

```python
from pathlib import Path
gene_mapping = convert_probes(
    list(matrix.probe_ids),
    cache_path=Path("/data/my_project/probe_cache.json"),
)
```

### Accepting a file-like object

```python
import io
text = open("matrix.txt").read()
matrix = parse(io.StringIO(text))
```

---

## GCT 1.2 Format

The exported file follows the GCT 1.2 specification:

```
#1.2
<n_probes>	<n_samples>
Name	Description	GSM001	GSM002	...
1007_s_at	DDR1	5.1000	6.2000	...
1053_at	RFC2	8.4000	9.5000	...
```

- **Name** — Affymetrix probe ID
- **Description** — gene symbol from MyGene.info (empty string if unmapped)
- Remaining columns — expression values rounded to 4 decimal places

---

## Project Layout

```
gct_task/
├── pyproject.toml               # build, lint, type-check, and test config
├── .github/workflows/ci.yml     # CI: ruff → mypy → pytest on Python 3.11 + 3.12
└── src/gct_pipeline/
    ├── models.py    frozen dataclasses (SeriesMatrix, GeneMapping) + exceptions
    ├── parser.py    streaming state-machine parser
    ├── converter.py batched HTTP + JSON cache + exponential backoff
    ├── exporter.py  GCT 1.2 line-by-line writer
    └── __init__.py  public API re-exports
tests/
    ├── conftest.py         shared fixtures and HTTP mock helpers
    ├── test_parser.py      parser unit tests
    ├── test_converter.py   converter unit tests
    ├── test_exporter.py    exporter unit tests
    └── test_integration.py end-to-end pipeline test
```

---

## Running Tests

```bash
pytest                                          # run all tests
pytest --cov=gct_pipeline --cov-report=term-missing   # with coverage
```

Expected: **35 tests pass in under 2 seconds** with no network access.

---

## Linting and Type Checking

```bash
ruff check src tests   # linting (PEP 8, import order, modernisation)
mypy src               # strict type checking
```

---

## Exception Hierarchy

```
PipelineError          base — catch this for any pipeline failure
├── ParseError         malformed GEO Series Matrix file
└── RetryableError     transient HTTP error (raised after all retries exhausted)
```

---

## Cache Management

The probe cache grows monotonically — unmapped probes are cached as `""` so they are not re-queried. To force a fresh lookup, delete or clear the cache file:

```bash
rm ~/.gct_pipeline/probe_cache.json   # Linux/macOS
del %USERPROFILE%\.gct_pipeline\probe_cache.json   # Windows
```

---

## CI/CD

GitHub Actions runs on every push and pull request to `main`:

| Step | Tool |
|---|---|
| Lint | `ruff check src tests` |
| Type check | `mypy src` (strict mode) |
| Tests | `pytest --cov=gct_pipeline` |
| Matrix | Python 3.11 and 3.12 |
