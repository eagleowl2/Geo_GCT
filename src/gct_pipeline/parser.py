"""Streaming parser for NCBI GEO Series Matrix files."""

from __future__ import annotations

from enum import Enum, auto
from pathlib import Path
from typing import IO

import numpy as np

from gct_pipeline.models import ParseError, SeriesMatrix


class _State(Enum):
    HEADER = auto()
    TABLE_HEADER = auto()
    TABLE_DATA = auto()
    DONE = auto()


def parse(source: str | Path | IO[str]) -> SeriesMatrix:
    """Parse a GEO Series Matrix file into an immutable SeriesMatrix.

    Args:
        source: File path (str/Path) or readable text stream (StringIO, open file).

    Returns:
        SeriesMatrix with metadata, sample IDs, probe IDs, and float32 expression array.

    Raises:
        ParseError: On malformed input (missing table markers, bad floats, etc.).
    """
    if isinstance(source, (str, Path)):
        with open(source, encoding="utf-8") as fh:
            return _parse_stream(fh)
    return _parse_stream(source)


def _parse_stream(fh: IO[str]) -> SeriesMatrix:
    metadata: dict[str, str] = {}
    sample_ids: tuple[str, ...] = ()
    sample_titles: tuple[str, ...] = ()
    probe_ids: list[str] = []
    rows: list[np.ndarray] = []
    state = _State.HEADER

    for lineno, raw in enumerate(fh, 1):
        line = raw.rstrip("\n\r")

        if state is _State.HEADER:
            if line.startswith("!series_matrix_table_begin"):
                if not sample_ids:
                    raise ParseError(
                        "!Sample_geo_accession not found before table begin"
                    )
                state = _State.TABLE_HEADER
                continue

            if line.startswith("!"):
                _parse_meta_line(line, metadata)
                if line.startswith("!Sample_geo_accession"):
                    sample_ids = _extract_sample_ids(line)
                elif line.startswith("!Sample_title"):
                    sample_titles = _extract_sample_ids(line)  # same tab-split logic

        elif state is _State.TABLE_HEADER:
            if not line or line.startswith("!"):
                continue
            cols = [c.strip('"') for c in line.split("\t")]
            # Drop trailing empty columns (artifact of trailing \t in some GEO files)
            while cols and not cols[-1]:
                cols.pop()
            if not cols or cols[0] != "ID_REF":
                got = cols[0] if cols else ""
                raise ParseError(
                    f"line {lineno}: expected table header starting with"
                    f" 'ID_REF', got {got!r}"
                )
            header_samples = tuple(cols[1:])
            if sample_ids and header_samples != sample_ids:
                # Allow header to be authoritative if metadata had none
                pass
            sample_ids = header_samples
            state = _State.TABLE_DATA

        elif state is _State.TABLE_DATA:
            if line.startswith("!series_matrix_table_end"):
                state = _State.DONE
                break
            if not line:
                continue
            cols = line.split("\t")
            probe_ids.append(cols[0].strip('"'))
            # Strip quotes; drop trailing empty columns (artifact of trailing \t)
            value_cols = [c.strip('"') for c in cols[1:]]
            while value_cols and not value_cols[-1]:
                value_cols.pop()
            try:
                values = np.array(value_cols, dtype=np.float32)
            except ValueError as exc:
                raise ParseError(f"line {lineno}: cannot parse floats: {exc}") from exc
            if len(values) != len(sample_ids):
                raise ParseError(
                    f"line {lineno}: got {len(values)} values,"
                    f" expected {len(sample_ids)}"
                )
            rows.append(values)

    if state is not _State.DONE:
        if state is _State.HEADER:
            raise ParseError("!series_matrix_table_begin not found")
        raise ParseError("!series_matrix_table_end not found")

    if not rows:
        expression = np.empty((0, len(sample_ids)), dtype=np.float32)
    else:
        expression = np.vstack(rows)

    return SeriesMatrix(
        metadata=metadata,
        sample_ids=sample_ids,
        sample_titles=sample_titles,
        probe_ids=tuple(probe_ids),
        expression=expression,
    )


def _parse_meta_line(line: str, metadata: dict[str, str]) -> None:
    """Split a GEO metadata line into the metadata dict.

    Handles two separator styles found in GEO Series Matrix files:
      - ``!Key = value``  (standard SOFT format, ` = ` with spaces)
      - ``!Key\tvalue``   (tab-separated, used by some fields)

    Only the first value is stored; multi-value tab-separated lines
    (e.g. !Sample_geo_accession) are handled separately by _extract_sample_ids.
    """
    stripped = line.lstrip("!")
    if " = " in stripped:
        key, _, raw_value = stripped.partition(" = ")
        # Take only the first tab-delimited value for the metadata dict
        value = raw_value.split("\t")[0].strip().strip('"')
        metadata[key.strip()] = value
    elif "\t" in stripped:
        key, _, raw_value = stripped.partition("\t")
        # Take only the first tab-delimited value (multi-sample lines
        # repeat the same value once per sample — keep just one).
        first_value = raw_value.split("\t")[0]
        metadata[key.strip()] = first_value.strip().strip('"')


def _extract_sample_ids(line: str) -> tuple[str, ...]:
    """Extract sample IDs from a tab-separated !Sample_geo_accession line."""
    parts = line.split("\t")
    # First part is '!Sample_geo_accession = GSM001' or just '!Sample_geo_accession'
    first = parts[0]
    if "=" in first:
        _, _, first_id = first.partition("=")
        first_id = first_id.strip().strip('"')
        ids = [first_id] + [p.strip().strip('"') for p in parts[1:]]
    else:
        ids = [p.strip().strip('"') for p in parts[1:]]
    return tuple(i for i in ids if i)
