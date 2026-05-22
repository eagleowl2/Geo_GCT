"""GCT 1.2 format writer for SeriesMatrix data."""

from __future__ import annotations

from pathlib import Path
from typing import IO

from gct_pipeline.models import SeriesMatrix


def export_gct(
    matrix: SeriesMatrix,
    gene_mapping: dict[str, str],
    dest: str | Path | IO[str],
) -> None:
    """Write *matrix* to GCT 1.2 format.

    Args:
        matrix: Parsed expression data.
        gene_mapping: probe_id → gene_symbol; missing probes get empty description.
        dest: Output file path (str/Path) or writable text stream.
    """
    if isinstance(dest, (str, Path)):
        with open(dest, "w", encoding="utf-8", newline="\n") as fh:
            _write_gct(matrix, gene_mapping, fh)
    else:
        _write_gct(matrix, gene_mapping, dest)


def _write_gct(
    matrix: SeriesMatrix,
    gene_mapping: dict[str, str],
    fh: IO[str],
) -> None:
    n_probes = len(matrix.probe_ids)
    n_samples = len(matrix.sample_ids)

    fh.write("#1.2\n")
    fh.write(f"{n_probes}\t{n_samples}\n")
    fh.write("Name\tDescription\t" + "\t".join(matrix.sample_ids) + "\n")

    for i, probe_id in enumerate(matrix.probe_ids):
        gene_symbol = gene_mapping.get(probe_id, "")
        values = "\t".join(f"{v:.4f}" for v in matrix.expression[i])
        fh.write(f"{probe_id}\t{gene_symbol}\t{values}\n")
