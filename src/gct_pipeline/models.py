from __future__ import annotations

from dataclasses import dataclass

import numpy as np


class PipelineError(Exception):
    """Base exception for all pipeline errors."""


class ParseError(PipelineError):
    """Raised when a GEO Series Matrix file cannot be parsed."""


class RetryableError(PipelineError):
    """Raised on transient HTTP errors (429, 5xx) — caller should retry."""


@dataclass(frozen=True)
class SeriesMatrix:
    """Parsed GEO Series Matrix: immutable, numpy-backed expression data."""

    metadata: dict[str, str]
    sample_ids: tuple[str, ...]
    probe_ids: tuple[str, ...]
    # shape: (n_probes, n_samples), dtype float32
    expression: np.ndarray

    def __post_init__(self) -> None:
        expected = (len(self.probe_ids), len(self.sample_ids))
        if self.expression.shape != expected:
            raise ValueError(
                f"expression shape {self.expression.shape} != expected {expected}"
            )


@dataclass(frozen=True)
class GeneMapping:
    """Probe ID → gene symbol mapping; empty string for unmapped probes."""

    mapping: dict[str, str]
