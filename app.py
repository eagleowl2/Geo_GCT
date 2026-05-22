"""Streamlit web UI for the GEO Series Matrix → GCT pipeline."""

from __future__ import annotations

import gzip
import io

import pandas as pd
import streamlit as st

from gct_pipeline import (
    ParseError,
    PipelineError,
    RetryableError,
    convert_probes,
    export_gct,
    parse,
)

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="GEO → GCT Pipeline",
    page_icon="🧬",
    layout="centered",
)

# ---------------------------------------------------------------------------
# Session-state initialisation
# ---------------------------------------------------------------------------

for _key in ("matrix", "gene_mapping", "gct_content", "_file_id"):
    if _key not in st.session_state:
        st.session_state[_key] = None


def _reset_state() -> None:
    st.session_state.matrix = None
    st.session_state.gene_mapping = None
    st.session_state.gct_content = None
    st.session_state["_file_id"] = None


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

st.title("GEO → GCT Pipeline")
st.caption(
    "Upload a GEO Series Matrix file (.txt or .gz), preview parsed data, "
    "convert Affymetrix probe IDs to gene symbols, and download the GCT 1.2 output."
)
st.divider()

# ---------------------------------------------------------------------------
# Step 1 — Upload
# ---------------------------------------------------------------------------

st.subheader("Step 1 — Upload Series Matrix file")

uploaded_file = st.file_uploader(
    "Choose a GEO Series Matrix file",
    type=["txt", "gz"],
    help="Accepts plain-text (.txt) or gzip-compressed (.gz) GEO Series Matrix files.",
)

# Detect file change or removal
if uploaded_file is None:
    if st.session_state["_file_id"] is not None:
        _reset_state()
else:
    if st.session_state["_file_id"] != uploaded_file.file_id:
        _reset_state()
        st.session_state["_file_id"] = uploaded_file.file_id

        # Parse immediately on upload
        raw = uploaded_file.read()
        text_stream: io.StringIO | io.TextIOWrapper | None = None
        try:
            if uploaded_file.name.endswith(".gz"):
                text_stream = io.TextIOWrapper(
                    gzip.open(io.BytesIO(raw)), encoding="utf-8"
                )
            else:
                text_stream = io.StringIO(raw.decode("utf-8"))
        except (OSError, UnicodeDecodeError) as exc:
            st.error(f"Could not read uploaded file: {exc}")
            _reset_state()

        if text_stream is not None:
            try:
                st.session_state.matrix = parse(text_stream)
            except ParseError as exc:
                st.error(f"Parse error: {exc}")
                _reset_state()
            except Exception as exc:
                st.error(f"Unexpected error during parsing: {exc}")
                _reset_state()

# ---------------------------------------------------------------------------
# Step 2 — Series Preview (only when matrix is available)
# ---------------------------------------------------------------------------

if st.session_state.matrix is not None:
    matrix = st.session_state.matrix
    st.divider()
    st.subheader("Step 2 — Series Preview")

    # --- top-level metrics ---
    accession = matrix.metadata.get("Series_geo_accession", "—")
    title = matrix.metadata.get("Series_title", "—")

    c1, c2, c3 = st.columns(3)
    c1.metric("Probes", f"{len(matrix.probe_ids):,}")
    c2.metric("Samples", f"{len(matrix.sample_ids):,}")
    c3.metric("Accession", accession)

    if title != "—":
        st.caption(f"**Series title:** {title}")

    # --- full metadata expander ---
    with st.expander("Full metadata", expanded=False):
        meta_df = pd.DataFrame(
            list(matrix.metadata.items()), columns=["Key", "Value"]
        )
        st.dataframe(meta_df, use_container_width=True, hide_index=True)

    # --- expression preview ---
    st.subheader("Expression matrix preview (first 50 rows)")
    n_preview = min(50, len(matrix.probe_ids))
    preview_df = pd.DataFrame(
        matrix.expression[:n_preview],
        index=list(matrix.probe_ids[:n_preview]),
        columns=list(matrix.sample_ids),
    )
    preview_df.index.name = "Probe ID"
    st.dataframe(preview_df, use_container_width=True)

    # ---------------------------------------------------------------------------
    # Step 3 — Convert & Export
    # ---------------------------------------------------------------------------

    st.divider()
    st.subheader("Step 3 — Convert probes & export GCT")

    # Show success + download if pipeline already ran this session
    if st.session_state.gct_content is not None:
        st.success("Pipeline complete — gene mapping ready.")
        st.download_button(
            label="Download GCT file",
            data=st.session_state.gct_content,
            file_name=f"{accession or 'output'}.gct",
            mime="text/plain",
        )
        if st.button("Re-run pipeline"):
            st.session_state.gene_mapping = None
            st.session_state.gct_content = None
            st.rerun()
    else:
        st.info(
            "Click **Run Pipeline** to query MyGene.info for gene symbols "
            "and generate the GCT file. Results are cached locally so repeat "
            "runs with the same probes are instant."
        )
        if st.button("Run Pipeline", type="primary"):
            # --- convert ---
            gene_mapping: dict[str, str] | None = None
            try:
                with st.spinner("Querying MyGene.info for gene symbols…"):
                    gene_mapping = convert_probes(list(matrix.probe_ids))
                st.session_state.gene_mapping = gene_mapping
            except RetryableError as exc:
                st.error(
                    f"Network error — MyGene.info could not be reached after "
                    f"all retries: {exc}\n\nCheck your internet connection and "
                    f"click Run Pipeline again."
                )
            except PipelineError as exc:
                st.error(f"Pipeline error during probe conversion: {exc}")
            except Exception as exc:
                st.error(f"Unexpected error during probe conversion: {exc}")

            # --- export (only if convert succeeded) ---
            if gene_mapping is not None:
                try:
                    buf = io.StringIO()
                    export_gct(matrix, gene_mapping, buf)
                    st.session_state.gct_content = buf.getvalue().encode("utf-8")
                    st.rerun()
                except PipelineError as exc:
                    st.error(f"Export error: {exc}")
                    st.session_state.gct_content = None
                except Exception as exc:
                    st.error(f"Unexpected error during export: {exc}")
                    st.session_state.gct_content = None
