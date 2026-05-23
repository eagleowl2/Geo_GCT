"""Streamlit web UI for the GEO Series Matrix → GCT pipeline."""

from __future__ import annotations

import gzip
import io
import json
import math
import sys
from pathlib import Path

# Ensure the local src/ tree is always used, even if an older installed
# version of gct-pipeline is cached in the Streamlit Cloud environment.
# We also pop any pre-imported gct_pipeline so Python re-resolves from src/.
_SRC_DIR = (Path(__file__).parent / "src").resolve()
sys.path.insert(0, str(_SRC_DIR))
for _stale in [m for m in sys.modules if m == "gct_pipeline" or m.startswith("gct_pipeline.")]:
    del sys.modules[_stale]

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

# Capture where gct_pipeline actually loaded from — surfaced in the sidebar
# diagnostics expander so we can verify the right code is running.
import gct_pipeline as _gp_mod
_GP_LOCATION = getattr(_gp_mod, "__file__", "unknown")

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="GEO → GCT Pipeline",
    page_icon="🧬",
    layout="wide",
)

# ---------------------------------------------------------------------------
# Session-state initialisation
# ---------------------------------------------------------------------------

for _key in ("matrix", "gene_mapping", "gct_content", "_file_id", "_platform"):
    if _key not in st.session_state:
        st.session_state[_key] = None


def _reset_pipeline() -> None:
    """Clear conversion results but keep the parsed matrix."""
    st.session_state.gene_mapping = None
    st.session_state.gct_content = None


def _reset_state() -> None:
    """Full reset — used when a new file is uploaded or on error."""
    st.session_state.matrix = None
    st.session_state.gene_mapping = None
    st.session_state.gct_content = None
    st.session_state["_file_id"] = None


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

_PLATFORMS: dict[str, str] = {
    "hgu133plus2":  "Human — HG-U133 Plus 2.0",
    "hgu133a":      "Human — HG-U133A",
    "hgu133b":      "Human — HG-U133B",
    "hgu95av2":     "Human — HG-U95Av2",
    "mouse430_2":   "Mouse — Mouse430 2.0",
    "mouse4302":    "Mouse — MOE430 2.0",
    "rat2302":      "Rat — RAE230 2.0",
}

_DEFAULT_CACHE = Path.home() / ".gct_pipeline" / "probe_cache.json"
_APP_VERSION = "0.5.0"

with st.sidebar:
    st.markdown("## 🧬 GEO → GCT Pipeline")
    st.markdown(
        "**Made by [Alexander Rusinov](https://github.com/eagleowl2)**  \n"
        "Bioinformatics data pipeline · GEO Series Matrix → GCT 1.2"
    )
    st.caption(f"v{_APP_VERSION}")
    st.caption(
        "Powered by [MyGene.info](https://mygene.info) · "
        "[Source on GitHub](https://github.com/eagleowl2/Geo_GCT)"
    )

    st.divider()

    # --- Platform selector ---
    st.markdown("### ⚙️ Settings")
    platform_key = st.selectbox(
        "Affymetrix platform",
        options=list(_PLATFORMS.keys()),
        format_func=lambda k: _PLATFORMS[k],
        help=(
            "Select the microarray platform used to generate your data. "
            "This is passed to MyGene.info as a scope hint to improve mapping accuracy."
        ),
    )

    # If platform changed mid-session, invalidate previous conversion result
    if st.session_state["_platform"] != platform_key:
        if st.session_state["_platform"] is not None:
            _reset_pipeline()
        st.session_state["_platform"] = platform_key

    st.divider()

    # --- Cache management ---
    st.markdown("### 🗄️ Probe cache")
    if _DEFAULT_CACHE.exists():
        try:
            _cached: dict[str, str] = json.loads(
                _DEFAULT_CACHE.read_text(encoding="utf-8")
            )
            n_cached = len(_cached)
        except (json.JSONDecodeError, OSError):
            n_cached = 0
        st.metric("Cached probes", f"{n_cached:,}")
        st.caption(f"`{_DEFAULT_CACHE}`")
        if st.button("🗑️ Clear cache", help="Remove all locally cached probe→gene mappings"):
            _DEFAULT_CACHE.unlink(missing_ok=True)
            _reset_pipeline()
            st.success("Cache cleared.")
            st.rerun()
    else:
        st.caption("Cache is empty — results are stored here after the first run.")

    st.divider()

    # --- About ---
    with st.expander("ℹ️ About this app"):
        st.markdown(
            """
            This pipeline converts raw microarray expression data from
            **NCBI GEO** into the **GCT 1.2** format used by tools such as
            GSEA and GenePattern.

            **Pipeline steps:**
            1. Parse the GEO Series Matrix file (streaming, numpy-backed)
            2. Map Affymetrix probe IDs → gene symbols via MyGene.info
               (batched 500/call, locally cached, exponential backoff retry)
            3. Export GCT 1.2 for download

            **Source code & docs:** [github.com/eagleowl2/Geo_GCT](https://github.com/eagleowl2/Geo_GCT)
            """
        )

    with st.expander("🔧 Runtime diagnostics"):
        st.code(
            f"app version : {_APP_VERSION}\n"
            f"python      : {sys.version.split()[0]}\n"
            f"gct_pipeline: {_GP_LOCATION}\n"
            f"src/ exists : {_SRC_DIR.exists()}\n"
            f"src path    : {_SRC_DIR}",
            language="text",
        )
        st.caption(
            "If `gct_pipeline` does not point to a path ending in "
            "`src/gct_pipeline/__init__.py`, an old installed copy is "
            "shadowing the source tree."
        )

# ---------------------------------------------------------------------------
# Main header
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
            with st.spinner("Parsing file…"):
                try:
                    st.session_state.matrix = parse(text_stream)
                except ParseError as exc:
                    st.error(f"Parse error: {exc}")
                    _reset_state()
                except Exception as exc:
                    st.error(f"Unexpected error during parsing: {exc}")
                    _reset_state()

# ---------------------------------------------------------------------------
# Step 2 — Series Preview
# ---------------------------------------------------------------------------

if st.session_state.matrix is not None:
    matrix = st.session_state.matrix
    st.divider()
    st.subheader("Step 2 — Series Preview")

    accession = matrix.metadata.get("Series_geo_accession", "—")
    title = matrix.metadata.get("Series_title", "—")
    # Real GEO Series Matrix files expose organism via Sample_organism_ch1
    # (per-sample). Series_sample_organism is rare — try both, plus taxid fallback.
    _TAXID_TO_ORGANISM = {"9606": "Homo sapiens", "10090": "Mus musculus", "10116": "Rattus norvegicus"}
    organism = (
        matrix.metadata.get("Sample_organism_ch1")
        or matrix.metadata.get("Series_sample_organism")
        or _TAXID_TO_ORGANISM.get(matrix.metadata.get("Series_sample_taxid", ""), "—")
    )

    # Top metrics row
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Probes", f"{len(matrix.probe_ids):,}")
    m2.metric("Samples", f"{len(matrix.sample_ids):,}")
    m3.metric("Accession", accession)
    m4.metric("Organism", organism)

    if title != "—":
        st.caption(f"**Series title:** {title}")

    # Full metadata expander
    with st.expander("📋 Full metadata", expanded=False):
        meta_df = pd.DataFrame(
            list(matrix.metadata.items()), columns=["Key", "Value"]
        )
        st.dataframe(meta_df, use_container_width=True, hide_index=True)

    # Expression preview
    st.subheader("Expression matrix preview (first 50 rows)")
    n_preview = min(50, len(matrix.probe_ids))
    has_titles = len(matrix.sample_titles) == len(matrix.sample_ids) and len(matrix.sample_ids) > 0

    if n_preview == 0:
        st.warning(
            "**No inline expression data found in this file.**  \n"
            "This GEO Series Matrix file contains sample metadata and "
            f"**{len(matrix.sample_ids):,} sample IDs** but no expression table rows.  \n\n"
            "This is common for datasets that store expression values in separate "
            "supplementary files (e.g. NanoString GeoMx, 10x Visium, bulk RNA-seq counts).  \n"
            f"Look for supplementary files on the "
            f"[GEO page for {accession}](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}).",
            icon="ℹ️",
        )
    else:
        # Column labels: prefer sample titles when available
        if has_titles:
            use_titles = st.toggle(
                "Show sample titles instead of GSM IDs",
                value=True,
                help="Switches column headers between human-readable titles "
                     "(e.g. 'BMDM, untreated, 1') and GEO accession IDs.",
            )
            col_labels = list(matrix.sample_titles if use_titles else matrix.sample_ids)
        else:
            col_labels = list(matrix.sample_ids)

        preview_df = pd.DataFrame(
            matrix.expression[:n_preview],
            index=list(matrix.probe_ids[:n_preview]),
            columns=col_labels,
        )
        preview_df.index.name = "Probe ID"
        st.dataframe(preview_df, use_container_width=True)

    # Sample title reference table (show for all files that have titles)
    if has_titles:
        with st.expander("🏷️ Sample ID ↔ title mapping", expanded=False):
            st.dataframe(
                pd.DataFrame(
                    {"GSM ID": list(matrix.sample_ids), "Title": list(matrix.sample_titles)}
                ),
                use_container_width=True,
                hide_index=True,
            )

    # ---------------------------------------------------------------------------
    # Step 3 — Convert & Export
    # ---------------------------------------------------------------------------

    st.divider()
    st.subheader("Step 3 — Convert probes & export GCT")

    if len(matrix.probe_ids) == 0:
        st.info(
            "**Run Pipeline is unavailable** — this file has no inline expression rows to convert.  \n"
            "This is common for spatial or sequencing datasets that store counts in separate "
            "supplementary files.  \n"
            f"Find them on the [GEO page for {accession}]"
            f"(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}).",
            icon="🚫",
        )
    elif st.session_state.gct_content is not None and st.session_state.gene_mapping is not None:
        gene_mapping = st.session_state.gene_mapping

        # Mapping summary
        mapped   = {k: v for k, v in gene_mapping.items() if v}
        unmapped = [k for k, v in gene_mapping.items() if not v]
        rate     = 100 * len(mapped) / len(gene_mapping) if gene_mapping else 0.0

        st.success(
            f"Pipeline complete — {len(mapped):,} / {len(gene_mapping):,} probes "
            f"mapped to gene symbols ({rate:.1f}%)."
        )

        s1, s2, s3 = st.columns(3)
        s1.metric("Mapped probes",   f"{len(mapped):,}")
        s2.metric("Unmapped probes", f"{len(unmapped):,}")
        s3.metric("Mapping rate",    f"{rate:.1f}%")

        # Gene mapping preview
        with st.expander("🔬 Gene mapping preview (first 50 probes)", expanded=False):
            preview_items = list(gene_mapping.items())[:50]
            mapping_df = pd.DataFrame(preview_items, columns=["Probe ID", "Gene Symbol"])
            mapping_df["Status"] = mapping_df["Gene Symbol"].apply(
                lambda s: "✅ Mapped" if s else "⚠️ Unmapped"
            )
            st.dataframe(mapping_df, use_container_width=True, hide_index=True)

        # Unmapped probes
        if unmapped:
            with st.expander(f"⚠️ Unmapped probes ({len(unmapped):,})", expanded=False):
                st.caption(
                    "These probe IDs were not found in the MyGene.info database "
                    "for the selected platform. They appear in the GCT file with "
                    "an empty Description field."
                )
                st.dataframe(
                    pd.DataFrame(unmapped, columns=["Probe ID"]),
                    use_container_width=True,
                    hide_index=True,
                )

        st.download_button(
            label="⬇️ Download GCT file",
            data=st.session_state.gct_content,
            file_name=f"{accession}.gct" if accession != "—" else "output.gct",
            mime="text/plain",
            type="primary",
        )

        if st.button("🔄 Re-run with current settings"):
            _reset_pipeline()
            st.rerun()

    else:
        # Calculate expected batch count for the info message
        n_probes = len(matrix.probe_ids)
        n_batches = max(1, math.ceil(n_probes / 500))
        platform_label = _PLATFORMS.get(platform_key, platform_key)

        st.info(
            f"Click **Run Pipeline** to map **{n_probes:,} probes** to gene symbols "
            f"using MyGene.info ({n_batches} batch{'es' if n_batches > 1 else ''} of ≤500 probes).  \n"
            f"Platform: **{platform_label}** · Results cached locally for instant repeat runs."
        )

        if st.button("▶️ Run Pipeline", type="primary"):
            gene_mapping_result: dict[str, str] | None = None
            try:
                with st.spinner(
                    f"Querying MyGene.info for {n_probes:,} probes "
                    f"({n_batches} batch{'es' if n_batches > 1 else ''})…"
                ):
                    gene_mapping_result = convert_probes(
                        list(matrix.probe_ids),
                        platform=platform_key,
                    )
                st.session_state.gene_mapping = gene_mapping_result
            except RetryableError as exc:
                st.error(
                    f"Network error — MyGene.info could not be reached after all retries: {exc}  \n"
                    "Check your internet connection and click Run Pipeline again."
                )
            except PipelineError as exc:
                st.error(f"Pipeline error during probe conversion: {exc}")
            except Exception as exc:
                st.error(f"Unexpected error during probe conversion: {exc}")

            if gene_mapping_result is not None:
                try:
                    buf = io.StringIO()
                    export_gct(matrix, gene_mapping_result, buf)
                    st.session_state.gct_content = buf.getvalue().encode("utf-8")
                    st.rerun()
                except PipelineError as exc:
                    st.error(f"Export error: {exc}")
                    st.session_state.gct_content = None
                except Exception as exc:
                    st.error(f"Unexpected error during export: {exc}")
                    st.session_state.gct_content = None

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------

st.divider()
st.caption(
    "🧬 GEO → GCT Pipeline · Made by **Alexander Rusinov** · "
    "[GitHub](https://github.com/eagleowl2/Geo_GCT) · "
    "Powered by [MyGene.info](https://mygene.info) & "
    "[Streamlit](https://streamlit.io)"
)
