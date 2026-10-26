
import pandas as pd
import streamlit as st

from src.id_resolver import parse_gene_ids, resolve_gene_ids

st.set_page_config(
    page_title="TP53 Locus-to-Function Explorer",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 TP53 Locus-to-Function Explorer")

st.markdown(
    """
    Explore human genes using Ensembl identifiers and build a foundation
    for integrating TP53 ChIP-seq and RNA-seq evidence.

    **Current stage:** Gene ID validation and annotation retrieval.
    """
)

st.divider()

# Sidebar
with st.sidebar:
    st.header("🔬 Workflow")
    st.markdown(
        """
        1. Enter Ensembl gene IDs
        2. Validate the IDs
        3. Retrieve gene annotations
        4. View and download results
        5. Future: integrate ChIP-seq evidence
        6. Future: integrate RNA-seq evidence
        """
    )

    st.divider()
    st.caption("Species: Homo sapiens")
    st.caption("Annotation source: Ensembl REST API")

# Input section
st.subheader("1. Enter your gene IDs")

st.write(
    "Enter one or more human Ensembl gene IDs, separated by new lines, "
    "spaces, commas, or semicolons."
)

example_ids = "ENSG00000141510\nENSG00000157764"

raw_ids = st.text_area(
    "Ensembl gene IDs",
    value=example_ids,
    height=120,
    help="Use human Ensembl gene IDs beginning with ENSG.",
)

st.caption(
    "Examples: ENSG00000141510 (TP53) and ENSG00000157764 (BRAF)."
)

run_button = st.button(
    "🚀 Run annotation workflow",
    type="primary",
    use_container_width=True,
)

if run_button:
    valid_ids, invalid_tokens = parse_gene_ids(raw_ids)

    if not valid_ids and not invalid_tokens:
        st.error("Please enter at least one Ensembl gene ID.")

    else:
        if invalid_tokens:
            st.warning(
                "These entries do not look like valid Ensembl gene IDs: "
                + ", ".join(invalid_tokens)
            )

        if valid_ids:
            with st.spinner(
                "Contacting Ensembl and retrieving gene annotations..."
            ):
                try:
                    results, unresolved = resolve_gene_ids(valid_ids)

                    st.session_state["annotation_results"] = results
                    st.session_state["annotation_errors"] = unresolved

                except RuntimeError as exc:
                    st.error(str(exc))

                except Exception as exc:
                    st.error(
                        "An unexpected error occurred while processing "
                        f"the IDs: {exc}"
                    )

# Results section
if "annotation_results" in st.session_state:
    results = st.session_state["annotation_results"]
    unresolved = st.session_state.get("annotation_errors", [])

    st.divider()
    st.subheader("2. Gene annotation results")

    col1, col2, col3 = st.columns(3)

    col1.metric("Resolved genes", len(results))
    col2.metric("Unresolved IDs", len(unresolved))

    assemblies = (
        sorted(results["assembly"].dropna().astype(str).unique())
        if not results.empty
        else []
    )

    col3.metric(
        "Genome assembly",
        ", ".join(assemblies) if assemblies else "Not available",
    )

    if not results.empty:
        st.success("Gene annotations retrieved successfully.")

        st.dataframe(
            results,
            use_container_width=True,
            hide_index=True,
        )

        csv_data = results.to_csv(index=False).encode("utf-8")

        st.download_button(
            label="⬇️ Download gene annotations (CSV)",
            data=csv_data,
            file_name="tp53_gene_annotations.csv",
            mime="text/csv",
            use_container_width=True,
        )

    else:
        st.info(
            "No gene IDs were resolved. Check the IDs or try again later."
        )

    if unresolved:
        st.subheader("3. IDs needing attention")

        unresolved_df = pd.DataFrame(unresolved)

        st.dataframe(
            unresolved_df,
            use_container_width=True,
            hide_index=True,
        )

        st.download_button(
            label="⬇️ Download unresolved IDs (CSV)",
            data=unresolved_df.to_csv(index=False).encode("utf-8"),
            file_name="unresolved_gene_ids.csv",
            mime="text/csv",
        )

    st.divider()
    st.subheader("🧪 Scientific interpretation")

    st.info(
        """
        Gene annotation provides information such as gene name, genomic
        location, and biotype. It does not establish that TP53 binds to a
        gene or changes its expression.

        Those questions will be investigated in later stages using
        documented public ChIP-seq and RNA-seq datasets.
        """
    )

# Project roadmap
with st.expander("📌 Planned future features"):
    st.markdown(
        """
        - **ChIP-seq:** identify genomic regions with TP53 binding evidence.
        - **RNA-seq:** examine gene expression changes under selected conditions.
        - **Motif analysis:** investigate potential TP53-binding motifs.
        - **Evidence integration:** prioritise candidate target genes.
        - **Gene reports:** export results with dataset references and methods.
        """
    )

st.caption(
    "Educational bioinformatics project. Interpret results in the context "
    "of the source datasets and their experimental conditions."
)
