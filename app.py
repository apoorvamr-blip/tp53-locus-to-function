
import pandas as pd
import streamlit as st

from src.id_resolver import parse_gene_ids, resolve_gene_ids
from src.analysis_tools import (
    read_uploaded_table,
    analyse_chipseq,
    analyse_rnaseq,
    analyse_motifs,
    integrate_evidence,
)

st.set_page_config(
    page_title="TP53 Locus-to-Function Explorer",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 TP53 Locus-to-Function Explorer")
st.write(
    "Explore gene annotations and analyse uploaded GEO-derived "
    "ChIP-seq, RNA-seq, and DNA sequence data."
)

st.info(
    "This version analyses processed files you upload. It does not "
    "download GEO datasets or process raw FASTQ files automatically."
)

with st.sidebar:
    st.header("Analysis workflow")
    st.markdown(
        """
        1. Enter human Ensembl gene IDs
        2. Retrieve gene annotations
        3. Upload ChIP-seq peaks
        4. Upload RNA-seq results
        5. Optionally scan FASTA sequences
        6. Integrate available evidence
        7. Download results
        """
    )
    st.caption("Species: Homo sapiens")

st.subheader("1. Select genes")

raw_ids = st.text_area(
    "Human Ensembl gene IDs",
    value="ENSG00000141510\nENSG00000157764",
    height=110,
    help="Separate IDs with commas, spaces, semicolons, or new lines.",
)

st.caption(
    "Example: ENSG00000141510 (TP53). "
    "Enter the genes you want to investigate."
)

st.subheader("2. Upload your GEO-derived data")

chip_file = st.file_uploader(
    "ChIP-seq peak file",
    type=["bed", "narrowpeak", "tsv", "txt"],
    help=(
        "Use a BED or narrowPeak file with chromosome, start, and end "
        "in the first three columns."
    ),
)

rna_file = st.file_uploader(
    "RNA-seq differential-expression table",
    type=["csv", "tsv", "txt"],
    help=(
        "Needs gene IDs, log2FoldChange (or equivalent), "
        "and padj (or equivalent adjusted p-value)."
    ),
)

fasta_file = st.file_uploader(
    "Optional DNA sequences in FASTA format",
    type=["fa", "fasta", "fna"],
    help=(
        "FASTA headers should contain Ensembl gene IDs. "
        "Use sequences relevant to the genes being studied."
    ),
)

promoter_bp = st.number_input(
    "Promoter window upstream of the transcription start site (bp)",
    min_value=100,
    max_value=10000,
    value=2000,
    step=100,
)

run = st.button(
    "🚀 Run analysis",
    type="primary",
    use_container_width=True,
)

if run:
    gene_ids, invalid_tokens = parse_gene_ids(raw_ids)

    if not gene_ids:
        st.error("Enter at least one valid-looking Ensembl gene ID.")

    elif not any([chip_file, rna_file, fasta_file]):
        st.error("Upload at least one analysis file.")

    else:
        if invalid_tokens:
            st.warning(
                "Ignored entries that do not look like Ensembl gene IDs: "
                + ", ".join(invalid_tokens)
            )

        try:
            with st.spinner("Retrieving Ensembl gene annotations..."):
                annotations, unresolved = resolve_gene_ids(gene_ids)

            if annotations.empty:
                st.error(
                    "None of the entered IDs could be resolved by Ensembl."
                )
                st.session_state.pop("analysis_output", None)
                st.stop()

            selected_ids = set(
                annotations["ensembl_gene_id"].astype(str).str.split(".").str[0]
            )

            chip_results = None
            rna_results = None
            motif_results = None
            integrated = None
            errors = []

            if chip_file is not None:
                with st.spinner("Analysing ChIP-seq peak overlaps..."):
                    peaks = read_uploaded_table(chip_file)
                    chip_results = analyse_chipseq(
                        peaks,
                        annotations,
                        promoter_bp=int(promoter_bp),
                    )

            if rna_file is not None:
                with st.spinner("Analysing RNA-seq results..."):
                    expression = read_uploaded_table(rna_file)
                    rna_results = analyse_rnaseq(expression)

                    # Keep only the genes selected by the user.
                    rna_results = rna_results[
                        rna_results["ensembl_gene_id"].isin(selected_ids)
                    ].copy()

                    if rna_results.empty:
                        errors.append(
                            "No RNA-seq rows matched the selected Ensembl "
                            "gene IDs. Check the identifier column and IDs."
                        )

            if fasta_file is not None:
                with st.spinner("Scanning FASTA sequences for motifs..."):
                    fasta_text = fasta_file.getvalue().decode(
                        "utf-8", errors="replace"
                    )
                    motif_results = analyse_motifs(fasta_text)

                    if not motif_results.empty:
                        motif_results = motif_results[
                            motif_results["ensembl_gene_id"].isin(selected_ids)
                        ].copy()

            if rna_results is not None and not rna_results.empty:
                integrated = integrate_evidence(
                    chip_results,
                    rna_results,
                    motif_results,
                )

                integrated = integrated.merge(
                    annotations[
                        [
                            "ensembl_gene_id",
                            "gene_name",
                            "description",
                            "chromosome",
                            "start",
                            "end",
                            "strand",
                            "assembly",
                        ]
                    ],
                    on="ensembl_gene_id",
                    how="left",
                    suffixes=("", "_annotation"),
                )

            st.session_state["analysis_output"] = {
                "annotations": annotations,
                "unresolved": unresolved,
                "chip": chip_results,
                "rna": rna_results,
                "motifs": motif_results,
                "integrated": integrated,
                "errors": errors,
            }

        except Exception as exc:
            st.error(f"Analysis failed: {exc}")
            st.info(
                "Check your file format, required columns, gene IDs, "
                "and whether the uploaded files are processed data."
            )

# Display saved results
output = st.session_state.get("analysis_output")

if output:
    annotations = output["annotations"]
    unresolved = output["unresolved"]
    chip_results = output["chip"]
    rna_results = output["rna"]
    motif_results = output["motifs"]
    integrated = output["integrated"]

    st.divider()
    st.header("3. Results")

    for message in output["errors"]:
        st.warning(message)

    metric_cols = st.columns(3)
    metric_cols[0].metric("Annotated genes", len(annotations))
    metric_cols[1].metric(
        "Genes with promoter peaks",
        (
            int((chip_results["promoter_peak_count"] > 0).sum())
            if chip_results is not None
            else 0
        ),
    )
    metric_cols[2].metric(
        "Significant RNA-seq results",
        (
            int(rna_results["significant"].sum())
            if rna_results is not None
            else 0
        ),
    )

    with st.expander("Gene annotations", expanded=True):
        st.dataframe(
            annotations,
            use_container_width=True,
            hide_index=True,
        )
        st.download_button(
            "Download annotations CSV",
            annotations.to_csv(index=False).encode("utf-8"),
            "gene_annotations.csv",
            "text/csv",
        )

    if chip_results is not None:
        with st.expander("ChIP-seq promoter overlap results", expanded=True):
            st.dataframe(
                chip_results,
                use_container_width=True,
                hide_index=True,
            )
            st.download_button(
                "Download ChIP-seq results CSV",
                chip_results.to_csv(index=False).encode("utf-8"),
                "chipseq_promoter_results.csv",
                "text/csv",
            )

    if rna_results is not None:
        with st.expander("RNA-seq results", expanded=True):
            st.dataframe(
                rna_results,
                use_container_width=True,
                hide_index=True,
            )
            st.download_button(
                "Download RNA-seq results CSV",
                rna_results.to_csv(index=False).encode("utf-8"),
                "rnaseq_results.csv",
                "text/csv",
            )

    if motif_results is not None:
        with st.expander("Motif scan results", expanded=True):
            if motif_results.empty:
                st.warning(
                    "No matching FASTA records were found for your selected "
                    "genes. Check the FASTA headers and Ensembl IDs."
                )
            else:
                st.dataframe(
                    motif_results,
                    use_container_width=True,
                    hide_index=True,
                )
                st.download_button(
                    "Download motif results CSV",
                    motif_results.to_csv(index=False).encode("utf-8"),
                    "motif_results.csv",
                    "text/csv",
                )

    if integrated is not None and not integrated.empty:
        st.subheader("4. Integrated evidence")
        st.write(
            "Genes are ranked by a simple tally of available evidence. "
            "This is exploratory prioritisation, not a validated predictor."
        )

        st.dataframe(
            integrated,
            use_container_width=True,
            hide_index=True,
        )

        st.download_button(
            "⬇️ Download integrated gene report (CSV)",
            integrated.to_csv(index=False).encode("utf-8"),
            "tp53_integrated_gene_report.csv",
            "text/csv",
            type="primary",
        )

    if unresolved:
        st.subheader("Unresolved gene IDs")
        unresolved_df = pd.DataFrame(unresolved)
        st.dataframe(
            unresolved_df,
            use_container_width=True,
            hide_index=True,
        )

    st.divider()
    st.subheader("5. Interpretation and limitations")

    st.markdown(
        """
        - A promoter peak overlap indicates genomic proximity, not proof
          that TP53 directly regulates the gene.
        - RNA-seq changes depend on the experimental comparison and
          do not by themselves establish TP53-specific effects.
        - Motif matches are preliminary sequence evidence, not confirmed
          binding sites.
        - Interpret integrated evidence only after checking dataset
          accession, genome assembly, cell type, treatment, controls,
          and experimental conditions.
        - Missing evidence must not be interpreted as evidence of absence.
        """
    )

    st.caption(
        "Educational research tool. Verify dataset provenance and methods "
        "before drawing biological conclusions."
    )

else:
    st.info(
        "Upload one or more processed files and run the workflow to "
        "view your results here."
    )
