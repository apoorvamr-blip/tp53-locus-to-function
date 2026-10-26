
"""Basic analysis functions for uploaded GEO-derived datasets."""

import re
import pandas as pd


def read_uploaded_table(uploaded_file):
    """Read a CSV or TSV file uploaded through Streamlit."""
    filename = uploaded_file.name.lower()

    if filename.endswith((".tsv", ".txt", ".bed", ".narrowpeak")):
        return pd.read_csv(
            uploaded_file,
            sep="\t",
            comment="#",
            header=None if filename.endswith((".bed", ".narrowpeak")) else "infer",
        )

    return pd.read_csv(uploaded_file)


def find_column(df, candidates):
    """Find a column using common alternative names."""
    lookup = {str(col).strip().lower(): col for col in df.columns}

    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]

    return None


def analyse_chipseq(peaks, annotations, promoter_bp=2000):
    """
    Count uploaded ChIP-seq peaks overlapping each selected gene's
    promoter region. BED coordinates are treated as zero-based,
    half-open intervals.
    """
    if peaks.shape[1] < 3:
        raise ValueError("The peak file needs at least chromosome, start, end.")

    peak_data = peaks.iloc[:, :3].copy()
    peak_data.columns = ["chromosome", "peak_start", "peak_end"]

    peak_data["chromosome"] = peak_data["chromosome"].astype(str)
    peak_data["peak_start"] = pd.to_numeric(
        peak_data["peak_start"], errors="coerce"
    )
    peak_data["peak_end"] = pd.to_numeric(
        peak_data["peak_end"], errors="coerce"
    )
    peak_data = peak_data.dropna(
        subset=["peak_start", "peak_end"]
    )

    results = []

    for _, gene in annotations.iterrows():
        chrom = str(gene.get("chromosome", ""))
        start = pd.to_numeric(gene.get("start"), errors="coerce")
        end = pd.to_numeric(gene.get("end"), errors="coerce")
        strand = gene.get("strand")

        if pd.isna(start) or pd.isna(end):
            continue

        # Ensembl coordinates are 1-based; BED intervals are 0-based.
        tss = int(start) - 1 if strand == 1 else int(end) - 1

        if strand == 1:
            promoter_start = max(0, tss - promoter_bp)
            promoter_end = tss + 1
        else:
            promoter_start = max(0, tss)
            promoter_end = tss + promoter_bp + 1

        matching = peak_data[
            (peak_data["chromosome"] == chrom)
            & (peak_data["peak_start"] < promoter_end)
            & (peak_data["peak_end"] > promoter_start)
        ]

        results.append({
            "ensembl_gene_id": gene.get("ensembl_gene_id"),
            "gene_name": gene.get("gene_name"),
            "promoter_peak_count": len(matching),
            "promoter_start_bed": promoter_start,
            "promoter_end_bed": promoter_end,
        })

    return pd.DataFrame(results)


def analyse_rnaseq(expression):
    """Standardise common differential-expression column names."""
    gene_col = find_column(
        expression,
        ["gene_id", "gene", "ensembl_gene_id", "Geneid", "ID"],
    )
    lfc_col = find_column(
        expression,
        ["log2FoldChange", "log2FC", "logFC", "log2_fold_change"],
    )
    padj_col = find_column(
        expression,
        ["padj", "FDR", "adj.P.Val", "adjusted_p_value", "qvalue"],
    )

    if not gene_col or not lfc_col or not padj_col:
        raise ValueError(
            "The RNA-seq table must contain gene IDs, a log2 fold-change "
            "column, and an adjusted p-value column (such as padj)."
        )

    result = expression[[gene_col, lfc_col, padj_col]].copy()
    result.columns = ["ensembl_gene_id", "log2FoldChange", "padj"]

    result["ensembl_gene_id"] = (
        result["ensembl_gene_id"].astype(str).str.strip()
    )
    result["ensembl_gene_id"] = result["ensembl_gene_id"].str.split(".").str[0]
    result["log2FoldChange"] = pd.to_numeric(
        result["log2FoldChange"], errors="coerce"
    )
    result["padj"] = pd.to_numeric(result["padj"], errors="coerce")

    result["significant"] = (
        (result["padj"] < 0.05)
        & (result["log2FoldChange"].abs() >= 1)
    )

    result["expression_direction"] = result["log2FoldChange"].map(
        lambda value: (
            "Upregulated" if value > 0
            else "Downregulated" if value < 0
            else "No change"
        ) if pd.notna(value) else "Unknown"
    )

    return result


def analyse_motifs(fasta_text):
    """
    Preliminary scan for a simplified TP53 response-element pattern.
    This is a consensus-pattern screen, not a full PWM analysis.
    FASTA headers should contain gene IDs beginning with ENSG.
    """
    # IUPAC consensus: R=A/G, W=A/T, Y=C/T.
    half_site = r"[AG][AG][AG]C[AT]W[AG][CT][CT]"

    # Two half-sites with 0-13 intervening nucleotides.
    motif = re.compile(half_site + r"[ACGTN]{0,13}" + half_site, re.I)
    reverse_motif = re.compile(
        half_site + r"[ACGTN]{0,13}" + half_site, re.I
    )

    records = re.findall(
        r">([^\n]+)\n([^>]*)",
        fasta_text,
        flags=re.MULTILINE,
    )

    results = []

    for header, sequence in records:
        match = re.search(r"ENSG\d+(?:\.\d+)?", header, re.I)
        if not match:
            continue

        gene_id = match.group(0).upper().split(".")[0]
        sequence = re.sub(r"\s+", "", sequence).upper()
        reverse_sequence = sequence.translate(
            str.maketrans("ACGTN", "TGCAN")
        )[::-1]

        forward_count = len(list(motif.finditer(sequence)))
        reverse_count = len(list(reverse_motif.finditer(reverse_sequence)))

        results.append({
            "ensembl_gene_id": gene_id,
            "sequence_length": len(sequence),
            "motif_count": forward_count + reverse_count,
            "motif_method": "Simplified TP53 consensus-pattern scan",
        })

    return pd.DataFrame(results)


def integrate_evidence(chip, rna, motifs=None):
    """Combine available evidence without claiming causal proof."""
    result = rna.copy()

    if chip is not None and not chip.empty:
        result = result.merge(
            chip,
            on="ensembl_gene_id",
            how="left",
        )
        result["promoter_peak_count"] = (
            result["promoter_peak_count"].fillna(0).astype(int)
        )
    else:
        result["promoter_peak_count"] = pd.NA

    if motifs is not None and not motifs.empty:
        result = result.merge(
            motifs[["ensembl_gene_id", "motif_count"]],
            on="ensembl_gene_id",
            how="left",
        )
    else:
        result["motif_count"] = pd.NA

    result["evidence_score"] = (
        result["significant"].fillna(False).astype(int)
        + (
            result["promoter_peak_count"].fillna(0).gt(0).astype(int)
            if result["promoter_peak_count"].notna().any()
            else 0
        )
        + (
            result["motif_count"].fillna(0).gt(0).astype(int)
            if result["motif_count"].notna().any()
            else 0
        )
    )

    result["interpretation"] = result["evidence_score"].map({
        0: "No detected evidence in the supplied analyses",
        1: "Evidence in one analysis",
        2: "Evidence in two analyses",
        3: "Evidence in all three analyses",
    })

    return result.sort_values(
        ["evidence_score", "padj"],
        ascending=[False, True],
        na_position="last",
    )
