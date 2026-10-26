
"""Resolve human Ensembl gene IDs using the official Ensembl REST API."""

import re
import time
import pandas as pd
import requests

ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/id"
ENSEMBL_GENE_ID_PATTERN = re.compile(r"^ENSG\d+(?:\.\d+)?$", re.IGNORECASE)
MAX_BATCH_SIZE = 500


def parse_gene_ids(raw_text: str) -> tuple[list[str], list[str]]:
    """Separate valid-looking Ensembl gene IDs from invalid input."""
    tokens = re.split(r"[\s,;]+", raw_text.strip()) if raw_text.strip() else []
    valid = []
    invalid = []
    seen = set()

    for token in tokens:
        if not token:
            continue

        token = token.strip()

        if ENSEMBL_GENE_ID_PATTERN.fullmatch(token):
            lookup_id = token.split(".", 1)[0].upper()

            if lookup_id not in seen:
                valid.append(token.upper())
                seen.add(lookup_id)
        else:
            invalid.append(token)

    return valid, invalid


def _lookup_batch(ids: list[str], session: requests.Session) -> dict:
    """Send a batch of gene IDs to Ensembl."""
    response = session.post(
        ENSEMBL_LOOKUP_URL,
        params={"species": "homo_sapiens"},
        json={"ids": ids},
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
        timeout=30,
    )

    response.raise_for_status()
    payload = response.json()

    if not isinstance(payload, dict):
        raise ValueError("Ensembl returned an unexpected response format.")

    return payload


def resolve_gene_ids(gene_ids: list[str]) -> tuple[pd.DataFrame, list[dict]]:
    """Resolve gene IDs and return annotations plus unresolved-ID errors."""
    unique_ids = []
    original_ids = {}
    errors = []

    for raw_id in gene_ids:
        original = str(raw_id).strip().upper()
        lookup_id = original.split(".", 1)[0]

        if not ENSEMBL_GENE_ID_PATTERN.fullmatch(original):
            errors.append({
                "input_id": original,
                "status": "Invalid Ensembl gene ID",
            })
            continue

        if lookup_id not in original_ids:
            unique_ids.append(lookup_id)
            original_ids[lookup_id] = original

    rows = []

    with requests.Session() as session:
        for start in range(0, len(unique_ids), MAX_BATCH_SIZE):
            batch = unique_ids[start:start + MAX_BATCH_SIZE]

            try:
                payload = _lookup_batch(batch, session)
            except requests.RequestException as exc:
                raise RuntimeError(
                    "Could not contact Ensembl. Check your internet "
                    f"connection and try again. Details: {exc}"
                ) from exc
            except ValueError as exc:
                raise RuntimeError(str(exc)) from exc

            for gene_id in batch:
                record = payload.get(gene_id)

                if not record:
                    errors.append({
                        "input_id": original_ids[gene_id],
                        "status": "Not found by Ensembl",
                    })
                    continue

                if record.get("object_type") != "Gene":
                    errors.append({
                        "input_id": original_ids[gene_id],
                        "status": "Resolved, but not a gene",
                    })
                    continue

                rows.append({
                    "input_id": original_ids[gene_id],
                    "ensembl_gene_id": record.get("id", gene_id),
                    "gene_name": record.get("display_name"),
                    "description": record.get("description"),
                    "biotype": record.get("biotype"),
                    "species": record.get("species"),
                    "chromosome": record.get("seq_region_name"),
                    "start": record.get("start"),
                    "end": record.get("end"),
                    "strand": record.get("strand"),
                    "assembly": record.get("assembly_name"),
                    "annotation_source": record.get("source"),
                    "annotation_status": "Resolved",
                })

            if start + MAX_BATCH_SIZE < len(unique_ids):
                time.sleep(0.1)

    columns = [
        "input_id",
        "ensembl_gene_id",
        "gene_name",
        "description",
        "biotype",
        "species",
        "chromosome",
        "start",
        "end",
        "strand",
        "assembly",
        "annotation_source",
        "annotation_status",
    ]

    return pd.DataFrame(rows, columns=columns), errors
