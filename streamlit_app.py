#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import os
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import SDWriter
from tqdm import tqdm


PAGE_TITLE = "SDF → CSV (SMILES)"
PAGE_ICON = "Santu Chall -> "

st.set_page_config(page_title=PAGE_TITLE, page_icon=PAGE_ICON, layout="wide")


def _initialize_defaults() -> None:
    """Populate Streamlit session defaults exactly once per session."""

    if "dataset_id_default" not in st.session_state:
        st.session_state["dataset_id_default"] = datetime.utcnow().strftime(
            "RUN-%Y%m%d-%H%M%S"
        )


def _build_alcoa_metadata(
    *,
    operator_name: str,
    contact_info: str,
    purpose: str,
    storage_plan: str,
    dataset_id: str,
    file_hash: str,
) -> Dict[str, str]:
    """Collect UI-provided metadata when ALCOA+ is enabled."""

    return {
        "operator": operator_name.strip(),
        "contact": contact_info.strip(),
        "purpose": purpose.strip(),
        "storage_plan": storage_plan.strip(),
        "dataset_id": dataset_id.strip(),
        "file_hash": file_hash,
        "processing_label": "streamlit_app_v1",
    }


def _build_alcoa_columns(
    *,
    run_timestamp_utc: str,
    source_file: str,
    metadata: Dict[str, str] | None,
) -> Dict[str, str]:
    """Return per-row ALCOA+ metadata columns."""

    metadata = metadata or {}
    dataset_identifier = metadata.get("dataset_id") or f"{source_file}::{run_timestamp_utc}"

    return {
        "alcoa_attributable_operator": metadata.get("operator", ""),
        "alcoa_legible_purpose": metadata.get("purpose", ""),
        "alcoa_contemporaneous_timestamp_utc": run_timestamp_utc,
        "alcoa_original_source_file": source_file,
        "alcoa_accurate_input_sha256": metadata.get("file_hash", ""),
        "alcoa_complete_dataset_id": dataset_identifier,
        "alcoa_consistent_processing_label": metadata.get(
            "processing_label", "sdf_2_smiles_v1"
        ),
        "alcoa_enduring_storage_plan": metadata.get("storage_plan", ""),
        "alcoa_available_contact": metadata.get("contact", ""),
    }


def count_molecules_in_sdf(sdf_path: str) -> int:
    """Count records by scanning for '$$$$' separators."""

    count = 0
    with open(sdf_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.strip() == "$$$$":
                count += 1
    return count


def process_sdf_records(
    *,
    sdf_path: str,
    bad_sdf_path: str | None = None,
    enforce_alcoa: bool = False,
    alcoa_metadata: Dict[str, str] | None = None,
    run_timestamp_utc: str | None = None,
    show_progress: bool = False,
) -> Tuple[List[Dict[str, str]], Dict[str, Any]]:
    """Convert SDF molecules into Python dict rows plus a summary payload."""

    if not os.path.isfile(sdf_path):
        raise FileNotFoundError(f"Input SDF not found: {sdf_path}")

    total_mols = count_molecules_in_sdf(sdf_path)
    run_timestamp_utc = run_timestamp_utc or (
        datetime.utcnow().isoformat(timespec="seconds") + "Z"
    )

    rows: List[Dict[str, str]] = []
    n_total = 0
    n_parsed = 0
    n_smiles_ok = 0
    n_parse_fail = 0
    n_smiles_fail = 0

    bad_writer = SDWriter(bad_sdf_path) if bad_sdf_path else None

    with open(sdf_path, "rb") as sdf_file:
        supplier = Chem.ForwardSDMolSupplier(
            sdf_file,
            sanitize=True,
            removeHs=False,
        )

        iterator = supplier
        if show_progress:
            iterator = tqdm(
                supplier,
                total=total_mols,
                desc="Converting SDF",
                unit="mol",
            )

        for idx, mol in enumerate(iterator, start=1):
            n_total += 1

            if mol is None:
                n_parse_fail += 1
                continue

            n_parsed += 1
            try:
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            except Exception:
                n_smiles_fail += 1
                if bad_writer:
                    bad_writer.write(mol)
                continue

            n_smiles_ok += 1
            props = mol.GetPropsAsDict(includePrivate=False, includeComputed=False)

            row: Dict[str, str] = {
                "record_index": str(idx),
                "smiles": smiles,
                "source_file": os.path.basename(sdf_path),
                "processing_timestamp_utc": run_timestamp_utc,
            }

            try:
                name = mol.GetProp("_Name")
            except KeyError:
                name = ""
            row["mol_name"] = name

            for key, value in props.items():
                key_str = str(key)
                if key_str in row:
                    key_str = f"prop_{key_str}"
                row[key_str] = str(value)

            if enforce_alcoa:
                row.update(
                    _build_alcoa_columns(
                        run_timestamp_utc=run_timestamp_utc,
                        source_file=os.path.basename(sdf_path),
                        metadata=alcoa_metadata,
                    )
                )

            rows.append(row)

    if bad_writer:
        bad_writer.close()

    summary: Dict[str, Any] = {
        "input_sdf": os.path.abspath(sdf_path),
        "run_timestamp_utc": run_timestamp_utc,
        "counts": {
            "total_records_seen": n_total,
            "total_records_expected_from_separators": total_mols,
            "parsed_ok": n_parsed,
            "smiles_converted_ok": n_smiles_ok,
            "parse_failures": n_parse_fail,
            "smiles_failures": n_smiles_fail,
        },
    }

    return rows, summary


def main() -> None:
    _initialize_defaults()

    st.title("ALCOA+ Ready Molecular Converter")
    st.caption(
        "Upload an SDF or SD file containing one or many small molecules and "
        "instantly generate a clean CSV with canonical SMILES plus the original "
        "SD tags."
    )

    info_col, checkbox_col = st.columns([3, 1])
    with info_col:
        st.markdown(
            "**Workflow**: upload → optional ALCOA+ annotation → convert → download `CSV`."
        )
    with checkbox_col:
        enable_alcoa = st.checkbox(
            "ALCOA+ mode",
            value=True,
            help="Adds attributable, contemporaneous, and enduring metadata columns",
        )

    uploaded_file = st.file_uploader(
        "Drop an SDF / SD file",
        type=["sdf", "sd"],
        help="Gzip/zip archives are not supported. Size limit is whatever your browser can handle.",
    )

    if uploaded_file:
        file_details = f"**Selected:** `{uploaded_file.name}` · {uploaded_file.size / 1024:.1f} kB"
        st.markdown(file_details)

    alcoa_details: Optional[Dict[str, str]] = None

    if enable_alcoa:
        st.success(
            "ALCOA+ tracking activated — every record will contain attributable, "
            "legible, contemporaneous, consistent, complete, enduring, and available metadata."
        )
        with st.expander("Provide ALCOA+ annotations", expanded=True):
            operator_name = st.text_input(
                "Operator / Analyst name",
                placeholder="Dr. Jane Doe",
            )
            contact_info = st.text_input(
                "Contact / Email",
                placeholder="jane.doe@lab.org",
            )
            purpose = st.text_area(
                "Purpose or study context",
                placeholder="Stability screening for candidate ligands",
            )
            dataset_id = st.text_input(
                "Dataset identifier",
                value=st.session_state["dataset_id_default"],
                help="Custom code tying this export back to your lab notebook or ELN entry.",
            )
            storage_plan = st.selectbox(
                "Planned long-term storage",
                (
                    "21 CFR Part 11 compliant document vault",
                    "Validated data lake",
                    "Regulated LIMS",
                    "Local secure drive (with routine backups)",
                ),
            )

        required_fields = [operator_name, contact_info, purpose]
        missing_required = any(not field.strip() for field in required_fields)
    else:
        missing_required = False

    st.divider()
    convert_btn = st.button(
        "Convert to CSV",
        type="primary",
        use_container_width=True,
        disabled=uploaded_file is None,
    )

    if convert_btn:
        if not uploaded_file:
            st.error("Please upload an SDF / SD file before converting.")
            return

        if enable_alcoa and missing_required:
            st.error("Operator name, contact, and purpose are required for ALCOA+ mode.")
            return

        file_bytes = uploaded_file.getvalue()
        if not file_bytes:
            st.error("The uploaded file seems empty.")
            return

        file_hash = hashlib.sha256(file_bytes).hexdigest()

        if enable_alcoa:
            alcoa_details = _build_alcoa_metadata(
                operator_name=operator_name,
                contact_info=contact_info,
                purpose=purpose,
                storage_plan=storage_plan,
                dataset_id=dataset_id,
                file_hash=file_hash,
            )

        suffix = Path(uploaded_file.name).suffix or ".sdf"
        with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
            tmp.write(file_bytes)
            tmp_path = tmp.name

        with st.spinner("Converting molecules with RDKit..."):
            try:
                rows, summary = process_sdf_records(
                    sdf_path=tmp_path,
                    bad_sdf_path=None,
                    enforce_alcoa=enable_alcoa,
                    alcoa_metadata=alcoa_details,
                    show_progress=False,
                )
            finally:
                os.unlink(tmp_path)

        if not rows:
            st.warning(
                "No molecules were successfully converted. Check the input file or "
                "try running it through the CLI for detailed logs."
            )
            return

        df = pd.DataFrame(rows)

        counts = summary["counts"]
        col_a, col_b, col_c, col_d = st.columns(4)
        col_a.metric("Total records (separators)", counts["total_records_expected_from_separators"])
        col_b.metric("Read successfully", counts["parsed_ok"])
        col_c.metric("SMILES generated", counts["smiles_converted_ok"])
        col_d.metric("Conversion issues", counts["parse_failures"] + counts["smiles_failures"])

        st.success(
            f"Completed conversion — {counts['smiles_converted_ok']} structures exported to CSV."
        )

        preview_rows = min(len(df), 500)
        st.markdown(f"**Preview ({preview_rows} of {len(df)} rows)**")
        st.dataframe(df.head(preview_rows), use_container_width=True, height=420)
        if len(df) > preview_rows:
            st.caption("Preview truncated for performance. Download for the full dataset.")

        csv_bytes = df.to_csv(index=False).encode("utf-8")
        default_name = f"{Path(uploaded_file.name).stem or 'converted'}_smiles.csv"
        st.download_button(
            "Download CSV",
            data=csv_bytes,
            file_name=default_name,
            mime="text/csv",
            use_container_width=True,
        )

        summary_payload = summary | {
            "output_rows": len(df),
            "alcoa_enabled": enable_alcoa,
        }

        with st.expander("Run summary / audit trail", expanded=False):
            st.json(summary_payload)
            st.caption(
                "Save this JSON alongside the CSV to maintain an at-a-glance audit trail."
            )


if __name__ == "__main__":
    main()
