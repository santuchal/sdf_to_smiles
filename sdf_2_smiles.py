#!/usr/bin/env python3


from __future__ import annotations

import argparse
import csv
import json
import os
from datetime import datetime
from typing import Any, Dict, List

from rdkit import Chem
from rdkit.Chem import SDWriter
from tqdm import tqdm


def _build_alcoa_columns(
    *,
    run_timestamp_utc: str,
    source_file: str,
    metadata: Dict[str, str] | None,
) -> Dict[str, str]:
    """Return a dictionary with the ALCOA+ metadata columns."""

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


def process_sdf_records(
    *,
    sdf_path: str,
    bad_sdf_path: str | None = None,
    enforce_alcoa: bool = False,
    alcoa_metadata: Dict[str, str] | None = None,
    run_timestamp_utc: str | None = None,
    show_progress: bool = True,
) -> tuple[List[Dict[str, str]], Dict[str, Any]]:

    if not os.path.isfile(sdf_path):
        raise FileNotFoundError(f"Input SDF not found: {sdf_path}")

    total_mols = count_molecules_in_sdf(sdf_path)

    run_timestamp_utc = run_timestamp_utc or (
        datetime.utcnow().isoformat(timespec="seconds") + "Z"
    )

    rows: list[dict[str, str]] = []
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

        mol_iterator = supplier
        if show_progress:
            mol_iterator = tqdm(
                supplier, total=total_mols, desc="Converting SDF", unit="mol"
            )

        for idx, mol in enumerate(mol_iterator, start=1):
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

            row: dict[str, str] = {
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


def count_molecules_in_sdf(sdf_path: str) -> int:
    """
    Count the number of molecules in an SDF by counting '$$$$' separators.

    This is I/O-bound but much cheaper than fully parsing with RDKit and
    gives us an accurate total for a proper progress bar.
    """
    count = 0
    with open(sdf_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.strip() == "$$$$":
                count += 1
    return count


def convert_sdf_to_smiles(
    sdf_path: str,
    out_csv: str,
    bad_sdf_path: str,
    summary_json: str | None = None,
) -> None:

    if not os.path.isfile(sdf_path):
        raise FileNotFoundError(f"Input SDF not found: {sdf_path}")

    rows, summary = process_sdf_records(
        sdf_path=sdf_path,
        bad_sdf_path=bad_sdf_path,
        enforce_alcoa=False,
        alcoa_metadata=None,
        show_progress=True,
    )

    if not rows:
        # Nothing converted successfully; still create an empty CSV with header
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "record_index",
                    "smiles",
                    "source_file",
                    "processing_timestamp_utc",
                    "mol_name",
                ]
            )
    else:
        # Build a stable, complete header from all keys (Complete & Consistent)
        fieldnames = sorted({key for row in rows for key in row.keys()})
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    summary.update(
        {
            "output_csv": os.path.abspath(out_csv),
            "bad_sdf": os.path.abspath(bad_sdf_path),
        }
    )

    # Optional JSON summary to document the run (Attributable / Available / Enduring)
    if summary_json:
        with open(summary_json, "w", encoding="utf-8") as jf:
            json.dump(summary, jf, indent=2)

    # Also print a quick human-readable summary
    print("\n=== Conversion summary ===")
    print(json.dumps(summary, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert an SDF with many small molecules into a CSV with SMILES, "
            "logging non-convertible records into bad_file.sdf."
        )
    )
    parser.add_argument(
        "sdf",
        help="Path to input .sdf file",
    )
    parser.add_argument(
        "--out-csv",
        dest="out_csv",
        help="Path to output .csv file (default: <input_basename>.csv)",
    )
    parser.add_argument(
        "--bad-sdf",
        dest="bad_sdf",
        help=(
            "Path to SDF containing molecules that failed SMILES conversion "
            "(default: bad_file.sdf in the same directory as the input)."
        ),
    )
    parser.add_argument(
        "--summary-json",
        dest="summary_json",
        help="Optional JSON file summarising the run for audit/traceability.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sdf_path = args.sdf

    if args.out_csv:
        out_csv = args.out_csv
    else:
        base = os.path.splitext(os.path.basename(sdf_path))[0]
        out_csv = os.path.join(os.path.dirname(sdf_path), base + ".csv")

    if args.bad_sdf:
        bad_sdf_path = args.bad_sdf
    else:
        bad_sdf_path = os.path.join(os.path.dirname(sdf_path), "bad_file.sdf")

    convert_sdf_to_smiles(
        sdf_path=sdf_path,
        out_csv=out_csv,
        bad_sdf_path=bad_sdf_path,
        summary_json=args.summary_json,
    )


if __name__ == "__main__":
    main()
