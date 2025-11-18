# SDF → SMILES Converter

Convert multi-record `.sdf` / `.sd` files into clean `.csv` tables containing canonical SMILES plus all SD tags. The new Streamlit UI wraps the original RDKit pipeline, adds instant previews, and optionally injects ALCOA+ metadata so that every record is attributable, legible, contemporaneous, original, accurate, complete, consistent, enduring, and available.

## Key Features
- Upload `.sdf` or `.sd` files with thousands of molecules and convert them entirely in-memory
- Canonical isomeric SMILES are calculated with RDKit and merged with every SD property tag
- Enable **ALCOA+ mode** (checkbox) to capture operator, purpose, storage plan, dataset ID, hashes, and timestamps on every record
- Metrics, interactive previews (first 500 rows), JSON audit trails, and a one-click CSV download button
- Legacy command-line interface (`sdf_2_smiles.py`) preserved for headless or scripted workflows

## Requirements
- Python 3.10+
- [RDKit](https://www.rdkit.org/docs/Install.html) (use `rdkit` on most platforms)
- Streamlit, pandas, tqdm

Install everything into a virtual environment:

```bash
pip install rdkit streamlit pandas tqdm
```

## Running the Streamlit App
1. Activate your environment and install the dependencies above.
2. From the project root run: `streamlit run streamlit_app.py`
3. Upload an `.sdf` / `.sd` file, keep the **ALCOA+ mode** box checked if you need compliance metadata, provide the required annotations, and click **Convert to CSV**.
4. Review the metrics/preview, then download the CSV and optional JSON summary for your audit trail.

### ALCOA+ Columns
When the checkbox is selected the app automatically adds the following columns to each row:

| Principle | Column(s) |
|-----------|-----------|
| Attributable | `alcoa_attributable_operator`, `alcoa_available_contact` |
| Legible | `alcoa_legible_purpose` |
| Contemporaneous | `alcoa_contemporaneous_timestamp_utc` |
| Original | `alcoa_original_source_file` |
| Accurate | `alcoa_accurate_input_sha256` |
| Complete | `alcoa_complete_dataset_id` |
| Consistent | `alcoa_consistent_processing_label` |
| Enduring | `alcoa_enduring_storage_plan` |

All timestamps use UTC ISO-8601, hashes are SHA-256 of the uploaded file, and dataset IDs default to a time-stamped token (editable in the UI). Save the JSON summary available in the **Run summary / audit trail** expander next to the CSV for full ALCOA+ coverage.

## Command-line Conversion (optional)
The original script still works for batch workflows:

```bash
python sdf_2_smiles.py input.sdf --out-csv molecules.csv --bad-sdf bad_file.sdf --summary-json run.json
```

This path provides the same conversion engine plus a JSON log but omits the Streamlit user interface.

## Troubleshooting
- **Large files**: RDKit streams molecules, so conversions only use modest memory. Still, huge uploads (>100 MB) may take time—keep the Streamlit session open until the spinner finishes.
- **Missing SMILES**: Molecules that cannot be parsed or converted are simply excluded (and counted in the metrics). Use the CLI with `--bad-sdf` to capture them for inspection.
- **ALCOA+ errors**: Operator name, contact details, and purpose are mandatory when the compliance checkbox is checked. Fill them in before converting.
- **`ModuleNotFoundError: No module named 'apt_pkg'`**: This means your Python runtime cannot load the Debian/Ubuntu `python-apt` bindings (common when using Python 3.11+).

Feel free to extend the UI layout, theme, or metadata prompts to match your lab’s SOPs.
