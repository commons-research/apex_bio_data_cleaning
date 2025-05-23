import os
from pathlib import Path

import pandas as pd


def is_quant_file_empty(filepath):
    """Return True if the file is empty (no rows), False otherwise."""
    try:
        df = pd.read_csv(filepath)
        return df.empty
    except pd.errors.EmptyDataError:
        return True
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return True


def check_column_exists(filepath, column_name):
    """Check if a specific column exists in the CSV file."""
    try:
        df = pd.read_csv(filepath, nrows=1)  # Only need headers
        return column_name in df.columns
    except Exception as e:
        print(f"Error checking columns in {filepath}: {e}")
        return False


def process_subdirs(root_dir):
    pos_summary = []
    neg_summary = []

    for subdir, _, _ in os.walk(root_dir):
        subdir_path = Path(subdir)
        if subdir_path.name.endswith("_pos") or subdir_path.name.endswith("_neg"):
            mode = "pos" if subdir_path.name.endswith("_pos") else "neg"
            quant_files = list(subdir_path.glob("*_quant_full.csv"))
            feature_list_files = list(subdir_path.glob("feature_list.csv"))

            if not quant_files or not feature_list_files:
                continue  # Skip if required files are missing

            quant_file = quant_files[0]
            feature_file = feature_list_files[0]

            is_empty = is_quant_file_empty(quant_file)
            has_iin_id = check_column_exists(quant_file, "ion_identities:iin_id")

            try:
                feature_df = pd.read_csv(feature_file)
                for _, row in feature_df.iterrows():
                    summary_row = {
                        "dir": str(subdir_path),
                        "name": row.get("name", ""),
                        "pubchem_id": row.get("pubchem_id", ""),
                        "smiles": row.get("smiles", ""),
                        "gnps_file_empty": is_empty,
                        "has_ion_identity_network": has_iin_id,
                    }
                    if mode == "pos":
                        pos_summary.append(summary_row)
                    else:
                        neg_summary.append(summary_row)
            except Exception as e:
                print(f"Error reading {feature_file}: {e}")

    return pos_summary, neg_summary


def write_summary(summary, filename):
    if summary:
        df = pd.DataFrame(summary)
        df.to_csv(filename, index=False)
        print(f"Wrote summary to {filename}")
    else:
        print(f"No data to write for {filename}")
