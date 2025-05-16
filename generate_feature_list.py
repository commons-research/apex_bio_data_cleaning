import os
import subprocess
from typing import List

import environments_utils
import pandas as pd
from tqdm import tqdm

from constants import LINUX_BIN_PATH, MACOS_BIN_PATH
from pos_batch_mode import generate_positive_batch_mode


def generate_and_save(df: pd.DataFrame, output_dir="data", ionisation="both"):

    if ionisation == "both":
        ionisations = ["pos", "neg"]
    elif ionisation == "pos":
        ionisations = ["pos"]
    elif ionisation == "neg":
        ionisations = ["neg"]
    else:
        raise ValueError("Ionisation must be 'pos', 'neg' or 'both'. ")

    unique_cas_numbers = df["apex_cas_number"].unique()

    for cas_number in tqdm(unique_cas_numbers, desc="Analyzing data", unit="molecule"):
        cas_number_for_dir = cas_number.replace("-", "_").split(",")[0]

        for ion in ionisations:
            file_name = cas_number_for_dir + "_" + ion
            path = os.path.join(output_dir, file_name)
            os.makedirs(path, exist_ok=True)
            create_feature_list(
                df,
                cas_number,
                path,
            )


def create_feature_list(df: pd.DataFrame, cas_number: str, path: str) -> None:
    row = df[df["apex_cas_number"] == cas_number]

    # Créer le DataFrame CSV
    df_csv = {}
    df_csv["name"] = (
        row["apex_item_name"].values
        if len(row) == 1
        else [name + f"_{i}" for i, name in enumerate(row["apex_item_name"].values)]
    )
    df_csv["neutral mass"] = row["super_parent_exact_mass"].values
    df_csv["formula"] = row["super_parent_molecular_formula"].values
    df_csv["rt"] = [0] * len(row)
    df_csv["pubchem_id"] = row["pubchem_cid"].values
    df_csv["smiles"] = row["super_parent_smiles"].values
    df_csv = pd.DataFrame(df_csv)
    csv_name = "feature_list.csv"
    df_csv.to_csv(os.path.join(path, csv_name), index=False)


def main():
    # Exemple d'utilisation
    df = pd.read_csv("apex_bio_cleaned.tsv", sep="\t")
    df_plate_3 = df[df.rack_number == 3]
    generate_and_save(df_plate_3, output_dir="data", ionisation="both")


def run_mzmine(working_dir: str):
    if environments_utils.is_macos():
        mzmine_bin_path = MACOS_BIN_PATH
    elif environments_utils.is_linux():
        mzmine_bin_path = LINUX_BIN_PATH
    else:
        raise ValueError(
            "For Windows, check the mzmine documentation for batch processing."
        )

    subprocess.run(
        [
            mzmine_bin_path,
            "-b",
            "batch_file.mzbatch",
        ],
        cwd=working_dir,
    )


def write_mzbatch_pos_file(parent_dir: str, output_file: str) -> None:
    batch_file = generate_positive_batch_mode(
        os.path.join(parent_dir, "TODO"),
        os.path.join(parent_dir, "feature_list.csv"),
        os.path.join(parent_dir, "TODO"),
        os.path.join(parent_dir, "TODO"),
    )
    with open(os.path.join(parent_dir, output_file), "w") as f:
        f.write(batch_file)


if __name__ == "__main__":
    main()
