import argparse
import os
import shutil

import pandas as pd
from dotenv import load_dotenv
from tqdm import tqdm

from src.mzmine_utils import run_mzmine, write_mzbatch_file

load_dotenv()
if not load_dotenv():
    raise ValueError(
        "There is no path to copy the mzml files. A zenodo repo to download the data will come soon."
    )


def run(df: pd.DataFrame, output_dir: str = "data", ionisation: str = "both") -> None:

    if ionisation == "both":
        ionisations = ["pos", "neg"]
    elif ionisation == "pos":
        ionisations = ["pos"]
    elif ionisation == "neg":
        ionisations = ["neg"]
    else:
        raise ValueError("Ionisation must be 'pos', 'neg' or 'both'. ")

    unique_pubchem_cid = df["pubchem_cid"].unique()

    for pubchem_id in tqdm(unique_pubchem_cid, desc="Analyzing data", unit="molecule"):
        pubchem_id = str(pubchem_id)

        for ion in ionisations:
            file_name = pubchem_id + "_" + ion
            path = os.path.join(output_dir, file_name)
            os.makedirs(path, exist_ok=True)
            create_feature_list(
                df,
                pubchem_id,
                path,
            )
            mzml_file_name = copy_mzml_into_dir(pubchem_id, ion, path)
            write_mzbatch_file(
                path,
                mzml_file_name,
                file_name + ".mzbatch",
                file_name,
                ionisation=ion,
            )
            run_mzmine(path, file_name + ".mzbatch")


def create_feature_list(df: pd.DataFrame, pubchem_id: str, path: str) -> None:
    """Takes as input the dataframe of all compounds in the ApexBio library
    and filters the rows where the pubchem_id is the one given in input (should be mostly one single row).
    Then it saves a csv that will be used by MZmine as the target feature list.
    """
    row = df[df["pubchem_cid"] == int(pubchem_id)]

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


def copy_mzml_into_dir(pubchem_id: str, ionisation_mode: str, output_dir: str) -> str:
    path_to_copy = os.getenv("MZML_DATA_DIR")
    files_in_dir = os.listdir(path_to_copy)
    matching_files = [
        file for file in files_in_dir if pubchem_id in file and ionisation_mode in file
    ]
    if len(matching_files) == 0:
        raise ValueError(f"No file found for {pubchem_id} and {ionisation_mode}")
    elif len(matching_files) > 1:
        # if we have multiple files, we take the one with the
        # biggest file size using the os.path.getsize function
        matching_files = sorted(
            matching_files,
            key=lambda x: os.path.getsize(os.path.join(path_to_copy, x)),
            reverse=True,
        )

    file_to_copy = os.path.join(path_to_copy, matching_files[0])
    output_path = os.path.join(os.getcwd(), output_dir, matching_files[0])
    shutil.copy(file_to_copy, output_path)
    return output_path


def main():
    arg_parser = argparse.ArgumentParser(
        description="Generate feature list and run mzmine for each molecule in the ApexBio library."
    )
    arg_parser.add_argument(
        "--output_dir",
        type=str,
        default="data",
        help="Directory to save the output files.",
    )
    arg_parser.add_argument(
        "--ionisation",
        type=str,
        default="both",
        choices=["pos", "neg", "both"],
        help="Ionisation mode to process: 'pos', 'neg' or 'both'. Default is 'both'.",
    )
    arg_parser.add_argument(
        "--plate_number",
        type=int,
        choices=[-1, 1, 2, 3, 4, 5, 6, 7],
        default=-1,
        help="Plate number to process (1-7). Default is -1, which processes all plates.",
    )
    df = pd.read_csv("apex_bio_cleaned.tsv", sep="\t")
    args = arg_parser.parse_args()
    if args.plate_number == -1:
        run(df, output_dir=args.output_dir, ionisation=args.ionisation)
        return
    if args.plate_number not in [1, 2, 3, 4, 5, 6, 7]:
        raise ValueError("Plate number must be between 1 and 7.")

    df_plate = df[df.rack_number == args.plate_number]
    run(df_plate, output_dir=args.output_dir, ionisation=args.ionisation)


if __name__ == "__main__":
    main()
