import os
from datetime import datetime

import pandas as pd
from tqdm import tqdm


def generate_and_save(df: pd.DataFrame, output_dir="data", ionisation="both"):

    unique_cas_numbers = df["apex_cas_number"].unique()

    for cas_number in tqdm(unique_cas_numbers):
        row = df[df["apex_cas_number"] == cas_number]
        cas_number = row["apex_cas_number"].values[0].replace("-", "_").split(",")[0]

        # Générer noms de fichiers
        filenames = []
        if ionisation == "pos":
            filenames = [f"{cas_number}_pos"]
        elif ionisation == "neg":
            filenames = [f"{cas_number}_neg"]
        elif ionisation == "both":
            filenames = [
                f"{cas_number}_pos",
                f"{cas_number}_neg",
            ]

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
        df_csv = pd.DataFrame(df_csv)
        csv_name = "feature_list.csv"

        # Créer les dossiers + sauvegarder CSV dans chacun
        for folder in filenames:
            path = os.path.join(output_dir, folder)
            os.makedirs(path, exist_ok=True)
            df_csv.to_csv(os.path.join(path, csv_name), index=False)


def main():
    # Exemple d'utilisation
    df = pd.read_csv("apex_bio_cleaned.tsv", sep="\t")
    generate_and_save(df, output_dir="data", ionisation="both")


if __name__ == "__main__":
    main()
