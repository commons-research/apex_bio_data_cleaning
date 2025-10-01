from typing import List

import pandas as pd
from cache_decorator import Cache
from tqdm import tqdm

from src.molecules import Molecule

DOUBLE_MOLECULES: List[str] = ["B1358", "B5054"]


@Cache()
def generate_molecules(df: pd.DataFrame) -> List[Molecule]:
    molecules: List[Molecule] = []
    for index, row in tqdm(df.iterrows(), total=len(df)):
        if row["CatalogNumber"] in DOUBLE_MOLECULES:
            molecules.append(
                Molecule(
                    apex_cas_number=row["CAS Number"],
                    apex_item_name=row["Item Name"],
                    apex_molecular_weight=row["M.w."],
                    apex_molecular_formula=row["Formula"],
                    apex_smiles=row["SMILES"].split(".")[0],
                    catalog_number=row["CatalogNumber"],
                    plate_location=row["Plate Location"],
                    rack_number_as_str=row["Rack Number"],
                )
            )

            molecules.append(
                Molecule(
                    apex_cas_number=row["CAS Number"],
                    apex_item_name=row["Item Name"],
                    apex_molecular_weight=row["M.w."],
                    apex_molecular_formula=row["Formula"],
                    apex_smiles=row["SMILES"].split(".")[1],
                    catalog_number=row["CatalogNumber"],
                    plate_location=row["Plate Location"],
                    rack_number_as_str=row["Rack Number"],
                )
            )
            continue

        mol = Molecule(
            apex_cas_number=row["CAS Number"],
            apex_item_name=row["Item Name"],
            apex_molecular_weight=row["M.w."],
            apex_molecular_formula=row["Formula"],
            apex_smiles=row["SMILES"],
            catalog_number=row["CatalogNumber"],
            plate_location=row["Plate Location"],
            rack_number_as_str=row["Rack Number"],
        )
        molecules.append(mol)
    return molecules


def save_cleaned_molecules(molecules: List[Molecule], output_file: str) -> None:
    """
    Save the cleaned molecules to a CSV file.
    """
    list_mol = [i.to_dataframe() for i in molecules]
    dataframe_mol = pd.concat(list_mol, ignore_index=True)

    if output_file.endswith(".csv"):
        dataframe_mol.to_csv(output_file, index=False)
    elif output_file.endswith(".tsv"):
        dataframe_mol.to_csv(output_file, sep="\t", index=False)
    else:
        raise ValueError("Output file must be a .csv or .tsv file.")


def main():
    excel_file = "excel_data/apex_bio_revised.xlsx"
    output_file = "apex_bio_cleaned.tsv"

    df = pd.read_excel(excel_file, sheet_name=1)

    # we want to save an unmodified version of the file to csv so that we can version control it
    df.to_csv("apex_bio_raw.tsv", index=False, sep="\t")
    molecules = generate_molecules(df)
    save_cleaned_molecules(molecules, output_file)
    return


if __name__ == "__main__":
    main()
