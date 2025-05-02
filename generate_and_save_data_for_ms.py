# Import some tools 
import pandas as pd
from molecules import Molecule
from tqdm.auto import tqdm
from typing import List
from datetime import datetime
from datetime import date 

# read excel
df = pd.read_excel("./apex_bio.xlsx", sheet_name="Chemical Data")

from cache_decorator import Cache

@Cache()
def generate_molecules(df): 
    molecules: List[Molecule] = []
    for index, row in tqdm(df.iterrows(), total=len(df)):
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
 
molecules = generate_molecules(df)

list_mol = [i.to_dataframe() for i in molecules]
dataframe_mol = pd.concat(list_mol, ignore_index=True)

from datetime import datetime
import os
import pandas as pd

path = "sub_directory_molecules"
os.makedirs(path, exist_ok=True)

def generate_and_save(df, initials="JDAN", ionisation="both"):
    current_time = datetime.now()
    date = f"{current_time.year:04}{current_time.month:02}{current_time.day:02}"

    for i in range(len(df)):
        row = df.iloc[[i]]
        cas_number = row["apex_cas_number"].values[0].replace("-", "_").split(",")[0]

        # Générer noms de fichiers
        filenames = []
        if ionisation == "pos":
            filenames = [f"{date}_{initials}_{cas_number}_pos"]
        elif ionisation == "neg":
            filenames = [f"{date}_{initials}_{cas_number}_neg"]
        elif ionisation == "both":
            filenames = [
                f"{date}_{initials}_{cas_number}_pos",
                f"{date}_{initials}_{cas_number}_neg"
            ]

        # Créer le DataFrame CSV
        df_csv = pd.DataFrame(columns=["name", "neutral mass", "formula", "rt"])
        df_csv.loc[0, "name"] = row["apex_item_name"].values[0]
        df_csv.loc[0, "neutral mass"] = row["super_parent_exact_mass"].values[0]
        df_csv.loc[0, "formula"] = row["super_parent_molecular_formula"].values[0]
        df_csv.loc[0, "rt"] = 0
        csv_name = f"features_molecule_{i+1}.csv"

        # Créer les dossiers + sauvegarder CSV dans chacun
        for folder in filenames:
            path = os.path.join("sub_directory_molecules", folder)
            os.makedirs(path, exist_ok=True)
            df_csv.to_csv(os.path.join(path, csv_name), index=False)

generate_and_save(dataframe_mol, initials="JDAN", ionisation="both")



