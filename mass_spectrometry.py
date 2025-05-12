# %%
import pandas as pd
from molecules import Molecule
from tqdm.auto import tqdm
from typing import List
from datetime import datetime
from datetime import date 

# %%
df = pd.read_excel("./apex_bio.xlsx", sheet_name="Chemical Data")

# %%
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

# %%
molecules = generate_molecules(df)

# %%
list_mol = [i.to_dataframe() for i in molecules]
dataframe_mol = pd.concat(list_mol, ignore_index=True)

# %%
from datetime import datetime
import os
import pandas as pd

path = "sub_directory_molecules"
os.makedirs(path, exist_ok=True)

def generate_and_save(df, initials="JDAN", ionisation="both"):
    current_time = datetime.now()
    date = f"{current_time.year:04}{current_time.month:02}{current_time.day:02}"

    for index, row in df.iterrows():  
        cas_number = row["apex_cas_number"].replace("-", "_").split(",")[0]

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
        df_csv.loc[0, "name"] = row["apex_item_name"]
        df_csv.loc[0, "neutral mass"] = row["super_parent_exact_mass"]
        df_csv.loc[0, "formula"] = row["super_parent_molecular_formula"]
        df_csv.loc[0, "rt"] = 0
        csv_name = f"features_molecule.csv"

        # Créer les dossiers + sauvegarder CSV dans chacun
        for folder in filenames:
            path = os.path.join("sub_directory_molecules", folder)
            os.makedirs(path, exist_ok=True)
            df_csv.to_csv(os.path.join(path, csv_name), index=False)

generate_and_save(dataframe_mol, initials="JDAN", ionisation="both")


# %%
plate_3 = dataframe_mol[dataframe_mol['rack_number'] == 3]

def generate_csv_mass_spec_pos(df, initials= "JDAN"): 
    current_time = datetime.now()
    date = f"{current_time.year:04}{current_time.month:02}{current_time.day:02}"

    csv_mass_spectrometry = pd.DataFrame(columns=["File Name", "Position"])
    for index, row in df.iterrows():  
        cas_number = row["apex_cas_number"].replace("-", "_").split(",")[0]
        location = row["plate_location"]
        location = location[0] + str(int(location[1:]))
        position = f"B:{location}"
        file_name = f"{date}_{initials}_{cas_number}_pos"

        csv_mass_spectrometry.loc[len(csv_mass_spectrometry)] = [file_name, position]

        csv_mass_spectrometry.to_csv("csv_ms_name_pos.csv", index=False)

generate_csv_mass_spec_pos(df=plate_3)


# %%
plate_3 = dataframe_mol[dataframe_mol['rack_number'] == 3]

def generate_csv_mass_spec_neg(df, initials= "JDAN"): 
    current_time = datetime.now()
    date = f"{current_time.year:04}{current_time.month:02}{current_time.day:02}"

    csv_mass_spectrometry = pd.DataFrame(columns=["File Name", "Position"])
    for index, row in df.iterrows():  
        cas_number = row["apex_cas_number"].replace("-", "_").split(",")[0]
        location = row["plate_location"]
        location = location[0] + str(int(location[1:]))
        position = f"B:{location}"
        file_name = f"{date}_{initials}_{cas_number}_neg"

        csv_mass_spectrometry.loc[len(csv_mass_spectrometry)] = [file_name, position]

        csv_mass_spectrometry.to_csv("csv_ms_name_neg.csv", index=False)

generate_csv_mass_spec_neg(df=plate_3)


