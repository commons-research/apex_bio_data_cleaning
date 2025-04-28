# import some stuff 
import pandas as pd
from molecules import Molecule
from tqdm.auto import tqdm
from typing import List
from datetime import datetime
from datetime import date 

# %%
df = pd.read_excel("./apex_bio.xlsx", sheet_name="Chemical Data")

# %%
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

# convert my list in dataframe 
list_mol = [i.to_dataframe() for i in molecules]
dataframe_mol = pd.concat(list_mol, ignore_index=True)

# define a function to convert my data in a df named "File Name"

def generate_file_name (df, initials, ionisation):
    current_time = datetime.now()
    date = f"{current_time.year:04}{current_time.month:02}{current_time.day:02}"
    pubchemcid = df["pubchem_cid"]
    if ionisation == "pos":
        return pd.DataFrame({"File Name":[f"{date}_{initials}_{pubchemcid}_pos"]})
    elif ionisation == "neg":
        return pd.DataFrame({"File Name":[f"{date}_{initials}_{pubchemcid}_neg"]})
    elif ionisation == "both": 
        negative_positive = [
            f"{date}_{initials}_{pubchemcid}_pos",
            f"{date}_{initials}_{pubchemcid}_neg"
        ]
        return pd.DataFrame({"File Name": negative_positive})

# Using my function  
complete_list = []
for index, row in dataframe_mol.iterrows():
    file_name = generate_file_name(df=row, initials="JDAN", ionisation="both")
    complete_list.append(file_name)

final_df = pd.concat(complete_list).reset_index(drop=True)



