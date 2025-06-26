import pandas as pd
from molecules import Molecule
from tqdm.auto import tqdm
from typing import List

df = pd.read_excel("./apex_bio.xlsx", sheet_name="Chemical Data")

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

#prendre tous les éléments de la liste "molécules" et 'ajouter à une nouvelle liste sous forme de dataframe 
#sauvegarder les molécules qui ont failed le cas number sous forme de fichier tsv 

molecules_dataframe = [i.to_dataframe() for i in molecules if i.failed_cas_number]
failed_cas_20 = pd.concat(molecules_dataframe, ignore_index=True)
type(failed_cas_20)

failed_cas_20.to_csv('failed_cas.tsv', sep='\t', index=False)

#sauvegarder les molécules qui ont failed le item name sous forme de fichier tsv 

molecules_dataframe_item_failed = [i.to_dataframe() for i in molecules if i.failed_item_name]
failed_iname = pd.concat(molecules_dataframe_item_failed, ignore_index=True)

failed_iname.to_csv("failed_iname.tsv", sep='\t', index=False)

#sauvegarder la(les) molécule.s qui a.ont failed le SMILES sous forme de fichier tsv 

molecules_dataframe_smiles_failed = [i.to_dataframe() for i in molecules if i.failed_smiles]
failed_smiles = pd.concat(molecules_dataframe_smiles_failed, ignore_index=True)
failed_smiles.to_csv("failed_smiles.tsv", sep="\t", index=False)

cleaned_molecules = [i.to_dataframe() for i in molecules]
cleaned_data_molecules = pd.concat(cleaned_molecules, ignore_index=True)
cleaned_data_molecules.to_csv("cleaned_molecules.tsv", sep="\t", index=False)