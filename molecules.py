from dataclasses import dataclass, field
from typing import List, Optional, Union

import pubchempy as pcp
from cache_decorator import Cache
from molvs import Standardizer
from rdkit import Chem
from rdkit.Chem import Mol


@dataclass
class Molecule:
    _apex_cas_number: str
    _apex_item_name: str
    _apex_molecular_weight: float
    _apex_molecular_formula: str
    _apex_smiles: Optional[str]
    _apex_cleaned_smiles: Optional[str] = field(init=False, default=None)
    _rdkit_molecule: Optional[Mol] = field(init=False, default=None)
    _pubchem_compound: Optional[List[pcp.Compound]] = field(init=False, default=None)

    def __post_init__(self):
        self._set_rdkit_molecule()
        self._fetch_pubchem_compound()

    @Cache()
    def _get_compound_from_pubchem(
        self, query: str, namespace: str = "name"
    ) -> Optional[List[pcp.Compound]]:
        data = pcp.get_compounds(query, namespace)
        return data

    def _set_rdkit_molecule(self):
        """
        Returns the RDKit molecule object.
        """
        if self._rdkit_molecule is None:
            if self._apex_smiles is not None:
                try:
                    self._rdkit_molecule = Chem.MolFromSmiles(self._apex_smiles)
                except:
                    print(f"Error parsing SMILES: {self._apex_smiles}")
                    self._rdkit_molecule = None

    def get_rdkit_molecule(self) -> Optional[Mol]:
        """
        Returns the RDKit molecule object.
        """
        if self._rdkit_molecule is None:
            self._set_rdkit_molecule()
        return self._rdkit_molecule

    def get_apex_molecular_formula_clean(self) -> str:
        """
        Returns the molecular formula without the salt of anything after the `.`
        """
        return self._apex_molecular_formula.split(".")[0].strip()

    def get_apex_molecular_weight(self) -> float:
        return self._apex_molecular_weight

    def has_compound(self):
        return self._pubchem_compound is not None

    def has_many_compounds(self):
        return len(self._pubchem_compound) > 1

    def _fetch_pubchem_compound(self):
        """
        Fetches the PubChem compound using the CAS number.
        """
        if self._pubchem_compound is not None:
            return

        # we need to iterate over the CAS number, the item name and the SMILES
        # to find the compound
        # first try the CAS number
        try:
            self._pubchem_compound = self._get_compound_from_pubchem(
                self._apex_cas_number, "name"
            )
        except KeyError:
            print(f"Error fetching compound by CAS number. Will try with Item Name")
            try:
                # if that fails, try the item name
                self._pubchem_compound = self._get_compound_from_pubchem(
                    self._apex_item_name, "name"
                )
            except KeyError:
                print(f"Error fetching compound by item name.")
                try:
                    # if that fails, try the SMILES
                    self._pubchem_compound = self._get_compound_from_pubchem(
                        self._apex_smiles, "smiles"
                    )
                except ValueError:
                    print(f"Error fetching compound by SMILES.")
                    self._pubchem_compound = None

    def get_pubchem_inchikey(self) -> Optional[Union[str, List[str]]]:
        """
        Returns the InChIKey of the Pubchem compound.
        """
        if self._pubchem_compound is None:
            return None
        if self.has_many_compounds():
            return [compound.inchikey for compound in self._pubchem_compound]
        else:
            return self._pubchem_compound[0].inchikey
