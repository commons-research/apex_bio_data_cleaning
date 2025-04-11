from dataclasses import dataclass, field
from typing import List, Optional, Union

import pubchempy as pcp
from cache_decorator import Cache
from dict_hash import Hashable, sha256
from molvs import Standardizer
from rdkit import Chem
from rdkit.Chem import Mol


@dataclass
class Molecule(Hashable):
    apex_cas_number: str
    apex_item_name: str
    apex_molecular_weight: float
    apex_molecular_formula: str
    apex_smiles: Optional[str]
    _apex_cleaned_smiles: Optional[str] = field(init=False, default=None)
    _rdkit_molecule: Optional[Mol] = field(init=False, default=None)
    _pubchem_compound: Optional[Union[List[pcp.Compound], pcp.Compound]] = field(
        init=False, default=None
    )

    def __post_init__(self):
        if isinstance(self.apex_smiles, int):
            self.apex_smiles = None
        if isinstance(self.apex_smiles, float):
            self.apex_smiles = None
        self._set_rdkit_molecule()
        self._fetch_pubchem_compound()
        self._clean_apex_smiles()

    def _clean_apex_smiles(self):
        """
        Cleans the SMILES string using the molvs library.
        """
        if self.apex_smiles and self._rdkit_molecule:
            standardizer = Standardizer()
            cleaned_mol = standardizer.super_parent(self._rdkit_molecule)
            self._apex_cleaned_smiles = Chem.MolToSmiles(cleaned_mol, canonical=True)
            return

        return

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
            if self.apex_smiles is not None and isinstance(self.apex_smiles, str):
                try:
                    self._rdkit_molecule = Chem.MolFromSmiles(self.apex_smiles)
                except:
                    print(f"Error parsing SMILES: {self.apex_smiles}")
                    self._rdkit_molecule = None

    def get_rdkit_molecule(self) -> Optional[Mol]:
        """
        Returns the RDKit molecule object.
        """
        return self._rdkit_molecule

    def get_apex_molecular_formula_clean(self) -> str:
        """
        Returns the molecular formula without the salt of anything after the `.`
        """
        return self.apex_molecular_formula.split(".")[0].strip()

    def get_apex_molecular_weight(self) -> float:
        return self.apex_molecular_weight

    def has_compound(self):
        return self._pubchem_compound is not None

    def has_many_compounds(self):
        if isinstance(self._pubchem_compound, pcp.Compound):
            return False
        if self._pubchem_compound is None:
            return False
        return len(self._pubchem_compound) > 1

    def get_clean_apex_smiles(self) -> Optional[str]:
        """
        Returns the cleaned SMILES string.
        """
        if self._apex_cleaned_smiles is None:
            return None
        return self._apex_cleaned_smiles

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
                self.apex_cas_number, "name"
            )
        except KeyError:
            print(f"Error fetching compound by CAS number. Will try with Item Name")
            try:
                # if that fails, try the item name
                self._pubchem_compound = self._get_compound_from_pubchem(
                    self.apex_cas_number, "name"
                )
            except KeyError:
                print(f"Error fetching compound by item name.")
                try:
                    # if that fails, try the SMILES
                    self._pubchem_compound = self._get_compound_from_pubchem(
                        self.apex_smiles, "smiles"
                    )
                except ValueError:
                    print(f"Error fetching compound by SMILES.")
                    self._pubchem_compound = None
                    return

        # now we only keep the molecules that have a different inchikey 2D
        compounds = []
        inchikeys = set()
        for compound in self._pubchem_compound:
            if compound.inchikey[:14] not in inchikeys:
                compounds.append(compound)
                inchikeys.add(compound.inchikey[:14])

        if len(compounds) == 1:
            self._pubchem_compound = compounds[0]
        else:
            self._pubchem_compound = compounds

    def get_pubchem_inchikey(self) -> Optional[Union[str, List[str]]]:
        """
        Returns the InChIKey of the Pubchem compound.
        """
        if self._pubchem_compound is None:
            return None
        if self.has_many_compounds():
            return [compound.inchikey for compound in self._pubchem_compound]
        else:
            return self._pubchem_compound.inchikey

    def consistent_hash(self, use_approximation=False):
        return sha256(
            {
                "apex_cas_number": self.apex_cas_number,
                "apex_item_name": self.apex_item_name,
            },
            use_approximation=use_approximation,
        )

    def matching_mass(self) -> bool:
        """
        Returns True if the mass of the molecule in the Apex Bio file
        is matching the mass of the Pubchem compound.
        """
        if self._pubchem_compound is None:
            return False
        if isinstance(self._pubchem_compound, pcp.Compound):
            return (
                abs(
                    float(self._pubchem_compound.molecular_weight)
                    - self.apex_molecular_weight
                )
                < 0.2
            )
        else:
            for compound in self._pubchem_compound:
                if (
                    abs(float(compound.molecular_weight) - self.apex_molecular_weight)
                    < 0.2
                ):
                    return True
        return False
