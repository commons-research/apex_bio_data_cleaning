import logging
from dataclasses import InitVar, dataclass, field
from typing import List, Optional, Union

import pandas as pd
import pubchempy as pcp
from cache_decorator import Cache
from dict_hash import Hashable, sha256
from molvs import Standardizer
from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem.Descriptors import ExactMolWt, MolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


@dataclass
class Molecule(Hashable):
    apex_cas_number: str
    failed_cas_number: bool = field(init=False, default=False)
    apex_item_name: str
    failed_item_name: bool = field(init=False, default=False)
    apex_molecular_weight: float
    apex_molecular_formula: str
    apex_smiles: Optional[str]
    failed_smiles: bool = field(init=False, default=False)
    super_parent_smiles: Optional[str] = field(init=False, default=None)
    rdkit_molecule: Optional[Mol] = field(init=False, default=None)
    super_parent_mol: Optional[Mol] = field(init=False, default=None)
    apex_valid_smiles: bool = field(init=False, default=False)
    pubchem_compound: Optional[Union[List[pcp.Compound], pcp.Compound]] = field(
        init=False, default=None
    )
    catalog_number: str
    plate_location: str
    rack_number_as_str: InitVar[str]
    rack_number: int = field(init=False)

    def __post_init__(self, rack_number_as_str: str):
        if isinstance(self.apex_smiles, int):
            self.apex_smiles = None
        if isinstance(self.apex_smiles, float):
            self.apex_smiles = None

        self._fetch_pubchem_compound()
        self._set_rdkit_molecule()
        self._clean_apex_smiles()

        self.rack_number = int(rack_number_as_str.split("-")[1])

    def _clean_apex_smiles(self):
        """
        Cleans the SMILES string using the molvs library.
        """
        if self.rdkit_molecule:
            standardizer = Standardizer()
            self.super_parent_mol = standardizer.super_parent(self.rdkit_molecule)
            self.super_parent_smiles = Chem.MolToSmiles(self.super_parent_mol)
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
        if not self.rdkit_molecule:
            try:
                self.rdkit_molecule = Chem.MolFromSmiles(self.apex_smiles)
                # since rdkit return None if the SMILES is not valid, we need to
                # check if the molecule is valid
                if self.rdkit_molecule:
                    self.apex_valid_smiles = True
                    return

                self.rdkit_molecule = Chem.MolFromSmiles(
                    self.pubchem_compound.canonical_smiles
                )

            except:
                try:
                    self.rdkit_molecule = Chem.MolFromSmiles(
                        self.pubchem_compound.canonical_smiles
                    )
                except:
                    self.rdkit_molecule = None

    def get_rdkit_molecule(self) -> Optional[Mol]:
        """
        Returns the RDKit molecule object.
        """
        return self.rdkit_molecule

    def get_apex_molecular_formula_clean(self) -> str:
        """
        Returns the molecular formula without the salt of anything after the `.`
        """
        return self.apex_molecular_formula.split(".")[0].strip()

    def get_apex_molecular_weight(self) -> float:
        return self.apex_molecular_weight

    def has_compound(self):
        """
        Returns True if at least one compound was found in PubChem.
        """
        return self.pubchem_compound is not None

    def has_many_compounds(self):
        if isinstance(self.pubchem_compound, pcp.Compound):
            return False
        if self.pubchem_compound is None:
            return False
        return len(self.pubchem_compound) > 1

    def get_super_parent_smiles(self) -> Optional[str]:
        """
        Returns the cleaned SMILES string.
        """
        return self.super_parent_smiles

    def _fetch_pubchem_compound(self):
        """
        Fetches the PubChem compound using the CAS number.
        """

        # we need to iterate over the CAS number, the item name and the SMILES
        # to find the compound
        # first try the CAS number
        self.pubchem_compound = self._get_compound_from_pubchem(
            self.apex_cas_number, "name"
        )
        if not self.pubchem_compound:
            self.failed_cas_number = True
            # if that fails, try the item name
            self.pubchem_compound = self._get_compound_from_pubchem(
                self.apex_item_name, "name"
            )

        if not self.pubchem_compound:
            self.failed_item_name = True
            # if that fails, try the SMILES
            if self.apex_smiles:
                self.pubchem_compound = self._get_compound_from_pubchem(
                    self.apex_smiles, "smiles"
                )
        if not self.pubchem_compound:
            self.failed_smiles = True
            print(f"Error fetching compound by SMILES.")
            self.pubchem_compound = None
            return

        # now we only keep the molecules that have a different inchikey 2D
        compounds = []
        inchikeys = set()
        for compound in self.pubchem_compound:
            if compound.inchikey[:14] not in inchikeys:
                compounds.append(compound)
                inchikeys.add(compound.inchikey[:14])

        if len(compounds) == 1:
            self.pubchem_compound = compounds[0]
        else:
            self.pubchem_compound = compounds

    def get_pubchem_inchikey(self) -> Optional[Union[str, List[str]]]:
        """
        Returns the InChIKey of the Pubchem compound.
        """
        if not self.pubchem_compound:
            return None
        if self.has_many_compounds():
            return [compound.inchikey for compound in self.pubchem_compound]
        else:
            return self.pubchem_compound.inchikey

    def matching_mass(self) -> bool:
        """
        Returns True if the mass of the molecule in the Apex Bio file
        is matching the mass of the Pubchem compound.
        """
        if not self.pubchem_compound:
            return False
        if isinstance(self.pubchem_compound, pcp.Compound):
            return (
                abs(
                    float(self.pubchem_compound.molecular_weight)
                    - self.apex_molecular_weight
                )
                < 0.2
            )
        else:
            for compound in self.pubchem_compound:
                if (
                    abs(float(compound.molecular_weight) - self.apex_molecular_weight)
                    < 0.2
                ):
                    return True
        return False

    def to_dataframe(self) -> pd.DataFrame:
        """
        Returns the molecule as a pandas dataframe with some metadata.
        """

        super_parent_exact_mass = (
            ExactMolWt(self.super_parent_mol) if self.super_parent_mol else None
        )
        super_parent_molecular_weight = (
            MolWt(self.super_parent_mol) if self.super_parent_mol else None
        )
        super_parent_formula = (
            CalcMolFormula(self.super_parent_mol) if self.super_parent_mol else None
        )

        data = {
            "catalog_number": self.catalog_number,
            "plate_location": self.plate_location,
            "rack_number": self.rack_number,
            "apex_cas_number": self.apex_cas_number,
            "apex_item_name": self.apex_item_name,
            "apex_molecular_weight": self.apex_molecular_weight,
            "apex_molecular_formula": self.apex_molecular_formula,
            "apex_smiles": self.apex_smiles,
            "failed_cas_number": self.failed_cas_number,
            "failed_item_name": self.failed_item_name,
            "failed_smiles": self.failed_smiles,
            "super_parent_smiles": self.super_parent_smiles,
            "super_parent_exact_mass": super_parent_exact_mass,
            "super_parent_molecular_weight": super_parent_molecular_weight,
            "super_parent_molecular_formula": super_parent_formula,
        }
        return pd.DataFrame(data, index=[0])

    def consistent_hash(self, use_approximation=True) -> str:
        return sha256(
            {
                "apex_cas_number": self.apex_cas_number,
                "apex_item_name": self.apex_item_name,
            },
            use_approximation=use_approximation,
        )
