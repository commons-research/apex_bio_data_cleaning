from dataclasses import dataclass
from typing import List, Optional

import pubchempy as pcp
from cache_decorator import Cache


@dataclass
class Molecule:
    _apex_cas_number: str
    _apex_item_name: str
    _apex_molecular_weight: float
    _apex_molecular_formula: str
    _apex_smiles: Optional[str]
    _pubchem_compound: Optional[List[pcp.Compound]] = None

    @Cache()
    def _get_compound_from_pubchem(
        self, query: str, namespace: str = "name"
    ) -> Optional[List[pcp.Compound]]:
        data = pcp.get_compounds(query, namespace)
        return data

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

    def fetch_pubchem_compound(self) -> "Molecule":
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
            return self
        except Exception as e:
            print(f"Error fetching compound by CAS number. Will try with Item Name")
            try:
                # if that fails, try the item name
                self._pubchem_compound = self._get_compound_from_pubchem(
                    self._apex_item_name, "name"
                )
                return self
            except Exception as e:
                print(f"Error fetching compound by item name: {e}")
                try:
                    # if that fails, try the SMILES
                    self._pubchem_compound = self._get_compound_from_pubchem(
                        self._apex_smiles, "smiles"
                    )
                    return self
                except Exception as e:
                    print(f"Error fetching compound by SMILES: {e}")
                    self._pubchem_compound = None
                    return self
