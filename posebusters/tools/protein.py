"""Protein related functions."""

from __future__ import annotations

import logging
from collections.abc import Iterable

from rdkit.Chem.rdchem import Atom, Mol

logger = logging.getLogger(__name__)
_inorganic_cofactor_elements = {
    "Li",
    "Be",
    "Na",
    "Mg",
    "Cl",
    "K",
    "Ca",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Br",
    "Rb",
    "Mo",
    "Cd",
}
_inorganic_cofactor_ccd_codes = {
    "FES",
    "MOS",
    "PO3",
    "PO4",
    "PPK",
    "SO3",
    "SO4",
    "VO4",
}
_water_ccd_codes = {
    "HOH",
}


def get_atom_type_mask(mol: Mol, ignore_types: Iterable[str]) -> list[bool]:
    """Get mask for atoms to keep."""
    ignore_types = set(ignore_types)
    unsupported = ignore_types - {"hydrogens", "protein", "organic_cofactors", "inorganic_cofactors", "waters"}
    if unsupported:
        raise ValueError(f"Ignore types {unsupported} not supported.")

    if mol.GetAtomWithIdx(0).GetPDBResidueInfo() is None:
        logger.warning("No PDB information found. Assuming organic molecule.")

    ignore_h = "hydrogens" in ignore_types
    ignore_protein = "protein" in ignore_types
    ignore_org_cof = "organic_cofactors" in ignore_types
    ignore_inorg_cof = "inorganic_cofactors" in ignore_types
    ignore_water = "waters" in ignore_types

    return [
        _keep_atom(a, ignore_h, ignore_protein, ignore_org_cof, ignore_inorg_cof, ignore_water) for a in mol.GetAtoms()
    ]


def _keep_atom(  # noqa: PLR0913, PLR0911
    atom: Atom, ignore_h: bool, ignore_protein: bool, ignore_org_cof: bool, ignore_inorg_cof: bool, ignore_water: bool
) -> bool:
    """Whether to keep atom for given ignore flags."""
    symbol = atom.GetSymbol()
    if ignore_h and symbol == "H":
        return False

    if ignore_inorg_cof and symbol in _inorganic_cofactor_elements:
        return False

    # if loaded from PDB file, we can use the residue names and the hetero flag
    info = atom.GetPDBResidueInfo()
    if info is None:
        if ignore_org_cof:
            return False
        return True

    is_hetero = info.GetIsHeteroAtom()
    if ignore_protein and not is_hetero:
        return False

    residue_name = info.GetResidueName()
    if ignore_water and residue_name in _water_ccd_codes:
        return False

    if ignore_inorg_cof and residue_name in _inorganic_cofactor_ccd_codes:
        return False

    return True
