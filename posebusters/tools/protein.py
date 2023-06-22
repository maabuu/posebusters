from rdkit.Chem.rdchem import Atom, Mol
from typing import Iterable

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
    "Cd",
    "I",
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


def get_mask(mol: Mol, ignore_h: bool, ignore_types: Iterable[str]) -> list[bool]:
    """Get mask for atoms to keep."""
    ignore_types = set(ignore_types)
    if ignore := ignore_types - {"protein", "organic_cofactors", "inorganic_cofactors"}:
        raise ValueError(f"Ignore types {ignore} not supported.")

    ignore_protein = "protein" in ignore_types
    ignore_org_cof = "organic_cofactors" in ignore_types
    ignore_inorg_cof = "inorganic_cofactors" in ignore_types

    return [_keep_atom(a, ignore_h, ignore_protein, ignore_org_cof, ignore_inorg_cof) for a in mol.GetAtoms()]


def _keep_atom(atom: Atom, ignore_h: bool, ignore_protein: bool, ignore_org_cof: bool, ignore_inorg_cof: bool) -> bool:
    """Whether to keep atom for given ignore flags."""
    if ignore_h and atom.GetSymbol() == "H":
        return False

    inorganic = atom.GetSymbol() in _inorganic_cofactor_elements
    if (inorganic and ignore_inorg_cof) or (not inorganic and ignore_org_cof):
        return False

    # if loaded from PDB, we can use the residue names and the hetero flag
    info = atom.GetPDBResidueInfo()
    if info is not None:
        is_hetero = info.GetIsHeteroAtom()
        if ignore_protein and not is_hetero:
            return False

        inorganic = info.GetResidueName() in _inorganic_cofactor_ccd_codes
        if (inorganic and ignore_inorg_cof) or (not inorganic and ignore_org_cof):
            return False

    return True
