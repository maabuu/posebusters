"""Module to check volume overlap between docked ligand and protein."""
from __future__ import annotations

import logging

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdShapeHelpers import ShapeProtrudeDist

from ..tools.molecules import delete_atoms

logger = logging.getLogger(__name__)
_inorganic_cofactors = {"Li", "Be", "Na", "Mg", "Cl", "K", "Ca", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Br", "Rb", "I"}


def check_volume_overlap(
    mol_pred: Mol,
    mol_cond: Mol,
    clash_cutoff: float = 0.05,
    vdw_scale: float = 0.8,
    ignore_hydrogens: bool = True,
    ignore_types: set[str] = set(),
) -> dict[str, dict]:
    """Check volume overlap between ligand and protein.

    Args:
        mol_pred: Predicted molecule (docked ligand) with one conformer.
        mol_cond: Conditioning molecule (protein) with one conformer.
        clash_cutoff: Threshold for how much volume overlap is allowed. This is the share of volume of `mol_pred`
            that overlaps with `mol_cond`. Defaults to 0.05.
        vdw_scale: Scaling factor for the van der Waals radii which define the volume around each atom. Defaults to 0.8.
        ignore_hydrogens: Whether to ignore hydrogens. Defaults to True.
        ignore_types: Which types of atoms to ignore. Defaults to {}. Possible values are "protein",
            "organic_cofactors", "inorganic_cofactors".

    Returns:
        PoseBusters results dictionary.
    """
    assert isinstance(mol_pred, Mol)
    assert isinstance(mol_cond, Mol)

    if ignore := set(ignore_types) - {"protein", "organic_cofactors", "inorganic_cofactors"}:
        raise ValueError(f"Ignore types {ignore} not supported.")
    indices = set()
    if "protein" in ignore_types:
        indices.update(a.GetIdx() for a in mol_cond.GetAtoms() if not a.GetPDBResidueInfo().GetIsHeteroAtom())
    if "organic_cofactors" in ignore_types:
        indices.update(a.GetIdx() for a in mol_cond.GetAtoms() if a.GetPDBResidueInfo().GetIsHeteroAtom())
    if "inorganic_cofactors" in ignore_types:
        indices.update(a.GetIdx() for a in mol_cond.GetAtoms() if a.GetSymbol() in _inorganic_cofactors)
    else:
        # hetero atoms include inorganic cofactors so remove again if not ignored
        [indices.discard(a.GetIdx()) for a in mol_cond.GetAtoms() if a.GetSymbol() in _inorganic_cofactors]
    mol_cond = delete_atoms(mol_cond, indices)

    overlap = 1 - ShapeProtrudeDist(mol_pred, mol_cond, vdwScale=vdw_scale, ignoreHs=ignore_hydrogens)

    results = {
        "volume_overlap": overlap,
        "no_volume_clash": overlap <= clash_cutoff,
    }

    return {"results": results}
