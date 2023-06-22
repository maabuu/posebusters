"""Module to check volume overlap between docked ligand and protein."""
from __future__ import annotations

import logging

import numpy as np
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdShapeHelpers import ShapeProtrudeDist

from ..tools.molecules import delete_atoms
from ..tools.protein import get_mask

logger = logging.getLogger(__name__)


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
        clash_cutoff: Threshold for how much volume overlap is allowed. This is the maximum share of volume of
            `mol_pred` allowed to overlap with `mol_cond`. Defaults to 0.05.
        vdw_scale: Scaling factor for the van der Waals radii which define the volume around each atom. Defaults to 0.8.
        ignore_hydrogens: Whether to ignore hydrogens. Defaults to True.
        ignore_types: Which types of atoms to ignore. Possible values are "protein", "organic_cofactors",
            "inorganic_cofactors". Defaults to none (empty set).

    Returns:
        PoseBusters results dictionary.
    """
    assert isinstance(mol_pred, Mol)
    assert isinstance(mol_cond, Mol)

    keep_mask = np.array(get_mask(mol_cond, ignore_hydrogens, ignore_types))
    if keep_mask.sum() == 0:
        return {"results": {"volume_overlap": np.nan, "no_volume_clash": True}}
    if keep_mask.sum() < len(keep_mask):
        mol_cond = delete_atoms(mol_cond, np.where(~keep_mask)[0].tolist())

    overlap = 1 - ShapeProtrudeDist(mol_pred, mol_cond, vdwScale=vdw_scale, ignoreHs=ignore_hydrogens)

    results = {
        "volume_overlap": overlap,
        "no_volume_clash": overlap <= clash_cutoff,
    }

    return {"results": results}
