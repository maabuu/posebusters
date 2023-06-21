"""Module to check volume overlap between docked ligand and protein."""
from __future__ import annotations

import logging

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdShapeHelpers import ShapeProtrudeDist

from ..tools.molecules import delete_hetatoms

logger = logging.getLogger(__name__)


def check_volume_overlap(
    mol_pred: Mol,
    mol_cond: Mol,
    clash_cutoff: float = 0.05,
    vdw_scale: float = 0.8,
    ignore_hydrogens: bool = True,
    ignore_hetatms: bool = False,
) -> dict[str, dict]:
    """Check volume overlap between ligand and protein.

    Args:
        mol_pred: Predicted molecule (docked ligand) with one conformer.
        mol_cond: Conditioning molecule (protein) with one conformer.
        clash_cutoff: Threshold for how much volume overlap is allowed. This is the share of volume of `mol_pred`
            that overlaps with `mol_cond`. Defaults to 0.05.
        vdw_scale: Scaling factor for the van der Waals radii which define the volume around each atom. Defaults to 0.8.
        ignore_hydrogens: Whether to ignore hydrogens. Defaults to True.
        ignore_hetatms: Whether to ignore heteroatoms in protein. Defaults to False.

    Returns:
        MolBusters results dictionary.
    """
    assert isinstance(mol_pred, Mol)
    assert isinstance(mol_cond, Mol)

    if ignore_hetatms:
        mol_cond = delete_hetatoms(mol_cond)

    overlap = 1 - ShapeProtrudeDist(mol_pred, mol_cond, vdwScale=vdw_scale, ignoreHs=ignore_hydrogens)

    results = {
        "volume_overlap": overlap,
        "no_volume_clash": overlap <= clash_cutoff,
    }

    return {"results": results}
