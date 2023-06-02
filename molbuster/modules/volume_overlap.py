"""Module to check volume overlap between docked ligand and protein."""
from __future__ import annotations

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdShapeHelpers import ShapeProtrudeDist


def check_volume_overlap(
    mol_pred: Mol,
    mol_cond: Mol,
    ignore_hydrogens: bool = True,
    clash_cutoff: float = 0.05,
    vdw_scale: float = 0.8,
) -> dict[str, dict]:
    """Check volume overlap between ligand and protein.

    Args:
        mol_pred: Predicted molecule (docked ligand) with one conformer.
        mol_cond: Conditioning molecule (protein) with one conformer.
        clash_cutoff: Threshold for how much volume overlap is allowed. This is the share of volume of `mol_pred`
            that overlaps with `mol_cond`. Defaults to 0.05.
        vdw_scale: Scaling factor for the van der Waals radii which define the volume around each atom. Defaults to 0.8.

    Returns:
        MolBuster results dictionary.
    """
    assert isinstance(mol_pred, Mol)
    assert isinstance(mol_cond, Mol)

    overlap = 1 - ShapeProtrudeDist(mol_pred, mol_cond, vdwScale=vdw_scale, ignoreHs=ignore_hydrogens)

    results = {
        "volume_overlap": overlap,
        "no_volume_clash": overlap <= clash_cutoff,
    }

    return {"results": results}
