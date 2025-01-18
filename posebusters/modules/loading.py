"""Module to check loading of protein and ligand files."""

from __future__ import annotations

from typing import Any

from rdkit.Chem.rdchem import Mol


def check_loading(mol_pred: Any = None, mol_true: Any = None, mol_cond: Any = None) -> dict[str, dict[str, bool]]:
    """Check that molecule files were loaded successfully.

    Args:
        mol_pred: Predicted molecule. Defaults to None.
        mol_true: Ground truth molecule. Defaults to None.
        mol_cond: Conditioning molecule. Defaults to None.

    Returns:
        PoseBusters results dictionary.
    """
    results = {
        "mol_pred_loaded": isinstance(mol_pred, Mol),
        "mol_true_loaded": isinstance(mol_true, Mol),
        "mol_cond_loaded": isinstance(mol_cond, Mol),
    }
    return {"results": results}
