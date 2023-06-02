"""Module to check loading of protein and ligand files."""
from __future__ import annotations

from typing import Any


def check_loading(mol_pred: Any = None, mol_true: Any = None, mol_cond: Any = None) -> dict[str, dict[str, bool]]:
    """Check that molecule files were loaded successfully.

    Args:
        mol_pred: Predicted molecule. Defaults to None.
        mol_true: Ground truth molecule. Defaults to None.
        mol_cond: Conditioning molecule. Defaults to None.

    Returns:
        MolBuster results dictionary.
    """
    results = {
        "mol_pred_loaded": mol_pred is not None,
        "mol_true_loaded": mol_true is not None,
        "mol_cond_loaded": mol_cond is not None,
    }
    return {"results": results}
