"""Module to check flatness of ligand substructures."""
from __future__ import annotations

from typing import Any, Iterable

import numpy as np
from numpy import ndarray as Array
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolFromSmarts

_flat = {
    "aromatic_5_membered_rings_sp2": "[ar5^2]1[ar5^2][ar5^2][ar5^2][ar5^2]1",
    "aromatic_6_membered_rings_sp2": "[ar6^2]1[ar6^2][ar6^2][ar6^2][ar6^2][ar6^2]1",
    "trigonal_planar_double_bonds": "[#6;X3;D2,D3;^2]=[#6;X3;D2,D3;^2]",
}


def check_flatness(
    mol_pred: Mol, threshold_flatness: float = 0.1, flat_systems: dict[str, str] = _flat
) -> dict[str, Any]:
    """Check whether substructures of molecule are as flat as chemistry predicts.

    Args:
        mol_pred: Predicted molecule (docked ligand) with exactly one conformer.
        threshold_flatness: _description_. Defaults to 0.1.
        flat_systems: _description_. Defaults to None.

    Returns:
        MolBuster results dictionary.
    """
    planar_groups = []
    types = []
    for flat_system, smarts in flat_systems.items():
        match = MolFromSmarts(smarts)
        atom_groups = list(mol_pred.GetSubstructMatches(match))
        planar_groups += atom_groups
        types += [flat_system] * len(atom_groups)

    # calculate distances to plane and check threshold
    coords = [_get_coords(mol_pred, group) for group in planar_groups]
    max_distances = [float(_get_distances_to_plane(X).max()) for X in coords]
    flatness_passes = [bool(d <= threshold_flatness) for d in max_distances]

    details = {
        "type": types,
        "planar_group": planar_groups,
        "max_distance": max_distances,
        "flatness_passes": flatness_passes,
    }

    results = {"flatness_passes": all(flatness_passes)}

    return {"results": results, "details": details}


def _get_distances_to_plane(X: Array) -> Array:
    """Get distances of points X to their common plane."""
    # center points X in R^(n x 3)
    X = X - X.mean(axis=0)
    # singular value decomposition
    _, _, V = np.linalg.svd(X)
    # last vector in V is normal vector to plane
    n = V[-1]
    # distances to plane are projections onto normal
    d = np.dot(X, n)
    return d


def _get_coords(mol: Mol, indices: Iterable[int]) -> Array:
    return np.array([mol.GetConformer().GetAtomPosition(i) for i in indices])
