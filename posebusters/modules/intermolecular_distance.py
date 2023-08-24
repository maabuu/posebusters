"""Module to check intermolecular distances between ligand and protein."""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from rdkit.Chem.rdchem import GetPeriodicTable, Mol

from ..tools.protein import get_atom_type_mask

_periodic_table = GetPeriodicTable()


def check_intermolecular_distance(
    mol_pred: Mol,
    mol_cond: Mol,
    radius_type: str = "vdw",
    radius_scale: float = 1.0,
    clash_cutoff: float = 0.75,
    ignore_types: set[str] = {"hydrogens"},
    max_distance: float = 5.0,
    search_distance: float = 6.0,
) -> dict[str, Any]:
    """Check that predicted molecule is not too close and not too far away from conditioning molecule.

    Args:
        mol_pred: Predicted molecule (docked ligand) with one conformer.
        mol_cond: Conditioning molecule (protein) with one conformer.
        radius_type: Type of atomic radius to use. Possible values are "vdw" (van der Waals) and "covalent".
            Defaults to "vdw".
        radius_scale: Scaling factor for the atomic radii. Defaults to 0.8.
        clash_cutoff: Threshold for how much the atoms may overlap before a clash is reported. Defaults
            to 0.05.
        ignore_types: Which types of atoms to ignore in mol_cond. Possible values to include are "hydrogens", "protein",
            "organic_cofactors", "inorganic_cofactors", "waters". Defaults to {"hydrogens"}.
        max_distance: Maximum distance (in Angstrom) predicted and conditioning molecule may be apart to be considered
            as valid. Defaults to 5.0.

    Returns:
        PoseBusters results dictionary.
    """
    coords_ligand = mol_pred.GetConformer().GetPositions()
    coords_protein = mol_cond.GetConformer().GetPositions()

    atoms_ligand = np.array([a.GetSymbol() for a in mol_pred.GetAtoms()])
    atoms_protein_all = np.array([a.GetSymbol() for a in mol_cond.GetAtoms()])

    mask = [a.GetSymbol() != "H" for a in mol_pred.GetAtoms()]
    coords_ligand = coords_ligand[mask, :]
    atoms_ligand = atoms_ligand[mask]

    if ignore_types:
        mask = get_atom_type_mask(mol_cond, ignore_types)
        coords_protein = coords_protein[mask, :]
        atoms_protein_all = atoms_protein_all[mask]

    # get radii
    if radius_type == "vdw":
        radius_ligand = np.array([_periodic_table.GetRvdw(a) for a in atoms_ligand])
        radius_protein_all = np.array([_periodic_table.GetRvdw(a) for a in atoms_protein_all])
    elif radius_type == "covalent":
        radius_ligand = np.array([_periodic_table.GetRcovalent(a) for a in atoms_ligand])
        radius_protein_all = np.array([_periodic_table.GetRcovalent(a) for a in atoms_protein_all])
    else:
        raise ValueError(f"Unknown radius type {radius_type}. Valid values are 'vdw' and 'covalent'.")

    distances_all = _pairwise_distance(coords_ligand, coords_protein)

    # select atoms that are close to ligand to check for clash
    mask_protein = distances_all.min(axis=0) <= search_distance
    distances = distances_all[:, mask_protein]
    radius_protein = radius_protein_all[mask_protein]
    atoms_protein = atoms_protein_all[mask_protein]

    radius_sum = radius_ligand[:, None] + radius_protein[None, :]
    relative_distance = distances / radius_sum
    violations = relative_distance < 1 / radius_scale

    if distances.size > 0:
        violations[np.unravel_index(distances.argmin(), distances.shape)] = True  # add smallest distances as info
        violations[np.unravel_index(relative_distance.argmin(), relative_distance.shape)] = True
    violation_ligand, violation_protein = np.where(violations)

    # collect details around those violations in a dataframe
    details = pd.DataFrame()
    details["ligand_atom_id"] = violation_ligand
    details["protein_atom_id"] = violation_protein
    details["ligand_element"] = [atoms_ligand[i] for i in violation_ligand]
    details["protein_element"] = [atoms_protein[i] for i in violation_protein]
    details["ligand_vdw"] = [radius_ligand[i] for i in violation_ligand]
    details["protein_vdw"] = [radius_protein[i] for i in violation_protein]
    details["sum_radii"] = details["ligand_vdw"] + details["protein_vdw"]
    details["distance"] = distances[violation_ligand, violation_protein]
    details["sum_radii_scaled"] = details["sum_radii"] * radius_scale
    details["relative_distance"] = details["distance"] / details["sum_radii_scaled"]
    details["clash"] = details["relative_distance"] < clash_cutoff

    results = {
        "smallest_distance": details["distance"].min(),
        "not_too_far_away": details["distance"].min() <= max_distance,
        "num_pairwise_clashes": details["clash"].sum(),
        "no_clashes": not details["clash"].any(),
    }

    i = np.argmin(details["relative_distance"]) if len(details) > 0 else None
    most_extreme = {"most_extreme_" + c: details.loc[i][str(c)] if i is not None else pd.NA for c in details.columns}
    results = {**results, **most_extreme}

    return {"results": results, "details": details}


def _pairwise_distance(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.linalg.norm(x[:, None, :] - y[None, :, :], axis=-1)
