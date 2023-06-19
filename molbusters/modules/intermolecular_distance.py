"""Module to check intermolecular distances between ligand and protein."""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from rdkit.Chem.rdchem import GetPeriodicTable, Mol

from ..tools.molecules import get_hbond_acceptors, get_hbond_donors

_periodic_table = GetPeriodicTable()


def check_intermolecular_distance(
    mol_pred: Mol,
    mol_cond: Mol,
    vdw_scale: float = 0.8,
    clash_cutoff: float = 0.05,
    ignore_hydrogens: bool = True,
    ignore_hydrogen_bonds: bool = False,
    max_distance: float = 5.0,
) -> dict[str, Any]:
    """Calculate pairwise intermolecular distances between ligand and protein atoms.

    Args:
        mol_pred: Predicted molecule (docked ligand) with one conformer.
        mol_cond: Conditioning molecule (protein) with one conformer.
        vdw_scale: Scaling factor for the van der Waals radii which define the volume around each atom. Defaults to 0.8.
        clash_cutoff: Threshold for how much scaled van der Waal radii may overlap before a clash is reported. Defaults
            to 0.05.
        ignore_hydrogens: Whether to ignore hydrogens. Defaults to True.
        ignore_hydrogen_bonds: Whether to ignore hydrogens which may be involved in hydrogen bonding. Defaults to False.
        max_distance: Maximum distance (in Angstrom) between ligand and protein to be considered as valid. Defaults to
            5.0.

    Returns:
        MolBusters results dictionary.
    """
    coords_ligand = mol_pred.GetConformer().GetPositions()
    coords_protein = mol_cond.GetConformer().GetPositions()

    atoms_ligand = _get_atoms(mol_pred)
    atoms_protein = _get_atoms(mol_cond)

    vdw_ligand = _get_vdw_radii(mol_pred)
    vdw_protein_all = _get_vdw_radii(mol_cond)

    if ignore_hydrogens:
        heavy_atoms_mask_ligand = np.asarray(atoms_ligand != "H")
        heavy_atoms_mask_protein = np.asarray(atoms_protein != "H")

        coords_ligand = coords_ligand[heavy_atoms_mask_ligand, :]
        coords_protein = coords_protein[heavy_atoms_mask_protein, :]

        atoms_ligand = atoms_ligand[heavy_atoms_mask_ligand]
        atoms_protein = atoms_protein[heavy_atoms_mask_protein]

        vdw_ligand = vdw_ligand[heavy_atoms_mask_ligand]
        vdw_protein_all = vdw_protein_all[heavy_atoms_mask_protein]

    elif ignore_hydrogen_bonds:
        raise NotImplementedError
        # hbonds = _is_hbond(mol_pred, mol_cond, zip(violation_ligand, violation_protein))
        # details["has_h"] = details.apply(_has_h, axis=1)
        # details["hbond_donor_acceptor"] = hbonds

    distances_all = _pairwise_distance(coords_ligand, coords_protein)

    # minimum distance
    smallest_distance = distances_all.min()

    # select atoms that are close to ligand to check for clash
    mask_protein = distances_all.min(axis=0) <= 4.5
    distances = distances_all[:, mask_protein]
    vdw_protein = vdw_protein_all[mask_protein]

    vdw_sum = vdw_ligand[:, None] + vdw_protein[None, :]
    violation_ligand, violation_protein = np.where(distances - vdw_sum < 0)

    # collect details around those possible violations
    details = pd.DataFrame()
    details["ligand_atom_id"] = violation_ligand
    details["protein_atom_id"] = violation_protein
    details["ligand_element"] = [atoms_ligand[i] for i in violation_ligand]
    details["protein_element"] = [atoms_protein[i] for i in violation_protein]
    details["ligand_vdw"] = [vdw_ligand[i] for i in violation_ligand]
    details["protein_vdw"] = [vdw_protein[i] for i in violation_protein]
    details["sum_vdw"] = details["ligand_vdw"] + details["protein_vdw"]
    details["distance"] = distances[violation_ligand, violation_protein]

    # identify clashes after scaling vdw and applying cutoff
    details["sum_vdw_scaled"] = details["sum_vdw"] * vdw_scale
    details["gap"] = details["distance"] - details["sum_vdw_scaled"]
    details["clash"] = details["gap"] < -clash_cutoff

    results = {
        "smallest_distance": smallest_distance,
        "not_too_far_away": max_distance >= smallest_distance,
        "num_pairwise_clashes": details["clash"].sum(),
        "min_gap": details["gap"].min(),
        "no_clashes": not details["clash"].any(),
    }

    return {"results": results, "details": details}


def _pairwise_distance(x, y):
    return np.linalg.norm(x[:, None, :] - y[None, :, :], axis=-1)


def _get_vdw_radii(mol: Mol) -> np.ndarray:
    return np.array([_periodic_table.GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()])


def _get_atoms(mol: Mol) -> np.ndarray:
    return np.array([a.GetSymbol() for a in mol.GetAtoms()])


def _has_h(row):
    if row["ligand_element"] == "H" or row["protein_element"] == "H":
        return True
    return False


def _is_hbond(ligand, protein, atom_pairs) -> list[bool]:
    ligand_hbond_accept = get_hbond_acceptors(ligand)
    ligand_hbond_donor = get_hbond_donors(ligand)
    protein_hbond_accept = get_hbond_acceptors(protein)
    protein_hbond_donor = get_hbond_donors(protein)

    out = []
    for ligand_i, protein_i in atom_pairs:
        if ligand_i in ligand_hbond_accept and protein_i in protein_hbond_donor:
            out.append(True)
        elif ligand_i in ligand_hbond_donor and protein_i in protein_hbond_accept:
            out.append(True)
        else:
            out.append(False)
    return out
