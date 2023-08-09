"""Module to check bond lengths, bond angles, and internal clash of ligand conformations."""
from __future__ import annotations

from copy import deepcopy
from logging import getLogger
from typing import Any, Iterable

import numpy as np
import pandas as pd
from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdDistGeom import GetMoleculeBoundsMatrix
from rdkit.Chem.rdmolops import SanitizeMol

from ..tools.molecules import remove_all_charges_and_hydrogens
from ..tools.stats import bape, bpe, pe

logger = getLogger(__name__)

col_lb = "lower_bound"
col_ub = "upper_bound"
col_pe = "percent_error"
col_bpe = "bound_percent_error"
col_bape = "bound_absolute_percent_error"

bound_matrix_params = {
    "set15bounds": True,
    "scaleVDW": True,
    "doTriangleSmoothing": True,
    "useMacrocycle14config": False,
}

_empty_results = {
    "results": {
        "number_bonds": np.nan,
        "shortest_bond_relative_length": np.nan,
        "longest_bond_relative_length": np.nan,
        "number_short_outlier_bonds": np.nan,
        "number_long_outlier_bonds": np.nan,
        "bond_lengths_within_bounds": np.nan,
        "number_angles": np.nan,
        "most_extreme_relative_angle": np.nan,
        "number_outlier_angles": np.nan,
        "bond_angles_within_bounds": np.nan,
        "number_noncov_pairs": np.nan,
        "shortest_noncovalent_relative_distance": np.nan,
        "number_clashes": np.nan,
        "no_internal_clash": np.nan,
    }
}


def check_geometry(
    mol_pred: Mol,
    threshold_bad_bond_length: float = 0.2,
    threshold_clash: float = 0.2,
    threshold_bad_angle: float = 0.2,
    bound_matrix_params: dict[str, Any] = bound_matrix_params,
    ignore_hydrogens: bool = True,
    sanitize: bool = True,
) -> dict[str, Any]:
    """Use RDKit distance geometry bounds to check the geometry of a molecule.

    Args:
        mol_pred: Predicted molecule (docked ligand) with exactly one conformer.
        threshold_bad_bond_length: Bond length threshold in relative percentage. 0.2 means that bonds may be up to 20%
            longer than DG bounds. Defaults to 0.2.
        threshold_clash: _description_. Defaults to 0.2.
        threshold_bad_angle: _description_. Defaults to 0.2.
        bound_matrix_params: _description_. Defaults to bound_matrix_params.
        ignore_hydrogens: Whether to ignore hydrogens. Defaults to True.
        sanitize: Sanitize molecule before running DG module (recommended). Defaults to True.

    Returns:
        PoseBusters results dictionary.
    """
    mol = deepcopy(mol_pred)
    assert mol.GetNumConformers() == 1, "Molecule must have exactly one conformer."

    # sanitize to ensure DG works or manually process molecule
    try:
        if sanitize:
            flags = SanitizeMol(mol)
            assert flags == 0, f"Sanitization failed with flags {flags}"
        else:
            # also removes stereochemistry information which we check elsewhere
            mol = remove_all_charges_and_hydrogens(mol)
    except:
        return _empty_results

    # get bonds and angles
    bond_set = sorted(_get_bond_atom_indices(mol))  # tuples
    angles = sorted(_get_angle_atom_indices(bond_set))  # triples
    angle_set = {(a[0], a[2]): a for a in angles}  # {tuples : triples}

    if len(bond_set) == 0:
        logger.warning("Molecule has no bonds.")

    # distance geometry bounds, lower triangle min distances, upper triangle max distances
    bounds = GetMoleculeBoundsMatrix(mol, **bound_matrix_params)

    # indices
    lower_triangle_idcs = np.tril_indices(mol.GetNumAtoms(), k=-1)
    upper_triangle_idcs = (lower_triangle_idcs[1], lower_triangle_idcs[0])

    # 1,2- distances
    df_12 = pd.DataFrame()
    df_12["atom_pair"] = list(zip(*upper_triangle_idcs))  # indices have i1 < i2
    df_12["atom_types"] = [
        "--".join(tuple(mol.GetAtomWithIdx(int(j)).GetSymbol() for j in i)) for i in df_12["atom_pair"]
    ]
    df_12["angle"] = df_12["atom_pair"].apply(lambda x: angle_set.get(x, None))
    df_12["has_hydrogen"] = [_has_hydrogen(mol, i) for i in df_12["atom_pair"]]
    df_12["is_bond"] = [i in bond_set for i in df_12["atom_pair"]]
    df_12["is_angle"] = df_12["angle"].apply(lambda x: x is not None)
    df_12[col_lb] = bounds[lower_triangle_idcs]
    df_12[col_ub] = bounds[upper_triangle_idcs]

    # add observed dimensions
    conformer = mol.GetConformer()
    conf_distances = _pairwise_distance(conformer.GetPositions())
    df_12["conf_id"] = conformer.GetId()
    df_12["distance"] = conf_distances[lower_triangle_idcs]

    if ignore_hydrogens:
        df_12 = df_12.loc[~df_12["has_hydrogen"], :]

    # calculate violations
    df_bonds = _bond_check(df_12)  # makes copy
    df_clash = _clash_check(df_12)  # makes copy
    df_angles = _angle_check(df_12)  # makes copy

    # bond statistics
    number_bonds = len(df_bonds)
    number_short_outlier_bonds = sum(df_bonds[col_pe] < -threshold_bad_bond_length)
    number_long_outlier_bonds = sum(df_bonds[col_pe] > threshold_bad_bond_length)
    number_valid_bonds = number_bonds - number_short_outlier_bonds - number_long_outlier_bonds

    shortest_bond_relative_length = (df_bonds["distance"] / df_bonds["lower_bound"]).min()
    longest_bond_relative_length = (df_bonds["distance"] / df_bonds["upper_bound"]).max()

    # angle statistics
    number_angles = len(df_angles)
    number_outlier_angles = sum(df_angles[col_bape] > threshold_bad_angle)
    number_valid_angles = number_angles - number_outlier_angles

    lb_extreme_angle = 2 - (df_angles["distance"] / df_angles["lower_bound"]).min()
    ub_extreme_angle = (df_angles["distance"] / df_angles["upper_bound"]).max()
    most_extreme_angle = max(lb_extreme_angle, ub_extreme_angle)

    # steric clash statistics
    number_noncov_pairs = len(df_clash)
    number_clashes = sum(df_clash[col_bpe] < -threshold_clash)
    number_valid_noncov_pairs = number_noncov_pairs - number_clashes

    shortest_noncovalent_distance = (df_clash["distance"] / df_clash["lower_bound"]).min()

    results = {
        "number_bonds": number_bonds,
        "shortest_bond_relative_length": shortest_bond_relative_length,
        "longest_bond_relative_length": longest_bond_relative_length,
        "number_short_outlier_bonds": number_short_outlier_bonds,
        "number_long_outlier_bonds": number_long_outlier_bonds,
        "bond_lengths_within_bounds": number_valid_bonds == number_bonds,
        "number_angles": number_angles,
        "most_extreme_relative_angle": most_extreme_angle,
        "number_outlier_angles": number_outlier_angles,
        "bond_angles_within_bounds": number_valid_angles == number_angles,
        "number_noncov_pairs": number_noncov_pairs,
        "shortest_noncovalent_relative_distance": shortest_noncovalent_distance,
        "number_clashes": number_clashes,
        "no_internal_clash": number_valid_noncov_pairs == number_noncov_pairs,
    }

    return {"results": results, "details": {"bonds": df_bonds, "clash": df_clash, "angles": df_angles}}


def _bond_check(df: pd.DataFrame) -> pd.DataFrame:
    # covalent bond length
    df = df.loc[df["is_bond"], :].copy()

    # bonds can be too short or too long
    df[col_pe] = df.apply(lambda x: bpe(*x[["distance", col_lb, col_ub]]), axis=1)

    return df


def _angle_check(df: pd.DataFrame) -> pd.DataFrame:
    # noncovalent bonds (1,2-distances informed by van der Waals radii)
    df = df[(~df["is_bond"]) & (df["is_angle"])].copy()

    # angles have no direction (we do not know if larger or bigger beyond bounds)
    df[col_bape] = df.apply(lambda x: bape(*x[["distance", col_lb, col_ub]]), axis=1)

    return df


def _clash_check(df: pd.DataFrame) -> pd.DataFrame:
    # noncovalent bonds (1,2-distances informed by van der Waals radii)
    df = df[(~df["is_bond"]) & (~df["is_angle"])].copy()

    # clash is only when lower bound is violated
    def _lb_pe(value, lower_bound):
        if value >= lower_bound:
            return 0.0
        return pe(value, lower_bound)

    df[col_bpe] = df.apply(lambda x: _lb_pe(*x[["distance", col_lb]]), axis=1)

    return df


def _get_bond_atom_indices(mol: Mol) -> list[tuple[int, int]]:
    bonds = []
    for bond in mol.GetBonds():
        bond_tuple = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        bond_tuple = _sort_bond(bond_tuple)
        bonds.append(bond_tuple)
    return bonds


def _get_angle_atom_indices(bonds: list[tuple[int, int]]) -> list[tuple[int, int, int]]:
    """Check all combinations of bonds to generate list of molecule angles."""
    angles = []
    bonds = list(bonds)
    for i in range(len(bonds)):
        for j in range(i + 1, len(bonds)):
            angle = _two_bonds_to_angle(bonds[i], bonds[j])
            if angle is not None:
                angles.append(angle)
    return angles


def _two_bonds_to_angle(bond1: tuple[int, int], bond2: tuple[int, int]) -> None | tuple[int, int, int]:
    set1 = set(bond1)
    set2 = set(bond2)
    all_atoms = set1 | set2
    # angle requires two bonds to share exactly one atom
    if len(all_atoms) != 3:
        return None
    # find shared atom
    shared_atom = set1 & set2
    other_atoms = all_atoms - shared_atom
    return (min(other_atoms), shared_atom.pop(), max(other_atoms))


def _sort_bond(bond: tuple[int, int]) -> tuple[int, int]:
    return (min(bond), max(bond))


def _get_torsions_atom_indices(mol: Mol) -> list[tuple[int, int, int, int]]:
    matches = mol.GetSubstructMatches(RotatableBondSmarts, uniquify=1)
    torsion_list = []
    for idx2, idx3 in matches:
        bond = mol.GetBondBetweenAtoms(idx2, idx3)
        atom2 = mol.GetAtomWithIdx(idx2)
        atom3 = mol.GetAtomWithIdx(idx3)
        # if (((atom2.GetHybridization() != HybridizationType.SP2)
        #     and (atom2.GetHybridization() != HybridizationType.SP3))
        #     or ((atom3.GetHybridization() != HybridizationType.SP2)
        #     and (atom3.GetHybridization() != HybridizationType.SP3))):
        #     continue
        torsion = None
        for b1 in atom2.GetBonds():
            if b1.GetIdx() == bond.GetIdx():
                continue
            idx1 = b1.GetOtherAtomIdx(idx2)
            for b2 in atom3.GetBonds():
                if (b2.GetIdx() == bond.GetIdx()) or (b2.GetIdx() == b1.GetIdx()):
                    continue
                idx4 = b2.GetOtherAtomIdx(idx3)
                # skip 3-membered rings
                if idx4 == idx1:
                    continue
                torsion = (idx1, idx2, idx3, idx4)
                break
            if torsion is not None:
                break
        if torsion is not None:
            torsion_list.append(torsion)
    return torsion_list


def _pairwise_distance(x):
    return np.linalg.norm(x[:, None, :] - x[None, :, :], axis=-1)


def _has_hydrogen(mol: Mol, idcs: Iterable[int]) -> bool:
    return any(_is_hydrogen(mol, idx) for idx in idcs)


def _is_hydrogen(mol: Mol, idx: int) -> bool:
    return mol.GetAtomWithIdx(int(idx)).GetAtomicNum() == 1
