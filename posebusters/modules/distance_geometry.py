"""Module to check bond lengths, bond angles, and internal clash of ligand conformations."""

from __future__ import annotations

from collections.abc import Iterable
from copy import deepcopy
from logging import getLogger
from typing import Any

import numpy as np
import pandas as pd
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdDistGeom import GetMoleculeBoundsMatrix
from rdkit.Chem.rdmolops import SanitizeMol

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

col_n_bonds = "number_bonds"
col_shortest_bond = "shortest_bond_relative_length"
col_longest_bond = "longest_bond_relative_length"
col_n_short_bonds = "number_short_outlier_bonds"
col_n_long_bonds = "number_long_outlier_bonds"
col_n_good_bonds = "number_valid_bonds"
col_bonds_result = "bond_lengths_within_bounds"
col_n_angles = "number_angles"
col_extremest_angle = "most_extreme_relative_angle"
col_n_bad_angles = "number_outlier_angles"
col_n_good_angles = "number_valid_angles"
col_angles_result = "bond_angles_within_bounds"
col_n_noncov = "number_noncov_pairs"
col_closest_noncov = "shortest_noncovalent_relative_distance"
col_n_clashes = "number_clashes"
col_n_good_noncov = "number_valid_noncov_pairs"
col_clash_result = "no_internal_clash"

_empty_results = {
    col_n_bonds: np.nan,
    col_shortest_bond: np.nan,
    col_longest_bond: np.nan,
    col_n_short_bonds: np.nan,
    col_n_long_bonds: np.nan,
    col_bonds_result: np.nan,
    col_n_angles: np.nan,
    col_extremest_angle: np.nan,
    col_n_bad_angles: np.nan,
    col_angles_result: np.nan,
    col_n_noncov: np.nan,
    col_closest_noncov: np.nan,
    col_n_clashes: np.nan,
    col_clash_result: np.nan,
}


def check_geometry(  # noqa: PLR0913, PLR0915
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
        mol_pred: Predicted molecule (docked ligand). Only the first conformer will be checked.
        threshold_bad_bond_length: Bond length threshold in relative percentage. 0.2 means that bonds may be up to 20%
            longer than DG bounds. Defaults to 0.2.
        threshold_clash: Threshold for how much overlap constitutes a clash. 0.2 means that the two atoms may be up to
            80% of the lower bound apart. Defaults to 0.2.
        threshold_bad_angle: Bond angle threshold in relative percentage. 0.2 means that bonds may be up to 20%
            longer than DG bounds. Defaults to 0.2.
        bound_matrix_params: Parameters passe to RDKit's GetMoleculeBoundsMatrix function.
        ignore_hydrogens: Whether to ignore hydrogens. Defaults to True.
        sanitize: Sanitize molecule before running DG module (recommended). Defaults to True.

    Returns:
        PoseBusters results dictionary.
    """
    mol = deepcopy(mol_pred)
    results = _empty_results.copy()

    if mol.GetNumConformers() == 0:
        logger.warning("Molecule does not have a conformer.")
        return {"results": results}

    if mol.GetNumAtoms() == 1:
        logger.warning(f"Molecule has only {mol.GetNumAtoms()} atoms.")
        results[col_angles_result] = True
        results[col_bonds_result] = True
        results[col_clash_result] = True
        return {"results": results}

    # sanitize to ensure DG works or manually process molecule
    try:
        if sanitize:
            flags = SanitizeMol(mol)
            assert flags == 0, f"Sanitization failed with flags {flags}"
    except Exception:
        return {"results": results}

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
    df_bonds = _bond_check(df_12)
    df_clash = _clash_check(df_12)
    df_angles = _angle_check(df_12)

    # bond statistics
    results[col_n_bonds] = len(df_bonds)
    results[col_n_short_bonds] = sum(df_bonds[col_pe] < -threshold_bad_bond_length)
    results[col_n_long_bonds] = sum(df_bonds[col_pe] > threshold_bad_bond_length)
    results[col_n_good_bonds] = results[col_n_bonds] - results[col_n_short_bonds] - results[col_n_long_bonds]
    results[col_bonds_result] = results[col_n_good_bonds] == results[col_n_bonds]
    results[col_shortest_bond] = (df_bonds["distance"] / df_bonds["lower_bound"]).min()
    results[col_longest_bond] = (df_bonds["distance"] / df_bonds["upper_bound"]).max()

    # angle statistics
    results[col_n_angles] = len(df_angles)
    results[col_n_bad_angles] = sum(df_angles[col_bape] > threshold_bad_angle)
    results[col_n_good_angles] = results[col_n_angles] - results[col_n_bad_angles]
    results[col_angles_result] = results[col_n_good_angles] == results[col_n_angles]
    lb_extreme_angle = 2 - (df_angles["distance"] / df_angles["lower_bound"]).min()
    ub_extreme_angle = (df_angles["distance"] / df_angles["upper_bound"]).max()
    results[col_extremest_angle] = max(lb_extreme_angle, ub_extreme_angle)

    # steric clash statistics
    results[col_n_noncov] = len(df_clash)
    results[col_n_clashes] = sum(df_clash[col_bpe] < -threshold_clash)
    results[col_n_good_noncov] = results[col_n_noncov] - results[col_n_clashes]
    results[col_clash_result] = results[col_n_good_noncov] == results[col_n_noncov]
    results[col_closest_noncov] = (df_clash["distance"] / df_clash["lower_bound"]).min()

    return {"results": results, "details": {"bonds": df_bonds, "clash": df_clash, "angles": df_angles}}


def _bond_check(df: pd.DataFrame) -> pd.DataFrame:
    # bonds can be too short or too long
    df = df[df["is_bond"]].copy()
    df[col_pe] = df.apply(lambda x: bpe(*x[["distance", col_lb, col_ub]]), axis=1)
    return df


def _angle_check(df: pd.DataFrame) -> pd.DataFrame:
    # angles have no direction (we do not know if larger or bigger beyond bounds)
    df = df[(~df["is_bond"]) & (df["is_angle"])].copy()
    df[col_bape] = df.apply(lambda x: bape(*x[["distance", col_lb, col_ub]]), axis=1)
    return df


def _clash_check(df: pd.DataFrame) -> pd.DataFrame:
    # clash is only when lower bound is violated
    df = df[(~df["is_bond"]) & (~df["is_angle"])].copy()

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
    # angle requires two bonds to share exactly one atom, that is we must have 3 atoms
    if len(all_atoms) != 3:  # noqa: PLR2004
        return None
    # find shared atom
    shared_atom = set1 & set2
    other_atoms = all_atoms - shared_atom
    return (min(other_atoms), shared_atom.pop(), max(other_atoms))


def _sort_bond(bond: tuple[int, int]) -> tuple[int, int]:
    return (min(bond), max(bond))


def _pairwise_distance(x):
    return np.linalg.norm(x[:, None, :] - x[None, :, :], axis=-1)


def _has_hydrogen(mol: Mol, idcs: Iterable[int]) -> bool:
    return any(_is_hydrogen(mol, idx) for idx in idcs)


def _is_hydrogen(mol: Mol, idx: int) -> bool:
    return mol.GetAtomWithIdx(int(idx)).GetAtomicNum() == 1
