"""Module to check bond lengths, bond angles, and internal clash of ligand conformations."""
from __future__ import annotations

from logging import getLogger
from typing import Any, Iterable

import numpy as np
import pandas as pd
from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdDistGeom import GetMoleculeBoundsMatrix

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


def check_geometry(
    mol_pred: Mol,
    threshold_bad_bond_length: float = 0.2,
    threshold_clash: float = 0.2,
    threshold_bad_angle: float = 0.2,
    bound_matrix_params: dict[str, Any] = bound_matrix_params,
    ignore_hydrogens: bool = True,
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

    Returns:
        MolBuster results dictionary.
    """
    molecule = mol_pred
    assert molecule.GetNumConformers() == 1

    # remove hydrogens, charges, stereochemistry tags
    molecule = remove_all_charges_and_hydrogens(molecule)

    # get bonds and angles
    bond_set = sorted(_get_bond_atom_indices(molecule))  # tuples
    angles = sorted(_get_angle_atom_indices(bond_set))  # triples
    angle_set = {(a[0], a[2]): a for a in angles}  # tuples : triples

    # distance geometry bounds, lower triangle min distances, upper triangle max distances
    bounds = GetMoleculeBoundsMatrix(molecule, **bound_matrix_params)

    # indices
    lower_triangle_idcs = np.tril_indices(molecule.GetNumAtoms(), k=-1)
    upper_triangle_idcs = (lower_triangle_idcs[1], lower_triangle_idcs[0])

    # 1,2- distances
    df_12 = pd.DataFrame()
    df_12["atom_pair"] = list(zip(*upper_triangle_idcs))  # indices have i1 < i2
    df_12["atom_types"] = [
        "--".join(tuple(molecule.GetAtomWithIdx(int(j)).GetSymbol() for j in i)) for i in df_12["atom_pair"]
    ]
    df_12["angle"] = df_12["atom_pair"].apply(lambda x: angle_set.get(x, None))
    df_12["has_hydrogen"] = [_has_hydrogen(molecule, i) for i in df_12["atom_pair"]]
    df_12["is_bond"] = [i in bond_set for i in df_12["atom_pair"]]
    df_12["is_angle"] = df_12["angle"].apply(lambda x: x is not None)
    df_12[col_lb] = bounds[lower_triangle_idcs]
    df_12[col_ub] = bounds[upper_triangle_idcs]

    # for internal clash - use 13 or 12 distances?
    # # 1,3- distances
    # df_13 = pd.DataFrame()
    # df_13["atom_triple"] = angles
    # df_13["has_hydrogen"] = [has_hydrogen(mol_pred, a) for a in angles]
    # df_13[col_lb] = [calc_angle_bounds(bounds, angle)[0] for angle in angles]
    # df_13[col_ub] = [calc_angle_bounds(bounds, angle)[1] for angle in angles]

    # calculate values and violations of bounds
    dfs_bonds, dfs_clash, dfs_angle = [], [], []
    conformer = molecule.GetConformer()
    # pairwise distances matrix (symmetric, zero diagonal)
    conf_distances = _pairwise_distance(conformer.GetPositions())
    df_12["conf_id"] = conformer.GetId()
    df_12["distance"] = conf_distances[lower_triangle_idcs]
    dfs_bonds.append(_bond_check(df_12, ignore_h=False))  # makes copy
    dfs_clash.append(_clash_check(df_12, ignore_h=False))  # makes copy
    dfs_angle.append(_angle_check(df_12, ignore_h=False))  # makes copy
    # conf_angles = [calc_angle(conf_distances, angle) for angle in angles]
    # df_13["true"] = true_flag
    # df_13["conf_id"] = conformer.GetId()
    # df_13["angle"] = conf_angles
    # dfs_angle.append(angle_check_old(df_13, ignore_h=False))  # makes copy

    df_bonds = pd.concat(dfs_bonds, ignore_index=True).reset_index(drop=True)
    df_clash = pd.concat(dfs_clash, ignore_index=True).reset_index(drop=True)
    df_angle = pd.concat(dfs_angle, ignore_index=True).reset_index(drop=True)

    sum = summarize_geometry(
        df_bonds,
        df_angle,
        df_clash,
        ignore_hydrogens=ignore_hydrogens,
        threshold_bad_bond_length=threshold_bad_bond_length,
        threshold_clash=threshold_clash,
        threshold_bad_angle=threshold_bad_angle,
    )

    bond_ratio = (sum["short_bonds"] + sum["long_bonds"]) / sum["number_bonds"]
    angle_ratio = sum["off_angles"] / sum["number_angles"] if sum["number_angles"] != 0.0 else np.nan
    clash_ratio = sum["number_clashes"] / sum["number_noncov_pairs"] if sum["number_noncov_pairs"] != 0.0 else np.nan

    passes_bonds = bond_ratio <= 0.0
    passes_angles = angle_ratio <= 0.0
    passes_clash = clash_ratio <= 0.0

    results = {
        "bond_lengths_within_bounds": passes_bonds,
        "bond_angles_within_bounds": passes_angles,
        "no_internal_clash": passes_clash,
    }

    return {"results": results, "details": {"bonds": df_bonds, "clash": df_clash, "angles": df_angle}}


def summarize_geometry(
    df_bonds: pd.DataFrame,
    df_angles: pd.DataFrame,
    df_clash: pd.DataFrame,
    threshold_bad_bond_length: float = 0.2,
    threshold_clash: float = 0.2,
    threshold_bad_angle: float = 0.2,
    ignore_hydrogens: bool = True,
) -> dict[str, float]:
    """Summarize the geometry checks of a molecule."""
    # remove hydrogen entries if not needed
    if ignore_hydrogens:
        df_bonds = df_bonds[~df_bonds["has_hydrogen"]]
        df_angles = df_angles[~df_angles["has_hydrogen"]]
        df_clash = df_clash[~df_clash["has_hydrogen"]]

    short_bonds = sum(df_bonds[col_pe] < -threshold_bad_bond_length)
    long_bonds = sum(df_bonds[col_pe] > threshold_bad_bond_length)
    number_clashes = sum(df_clash[col_bpe].abs() > threshold_clash)
    off_angles = sum(df_angles[col_bape] > threshold_bad_angle)

    return {
        "number_bonds": len(df_bonds),
        "good_bonds": len(df_bonds) - short_bonds - long_bonds,
        "short_bonds": short_bonds,
        "long_bonds": long_bonds,
        "number_noncov_pairs": len(df_clash),
        "good_noncovalent_pairs": len(df_clash) - number_clashes,
        "number_clashes": number_clashes,
        "number_angles": len(df_angles),
        "good_angles": len(df_angles) - off_angles,
        "off_angles": off_angles,
    }


def _bond_check(df: pd.DataFrame, ignore_h: bool = False) -> pd.DataFrame:
    # covalent bond length
    df = df[df["is_bond"]]
    if ignore_h:
        df = df[~df["has_hydrogen"]]
    df = df.copy()

    # bonds can be too short or too long
    df[col_pe] = df.apply(lambda x: bpe(*x[["distance", col_lb, col_ub]]), axis=1)

    return df


def _angle_check(df: pd.DataFrame, ignore_h: bool = False) -> pd.DataFrame:
    # noncovalent bonds (1,2-distances informed by van der Waals radii)
    df = df[(~df["is_bond"]) & (df["is_angle"])]
    if ignore_h:
        df = df[~df["has_hydrogen"]]
    df = df.copy()

    # angles have no direction (we do not know if larger or bigger beyond bounds)
    df[col_bape] = df.apply(lambda x: bape(*x[["distance", col_lb, col_ub]]), axis=1)

    return df


def _clash_check(df: pd.DataFrame, ignore_h: bool = False) -> pd.DataFrame:
    # noncovalent bonds (1,2-distances informed by van der Waals radii)
    df = df[(~df["is_bond"]) & (~df["is_angle"])]
    if ignore_h:
        df = df[~df["has_hydrogen"]]
    df = df.copy()

    # clash is when lower bound is violated
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
