"""Provides functions for manipulating molecules."""

from __future__ import annotations

from collections.abc import Iterable
from copy import deepcopy
from logging import getLogger

import numpy as np
from rdkit import RDLogger
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
from rdkit.Chem.Lipinski import HAcceptorSmarts, HDonorSmarts
from rdkit.Chem.rdchem import AtomValenceException, Bond, Conformer, GetPeriodicTable, Mol, RWMol
from rdkit.Chem.rdMolAlign import GetBestAlignmentTransform
from rdkit.Chem.rdmolfiles import MolFromSmarts
from rdkit.Chem.rdmolops import AddHs, RemoveHs, RemoveStereochemistry, RenumberAtoms, SanitizeMol
from rdkit.Chem.rdMolTransforms import TransformConformer

logger = getLogger(__name__)


def copy_pos_to_template_mol(mol_ref: Mol, mol_pos: Mol) -> Mol:
    """Copy positions from mol_pos to mol_ref if possible."""
    mol_new = deepcopy(mol_ref)
    mol_new.RemoveAllConformers()
    RemoveStereochemistry(mol_new)

    RDLogger.DisableLog("rdApp.")
    mol_pos = AssignBondOrdersFromTemplate(refmol=mol_new, mol=mol_pos)
    RDLogger.EnableLog("rdApp.")
    assert mol_pos is not None

    # get atom mapping between molecules
    matches = mol_pos.GetSubstructMatches(mol_new, uniquify=True)
    assert len(matches) >= 1, "No matches found between molecules."

    # get positions and reorder to match atom mapping
    assert mol_pos.GetNumConformers() >= 1, "No conformers found in mol_pos."
    pos_raw = mol_pos.GetConformer().GetPositions()
    pos_new = pos_raw[matches[0], :]

    # create new conformer
    conformer = Conformer(mol_new.GetNumAtoms())
    for i, pos in enumerate(pos_new):
        conformer.SetAtomPosition(i, pos)

    # add conformer to mol
    mol_new.AddConformer(conformer)
    assert mol_new.GetNumConformers() == 1

    SanitizeMol(mol_new)
    return mol_new


def neutralize_atoms(mol: Mol) -> Mol:
    """Add and remove hydrogens to neutralize charges ignoring overall charge."""
    # https://www.rdkit.org/docs/Cookbook.html#neutralizing-charged-molecules
    # stronger than rdkit.Chem.MolStandardize.rdMolStandardize.Uncharger
    try:
        pattern = MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
        at_matches = mol.GetSubstructMatches(pattern)
        at_matches_list = [y[0] for y in at_matches]
        if len(at_matches_list) > 0:
            for at_idx in at_matches_list:
                atom = mol.GetAtomWithIdx(at_idx)
                chg = atom.GetFormalCharge()
                hcount = atom.GetTotalNumHs()
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(hcount - chg)
                atom.UpdatePropertyCache()
    except AtomValenceException:
        logger.warning("AtomValenceException raised while neutralizing molecule. Continuing with original molecule.")
        return mol
    return mol


def remove_all_charges_and_hydrogens(mol: Mol) -> Mol:
    """Remove all charges and hydrogens from molecule."""
    try:
        # rdkit keeps hydrogens that define sterochemistry so remove stereo first
        RemoveStereochemistry(mol)
        mol = RemoveHs(mol)
        # remove charges
        mol = neutralize_atoms(mol)
    except AtomValenceException:
        # from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger
        # mol = Uncharger().uncharge(mol)
        logger.warning("AtomValenceException raised while neutralizing molecule. Continuing with original molecule.")
        return mol
    return mol


def add_stereo_hydrogens(mol: Mol) -> Mol:
    """Add all hydrogens but those on primary ketimine."""
    exclude = {match[1] for match in mol.GetSubstructMatches(MolFromSmarts("[CX3]=[NH1]"))}
    atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1 if a.GetIdx() not in exclude]
    mol = AddHs(mol, onlyOnAtoms=atoms, addCoords=True)
    return mol


def remove_isotopic_info(mol: Mol) -> Mol:
    """Remove isotopic information from molecule."""
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    return mol


def get_hbond_acceptors(mol: Mol) -> set[int]:
    """Return indices of atoms that can act as hydrogen bond acceptors."""
    return {s[0] for s in mol.GetSubstructMatches(HAcceptorSmarts, uniquify=1)}


def get_hbond_donors(mol: Mol) -> set[int]:
    """Return indices of atoms that can act as hydrogen bond donors."""
    return {s[0] for s in mol.GetSubstructMatches(HDonorSmarts, uniquify=1)}


def delete_atoms(mol: Mol, indices: Iterable[int]) -> Mol:
    """Delete atoms from molecule.

    Args:
        mol: Molecule to delete atoms from.

    Returns:
        Molecule without atoms.
    """
    # delete in reverse order to avoid reindexing issues
    indices = sorted(indices, reverse=True)
    if len(indices) == 0:
        return mol
    mol = RWMol(mol)
    for index in indices:
        mol.RemoveAtom(index)
    return Mol(mol)


def _align_and_renumber(mol_true: Mol, mol_pred: Mol) -> Mol:
    """Align mol_pred to mol_true and renumbers atoms in mol_true to match mol_pred."""
    alignment = GetBestAlignmentTransform(mol_true, mol_pred, symmetrizeConjugatedTerminalGroups=False)
    # rotate and translate
    TransformConformer(mol_true.GetConformer(0), alignment[1])
    # renumber atoms
    map = [m[0] for m in sorted(alignment[2], key=lambda x: x[1])]  # reverse map
    mol_true = RenumberAtoms(mol_true, map)
    return mol_true


def _renumber_mol1_to_match_mol2(mol1: Mol, mol2: Mol) -> Mol:
    """Renumber atoms in mol1 to match mol2 (in terms of RMSD)."""
    rmsd, transform, atom_map = GetBestAlignmentTransform(mol1, mol2, reflect=False)
    is_identical_idx = [p[0] == p[1] for p in atom_map]
    n_atoms = mol1.GetNumAtoms()
    n_different = n_atoms - sum(is_identical_idx)
    if all(is_identical_idx):
        logger.info("All %d atoms are already in the same order", n_atoms)
    else:
        logger.info("Swapping %d out of %d indices", n_different, n_atoms)
        mol1 = RenumberAtoms(mol1, [m[1] for m in atom_map])
    return mol1


_periodic_table = GetPeriodicTable()


def _get_atomic_number(atomic_symbol: str):
    try:
        symbol = atomic_symbol[0].upper() + atomic_symbol[1:].lower()
        if symbol == "D":
            symbol = "H"
        return _periodic_table.GetAtomicNumber(symbol)
    except Exception:
        logger.error("Unknown atomic symbol: %s", atomic_symbol)
        return atomic_symbol


def _get_bond_length(mol: Mol, bond: Bond) -> float:
    b1, b2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    conf = mol.GetConformer()
    length = conf.GetAtomPosition(b1).Distance(conf.GetAtomPosition(b2))
    return float(length)


def _calc_angle(distances: np.ndarray, angle: tuple[int, int, int]) -> float:
    # calc alpha (opposite angle of side a in triangle with side lengths a, b, c)
    a = distances[angle[0], angle[2]]
    b = distances[angle[0], angle[1]]
    c = distances[angle[1], angle[2]]
    return float(np.degrees(np.arccos((b**2 + c**2 - a**2) / (2 * b * c))))
