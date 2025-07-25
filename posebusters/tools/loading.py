"""Provides functions for loading molecules from files."""

from __future__ import annotations

import logging
from collections.abc import Generator, Iterable
from pathlib import Path

from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import (
    MolFromMol2File,
    MolFromMolBlock,
    MolFromMolFile,
    MolFromPDBFile,
    MolFromSmiles,
    SDMolSupplier,
)
from rdkit.Chem.rdmolops import AddHs, AssignStereochemistryFrom3D, Cleanup, SanitizeMol
from rdkit.rdBase import LogToPythonLogger

from ..tools.logging import CaptureLogger

LogToPythonLogger()
logger = logging.getLogger(__name__)


def safe_load_mol(path: Path | Mol, load_all: bool = False, **load_params) -> Mol | None:
    """Load one molecule from a file, optionally adding hydrogens and assigning bond orders.

    Args:
        path: Path to file containing molecule.

    Returns:
        Molecule object or None if loading failed.
    """
    if isinstance(path, Mol):
        return path

    path = Path(path)
    if not path.exists():
        logger.warning("File %s does not exist", path)
        return None

    try:
        with CaptureLogger():
            mol = _load_mol(path, load_all=load_all, **load_params)
        return mol
    except Exception as exception:
        logger.warning("Could not load molecule from %s with error: %s", path, exception)
    return None


def safe_supply_mols(
    path: Path, load_all=True, sanitize=True, indices: Iterable[int] | None = None, **load_params
) -> Generator[Mol | None, None, None]:
    """Supply molecules from a file, optionally adding hydrogens and assigning bond orders.

    Args:
        path: Path to file containing molecules.

    Returns:
        Molecule object or None if loading failed.
    """
    if isinstance(path, Mol):
        yield path
        return None

    path = Path(path)
    if not path.exists():
        logger.warning("File %s not found.", path)
        return None
    if path.suffix == ".sdf":
        pass
    elif path.suffix in {".mol", ".mol2"}:
        yield safe_load_mol(path, sanitize=True, **load_params)
        return None
    else:
        raise ValueError(f"Molecule file {path} has unknown format. Only .sdf, .mol and .mol2 are supported.")

    with SDMolSupplier(str(path), sanitize=False, strictParsing=True) as supplier:
        for index in range(len(supplier)) if indices is None else indices:
            mol_clean = _process_mol(supplier[index], sanitize=sanitize, **load_params)
            if mol_clean is not None:
                mol_clean.SetProp("_Path", str(path))
                mol_clean.SetProp("_Index", str(index))
            yield mol_clean


def get_num_mols(path: Path | Mol) -> int:
    """Get number of molecules in a file. Only supports multiple molecules in SDF files."""

    if isinstance(path, Mol):
        return 1

    path = Path(path)
    if not path.exists():
        return 0
    if path.suffix == ".sdf":
        with SDMolSupplier(str(path), sanitize=False, strictParsing=False) as supplier:
            return len(supplier)
    return 1


def _load_mol(  # noqa: PLR0913
    path: Path,
    load_all=False,
    sanitize=False,
    removeHs=False,
    strictParsing=False,
    proximityBonding=False,
    cleanupSubstructures=True,
    index: int | None = None,
    **params,
) -> Mol | None:
    """Load molecule(s) from a file, picking the right RDKit function."""

    if path.suffix == ".sdf":
        if load_all:
            mol = _load_and_combine_mols(path, sanitize=False, removeHs=removeHs, strictParsing=strictParsing)
        elif index is not None:
            with SDMolSupplier(str(path), sanitize=False, removeHs=removeHs, strictParsing=strictParsing) as supplier:
                mol = supplier[index]
        else:
            mol = MolFromMolFile(str(path), sanitize=False, removeHs=removeHs, strictParsing=strictParsing)
    elif path.suffix == ".mol2":
        # MolFromMol2File only loads only first molecule from mol2 file
        with open(path) as file:
            if load_all and sum(ln.strip().startswith("@<TRIPOS>MOLECULE") for ln in file.readlines()) > 1:
                logger.error("Cannot load multiple molecules from mol2 file, only loading first.")
        mol = MolFromMol2File(str(path), sanitize=False, removeHs=removeHs, cleanupSubstructures=cleanupSubstructures)
    elif path.suffix == ".pdb":
        mol = MolFromPDBFile(str(path), sanitize=False, removeHs=removeHs, proximityBonding=proximityBonding)
    elif path.suffix == ".mol":
        # .mol files only contain one molecule
        with open(path) as file:
            block = "".join(file.readlines()).strip() + "\nM  END"
            mol = MolFromMolBlock(block, sanitize=False, removeHs=removeHs, strictParsing=strictParsing)
    else:
        raise ValueError(f"Unknown file type {path.suffix}")

    if mol is not None:
        mol.SetProp("_Path", str(path))

    mol = _process_mol(mol, sanitize=sanitize, **params)

    return mol


def _load_and_combine_mols(path: Path, sanitize=True, removeHs=True, strictParsing=True) -> Mol | None:
    """Load mols from SDF file and combine into one molecule."""
    with SDMolSupplier(str(path), sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing) as supplier:
        # warning: combines molecules without checking identity or atom order
        mol = next(supplier)
        while mol is None and supplier.atEnd() is False:
            mol = next(supplier)
        for mol_next in supplier:
            if mol_next is not None:
                mol.AddConformer(mol_next.GetConformer(), assignId=True)
        return mol
    return None


def _process_mol(  # noqa: PLR0913
    mol: Mol | None,
    smiles: str | None = None,
    cleanup=False,
    sanitize=True,
    add_hs=False,
    assign_stereo=True,
) -> Mol:
    """Process a molecule, optionally adding hydrogens and assigning bond orders."""
    if mol is None:
        raise ValueError("Could not load molecule.")
    if smiles is not None:
        mol = _assign_bond_order(mol, smiles)
    if cleanup:
        mol = _cleanup(mol)
    if sanitize:
        mol = _sanitize(mol)
    if add_hs:
        mol = _add_hydrogens(mol)
    if assign_stereo:
        mol = _assign_stereo(mol)
    return mol


def _assign_bond_order(mol: Mol, smiles) -> Mol:
    template = MolFromSmiles(smiles)
    mol = AssignBondOrdersFromTemplate(template, mol)
    if mol is None:
        raise ValueError("Could not assign bond orders to molecule.")
    return mol


def _cleanup(mol: Mol) -> Mol:
    Cleanup(mol)
    if mol is None:
        raise ValueError("Could not cleanup molecule.")
    return mol


def _sanitize(mol: Mol) -> Mol:
    flags = SanitizeMol(mol)
    if mol is None or flags != 0:
        raise ValueError("Could not sanitize molecule.")
    return mol


def _add_hydrogens(mol: Mol) -> Mol:
    mol = AddHs(mol, addCoords=True)
    if mol is None:
        raise ValueError("Could not add hydrogens to molecule.")
    return mol


def _assign_stereo(mol: Mol) -> Mol:
    AssignStereochemistryFrom3D(mol)
    if mol is None:
        raise ValueError("Could not assign stereo to molecule.")
    return mol
