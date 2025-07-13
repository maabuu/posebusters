"""Module to check chemistry of molecules."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

import pandas as pd
from rdkit import RDLogger
from rdkit.Chem.inchi import InchiReadWriteError
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import DetectChemistryProblems, GetMolFrags, SanitizeFlags

from ..tools.inchi import get_inchi


def check_chemistry_using_rdkit(mol_pred: Mol) -> dict[str, dict[str, bool] | dict[str, bool | str] | pd.DataFrame]:
    """Check sanity of molecule using RDKit sanitization rules."""
    assert isinstance(mol_pred, Mol)
    mol = deepcopy(mol_pred)

    RDLogger.DisableLog("rdApp.*")  # type: ignore[attr-defined]
    errors = DetectChemistryProblems(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL)
    RDLogger.EnableLog("rdApp.*")  # type: ignore[attr-defined]

    messages = [error.Message() for error in errors]
    atom_indices = [e.GetAtomIndices() if hasattr(e, "GetAtomIndices") else e.GetAtomIdx() for e in errors]
    types = [error.GetType() for error in errors]

    details = pd.DataFrame(dict(messages=messages, atom_indices=atom_indices, types=types))

    results = {
        "passes_valence_checks": "AtomValenceException" not in types,
        "passes_kekulization": "AtomKekulizeException" not in types,
        "passes_rdkit_sanity_checks": len(errors) == 0,
    }
    return {"results": results, "details": details}


def check_chemistry_using_inchi(mol_pred: Mol) -> dict[str, dict[str, str | bool | None]]:
    """Check sanity of a molecule using InChI rules."""
    assert isinstance(mol_pred, Mol)
    mol = deepcopy(mol_pred)

    passes = False

    try:
        inchi = get_inchi(mol, inchi_strict=False)
        passes = True
        errors = ""
    except InchiReadWriteError as e:
        errors = str(e.args[1])
        inchi = None
    except Exception as e:
        errors = str(e)
        inchi = None

    return {"results": {"inchi_convertible": passes}, "details": dict(errors=errors, inchi=inchi)}


def check_all_atoms_connected(mol_pred: Mol) -> dict[str, dict[str, bool] | dict[str, int]]:
    """Check if all atoms in a molecule are connected."""
    assert isinstance(mol_pred, Mol)
    mol = deepcopy(mol_pred)

    num_frags = len(GetMolFrags(mol, asMols=False, sanitizeFrags=False))
    results = {"all_atoms_connected": num_frags <= 1}
    details = {"number_fragments": num_frags}

    return {"results": results, "details": details}


def check_chemistry(mol_pred: Mol) -> dict[str, Any]:
    """Check sanity using RDKit and connectedness.

    Notes:
    - Only for backward compatibility. Use `check_chemistry_using_rdkit` and `check_all_atoms_connected` separately.
    """

    rdkit_sanity = check_chemistry_using_rdkit(mol_pred)
    connectedness = check_all_atoms_connected(mol_pred)
    return {
        "results": rdkit_sanity["results"] | connectedness["results"],
        "details": rdkit_sanity["details"],
    }
